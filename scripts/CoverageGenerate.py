import sys, getopt
import os
import pandas as p
import numpy as np
import scipy.stats as ss
import scipy as sp
import scipy.misc as spm
import math
import argparse
import cPickle
import logging
import json
import glob
import re
import subprocess
import gzip
import atexit

from pprint import pprint
from Bio import SeqIO
from random import shuffle
from collections import defaultdict
from numpy.random import RandomState

TMP_INPUT = 'seq.tmp'
TMP_OUTPUT = 'reads.tmp'

# low-high seeds, giving 5M values
LOW_SEED_VALUE = 1000000
HIGH_SEED_VALUE = 6000000

class Species_Profiles():
    """Generates species proportions from log-normal"""

    def __init__(self, randomState, speciesMap, no_Species,no_Samples, speciesDistParams, readParams):

        self.randomState = randomState

        self.mean_log_Mean = speciesDistParams['mean_log_Mean']
    
        self.sd_log_Mean = speciesDistParams['sd_log_Mean']
    
        self.k_log_Sd = speciesDistParams['k_log_Sd']
    
        self.theta_log_Sd = speciesDistParams['theta_log_Sd']

        self.beta = speciesDistParams['beta']
        
        self.alpha = speciesDistParams['alpha']
        
        self.no_Species = no_Species
        
        self.no_Samples = no_Samples
        
        self.no_Reads = readParams['no_Reads']
        
        self.read_length = readParams['length']
        
        self.speciesMap = speciesMap
        
        self.freqArray = np.zeros((no_Species,no_Samples))
        
    def generateProfiles(self):

        self.log_Means = self.randomState.normal(self.mean_log_Mean, self.sd_log_Mean,self.no_Species)
        
        self.log_Sds = self.randomState.gamma(self.k_log_Sd,self.theta_log_Sd,self.no_Species)

        self.betaFilter = self.randomState.beta(self.beta,self.beta, self.no_Species) 
        
        #generate log abundances
        logAbundArray = np.zeros((self.no_Species,self.no_Samples))
        speciesFilter = np.zeros((self.no_Species,self.no_Samples),dtype=np.int)
            
        for s in range(self.no_Species):
            logAbundArray[s,:] = self.randomState.lognormal(self.log_Means[s],self.log_Sds[s],self.no_Samples)
            speciesFilter[s,:] = np.random.binomial(1, self.betaFilter[s], size=self.no_Samples)
    
        #logAbundArray = logAbundArray*speciesFilter
        logAbundArray = logAbundArray.transpose()
        row_sums = logAbundArray.sum(axis=1)
        self.freqArray = logAbundArray / row_sums[:, np.newaxis]
        
    def addStrainCoverages(self, speciesStrainSelect):
    
        for speciesId, strain_map in speciesStrainSelect.iteritems():
            s = self.speciesMap[speciesId]
            
            sFreq = self.freqArray[:,s]
        
            nStrains = len(strain_map)
             
            strainFreq = np.zeros((nStrains,self.no_Samples))
            for s in range(self.no_Samples):
                strainFreq[:,s] = self.randomState.dirichlet([self.alpha]*nStrains)
            
            idx = 0
            for strainId, [strainDir,nSeqs,all_records] in strain_map.iteritems():
                
                totalLen = sum([len(rec.seq) for rec in all_records])
                sstrainFreq = strainFreq[idx,:]*sFreq
                factor = (float(self.no_Reads)*float(self.read_length))/float(totalLen)
                strainCov = sstrainFreq*factor
                
                strain_map[strainId].append(totalLen)
                strain_map[strainId].append(strainCov)
                strain_map[strainId].append(sstrainFreq)
                                
                idx = idx + 1
        

def main(argv):
    
    parser = argparse.ArgumentParser()
    
    parser = argparse.ArgumentParser(description='Simulate a metagenomic data set')
    
    parser.add_argument("config_file", help="json config file")
    
    parser.add_argument('output_dir', metavar='DIR', help='Output directory')
    
    parser.add_argument('--log', default='MetagenomeSim.log', type=argparse.FileType('w'), help='Log file name')

    parser.add_argument('--coverage-out', metavar='FILE', default='coverage.tsv',
                        help='Output file for simulated genome coverage table', required=False)
    
    parser.add_argument('--select-out', metavar='FILE', default='select.tsv',
                        help='Selected strain for each species', required=False)
    
    parser.add_argument('-s','--random_seed',default=23724839, type=int, 
        help=("specifies seed for numpy random number generator defaults to 23724839"))
    
    args = parser.parse_args()
    
    @atexit.register
    def close_cov():
        coverage_file.close()
    
    #create new random state
    prng = RandomState(args.random_seed)

    r1_tmp = os.path.join(args.output_dir, '{0}1.fq'.format(TMP_OUTPUT))
    r2_tmp = os.path.join(args.output_dir, '{0}2.fq'.format(TMP_OUTPUT))
    
    #import ipdb; ipdb.set_trace()
    
    with open(args.config_file) as config_file:    
        config = json.load(config_file)

    #pprint(config)
    
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    select_file = open(os.path.join(args.output_dir, args.select_out), 'w')
    
    NSpecies = len(config['species'])
    speciesStrainSelect = {}

    speciesMap = {}
    strainFiles = {}
    s = 0
    #loop species
    for species in config['species']:
        #open scg file
    
        species_dir = species['dir']
    
        m = re.search(r"Strain_(\d+)", species_dir)

        species_id = m.group(1)
        
        speciesMap[species_id] = s
        s = s + 1
        
        hg_ident_file = species_dir + "/IdentHG.csv"
        scg_file = species_dir + "/strain_map_scg.csv"
        strain_map_file = species_dir + "/strain_map.csv"
    
        compStrains = defaultdict(dict)
        
        for line in open(hg_ident_file):
            
            line = line.rstrip()
            
            (taxaA, taxaB, identPH, identG) = line.split(",")  
    
            compStrains[taxaA[2:]][taxaB[2:]] = 100. - float(identPH)
            compStrains[taxaB[2:]][taxaA[2:]] = 100. - float(identPH)
    
        scg = p.read_csv(scg_file, header=0, index_col=0)
        scg_sum = (scg == 1).sum(axis=1)
        
        #loop through strains in species
        strain_map = {}
        seq_no = {}
        
        for line in open(strain_map_file):
            line = line.rstrip()
        
            (strainDir, taxaId) = line.split(",")
                   
            nSeqs = 0
            all_records = []
            for fasta_file in glob.glob(species_dir + "/" + strainDir+"/*.fna"):
                seq_records = list(SeqIO.parse(fasta_file, "fasta"))
                seq_index = SeqIO.index(fasta_file, 'fasta')
                all_records = all_records + seq_records
                nSeqs += len(seq_records)
            
            strain_map[taxaId] = [strainDir,nSeqs,all_records]
            
        strain_filter = []
        for taxaId, [strainDir,nSeqs,all_records] in strain_map.iteritems():
            if scg_sum[strainDir] >= config['min_scg'] and nSeqs < config['max_seqs']:
                strain_filter.append(taxaId) 
    #            print "Filtered: " + taxaId + ",Dir: " + strainDir + ",Scgs: " + str(scg_sum[strainDir]) + ",Seqs: " + str(nSeqs)
        #    else:
     #           print "Removed: " + taxaId + ",Dir: " + strainDir + ",Scgs: " + str(scg_sum[strainDir]) + ",Seqs: " + str(nSeqs)    
        strain_select = []
 
        
        #select strains in a rather stupid iterative way
        print "No. strains for :" + species_dir + " = " + str(len(strain_filter)) + " require: " + species['nStrains'] 
        nAttempts = 0
        if int(species['nStrains']) > 1:

            while nAttempts < 10000:
                strain_select = prng.choice(strain_filter, size=int(species['nStrains']),replace=False)
            
                distH = []
                
                strain_select2 = strain_select[:]
                for x in strain_select:
                    for y in strain_select2:

                        if x != y:
                            distH.append(compStrains[x][y])
            
                minH = min(distH)

                maxH = max(distH)

                if minH > float(config['min_div']) and maxH < float(config['max_div']):
                    break 
            
                nAttempts+=1

            if nAttempts == 10000:
        
                if len(strain_filter) > 0:
                    strain_select = [strain_filter.pop(0)]
                else:
                    print "No strains for :" + species_dir
                    sys.exit()
            
                prng.shuffle(strain_filter)
    
                for strain in strain_filter:
            
                    if len(strain_select) == int(species['nStrains']):
                        break 
       
                    distH = [compStrains[strain][x] for x in strain_select]
            
                    minH = min(distH)
            
                    maxH = max(distH)
            
                    if minH > float(config['min_div']) and maxH < float(config['max_div']):
                        strain_select.append(strain)
        else:
            prng.shuffle(strain_filter)
            strain_select = [strain_filter[0]] 
        nSelected =  len(strain_select)
        
        speciesStrainSelect[species_id] = dict((k, strain_map[k]) for k in strain_select)
        
        for strainId, [strainDir,nSeqs,all_records] in speciesStrainSelect[species_id].iteritems():
            select_file.write('{0}\t{1}\t{2}\t{3}\n'.format(species_id, strainId, nSeqs, strainDir))
            seq_tmp = os.path.join(args.output_dir, str(strainId) + TMP_INPUT)
            SeqIO.write(all_records, seq_tmp, 'fasta')
            strainFiles[strainId] = seq_tmp
        
        if nSelected < int(species['nStrains']):
            print species_dir + "," + str(nSelected) + "," + species['nStrains']   
    
    select_file.close()
    
       #now generate species profiles
    speciesProfiles = Species_Profiles(prng,speciesMap,config['no_Species'],config['no_Samples'],config['species_dist_Params'],config['reads'])
    speciesProfiles.generateProfiles()
    speciesProfiles.addStrainCoverages(speciesStrainSelect)
    
    #now can generate actual reads one sample at a time
    
    NSamples = int(config['no_Samples'])
    insert_len = int(config['reads']['insert_length'])
    insert_sd = int(config['reads']['insert_sd'])
    read_len = int(config['reads']['length'])
    
    child_seeds = prng.randint(LOW_SEED_VALUE, HIGH_SEED_VALUE, NSamples).tolist()

    coverage_file = open(os.path.join(args.output_dir, args.coverage_out), 'w')
    
    for n in range(NSamples):
        

        #loop species
        for speciesId, strain_map in speciesStrainSelect.iteritems():
            s = speciesMap[speciesId]
            
            for strainId, [strainDir,nSeqs,all_records,totalLen, strain_Cov, strain_Freq] in strain_map.iteritems():
    
                ssCov = strain_Cov[n]
                
                ssFreq = strain_Freq[n]
                
                coverage_file.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(n, speciesId, strainId, ssCov, ssFreq))

                
    
    
if __name__ == "__main__":
    main(sys.argv[1:])
