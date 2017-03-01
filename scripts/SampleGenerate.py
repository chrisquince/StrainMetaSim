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


def main(argv):
    
    parser = argparse.ArgumentParser()
    
    parser = argparse.ArgumentParser(description='Simulate a metagenomic data set')
    
    parser.add_argument("config_file", help="json config file")
    
    parser.add_argument('output_dir', metavar='DIR', help='Output directory')

    parser.add_argument("coverage_file", help="strain coverages")
    
    parser.add_argument('-n', '--output-name', metavar='PATH', help='Output file base name', required=True)

    parser.add_argument('-s', '--sample-no', metavar='PATH', help='Sample number to generate', required=True)

    parser.add_argument('--art-path', default='art_illumina', help='Path to ART executable [default: art_illumina]')

    parser.add_argument('--log', default='MetagenomeSim.log', type=argparse.FileType('w'), help='Log file name')
    
    parser.add_argument('-r','--random_seed',default=23724839, type=int, 
        help=("specifies seed for numpy random number generator defaults to 23724839"))
    
    args = parser.parse_args()
    
    #create new random state
    prng = RandomState(args.random_seed)
        
    base_name = os.path.join(args.output_dir, args.output_name)

    r1_tmp = os.path.join(args.output_dir, '{0}1.fq'.format(TMP_OUTPUT))
    r2_tmp = os.path.join(args.output_dir, '{0}2.fq'.format(TMP_OUTPUT))
    
    #import ipdb; ipdb.set_trace()
    
    with open(args.config_file) as config_file:    
        config = json.load(config_file)

    #pprint(config)

    strainCov = defaultdict(dict)
 
    for line in open(args.coverage_file):
        line = line.rstrip()
        (sample_no, species, strain, cov, freq) = line.split("\t")
        strainCov[int(sample_no)][strain] = float(cov)
          
    
    #now can generate actual reads one sample at a time    
    NSamples = int(config['no_Samples'])
    insert_len = int(config['reads']['insert_length'])
    insert_sd = int(config['reads']['insert_sd'])
    read_len = int(config['reads']['length'])
    
    child_seeds = prng.randint(LOW_SEED_VALUE, HIGH_SEED_VALUE, NSamples).tolist()
    
    n = int(args.sample_no)
    
    r1_final = '{0}.{1}.r1.fq.gz'.format(base_name, n)
    r2_final = '{0}.{1}.r2.fq.gz'.format(base_name, n)
        
    r1_tmp = os.path.join(args.output_dir, '{0}_{1}1.fq'.format(TMP_OUTPUT, n))
    r2_tmp = os.path.join(args.output_dir, '{0}_{1}2.fq'.format(TMP_OUTPUT, n))
    out_tmp = os.path.join(args.output_dir, '{0}_{1}'.format(TMP_OUTPUT, n))
    output_R1 = gzip.open(r1_final, 'wb')
    output_R2 = gzip.open(r2_final, 'wb')

    #loop strains
    for strainId, ssCov in strainCov[n].iteritems():
        print '\tRequesting {0:.4f} coverage for {1}'.format(ssCov, strainId)
        totalLen = 0    
        seq_tmp = os.path.join(args.output_dir, str(strainId) + TMP_INPUT)
            
        FastaFile = open(seq_tmp, 'rU')
        for rec in SeqIO.parse(FastaFile, 'fasta'):
            totalLen +=len(rec)
        FastaFile.close()

        # iteration target for ART
        try:
            subprocess.check_call([args.art_path,
                                           '-p',   # paired-end sequencing
                                           '-na',  # no alignment file
                                           '-rs', str(child_seeds[n]),
                                           '-m', str(insert_len),
                                           '-s', str(insert_sd),
                                           '-l', str(read_len),
                                           '-f', str(ssCov),
                                           '-i', seq_tmp,
                                           '-o', out_tmp],
                                          stdout=args.log, stderr=args.log)

        except OSError as e:
            print "There was an error executing \"art_illumina\"."
            print "Check that it is either on your PATH or specify it at runtime."
            raise e
        except subprocess.CalledProcessError as e:
            print e
            raise e
                    
        # count generated reads
        r1_n = 0
        for seq in SeqIO.parse(r1_tmp, 'fastq'):
            r1_n += 1

        r2_n = 0
        for seq in SeqIO.parse(r2_tmp, 'fastq'):
            r2_n += 1

        assert r1_n == r2_n, 'Error: failed to generate an equal number of fwd and rev reads'

        effective_cov = read_len * (r1_n + r2_n) / float(totalLen)
        print '\tGenerated {0} pairs for {1}, {2:.3f} coverage'.format(r1_n, strainId, effective_cov)

        if r1_n != r2_n:
            print 'Error: paired-end counts do not match {0} vs {1}'.format(r1_n, r2_n)
            sys.exit(1)

        with open(r1_tmp) as f:
            output_R1.write(f.read())
                
        with open(r2_tmp) as f:
            output_R2.write(f.read())

        os.remove(r1_tmp)
        os.remove(r2_tmp)
    
    
if __name__ == "__main__":
    main(sys.argv[1:])
