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

from operator import mul, div, eq, ne, add, ge, le, itemgetter
from itertools import izip
from numpy import array, log, exp
from scipy.special import gammaln
from scipy.optimize import minimize_scalar

from scipy.stats import chi2

from sklearn.metrics import roc_curve, auc, accuracy_score
from collections import defaultdict
from collections import Counter
from collections import Set

def variableTau(tau):
    """Calculates positions with variable bases"""
    N = tau.shape[0]
    G = tau.shape[1]
    variable_tau = np.zeros((N), dtype=bool)
    for v in range(N):
        diff = False
        id0 = np.where(tau[v,0,:] == 1)[0]
        for g in range(1,G):
            idg = np.where(tau[v,g,:] == 1)[0]
            if(idg[0] != id0[0]):
                diff = True 
            
            variable_tau[v] = diff
            
    return variable_tau

def compSND(tau1,tau2):
    G1 = tau1.shape[1]
    G2 = tau2.shape[1]
        
    snd = np.zeros((G1,G2),dtype=np.int)
    N = tau1.shape[0]
    for g in range(G1):
        #snd[g,g] = 0
        for h in range(G2):
            overlap = 0.0;
            for v in range(N):
                idg = np.where(tau1[v,g,:] == 1)[0]
                idh = np.where(tau2[v,h,:] == 1)[0]
                if(idg[0] == idh[0]):
                    overlap += 1 
                
            snd[g,h] = N - overlap
                
    return snd

def main(argv):

    #import ipdb; ipdb.set_trace()

    parser = argparse.ArgumentParser()
    parser.add_argument("sel_var_file", help="predicted variant positions")
        
    parser.add_argument("tau_file", help="known variants")
    
    parser.add_argument("--gene_info", help="output per gene",
                                action="store_true")

    args = parser.parse_args()
    
    sel_var_file = args.sel_var_file
    tau_file = args.tau_file

    sel_var = p.read_csv(sel_var_file, header=0, index_col=0)
    
    tau = p.read_csv(tau_file, header=0, index_col=0)

    idx = 0
    taudict = defaultdict(dict)
    genes = set()
    gene_true = Counter()
    for index, row in tau.iterrows():
        posn = int(row['Position'])
        taudict[index][posn] = idx
        gene_true[index] += 1
        idx = idx + 1
        genes.add(index)
    pindex = 0
    comp = 0
    tp = 0
    fp = 0
    gene_pred = Counter()
    gene_correct = Counter()

    for index, row in sel_var.iterrows():
        pposn = int(row['Position'])
        comp = comp + 1
        gene_pred[index] += 1
        genes.add(index)
        if pposn in taudict[index]: 
            tp = tp + 1
            gene_correct[index] += 1

    if args.gene_info:
        print "Gene,Var,PredVar,Recall,Prec"
        for gene in genes:
            if gene_pred[gene] > 0:
                precn = float(gene_correct[gene])/float(gene_pred[gene])
            else:
                precn = "na"
            if gene_true[gene] > 0:
                recall = float(gene_correct[gene])/float(gene_true[gene])
            else:
                recall = "na"
            print gene + "," + str(gene_true[gene]) + "," + str(gene_pred[gene]) + "," + str(recall) + "," + str(precn)

    if comp > 0:
        precn = float(tp)/float(comp)
    else:
        precn = "na"
    
    if idx > 0:
        recall = float(tp)/float(idx)
    else:
        recall = "na"


    print str(tp) + "," + str(idx) + "," + str(comp) + "," + str(recall) + "," + str(precn)
if __name__ == "__main__":
    main(sys.argv[1:])
