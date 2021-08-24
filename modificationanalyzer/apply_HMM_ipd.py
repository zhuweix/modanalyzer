import pysam
import argparse
import os
import itertools
import pickle
import numpy as np
from pomegranate import HiddenMarkovModel, DiscreteDistribution, State

def cigartuple_position(cigartuple):

    cigarlist = []
    for tp, num in cigartuple:
        if tp == 0: # M, unmethylated
            cigarlist.extend(['u'] * (num - 1))
        elif tp == 1: # I, methylated
            cigarlist.append('m')

    return cigarlist

def cigarMask(cigarlist, chrom, readstart, motifdict):

	maskedlist = [value if (index + readstart) in motifdict[chrom] else np.nan for index,value in enumerate(cigarlist)]

	return maskedlist

def HMMSamWriter(bamfile, nucmodel, motifpositiondict, outfile):
    tp_map = {
        'L': 'M',
        'N': 'D'
    }
    with pysam.AlignmentFile(bamfile, 'rb') as bam:
        with pysam.AlignmentFile(outfile, 'wb', header=bam.header) as outbam:
            for read in bam:

                cigartuple = read.cigartuples
                start = read.reference_start + 1 # 0-based to 1-based
                chrom = read.reference_name
                methlist = cigartuple_position(cigartuple)
                methlist =  cigarMask(methlist, chrom, start, motifpositiondict)

                logp, path = nucmodel.viterbi(methlist)

                statepath = [state.name[0] for idx, state in path[1::]]

                cigar = ''.join([str(len(list(g)))+tp_map[key] for key,g in (itertools.groupby(statepath))])
                
                read.cigarstring = cigar
                outbam.write(read)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bamfile')
    parser.add_argument('-o', '--output')        
    parser.add_argument('-t', '--trainedmodel')    
    parser.add_argument('-m', '--motifpositionfile')
    args = parser.parse_args()


    modeldecoy = HiddenMarkovModel('Nuc Finder')
    nucmodel = modeldecoy.from_json(args.trainedmodel)
    
    with open(args.motifpositionfile, 'rb') as filep:
        motifpositiondict = pickle.load(filep)
    
    HMMSamWriter(
        bamfile=args.bamfile, 
        nucmodel=nucmodel, 
        motifpositiondict=motifpositiondict, 
        outfile=args.output)