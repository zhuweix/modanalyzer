import pysam
import argparse
import os
import numpy as np
import pickle
from pomegranate import HiddenMarkovModel, DiscreteDistribution, State

def nucHMM():

    model = HiddenMarkovModel('Nuc Finder')

    #these are the frequencies of methylation vs unmethylation for the windows
    d1 = DiscreteDistribution({'m': 0.60,'u': 0.40})
    d2 = DiscreteDistribution({'m': 0.01,'u': 0.99})


    #the linker is highly mehtylated d1
    l0 = State( d1, name='L0' )
    l1 = State( d1, name='L1' )
    l2 = State( d1, name='L2' )#can loop as many times as it wants
    l3 = State( d1, name='L3' )
    l4 = State( d1, name='L4' )
    l5 = State( d1, name='L5' )



    #each of these represents a 10bp window inside the nucleosome each with low methylation
    n0 = State( d2, name='N0')
    n1 = State( d2, name='N1')
    n2 = State( d2, name='N2')#can loop as many times as it wants
    n3 = State( d2, name='N3')
    n4 = State( d2, name='N4')
    n5 = State( d2, name='N5')


    #add states to model
    model.add_states( [l0,l1,l2,l3,l4,l5,n0,n1,n2,n3,n4,n5] )

    model.add_transitions( model.start, [l0,l1,l2,l3,l4,l5,n0,n1,n2,n3,n4,n5], [0.02,0.02,0.02,0.02,0.02,0.02,0.18,0.18,0.18,0.18,0.18,0.18] )

    model.add_transition(  l0, l1, 0.99)

    model.add_transition(  l1, l2, 0.99)
    model.add_transition(  l2, l2, 0.1)
    model.add_transition(  l2, l3, 0.9)
    model.add_transition(  l3, l4, 0.99)
    model.add_transition(  l4, l5, 0.99)
    model.add_transition(  l5, n0, 0.99)

    model.add_transition(  n0, n1, 0.99)
    model.add_transition(  n1, n2, 0.99)
    model.add_transition(  n2, n2, 0.97)
    model.add_transition(  n2, n3, 0.03)
    model.add_transition(  n3, n4, 0.99)
    model.add_transition(  n4, n5, 0.99)

    model.add_transition(  n5, l0, 0.99)

    model.bake()

    return model


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


def nucHMMLearn(nucmodel, bamfile, chromtolearnon, motifpositiondict, output, threads, max_iterations=600):
    print('Getting read names')

    learnstrings = []
    with pysam.AlignmentFile(bamfile, 'rb') as bam:
        for read in bam:
            cigarstring = read.cigarstring
            cigartuple = read.cigartuples
            start = read.reference_start + 1 # 0-based to 1-based
            chrom = read.reference_name
            if chrom != chromtolearnon:
                continue
            cigarlist = cigartuple_position(cigartuple)
            cigarlist =  cigarMask(
                cigarlist,
                chrom,
                start,
                motifpositiondict)

            learnstrings.append(cigarlist)

    print('Number of reads:',len(learnstrings))
    print('Fitting model')
    nucmodel.fit(learnstrings, verbose=True, n_jobs=threads, max_iterations=max_iterations)

    with open(output, 'w') as filep:
        filep.write(nucmodel.to_json())
    return nucmodel


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bamfile')
    parser.add_argument('-o', '--outputmodel')        
    parser.add_argument('-m', '--motifpositionfile')
    parser.add_argument('-t', '--threads')
    parser.add_argument('-c','--chromtolearnon')
    parser.add_argument('-i', '--max_iterations', default=600)
    args = parser.parse_args()

    nucmodel = nucHMM()
    
    with open(args.motifpositionfile, 'rb') as filep:
        motifpositiondict = pickle.load(filep)
        
    nucmodel = nucHMMLearn(
        nucmodel=nucmodel,
        bamfile=args.bamfile,
        chromtolearnon=args.chromtolearnon,
        motifpositiondict=motifpositiondict,
        output=args.outputmodel,
        threads=int(args.threads),
        max_iterations=int(args.max_iterations))