import pysam
import hashlib
from collections import Counter
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
import re
import multiprocessing
import numpy as np
import pickle

def get_ZMWs(bamfile:str, coveragecutoff:int, outfile=None):
    '''Obtain the name of ZMW with coverage >= coveragecutoff
    
    Args:
        bamfile (str): name of the bamfile
        coveragecutoff (int): minimal coverage of subreads per ZMW
        outfile (str): name of the output file with the names of filtered ZMWs (None: no output)
    Returns:
        list: The name of filterd ZMWs
    '''
    print ('Analyzing ZMW read from bam files')
    ZMWholes = []
    with pysam.AlignmentFile(bamfile, "rb", check_sq=False) as bam:
        for read in bam:
            try:
                holeNumber = read.query_name.split('/')[1] # Get ZMW name: Name of ZMW read: {movieName}/{ZMWNumber}/{qStart}_{qEnd}
                if holeNumber:
                    ZMWholes.append(holeNumber)
            except:
                continue

    print ('Counting ZMW holes with >=%d coverage (subreads)' % coveragecutoff)
    ZMWholes = Counter(ZMWholes)
    filter_ZMWholes = []
    for hole, count in ZMWholes.items():
        if count >= coveragecutoff:
            filter_ZMWholes.append(hole)
    print ('%d ZMWs counted' % len(filter_ZMWholes))
    if not (outfile is None):
        with open(outfile, 'w') as filep:
            filep.write('\n'.join(filter_ZMWholes))
            print('Filtered ZMWs are stores in {}'.format(outfile))
    return filter_ZMWholes

def corrret_read_tag(bamfile: str, outfile: str):
    '''Correct the RGID of the bamfile.

    According to Current PacBio Bam file specification,
    the @RG ID tag of the read should be:
            RGID_STRING := md5(movieName + "//" + readType))[:8]
        where:
            movieName := @RG PU tag
            readType := @RG DS tag
    Note: For the purpose of PacBio IPDSummary to predict modification, 
    only the first @RG in the header is kept.

    Args:
        bamfile (str): name of the input bam file
        outfile (str): name of the output bam file
    Returns:
        None. The output is stored in the outfile

    '''
    with pysam.AlignmentFile(bamfile, 'rb', check_sq=False) as inbam:
        # correct header
        header = dict(inbam.header)
        moviename = str(header['RG'][0]['PU'])
        readtype = str(header['RG'][0]['DS'])   
        rgid =  str(hashlib.md5((moviename + '//' + readtype).encode()).hexdigest()[:8])
        header['RG'][0]['ID'] = rgid
        header['RG'] = [header['RG'][0]]

        with pysam.AlignmentFile(outfile, 'wb', header=header) as outbam:
            for read in inbam:
                read.set_tag('RG', rgid, replace=True)              
                outbam.write(read)

def split_ZMWs(bamfile: str, filter_ZMWs: list, number_splits: int, prefix: str):
    '''Split PacBio subread bamfile using filtered ZMW names
    
    Each split contains the same number of ZMWs.
    Args:
        bamfile (str): input PacBio bam file
        filter_ZMWs (list): list of ZMW names
        number_splits (int); number of the split
        prefix (str): prefix of the split file
    Returns:
        None. The output is stored in {prefix}_{No}.bam

    '''
    for i, zmws in enumerate(np.array_split(filter_ZMWs, number_splits)):
        outfile = '{}_{}.bam'.format(prefix, i)
        zmw_dict = {z:True for z in zmws}
        with pysam.AlignmentFile(bamfile, 'rb', check_sq=False) as inbam:
            with pysam.AlignmentFile(outfile, 'wb', template=inbam) as outbam:
                for read in inbam:
                    name = read.query_name
                    zmw = name.split('/')[1]
                    if zmw in zmw_dict:
                        outbam.write(read)

def makeReferenceDict(reffile: str):
    '''Generate Reference Dictionary using SeqIO.'''

    referenceDict = SeqIO.to_dict(SeqIO.parse(reffile, "fasta"))
    return referenceDict

def makeMotifPositionDict(refdict: dict, motif: str,motifmodposition: int):
    '''Generate Location of the Interested Nucleotide within the Motif
    
    Args:
        refdict (dict): Sequence Dictionary from SeqIO
        motif (str): Sequence of the Motif
        motifmodposition (int): 1-based relative location of the nucleotide in the motif
    
    Returns:
        dictionary: {sequence name: set of nucleotide locations (1-based)}
    '''
    motifpositiondict = {}
    for chro in refdict.keys():
        sequencestring = str(refdict[chro].seq.upper())
        forwardpositions = [m.start()+motifmodposition for m in re.finditer(motif,sequencestring)]
        sequencestringrev = str(refdict[chro].seq.complement().upper())
        reversepositions = [m.end()-motifmodposition+1 for m in re.finditer(motif[::-1],sequencestringrev)]
        motifpositiondict[chro] = set(forwardpositions + reversepositions)

    return motifpositiondict                        


def ipd_motif_finder(gff: str, motifpositionlist: set, avecov: int, mod: str, start: int, end:int):
    '''Filter the modified Motif Location from PacBio IPDSummary
    
    Args:
        gff (str): name of the gff file from PacBio IPDSummary
        motifpositionlist (set): 1-based modifided nucleotide location
        avecov (int): average coverage for Subreads per ZMW
        mod (str): Type of the modification predicted by PacBio IPDSummary

    Returns:
        list: sorted (1-based) positions of the modified nucleotides
    '''
    scorecutoffdict = {
                        7:12,8:12,9:12,10:12,11:12,12:12,13:12,14:12,15:12,16:12,17:12,18:12,19:12,
                        20:12,21:12,22:12,23:12,24:12,25:12,26:13,27:13,28:13,29:13,30:13,31:13,
                        32:13,33:13,34:13,35:13,36:13,37:13,38:13,39:13,40:13,41:13,42:13,43:13,
                        44:13,45:13,46:14,47:14,48:14,49:14,50:14,51:14,52:14,53:14,54:14,55:14,
                        56:14,57:14,58:14,59:14,60:14,61:14,62:14,63:14,64:14,65:14,66:15,67:15,
                        68:15,69:15,70:15,71:15,72:15,73:15,74:15,75:15,76:15,77:15,78:15,79:15,
                        80:15,81:15,82:15,83:15,84:15,85:15
                        }

    if avecov > 85:
        pvalue = 16
    else:
        pvalue = scorecutoffdict[avecov]

    modpositions = {}
    with open(gff) as filep:
        for line in filep:
            if '##' not in line:
                line = line.strip().split('\t')
                # if line[2] != mod: # Include modificaiton as well
                #     continue
                pos = int(line[3])
                if start <= pos <= end and pos in motifpositionlist and int(line[5]) >= pvalue:
                    modpositions[int(line[3])] = True
    modpositions = sorted(list(modpositions.keys()))
    return modpositions


def cigar_writer(modpositions: list, start: int,stop: int):
    '''Generate the Cigar string for the pseudo read of each ZMW.
    
    Args:
        modpositions (list): sorted list of identified motif locations (1-based)
        start (int): start position for read alignment
        stop (int): end position for read alignment

    Returns:
        str: Cigar string I for modification M for match (no modification)
    '''

    cigarstring = ''
    lastposition = start-1
    for position in modpositions:
        cigarstring += (str(position-lastposition)+'M')
        lastposition = position
        cigarstring += '1I'
         
    cigarstring += (str(stop-lastposition)+'M')
    return cigarstring

def generate_ipd_zmw_reads(bamfile: str, outfile: str, motifpositiondict: dict, 
                           reference: str, threads: int, coveragecutoff: int):
    '''Generate pseudo reads for each ZMW with modification predicted by IPDSummary.

    Args:
        bamfile (str): input aligned Bam file
        outfile (str): output Bam file
        motifpositiondict (dict): Nucleotide locations for modification analysis. Generated by makeMotifPositionDict
        reference (str): name of the reference file (should be indexed for IPDSummary)
        threads (int): number of threads for IPDSummary
        coveragecutoff (int): minimal coverage of ZMW

    Returns:
        None. Stored in outfile

    '''
    zmw_dict = {}
    final_bam = []
    with pysam.AlignmentFile(bamfile, 'rb', check_sq=False) as inbam:
        header = dict(inbam.header)
        # Edit header for new Bam file
        # label unsorted
        header['HD']['SO'] = 'unsorted'
        # add comment
        header['CO'] = ['SourceBamFile:%s' % bamfile]
        for read in inbam:
            zmw = read.query_name.split('/')[1]
            zmw_dict.setdefault(zmw, [])
            zmw_dict[zmw].append(read)
    # analysis for each zmw
    tmp_bam = '{}_tmp.bam'.format(bamfile[:-4])
    tmp_gff = '{}_tmp.gff'.format(bamfile[:-4])
    for count, (zmw, read_list) in enumerate(zmw_dict.items()):
        with pysam.AlignmentFile(tmp_bam, 'wb', header=header) as outbam:
            for read in read_list:
                outbam.write(read)
            read_header = read.header
        # count read depth
        cmd = 'samtools depth {}'.format(tmp_bam)
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True, text=True)
        depths = filter(lambda x:len(x)>0,(tuple(line.strip().split('\t')) for line in proc.stdout))
        depths = [d for d in depths if int(d[2])>=coveragecutoff]  
        if len(depths)>0:
            chrom = depths[0][0]
            start = int(depths[0][1])
            stop = int(depths[-1][1])
            avecov = int(round(sum([int(d[2]) for d in depths])/len(depths)))
        else:
            print('Error counting depth')
            return
            # continue
        subprocess.call('samtools index {}'.format(tmp_bam), shell=True)
        # calulcate p-value for motif position using ipdSummary from PacBio
        cmd = ('ipdSummary {} --reference {} '
               '--gff {} --identify m6A -w {}:{}-{} '
               '--quiet -j {} --identifyMinCov 3 --pvalue 0.9'.format(
                   tmp_bam, reference, tmp_gff, chrom,start,stop,threads)) #--identify 6mA
        proc = subprocess.run(cmd, shell=True, capture_output=True, check=True)
        # filter for high-quality motif position for each zmw
        motifpositionlist = motifpositiondict[chrom]
        modpositions = ipd_motif_finder(gff=tmp_gff,
                                         motifpositionlist=motifpositionlist,
                                         avecov=avecov,
                                         motif='m6A',
                                         start=start,
                                         end=stop)
        
        # generate cigar
        cigar = cigarWriter(modpositions,start,stop)
        # generate synthesized read for mod analysis
        new_read = pysam.AlignedSegment(header=read_header)
        new_read.query_name = zmw
        new_read.flag = 0           
        new_read.cigarstring = cigar
        new_read.mapping_quality = 255
        new_read.reference_name = chrom
        new_read.reference_start = start - 1 # pysam using 0-based coordinates
        new_read.set_tags([('cv', avecov, 'i'), ('ln', stop-start,'i')])
        final_bam.append(new_read)
    with pysam.AlignmentFile(outfile, 'wb', header=header) as outbam:
        for read in final_bam:
            outbam.write(read)    


def _worker_init(func):
    global _func
    _func = func


def _worker(x):
    return _func(x)

def _xmap(func, iterable, processes=None):
    with multiprocessing.Pool(processes, initializer=_worker_init, initargs=(func,)) as p:
        return p.map(_worker, iterable)  


def gen_ipd_zmw_multi(bam_prefix: str, out_prefix: str, idx: int, reference:str, motifpositiondict, coveragecutoff=9):
    '''Wrapper for multiprocess analysis of IPD prediction.
    '''

    print('Start working on {}_{}.bam'.format(bam_prefix, idx))
    generate_ipd_zmw_reads(bamfile='{}_{}.bam'.format(bam_prefix, idx), 
                           outfile='{}_{}.bam'.format(out_prefix, idx), 
                           motifpositiondict=motifpositiondict, 
                           reference=reference, 
                           threads=1,
                           coveragecutoff=coveragecutoff)    
    print('\nWriting {}_{}.bam'.format(out_prefix, idx))            


def run_multiprocess_ipd_analysis(bam_prefix: str, out_prefix: str, total_split: int, thread: int,
                                  reference: str, motifpositiondict: dict, coveragecutoff: int):
    '''
    Analyzer nucleotide modification using multiprocess
    '''
    _xmap(lambda i: gen_ipd_zmw_multi(
            bam_prefix=bam_prefix, 
            out_prefix=out_prefix,
            idx=i, 
            reference=reference, 
            motifpositiondict=motifpositiondict, 
            coveragecutoff=coveragecutoff), range(total_split), total_split)


def filter_low_coverage_mod(bamfile: str, reference: str, outfile: str, percentage_cutoff: int, motif: str, motifmodposition: int):

    '''Filter the locations with low percentage of modified nucleotides.

    Args:
        bamfile (str): input bam file
        reference (str): refenece file
        outfile (str): output file with all the percentage info
        filterfile (str): output file with the low percentage locations
        percentage_cutoff (int/float): lowcations with percentage < percentage_cutoff are stored in filterfile
        motif (str): motif to analyze modified nucleotide
        motifmodposition (int): 1-based position of modified nucleotide within the motif

    Returns:
        None. The results are stored in outfile and filterfile
    '''
    
    # preprocess
    refdict = makeReferenceDict(reference)
    motifdict = makeMotifPositionDict(
        refdict=refdict, 
        motif=motif, 
        motifmodposition=motifmodposition)

    motifcoverage = {}

    for chrom, location in motifdict.items():
        tmp = {p: [0,0] for p in location}
        motifcoverage[chrom] = tmp
    
    # Count modified nucleotide and total nucleotide in the reads

    with pysam.AlignmentFile(bamfile, 'rb', check_sq=False) as bam:
        count_issue = 0
        count_all = 0
        for read in bam:
            start = read.reference_start+1 # 0-based to 1-based
            end = read.reference_end # close location
            chrom = read.reference_name
            cigar_tuples = read.cigartuples # (Type, number) Type: 0 for M, 1 for I
            # first add total count
            test_pos = []
            for p in range(start, end + 1):
                if p in motifcoverage[chrom]:
                    motifcoverage[chrom][p][0] += 1
            # second add the modified count
            position = start - 1

            if cigar_tuples:
                for tp, num in cigar_tuples:
                    if tp == 1:
                        count_all += 1
                        position += num
                        if position not in motifcoverage[chrom]:
                            count_issue += 1
                        else:
                            motifcoverage[chrom][position][1] += 1
                    else:
                        position += num - 1
        # generate low cov list:
        all_percent = []
        low_percent = []
        for chrom, cdict in sorted(motifcoverage.items()):
            for pos, (total, mod) in sorted(cdict.items()):
                if total == 0:
                    continue
                percent = mod / total * 100
                all_percent.append([chrom, pos, percent])
                if percent < percentage_cutoff:
                    low_percent.append([chrom, pos, percent])
        with open(outfile, 'w') as filep:
            filep.write(''.join(map(lambda x: '{}\t{}\t{:.2f}\n'.format(x[0], x[1], x[2]), all_percent)))
        with open(filterfile, 'w') as filep:
            filep.write(''.join(map(lambda x: '{}\t{}\t{:.2f}\n'.format(x[0], x[1], x[2]), low_percent)))       
    if count_issue != 0:
        print('Error reading bam file! {} of {} detected modified As are not located in reference')   
        
        
def calculate_percent_m6A(bamfile: str, reference: str, bg='', data=''):
    '''Calculate the perbase methylation rate  '''
    # prepare A locations for m6A
    refdict = makeReferenceDict(reference)
    motifmodposition = makeMotifPositionDict(
        refdict=refdict, 
        motif='A',
        motifmodposition=1)    
    perbase_count = {}
    err = 0
    with pysam.AlignmentFile(bamfile, 'rb') as bam:
        for read in bam:
            chrom = read.reference_name
            perbase_count.setdefault(chrom, {})
            start = read.reference_start
            cigartuple = read.cigartuples
            pos = start - 1
            for tp, num in cigartuple:
                if tp == 0:
                    for i in range(num-1):
                        pos += 1
                        perbase_count[chrom].setdefault(pos, [0,0]) 
                        perbase_count[chrom][pos][1] += 1                    
                elif tp == 1: # Modified
                    pos += 1
                    if pos not in motifmodposition[chrom]:
                        err += 1
                    perbase_count[chrom].setdefault(pos, [0,0])                    
                    perbase_count[chrom][pos][0] += 1             
                    perbase_count[chrom][pos][1] += 1
    if err !=0:
        print('Error: {} m6A coordinates didnot match'.format(err))
    else:
        print('Analysis of {} is finished. Writing output...'.format(bamfile))
    if data:
        with open(data, 'wb') as filep:
            pickle.dump(perbase_count, filep, protocol=pickle.HIGHEST_PROTOCOL)

    # Generate bedgraph file
    if bg:
        content = ['track type=bedGraph']
        for chrom, cdict in perbase_count.items():
            if not cdict:
                continue
            for pos, (ca, ct) in sorted(cdict.items()):
                if ca == 0:
                    continue
                content.append('{}\t{}\t{}\t{:.1f}\n'.format(chrom, pos-1, pos,ca/ct*100))
        with open(bg, 'w') as filep:
            filep.write('\n'.join(content))            

