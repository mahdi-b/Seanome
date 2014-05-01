#!/usr/bin/python
import os
import sys
import argparse
import logging



from Bio import SeqIO
from multiprocessing import Pool

from classes.ProgramRunner import *

version = "alpha 0.01"

FORMAT = "%(asctime)-15s  %(message)s"
logging.basicConfig(format=FORMAT, level=logging.DEBUG)

def maskFile(samFile, inFileName, outFileName):
    samfile = pysam.Samfile(samFile)
    alis = samfile.fetch()
    mySeq =  SeqIO.read(inFileName, 'fasta')
    mseq = mySeq.seq.tomutable()
    for ali in alis:
        tags =  dict(ali.tags)
        # if alignment hits with at least the minimum similarity, then mask the hit region                                                                                          
        if tags['XM'] and (1 - float(tags['XM'])/int(ali.qlen) >= args.min_similarity):
            # TODO: collect the intervals and do it at once since intervals overalap
            maskRange = (min(ali.positions), max(ali.positions))
            print "masking region %s to %s" % (min(ali.positions), max(ali.positions))
            mseq[maskRange[0]: maskRange[1]+1] = "N" * (maskRange[1] - maskRange[0] +1)
            mySeq = SeqRecord(mseq, id=mySeq.id)
    SeqIO.write(mySeq, open(outFileName, 'w'), 'fasta')


def maskGenome(args):
    print args.input
    name =  os.path.splitext(os.path.basename(args.input))[0]    
    # Can eventually be parallelized using the pool of threads
    logging.debug("Starting the masking of input %s" % args.input)
    # run 
    logging.debug("Computing frequencies for finding repeats" )
    prog = ProgramRunner("build_lmer_table", [args.input, name+".freqs" ] )

    prog.run()
    logging.debug("Finding repeats in the %s" % args.input)
    prog = ProgramRunner("repeatScout", [args.input, name+"_repeats.fa", name+".freqs" ] )

    prog.run()
    # Filter the repeats that do not pass the requirements
    # for now, this only consists in dropping short reads < N
    logging.debug("Filetering repeats")
    repeats = SeqIO.parse(name+"_repeats.fa", 'fasta')
    filterReads=[]
    for read in repeats:
        if len(read.seq) >= args.min_length:
            filterReads.append(read)
    SeqIO.write(filterReads, open(name+"_repeats_filtered.fa", 'w'), 'fasta')
    if len(filterReads) >= 1:
        prog = ProgramRunner("bowtie-build", [args.input, name])
        prog.run()
        prog = ProgramRunner("bowtie-align", [name, name+"_repeats_filtered.fa", name+"_repeats_filtered.sam"])
        prog.run()
        # mask the regions that align to a repeat
        samfile = pysam.Samfile(X)
        alis = samfile.fetch()
        mySeq =  SeqIO.read(, 'fasta')
        mseq = mySeq.seq.tomutable()






    else:
        logging.debug("No repeats found is %s" % args.input)
    logging.debug("Done masking input %s" % args.input)





def main(argv):
    parser = argparse.ArgumentParser(description="Seanome description", epilog="Seanome long text description")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s '+version)
    subparsers = parser.add_subparsers(dest='action', help='Available commands')


    # Mask Genome
    parser_mask = subparsers.add_parser('mask')
    parser_mask.add_argument('-i', '--input',  required=True, help=" Input file to maks")
    parser_mask.add_argument('-o', '--output', required=True, help=" Masked output file")
    parser_mask.add_argument('-l', '--min_length', default=80, help=" Minimum alignment length")
    parser_mask.add_argument('-s', '--min_similarity', default=86, help=" Minimum alignment similarity")

    parser_mask.set_defaults(func=maskGenome)

    args = parser.parse_args()
    logging.debug("Initial ARGS are:")
    logging.debug(args)
    args.func(args)



if __name__ == "__main__":
    main(sys.argv)
