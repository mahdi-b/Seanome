#!/usr/bin/python
import os
import sys
import argparse
import logging
import pysam
import tempfile

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from multiprocessing import Pool

from classes.ProgramRunner import *
from classes.NHMMER_TO_ALI import *



def runInstance(myInstance):
   dryRunInstance(myInstance)
   myInstance.run()

def dryRunInstance(myInstance):
   logging.warning(myInstance.dryRun())


def makeDirOrdie(dirPath):
   if not os.path.isdir(dirPath):
      os.makedirs(dirPath)
   else:
       # TODO; Uncomment after testing done
      #logging.error("Split fasta directory %s already exists " % hmmerOutputDir)
      #sys.exit()
      pass
   return dirPath


def maskFile(args,samFile):
   '''
   Used in maskGenome
   '''
   samfile = pysam.Samfile(samFile)
   alis = samfile.fetch()
   mySeq =  SeqIO.read(args.input, 'fasta')
   mseq = mySeq.seq.tomutable()
   for ali in alis:
      tags =  dict(ali.tags)
        # if alignment hits with at least the minimum similarity, then mask the hit region                                                                                          
      if tags['XM'] and (1 - float(tags['XM'])/int(ali.qlen) >= args.min_similarity):
         # TODO: collect the intervals and do it at once since intervals overalap
         maskRange = (min(ali.positions), max(ali.positions))
         mseq[maskRange[0]: maskRange[1]+1] = "N" * (maskRange[1] - maskRange[0] +1)
         mySeq = SeqRecord(mseq, id=mySeq.id)
   SeqIO.write(mySeq, open(args.output, 'w'), 'fasta')


def maskGenome(args, pool):
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
        maskFile(args, name+"_repeats_filtered.sam")
    else:
        logging.debug("No repeats found is %s" % args.input)
    logging.debug("Done masking input %s" % args.input)


def generateAlis(args, lastzOutFile):
    '''
    Used in find_csr
    '''
    inFile = open(lastzOutFile, 'r')
    outputDir= args.output_dir
    outFileNum=0
    for line in inFile:
        line = line.rstrip()
        data = line.split()
        aliLength = min(len(data[2]),len(data[3]))
        simRatio = float(sum([data[2][i]==data[3][i] for i in range(len(data[2]))]))/ len(data[2])
        if aliLength > args.min_csr_len and simRatio > args.min_csr_sim:
            outFile = open(outputDir+"/outFile_%s_%s" % (str(outFileNum), str(aliLength)), 'w')
            outFile.write( ">%s_%s_%s\n%s\n>%s_%s_%s\n%s\n" % (data[0],data[7],data[8],data[2],data[1],data[9], data[10],data[3],))
            outFile.close()
            outFileNum+=1
    inFile.close()
    


def find_seed_csr(args, pool):
    outName = os.path.splitext(os.path.basename(args.input_1))[0]+"_"+os.path.splitext(os.path.basename(args.input_2))[0]+".lastz"
    logging.debug("Finding seed common shared regions between %s and %s: \n Starting... " % (args.input_1, args.input_2))
    prog = ProgramRunner("lastz", [args.input_1, args.input_2, args.min_csr_sim, outName])
    prog.run()
    generateAlis(args, outName)
    logging.debug("Finished finding seed common shared regions between %s and %s" % (args.input_1, args.input_2))
    
def find_csr(args, pool):
   # TODOL
   # validate the params
   # throw error if program returns error

   tmpDir = tempfile.mkdtemp(prefix="seanome_", dir=".")
   nhmmer_dir = os.path.join(tmpDir, "nhmmer_dir")
   os.mkdir(nhmmer_dir)
   aliFiles = os.listdir(args.alignments_dir)
   logging.debug("Finding %s alignments from  %s in ref %s: \n" % (len(aliFiles), args.alignments_dir, args.genome))
   pool.map(runInstance, [ProgramRunner("nhmmer", [os.path.join(args.alignments_dir, x), args.genome , os.path.join(nhmmer_dir,x.split(".")[0])] ) for x in aliFiles])

   logging.debug("Parsing outputs")
   nhmmerToAli = NHMMER_TO_ALI( args, args.alignments_dir, nhmmer_dir, args.genome, args.fasta_output)
   nhmmerToAli.run()

   logging.debug("Aligning hits")
   #Run mutiple sequence alignment 
   pool.map(runInstance, [ProgramRunner("muscle", [os.path.join(args.fasta_output, x), os.path.join(args.alis_output, x)] ) for x in aliFiles])

   logging.debug("Finished finding common shared regions between %s and %s" % (args.alignments_dir, args.genome))

def consensus(args, pool):
   aliFiles = os.listdir(args.alis_dir)
   pool.map(runInstance, [ProgramRunner("addConsensus", [os.path.join(args.alis_dir, x), os.path.join(args.cons_output, x)] ) for x in aliFiles])

def main(argv):
    parser = argparse.ArgumentParser(description="Seanome description", epilog="Seanome long text description")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s '+version)
    parser.add_argument('-t', '--threads', type=int, default = 1)
    subparsers = parser.add_subparsers(dest='action', help='Available commands')

    # Mask Genome
    parser_mask = subparsers.add_parser('mask')
    parser_mask.add_argument('-i', '--input',  required=True, help=" Input file to maks")
    parser_mask.add_argument('-o', '--output', required=True, help=" Masked output file")
    parser_mask.add_argument('-l', '--min_length', default=80, help=" Minimum alignment length")
    parser_mask.add_argument('-s', '--min_similarity', default=0.86, help=" Minimum alignment similarity")
    parser_mask.set_defaults(func=maskGenome)

    # seed_csr
    parser_mask = subparsers.add_parser('seed_csr')
    parser_mask.add_argument('-i1', '--input_1',  required=True, help="First input contig")
    parser_mask.add_argument('-i2', '--input_2', required=True, help="Second input contig")
    parser_mask.add_argument('-o', '--output_dir', type=makeDirOrdie, required=True, help=" Minimum common shared region similarity")
    parser_mask.add_argument('-l', '--min_csr_len', default=150, help=" Minimum common shared region length")
    parser_mask.add_argument('-s', '--min_csr_sim', default=0.88, help=" Minimum common shared region similarity")
    parser_mask.set_defaults(func=find_seed_csr)

    # seed_csr
    parser_mask = subparsers.add_parser('find_csr')
    parser_mask.add_argument('-a', '--alignments_dir',  required=True, help="Inut alignments directory")
    parser_mask.add_argument('-g', '--genome', required=True, help="Reference genome")
    parser_mask.add_argument('-f', '--fasta_output', required=True, type=makeDirOrdie, help="fasata sequences output")
    parser_mask.add_argument('-o', '--alis_output', required=True, type=makeDirOrdie, help="fasata sequences output")
    parser_mask.add_argument('-l', '--min_csr_len', default=150, help=" Minimum common shared region length")
    parser_mask.add_argument('-s', '--min_csr_sim', default=0.88, help=" Minimum common shared region similarity")
    parser_mask.set_defaults(func=find_csr)

    # consensus
    parser_mask = subparsers.add_parser('consensus')
    parser_mask.add_argument('-a', '--alis_dir',  required=True, help="Inut alignments directory")
    parser_mask.add_argument('-c', '--cons_output', required=True, type=makeDirOrdie, help="Output Consensus Directory")
    parser_mask.set_defaults(func=consensus)

    # Parse arguments
    args = parser.parse_args()
    pool = Pool(processes=args.threads)
    logging.debug("Initial ARGS are:")
    logging.debug(args)
    args.func(args, pool)

if __name__ == "__main__":
   version = "alpha 0.01"
   FORMAT = "%(asctime)-15s  %(message)s"
   logging.basicConfig(format=FORMAT, level=logging.DEBUG)
   main(sys.argv)
