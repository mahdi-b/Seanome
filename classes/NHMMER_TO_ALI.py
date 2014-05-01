import re
import os
from Bio import SeqIO


class NHMMER_TO_ALI():
    def __init__(self, args, initialAliDir, hmmDir, genomeFile, outFastaDir):
        self.args = args
        self.hmmDir = hmmDir
        self.outFastaDir = outFastaDir
        self.genomeFile = genomeFile
        self.initialAliDir = initialAliDir
        self.hitSeq = None
        self.hitCoords = None
        self.strand = None


    def __ungapFastaFile__(self,inFile, outFile):
        seqs = []
        for seq in SeqIO.parse(inFile, 'fasta'):
            seq.seq = seq.seq.ungap("-")
            seqs.append(seq)
        SeqIO.write(seqs, open(outFile, 'w'), 'fasta')


    # Reads the best hit fromt the nhmmer output
    # an returns it as a SeqIO.Seq object
    def __getHitFromHmm__(self, hmmFile):
        lines =  open(os.path.join(self.hmmDir, hmmFile)).readlines()

        hit_re = re.compile(r'score\s*bias\s*Evalue\s*hmmfrom\s*hmm\s*to\s*alifrom')
        length=0 # NOT USED FOR NOW
        query_length = re.compile("Query:\s*ali\s*\[M=(\d+)")
        for i in range(0,len(lines)):
            if query_length.search(lines[i]):
                length = query_length.search(lines[i]).groups()[0]
            elif hit_re.search(lines[i]):
                i +=2
                data = lines[i].split()
                qStart, qEnd, self.hStart, self.hEnd, acc = [int(data[i]) for i in [4,5,10,11]]+[data[14]]
                print lines[i]
                print qStart, qEnd, self.hStart, self.hEnd, acc,"\n"
                seqLen = 0
                if  self.hEnd > self.hStart:
                    self.strand="+"
                    seqLen = self.hEnd - self.hStart + 1
                else:
                    seqLen = self.hStart - self.hEnd + 1
                    self.strand="-"
                if int(seqLen) > self.args.min_csr_len and float(acc) > self.args.min_csr_sim:
                    if self.strand=="+":
                        self.hitSeq = self.__getSubSeqFromFasta__(self.hStart-1, self.hEnd) # nhmmer alignment is 1 based
                    else:
                        self.hitSeq = self.__getSubSeqFromFasta__(self.hEnd-1, self.hStart, True) # nhmmer alignment is 1 based
                    return True
                else:
                    print "Parameters not met for alignment %s " % hmmFile

    def __getSubSeqFromFasta__(self, start, end, reverse=False):
        seq = SeqIO.read(self.genomeFile, 'fasta')
        if reverse:
            seq.id = "%s_%s_%s_R" % (seq.id, start+1, end)
            seq.seq = seq[start:end].reverse_complement().seq
        else:
            seq.id = "%s_%s_%s" % (seq.id, start+1, end)
            seq.seq = seq[start:end].seq
        return seq

    def run(self):
        # get self.hitSeq = None, self.hitCoords = None, self.strand = None
        for output in os.listdir(self.hmmDir):
            self.__getHitFromHmm__(output)
            # write the ali (gapped) sequences found previously to a new (ungapped) file
            # we delay the writing in case seq is not found in hmmer, therefore we don't print a file
            self.__ungapFastaFile__(os.path.join(self.initialAliDir, output), os.path.join(self.outFastaDir,output))
            if self.hitSeq:
                SeqIO.write(self.hitSeq, open(os.path.join(self.outFastaDir,output), 'a'), 'fasta')

    def dryRun(self):
        pass
    
