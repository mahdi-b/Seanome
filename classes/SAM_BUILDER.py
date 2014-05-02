import re
import os
import copy
import pysam
import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq



class SAM_BUILDER():
    
    cigar_re = re.compile(r'([0-9]+)([M=XID])')

    def __init__(self, args, nopile=True, makeSmall = False):
        self.args = args
        self.nopile = nopile
        self.makeSmall = makeSmall


    def run(self):
        outfile = None
        samout= os.path.join(self.args.output_dir, os.path.basename(self.args.consensus_ali))
        refs = [ (r.id, len(r.seq) ) for r in  SeqIO.parse(self.args.consensus_ali, "fasta")]
        header = dict(HD = dict(VN = '1.0'), SQ = [ {'LN': refs[0][1], 'SN': refs[0][0] }] )
        if samout:
            outfile = pysam.Samfile( samout, "wh", header = header )
        refcnt = 0
        rcnt = 0
        samdat = None
        for seq in SeqIO.parse(self.args.consensus_ali, "fasta"):
            if seq.id.startswith('Con'):
                rcnt += 1
                continue
            small = None
            if self.makeSmall:
                small = pysam.Samfile( "%s.sam"%(refs[rcnt][0]), "wh", header = header )

            info = seq.id.split("_")
            rev = False
            if len(info) == 4 and info[3] == 'R':
                rev = True
            outfname = "%s.pile"%(info[0])
            out = None
            if not self.nopile:
                out = open(outfname, "w")

            samid = os.path.join(self.args.split_sam_dir,"%s.bam" % info[0] ) 
            samfile = pysam.Samfile(samid, "rb" )
            myReads = samfile.fetch( info[0], int(info[1]) , int(info[2]) )


            gaps = self.__accumulateGaps__(str(seq.seq))
            reglen = 0
            # pysam converts to python convention for pos.. 0 based
            # convert the incoming values from the fasta file to 0 base as well
            if not rev:
                region = (int(info[1]) - 1 , int(info[2]), )
            else:
                region = (int(info[1]) , int(info[2]) , )
                reglen = int(info[2]) - int(info[1])

            for read in (m for m in myReads if not m.is_unmapped):
            # trim everything to its proper size
                seqStr, qual = self.__addSeqGaps__(read)
                seqStr = seqStr[region[0] : region[1]]
                qual = qual[region[0] : region[1]].strip()
                if rev:
                    space = seqStr.count(" ")
                    seqStr = str(Seq(seqStr).reverse_complement())
                    # subtract leading spaces from region length
                    # this will give us how far of the right side we should be (minus the seq itself)
                    end = reglen - space
                    # subtract the sequence length to obtain the start position
                    start = end - len(seqStr)
                    # add the spaces back to the start of the now stripped sequence & reverse the quality
                    seqStr = " " * (start + space) + seqStr.strip()
                    qual = qual[::-1]
                    # as a final step, apply spaces or gaps based on the reference
                seqStr = self.__applyGaps__(gaps, seqStr)
                if out:
                    print >> out,  seqStr
                if samout or self.makeSmall:
                    samdat = self.__generateSamInfo__(read, seqStr, 0, qual)
                if samout:
                    outfile.write(samdat)
                if self.makeSmall:
                    small.write(samdat)

            if self.makeSmall:
                small.close()
            rcnt += 1
            if out:
                out.close()
        if outfile:
            outfile.close()



    def dryRun(self):
        pass

    def __addSeqGaps__(self, read):
        myStr = ""
        myQual = ""
        pos = 0
        for m in self.cigar_re.finditer(read.cigarstring):
            if m.group(2)== "D":
                myStr += '-' * int(m.group(1))
                myQual += ' ' * int(m.group(1))
            else:
                myStr += read.query[pos:pos+int(m.group(1))]
                myQual += read.qual[pos:pos+int(m.group(1))]
                pos += int(m.group(1))
        return (' ' * read.pos + myStr , " " * read.pos + myQual)

    
    def __makeCigar__(self, seq):
        c = ""
        cnt = 0
        cig = []
        for j in ("M" if (i != '-') else 'D' for i in seq):
            if c != j:
                if cnt != 0:
                    cig.append("%s%s"%(cnt, c))
                cnt = 1
                c = j
            else:
                cnt += 1
        if cnt != 0:
            cig.append("%s%s"%(cnt, c))
        return "".join(cig)



    def __accumulateGaps__(self, seq):
        gaps = []
        start = None
        for idx, i in enumerate(seq):
            if i == '-':
                if start == None:
                    start = idx
                else:
                    if start != None:
                        gaps.append( (start, idx, idx - start, ) )
                    start = None
        if start != None:
            gaps.append( (start, len(seq), len(seq) - start, ) )
        return gaps


    def __generateSamInfo__(self, sinfo, seqstr, refname, qual):
        spaces = seqstr.count(" ")
        a = pysam.AlignedRead()
        a.rname = refname
        a.qname = sinfo.qname
        a.seq= seqstr.strip().replace("-","")
        a.flag = sinfo.flag & 16
        a.pos =  spaces
        a.mapq = sinfo.mapq
        a.cigarstring =  self.__makeCigar__(seqstr.strip())
        a.rnext = -1
        a.pnext= -1 
        a.tlen = 0
        a.qual = qual.strip().replace(" ","")
        return a

    def __applyGaps__(seld, gaps, seqStr):
        for p in gaps:
            if p[0] > len(seqStr):
                pass
            elif p[0] == len(seqStr):
                seqStr = seqStr + '-' * p[-1]
            else:
                c = '-'
                if seqStr[p[0]] == ' ' or p[0] == 0:
                    c = ' '
                seqStr = seqStr[:p[0]] + c * p[-1] + seqStr[p[0]:]
        return seqStr
