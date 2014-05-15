import os
import copy

import pysam
from Bio import SeqIO
from Bio.Seq import Seq
import re
import sys


class SAM_BUILDER():
    cigar_re = re.compile(r'([0-9]+)([M=XID])')
    soft_clip = re.compile(r'(^[0-9]+)([S])')

    def __init__(self, args, nopile=True, makeSmall = False):
        self.args = args
        self.nopile = nopile
        self.makeSmall = makeSmall


    def getShiftedPosRef(self, seq, startInRef):
        pos = 0
        nbMatches = 0
        if startInRef == 0:
            return pos
        while nbMatches < startInRef:
            if seq[pos] != "-":
                nbMatches += 1
                pos += 1
            else:
                pos += 1
        return pos


    # def getShiftedPosRead(self, cigar, start):
    #     pos = 0
    #     shiftedReadStart=start
    #     while pos <= start + 1:
    #         if cigar[pos] == "I" :
    #             shiftedReadStart += 1
    #         # DSL - added since positions in tablet for reads that fit this case were off.
    #         elif cigar[pos] == 'D':
    #             shiftedReadStart -= 1            
    #         pos += 1
    #     return shiftedReadStart

    def getShiftedPosRead(self, cigar, start):
        pos = 0
        shiftedReadStart=start
        nbMatches = 0
        while nbMatches < start:
            if cigar[pos] == "I" :
                shiftedReadStart += 1
                # DSL - added since positions in tablet for reads that fit this case were off
            elif cigar[pos] == 'D':
                shiftedReadStart -= 1
            else:
                nbMatches+=1
            pos += 1
        return shiftedReadStart



    def expandCigar(self, cigar):
        cigstr = ""
        for m in self.cigar_re.finditer(cigar):
            cigstr +=  m.group(2) * int(m.group(1))
        return cigstr


    def compressCigar(self,cigar):
        cigar = cigar.rstrip("D")
        c = None
        cnt = 0
        ciglst = []
        for b in cigar:
            if b != c:
                if c!= None:
                    ciglst.append("%s%s"%(cnt, c)) 
                cnt = 1
                c = b
            else:
                cnt += 1
        if c!= None:
            ciglst.append("%s%s"%(cnt, c))     
        return "".join(ciglst) 


    def __getTag__(self, sinfo, tag):
        try:
            return sinfo.opt(tag)
        except KeyError:
            return None

    def generateSamInfo(self, sinfo, seq, cig, startpos, refname, qual):
        a = pysam.AlignedRead()
        a.rname = refname
        a.qname = sinfo.qname
        a.seq= seq
        a.flag = 0 #sinfo.flag & 16
        a.pos =  startpos
        a.mapq = sinfo.mapq
        a.cigarstring =  cig
        a.rnext = -1
        a.pnext= -1 #a.pos
        #a.isize = sinfo.isize
        a.tlen = 0
        a.qual = qual.strip().replace(" ","")
        tags=  []
        tmp = self.__getTag__(sinfo, "RG")
        if tmp:
            tags.append( ("RG", tmp,) )

        tmp = self.__getTag__(sinfo, "X0")
        if tmp:
            tags.append( ("X0", tmp,) )
        tmp = self.__getTag__(sinfo, "AS")
        if tmp:
            tags.append( ("AS", tmp,) )
        tmp = self.__getTag__(sinfo, "XS")
        if tmp:
            tags.append( ("XS", tmp,) )
        tmp = self.__getTag__(sinfo, "YS")
        if tmp:
            tags.append( ("YS", tmp,) )
        tmp = self.__getTag__(sinfo, "YT")
        if tmp:
            tags.append( ("YT", tmp,) )
        a.tags = tags
        #a.tags = sinfo.tags
        return a


    def parseMSA(self, ffile):
        msaRefs = dict()
        for seq in  SeqIO.parse(ffile, 'fasta'):
            if not seq.id.startswith("Con"):
                nameInfo = seq.id.split("_")
                msaRefs[nameInfo[0]]=[int(nameInfo[1])-1, int(nameInfo[2])-1, str(seq.seq)]
                if len(nameInfo)==4:
                    msaRefs[nameInfo[0]].append("Reversed")
        return msaRefs

    
    def computeStartPositions(self, read, refStart, refEnd, isReverse):
        if not isReverse: # this means ref is forward
            #if read starts Before msaRef
            if read.pos < refStart:
                startInRef = 0
                startInRead = refStart - read.pos
            else:
                startInRead = 0
                startInRef = read.pos - refStart
        else:
            endbase = read.aend - 1
            if endbase < refEnd:
                startInRef  =  refEnd - endbase
                startInRead =  0
            else:
                startInRef  =  0 
                startInRead =  endbase - refEnd
        return startInRef, startInRead


    def processSingleRef(self, samfile, name, refStart, refEnd, refSeq, isReversed, combinedOut, header):
        if self.makeSmall:
            small = pysam.Samfile( "%s.sam"%(name), "wh", header = header )
        else:
            small = None
        myReads = samfile.fetch(name, refStart, refEnd)

        for read in (q for q in myReads if not q.is_unmapped):
            startInRef, startInRead = self.computeStartPositions(read, refStart, refEnd, isReversed)
            shiftedStartRef = self.getShiftedPosRef(refSeq, startInRef)
            expandedCigar  = self.expandCigar(read.cigarstring)
            shiftedStartRead = self.getShiftedPosRead(expandedCigar, startInRead)

            if not isReversed:
                readString = read.query[shiftedStartRead:]
                readCigar = expandedCigar[shiftedStartRead:]
                quality = read.qqual[shiftedStartRead:]
            else:
                readString = str(Seq(read.query).reverse_complement()[shiftedStartRead:])
                quality = read.qqual[::-1][shiftedStartRead:]
                readCigar = expandedCigar[::-1][shiftedStartRead:]


            refString = refSeq[shiftedStartRef:]
            #DEBUGGING
            print "read is %s" % read.qname
            print "orignal Cigar is %s " % read.cigarstring
            print "original position %s " % read.pos
            print "original read str %s " % read.query
            
            print "startInRef is %s" % (startInRef)
            print "shiftedStartRef is %s"  % shiftedStartRef

            print "startInRead is %s"  % startInRead
            print "shiftedStartRead is %s" % (shiftedStartRead)

            print readCigar
            print refString
            print readString
            # inspecting columns and inserting gaps or deleting Inserts from read as necessary
            posRef = 0
            posRead = 0
            posCigar = 0
            newReadString = ""
            newCigarString = ""
            tempString = "" # Just visualization purposes for now
            newQualString = ""

            while posRef < len(refString) and posRead < len(readString):
                #print refString[:posRef+1]
                #print readString[:posRead+1]
                #print readCigar[:posCigar+1]
                if refString[posRef] == "-":
                    if readCigar[posCigar] == "I":
                        newReadString += readString[posRead]
                        tempString += readString[posRead]
                        newQualString += quality[posRead]
                        newCigarString += "M"
                        posRef += 1
                        posRead += 1
                        posCigar += 1
                    elif readCigar[posCigar] == "D":
                        newCigarString += readCigar[posCigar]
                        posRef += 1
                        #posCigar += 1
                    else:
                        tempString += "-"
                        newCigarString += "D"
                        posRef += 1                    
                elif readCigar[posCigar] == "I": # we know ref has a non-gap
                    posRead += 1
                    posCigar += 1
                    #print "I encountered"
                elif readCigar[posCigar] == "D":
                    tempString += "-"
                    newCigarString += readCigar[posCigar]
                    posRef += 1
                    posCigar += 1
                else:
                    newReadString += readString[posRead]
                    newQualString += quality[posRead]
                    tempString += readString[posRead]
                    newCigarString += readCigar[posCigar]
                    posRef += 1
                    posRead += 1
                    posCigar += 1
                #print tempString
                #print ""

            compactCigar = self.compressCigar(newCigarString)
            samData = self.generateSamInfo(read, newReadString, compactCigar, shiftedStartRef, 0, newQualString)
            #DEBUGGING
            print "refStr : %s" % refString
            print "tempStr: %s" % tempString
            print "cigrStr: %s" % newCigarString
            print "cmpCigr: %s\n" % compactCigar
            combinedOut.write(samData) 
            if small:
                small.write(samData) 
        if small:
            small.close()


    def run(self):
        msaInfo = self.parseMSA(self.args.consensus_ali)
        samout = os.path.join(self.args.output_dir, os.path.basename(self.args.consensus_ali))

        refs = [ (r.id, len(r.seq) ) for r in  SeqIO.parse(self.args.consensus_ali, "fasta")]
        header = dict(HD = dict(VN = '1.0'), SQ = [ {'LN': refs[0][1], 'SN': refs[0][0] }] )
        outfile = pysam.Samfile(samout, "wh", header = header )
        
        for refName, refVals in msaInfo.iteritems():
            isReversed = True
            if len(refVals) == 3:
                isReversed = False
            samid = os.path.join(self.args.split_sam_dir,"%s_sorted_WRG.bam" % refName )
            #samid = "./concat_contigs/%s.bam"%(refName)
            samfile = pysam.Samfile(samid, "rb" )        
            self.processSingleRef(samfile, refName, refVals[0], refVals[1], refVals[2], isReversed, outfile, header)
            samfile.close()
        outfile.close()


    def dryRun(self):
        pass

