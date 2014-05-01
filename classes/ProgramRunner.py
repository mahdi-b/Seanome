import os
# RUNNABLE OBJECTS that we use in pool map by calling run instance on 


# Paths to Binaries
PROGS_ROOT_repeatscout = "/home/mahdi//programs/RepeatScout-1/"
PROGS_ROOT_bowtie2 = "~/programs/bowtie2-2.2.1/"
PROGS_ROOT_lastz = "/home/mahdi/programs/lastz-distrib-1.02.00/src/"
PROGS_ROOT_hmmer = "/home/mahdi/programs/hmmer-3.1b1-linux-intel-x86_64/binaries/"
PROGS_ROOT_muscle = ""


LASTZ_OUTPUT_FORMAT = "general:name1,name2,text1,text2,score,strand1,strand2,start1,end1,start2,end2"



class ProgramRunner():
    commands = {
        "build_lmer_table" : PROGS_ROOT_repeatscout+"build_lmer_table -l 12 -sequence %s -freq %s",
        "repeatScout" : PROGS_ROOT_repeatscout+"RepeatScout -sequence %s -output %s -freq %s -l 12",
        "bowtie-build" : PROGS_ROOT_bowtie2+"bowtie2-build %s %s",
        "bowtie-align" : PROGS_ROOT_bowtie2+"bowtie2 -f --score-min L,-0.6,-20  --all -N 1 -L 15    --end-to-end -x %s -U %s -S %s",
        "lastz" : PROGS_ROOT_lastz+"lastz %s %s --strand=both --identity=%s --format="+LASTZ_OUTPUT_FORMAT+" --masking=1  > %s",
        "nhmmer" : PROGS_ROOT_hmmer+"nhmmer  %s %s > %s",
        "muscle" : PROGS_ROOT_muscle+"muscle -in %s  -out %s",
        }

    def __init__(self, program, params):
        self.program = program
        self.command = self.commands[program] % tuple(params)

    def run(self):
        os.system(self.command)

    def dryRun(self):
        return self.command
    

        
