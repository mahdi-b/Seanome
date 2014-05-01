import os
# RUNNABLE OBJECTS that we use in pool map by calling run instance on 

PROGS_ROOT_repeatscout="/home/mahdi//programs/RepeatScout-1/"
PROGS_ROOT_bowtie2="~/programs/bowtie2-2.2.1/"


class ProgramRunner():
    commands = {
        "build_lmer_table" : PROGS_ROOT_repeatscout+"build_lmer_table -l 12 -sequence %s -freq %s",
        "repeatScout" : PROGS_ROOT_repeatscout+"RepeatScout -sequence %s -output %s -freq %s -l 12",
        "bowtie-build" : PROGS_ROOT_bowtie2+"bowtie2-build %s %s",
        "bowtie-align" : PROGS_ROOT_bowtie2+"bowtie2 -f --score-min L,-0.6,-20  --all -N 1 -L 15    --end-to-end -x %s -U %s -S %s",

        }

    def __init__(self, program, params):
        self.program = program
        self.command = self.commands[program] % tuple(params)

    def run(self):
        os.system(self.command)

    def dryRun(self):
        return self.command
    

        
