# -*- coding: utf-8 -*-

import os, sys
from glob import glob
import subprocess as sp


class AutoRunMD :

    def __init__(self, topFile, taskName, grompp, mdrun, verbose=True, qsub=False):
        self.top = topFile

        self.PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))

    def run_suprocess(self, cmd):
        job = sp.Popen(cmd, shell=True)
        job.communicate()

        return self

    def generate_top(self, inpdb, ff="amber99sb-ildn", water='tip3p',
                     outgro="output", top="topol", ignh=True):

        cmd = "gmx pdb2gmx -f % -o %s -p %s -ff %s -water %s " %\
              (inpdb, outgro, top, ff, water)
        if ignh: cmd += "-ignh "

        self.run_suprocess(cmd)

        return self

    def modify_mdp(self, inmdp, outmdp, parameters):

        tofile = open(outmdp, 'wb')
        with open(inmdp) as lines :
            for s in lines :
                if len(s.split()) > 0 and s[0] != ";" :
                    if s.split()[0] in parameters.keys() \
                            and len(parameters[s.split()[0]]) > 0 :
                        tofile.write("%s    = ")
                        for item in parameters[s.split()[0]] :
                            tofile.write("%s "%str(item))
                        tofile.write(" \n")
                    else:
                        tofile.write(s)
                else :
                    tofile.write(s)

        return self

    def add_box(self, ingro, outgro, ):
        cmd = "editconf -f %s -o %s -c -d 1.2 -bt cubic" % (ingro, outgro)

        self.run_suprocess(cmd)

        return self

    def add_solvent(self, ingro, outgro, intop, spc="spc903.gro"):

        if not os.path.exists(spc):
            spc = os.path.join(self.PROJECT_ROOT, "data/spc903.gro")

        cmd = "gmx solvate -cp %s -cs %s -o %s -p %s " % (ingro, spc, outgro, intop)
        self.run_suprocess(cmd)

        return self

    def add_ions(self, emmdp, ingro, intop, outgro, ion_conc=0.15):
        if not os.path.exists(emmdp):
            emmdp = os.path.join(self.PROJECT_ROOT, "data/em_sol.mdp")

        cmd1 = "gmx grompp -f %s -c %s -p %s -o %s -maxwarn 100" % (emmdp, ingro, intop, "addion")
        self.run_suprocess(cmd1)

        cmd2 = "echo \"13 0 0 \" | gmx genion -s %s -p %s -o %s -neutral -conc %.2f " %\
               ("addion", intop, outgro, ion_conc)
        self.run_suprocess(cmd2)

        return self

    def minimize(self, ingro, outgro, emmdp, intop, nt=4):
        if not os.path.exists(emmdp):
            emmdp = os.path.join(self.PROJECT_ROOT, "data/em_sol.mdp")

        cmd1 = "gmx grompp -f %s -c %s -p %s -o %s -maxwarn 100" % (emmdp, ingro, intop, outgro)
        self.run_suprocess(cmd1)

        cmd2 = "gmx mdrun -deffnm %s -nt %d -v -gpu_id 0" % (outgro, nt)
        self.run_suprocess(cmd2)

        return self

    def md(self, ingro, outgro, nptmdp, intop, nt=4, restraints=False):

        if not restraints:
            mdp = nptmdp
            if not os.path.exists(nptmdp):
                mdp = os.path.join(self.PROJECT_ROOT, "data/npt.mdp")
        else:
            mdp = "npt_rest.mdp"
            if not os.path.exists(mdp):
                mdp = os.path.join(self.PROJECT_ROOT, "data/npt_rest.mdp")

        cmd1 = "gmx grompp -f %s -c %s -p %s -o %s -maxwarn 100" % (mdp, ingro, intop, outgro)
        self.run_suprocess(cmd1)

        cmd2 = "gmx mdrun -deffnm %s -nt %d -v -gpu_id 0" % (outgro, nt)
        self.run_suprocess(cmd2)

        return self

