#!/usr/bin/env python

import os, sys
import subprocess as sp
import argparse
from argparse import RawDescriptionHelpFormatter


class AutoRunMD :

    def __init__(self):
        self.top = ""

        self.PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))

    def run_suprocess(self, cmd):
        job = sp.Popen(cmd, shell=True)
        job.communicate()

        return self

    def generate_top(self, inpdb, ff="amber99sb-ildn", water='tip3p',
                     outgro="output", top="topol", ignh=True):

        cmd = "gmx pdb2gmx -f %s -o %s -p %s -ff %s -water %s " %\
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
        cmd = "gmx editconf -f %s -o %s -c -d 1.2 -bt cubic" % (ingro, outgro)

        self.run_suprocess(cmd)

        return self

    def add_solvent(self, ingro, outgro, intop="topol", spc="spc903.gro"):

        if not os.path.exists(spc):
            spc = os.path.join(self.PROJECT_ROOT, "data/spc903.gro")

        cmd = "gmx solvate -cp %s -cs %s -o %s -p %s " % (ingro, spc, outgro, intop)
        self.run_suprocess(cmd)

        return self

    def add_ions(self, ingro, outgro, emmdp="em_sol.mdp", intop="topol", ion_conc=0.15):
        if not os.path.exists(emmdp):
            emmdp = os.path.join(self.PROJECT_ROOT, "data/em_sol.mdp")

        cmd1 = "gmx grompp -f %s -c %s -p %s -o %s -maxwarn 100" % (emmdp, ingro, intop, "addion")
        self.run_suprocess(cmd1)

        cmd2 = "echo \"13 0 0 \" | gmx genion -s %s -p %s -o %s -neutral -conc %.2f " %\
               ("addion", intop, outgro, ion_conc)
        self.run_suprocess(cmd2)

        self.top = intop

        return self

    def minimize(self, ingro, outgro, emmdp="em_sol.mdp", intop="topol", nt=4):
        if not os.path.exists(emmdp):
            emmdp = os.path.join(self.PROJECT_ROOT, "data/em_sol.mdp")

        cmd1 = "gmx grompp -f %s -c %s -p %s -o %s -maxwarn 100" % (emmdp, ingro, intop, outgro)
        self.run_suprocess(cmd1)

        cmd2 = "gmx mdrun -deffnm %s -nt %d -v -gpu_id 0" % (outgro, nt)
        self.run_suprocess(cmd2)

        return self

    def md(self, ingro, outgro, nptmdp="npt.mdp", intop="topol", nt=4, restraints=False):

        if not restraints:
            mdp = nptmdp
            if not os.path.exists(nptmdp):
                mdp = os.path.join(self.PROJECT_ROOT+"/data/", nptmdp)
        else:
            mdp = "npt_rest.mdp"
            if not os.path.exists(mdp):
                mdp = os.path.join(self.PROJECT_ROOT+"/data/", mdp)

        cmd1 = "gmx grompp -f %s -c %s -p %s -o %s -maxwarn 100" % (mdp, ingro, intop, outgro)
        self.run_suprocess(cmd1)

        cmd2 = "gmx mdrun -deffnm %s -nt %d -v -gpu_id 1" % (outgro, nt)
        self.run_suprocess(cmd2)

        print("MD Simulation completed. ")

        return self

    def run_app(self, inpdb, outname, mode="solvated", production_run=True):

        if mode == "solvated":
            self.generate_top(inpdb, outgro=outname, top="topol")
            self.add_box(ingro=outname, outgro="box_"+outname)
            self.add_solvent(ingro="box_"+outname, outgro="wat_"+outname, )
            self.add_ions(ingro="wat_"+outname, outgro="ion_"+outname)
            self.minimize(ingro="ion_"+outname, outgro="em_"+outname)

            if production_run:
                self.md("em_"+outname, outgro="npt_"+outname, nptmdp="npt.mdp")

        elif mode == "gbsa":
            self.generate_top(inpdb, outgro=outname, top="topol")
            self.add_box(ingro=outname, outgro="box_"+outname)
            self.minimize(ingro="box_"+outname, outgro="em_"+outname, emmdp="em_sol.mdp")

            if production_run:
                mdp = os.path.join(self.PROJECT_ROOT, "data/gbsa.mdp")
                self.md("em_"+outname, outgro="npt_"+outname, nptmdp=mdp)

        return self

if __name__ == "__main__":

    d = """
    Run gromacs simulation in one-liner.
    """
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-f", type=str, default="input.pdb",
                        help="Input, str. The input pdb file. Default is input.pdb")
    parser.add_argument("-o", type=str, default="protein",
                        help="Input, str, optional. The output file name pattern. \n"
                             "Default is protein ")
    parser.add_argument("-gbsa", type=lambda x: (str(x).lower() == "true"), default=False,
                        help="Input, bool, optional. Run GBSA implict water simulation. \n"
                             "Default is False. ")
    parser.add_argument("-product", type=lambda x: (str(x).lower() == "true"), default=True,
                        help="Input, bool, optional. Run production MD simulations. \n"
                             "Default is True. ")

    args = parser.parse_args()

    app = AutoRunMD()

    if args.gbsa:
        mode = "gbsa"
    else:
        mode = "solvated"

    app.run_app(inpdb=args.f, outname=args.o, gbsa=mode, production_run=args.product)

