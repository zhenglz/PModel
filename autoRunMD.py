from simtk.openmm.app import *
from simtk.openmm import *
from simtk import openmm as mm
from simtk import unit as u
from datetime import datetime
import pandas as pd
import argparse
import os
from argparse import RawDescriptionHelpFormatter

class ParseConstraints:

    def __init__(self, cst_type="dihedral", filen="restraint.txt"):
        self.cst_type = cst_type
        self.cst_file = filen
        self.cst_df = pd.DataFrame()
        self.forces = None

    def read_restraints(self, usecols):
        dat = []
        with open(self.cst_file) as lines:
            for s in lines:
                l = []
                for i in usecols:
                    l.append(s.split()[i])
                dat.append(l)

        self.cst_df = pd.DataFrame(dat)

        return self

    def force_type(self):
        if self.cst_type == "dihedral":
            self.forces = mm.CustomTorsionForce('0.5*k*(theta-theta0)^2')
            self.forces.addPerTorsionParameter('theta0',)
            self.forces.addPerTorsionParameter('k',)

        elif self.cst_type == "distance":
            self.forces = mm.CustomBondForce('0.5*k*(r-r0)^2')
            self.forces.addPerTorsionParameter('r0',)
            self.forces.addPerTorsionParameter('k',)

        return self

    def add_forces(self):
        if self.cst_type == "dihedral":
            for items in self.cst_df.values:
                atom_index_i = int(items[0])
                atom_index_j = int(items[1])
                atom_index_k = int(items[2])
                atom_index_l = int(items[3])

                r0 = 3.1415 * float(items[-2]) / 180.0
                k = float(items[-1])

                self.forces.addTorsion(atom_index_i, atom_index_j,
                                             atom_index_k, atom_index_l,
                                             [r0, k])
        if self.cst_type == "distance":
            for items in self.cst_df.values:
                atom_index_i = int(items[0])
                atom_index_j = int(items[1])
                #atom_index_k = int(items[2])
                #atom_index_l = int(items[3])

                r0 = float(items[-3])
                k = float(items[-1])

                self.forces.addBond(atom_index_i, atom_index_j, [r0, k])

        return self

def dihedral_restraints(rst="restraints.txt", use_cols=[0, 1, 2, 3, -3, -1], rst_type="dihedral"):

    restraints = ParseConstraints(rst_type, rst)
    restraints.read_restraints(use_cols)
    restraints.force_type()
    restraints.add_forces()

    return restraints.forces


def process_pdb(inp, ):
    # read in a pdb file
    pdb = PDBFile(inp)

    return pdb

def prepare_system(pdb, gbsa=False, add_forces=None, add_hydrogens=False):

    # load force field
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    if gbsa:
        pdb = Modeller(pdb.topology, pdb.positions)

        if add_hydrogens:
            pdb.addHydrogens(forcefield)

        system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
                                         nonbondedCutoff=1*u.nanometer, constraints=HBonds,
                                         implicitSolvent=app.GBn2,
                                         implicitSolventSaltConc=0.1*u.moles/u.liter,)
        if add_forces is not None:
            system.addForce(add_forces)

        return system, pdb

    else:
        modeller = Modeller(pdb.topology, pdb.positions)

        # add hydrogen atoms when necessary
        if add_hydrogens:
            modeller.addHydrogens(forcefield)

        modeller.addSolvent(forcefield, padding=1.0*u.nanometers, model="tip3p",
                            ionicStrength=0.15*u.molar, negativeIon="Cl-", positiveIon="Na+")

        # prepare the simulation system by provide topology and algorithms for constraints
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,
                                         nonbondedCutoff=1*u.nanometer,
                                         constraints=HBonds)

        if add_forces is not None:
            system.addForce(add_forces)

        return system, modeller

def run_NPT(pdb, system, out, log, nsteps, temperature, gbsa=True):

    if not gbsa:
        system.addForce(MonteCarloBarostat(1*u.bar, 300*u.kelvin))

    # Temperature coupling
    system.addForce(AndersenThermostat(300*u.kelvin, 1/u.picosecond))


    # setup MD integrator
    integrator = LangevinIntegrator(temperature * u.kelvin, 1/u.picosecond, 0.002*u.picoseconds)

    # use cuda device
    #platform = mm.Platform.getPlatformByName('CUDA:0')

    # create a simulation object, provide the topology and other information
    simulation = Simulation(pdb.topology, system, integrator) #, platform)

    # provide the initial simulation start point
    simulation.context.setPositions(pdb.positions)

    # run energy minimization
    simulation.minimizeEnergy()

    # output information
    simulation.reporters.append(PDBReporter(out, 5000))
    simulation.reporters.append(StateDataReporter(log, 500, step=True,
                                potentialEnergy=True, temperature=True,
                                density=True, progress=True,
                                remainingTime=True, speed=True,
                                totalSteps=nsteps))
    simulation.reporters.append(CheckpointReporter('checkpoint.chk', 50000))

    simulation.step(nsteps)
    #simulation.saveState("state.xml")

    print("Simulation completed!")

if __name__ == "__main__":
    d = """A simple peptide simulation function with OpenMM.
    Examples:
        
        >>> # run a quick MD simulation
        >>> python autoRunMD.py --pdbin input.pdb -nsteps 10000
        >>> # run MD simulation with GBSA method
        >>> python autoRunMD.py --pdbin input.pdb --nsteps 10000 --gbsa True
        >>> # run MD simulation with dihedral restraints
        >>> python autoRunMD.py --pdbin input.pdb --nsteps 100000 --gbsa True --diherst restraints.txt
    """

    parser = argparse.ArgumentParser(description=d,
                                     formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("--pdbin", default="input.pdb", type=str,
                        help="Input, string. A pdb file for simulation. \n"
                             "Default is input.pdb")
    parser.add_argument("--pdbout", default="output.pdb", type=str,
                        help="Output, string. A pdb file for output trajectory. \n"
                             "Default is output.pdb")
    parser.add_argument("--logfile", default="output.log", type=str,
                        help="Output, string. A log file containing the simulation information. \n"
                             "Default is output.log")
    parser.add_argument("--nsteps", default=10000000, type=int,
                        help="Input, int, optional. Number of simulation steps. A step is 0.002 ps.\n"
                             "Default is 10000000, which makes 20 ns.")
    parser.add_argument("--gbsa", type=lambda x: (str(x).lower() == "true"), default=False,
                        help="Input, bool, optional. Whether run implict GBSA simulation.\n"
                             "Default is True. ")
    parser.add_argument("--temperature", default=300, type=int,
                        help="Input, int, optional. The simulation temperature, default is 300 K. ")
    parser.add_argument("--diherst", default=None, type=str,
                        help="Input, str, optional. Apply dihedral restraints to the system. ")

    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()

    now = datetime.now()

    # prepare the pdb file
    pdb = process_pdb(args.pdbin)

    if args.diherst is not None and os.path.exists(args.diherst):
        forces = dihedral_restraints(args.diherst, rst_type="dihedral")
    else:
        forces = None

    # prepare the system
    system, pdb = prepare_system(pdb, gbsa=args.gbsa, add_forces=forces)

    # run NVT simulations
    run_NPT(pdb, system=system, out=args.pdbout, log=args.logfile,
            nsteps=args.nsteps, temperature=args.temperature, gbsa=args.gbsa)

    time_used = datetime.now() - now
    print("Total Time Usage: ", time_used)
