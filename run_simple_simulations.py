from simtk.openmm.app import *
from simtk.openmm import *
from simtk import openmm as mm
from simtk import unit as u
#from sys import stdout
from datetime import datetime
import argparse
from argparse import RawDescriptionHelpFormatter


def dihedral_restraints(rst="restraints.txt", use_cols=[0, 1, 2, 3, -2, -1]):
    pass


def process_pdb(inp, addH=False, ):
    # read in a pdb file
    pdb = PDBFile(inp)

    return pdb

def prepare_system(pdb, gbsa=False):

    if gbsa:

        forcefield = ForceField('amber14-all.xml') #, 'amber14/tip3pfb.xml')

        system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
                                         nonbondedCutoff=1*u.nanometer, constraints=HBonds,
                                         implicitSolvent=app.GBn2,
                                         implicitSolventSaltConc=0.1*u.moles/u.liter,)

    else:
        # load force field
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

        # prepare the simulation system by provide topology and algorithms for constraints
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
                                         nonbondedCutoff=1*u.nanometer,
                                         constraints=HBonds)

    return system

def run_NVT(pdb, system, out, log, nsteps, temperature, ):

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
    simulation.reporters.append(PDBReporter(out, 1000))
    simulation.reporters.append(StateDataReporter(log, 1000, step=False,
                                potentialEnergy=True, temperature=True, progress=False))

    simulation.step(nsteps)

    print("Simulation completed!")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="A simple peptide simulation function. ")
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

    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()

    now = datetime.now()

    pdb = process_pdb(args.pdbin)

    system = prepare_system(pdb, gbsa=args.gbsa)

    run_NVT(pdb, system=system, out=args.pdbout, log=args.logfile,
            nsteps=args.nsteps, temperature=args.temperature)

    time_used = datetime.now() - now
    print("Total Time Usage: ", time_used)
