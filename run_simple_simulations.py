from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
#from sys import stdout

# read in a pdb file
pdb = PDBFile("input.pdb")

# load force field
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

# prepare the simulation system
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)

# setup MD integrator
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

# use cuda device
platform = mm.Platform.getPlatformByName('CUDA')

simulation = Simulation(pdb.topology, system, integrator, platform)

simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()

# output information
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter("output.log", 1000, step=True,
        potentialEnergy=True, temperature=True, progress=True))


simulation.step(10000)
