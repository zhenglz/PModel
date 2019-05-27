import simtk.openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import sys

pdb = PDBFile('input.pdb')

forcefield = ForceField('amber14-all.xml') #, 'amber14/tip3pfb.xml')

system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds,
        implicitSolvent=app.GBn2,
        implicitSolventSaltConc=0.1*moles/liter,)

# constraint the system with specific atoms
flat_bottom_force = mm.CustomTorsionForce("0.5*k*(theta-theta0)^2")
flat_bottom_force.addPerTorsionParameter('theta',)
flat_bottom_force.addPerTorsionParameter('k',)

with open('restraints.txt') as input_file:
    for line in input_file:
        columns = line.split()
        atom_index_i = int(columns[0])
        atom_index_j = int(columns[1])
        atom_index_k = int(columns[2])
        atom_index_l = int(columns[3])

        r0 = float(columns[-3])
        k = float(columns[-1])

        flat_bottom_force.addTorsion(atom_index_i, atom_index_j, atom_index_k, atom_index_l, [r0, k])

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

simulation = Simulation(pdb.topology, system, integrator)

simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output.pdb', 1000))

simulation.reporters.append(StateDataReporter("output.log", 1000, step=True,
        potentialEnergy=True, temperature=True))

simulation.step(10000000)
