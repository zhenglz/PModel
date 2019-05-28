from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
#from sys import stdout
import simtk.openmm as mm
import pandas as pd


class ParseConstraints:

    def __init__(self, cst_type="dihedral", filen="restraint.txt"):
        self.cst_type = cst_type
        self.cst_file = filen
        self.cst_df = pd.DataFrame()
        self.forces = None

    def read_restraints(self):
        with open(self.cst_file) as lines:

            self.cst_df = pd.DataFrame([x.split() for x in lines])

        return self

    def force_type(self):
        if self.cst_type == "dihedral":
            self.forces = mm.CustomTorsionForce('0.5*k*(theta-theta0)^2')
            self.forces.addPerTorsionParameter('theta',)
            self.forces.addPerTorsionParameter('k',)

        elif self.cst_type == "distance":
            self.forces = mm.CustomBondForce('0.5*k*(r-r0)^2')
            self.forces.addPerTorsionParameter('theta',)
            self.forces.addPerTorsionParameter('k',)

        else :

            self.forces == mm.CustomAngleForce('0.5*k*(d-d0)^2')

        return self

    def add_forces(self):
        if self.cst_type == "dihedral":
            for items in self.cst_df.values:
                atom_index_i = int(items[0])
                atom_index_j = int(items[1])
                atom_index_k = int(items[2])
                atom_index_l = int(items[3])

                r0 = float(items[-3])
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


if __name__ == "__main__":
    # read in a pdb file
    pdb = PDBFile("input.pdb")

    # load force field
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    # prepare the simulation system
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
            nonbondedCutoff=1*nanometer, constraints=HBonds)


    parsecont = ParseConstraints(filen='restraints.txt', cst_type='dihedral')
    parsecont.read_restraints()
    parsecont.force_type()
    parsecont.add_forces()

    # This line assumes that an openmm System object
    # has been created.
    # system = ...
    system.addForce(parsecont.forces)

    # setup MD integrator
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

    # use cuda device
    #platform = mm.Platform.getPlatformByName('CUDA')

    simulation = Simulation(pdb.topology, system, integrator)#, platform)

    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy()

    # output information
    simulation.reporters.append(PDBReporter('output.pdb', 1000))
    simulation.reporters.append(StateDataReporter("output.log", 1000, step=True,
            potentialEnergy=True, temperature=True,))


    simulation.step(10000)
