# PModel
## Protein/Peptide modeling with full-atomic simulations in explict/implict environment

# Author Information:
#### Zheng Liangzhen, zhenglz@outlook.com

# 1. Installation
Provide that you are using python 3.6, you should have the following python libraries:

    numpy
    pandas
    openmm

## To install openmm, use the following command:

    >>> conda install -c omnia -c conda-forge openmm

# 2. Usage
You could refer to the following examples to run the simulations:

    >>> # run a quick MD simulation
    >>> python autoRunMD.py --pdbin input.pdb -nsteps 10000
    >>> # run MD simulation with GBSA method
    >>> python autoRunMD.py --pdbin input.pdb --nsteps 10000 --gbsa True
    >>> # run MD simulation with dihedral restraints
    >>> python autoRunMD.py --pdbin input.pdb --nsteps 100000 --gbsa True --diherst restraints.txt

# 3. Bug and issues
Please report the bugs and issues if you find any. It is also welcomed to discuss about the codes.

