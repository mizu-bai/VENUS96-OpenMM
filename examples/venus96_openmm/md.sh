#!/bin/bash

# load gmx
source /opt/gromacs-2019.6/bin/GMXRC

# run md
../../src/venus96_openmm.x < ch4.dt5 | tee ch4.out
