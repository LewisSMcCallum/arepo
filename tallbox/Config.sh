#!/bin/bash            # this line only there to enable syntax highlighting in this file

## examples/Noh_3d/Config.sh
## config file for 3d Noh probelm

#--------------------------------------- Basic operation mode of code


LONG_Z=6.0
REFLECTIVE_X=2
REFLECTIVE_Y=2
REFLECTIVE_Z=2                           # in-/outflow boundary conditions in z direction

#--------------------------------------- Static Isothermal Sphere Potential
#EXTERNALGRAVITY
#STATICISO                     # static gravitational isothermal sphere potential
#ISO_M200=10000                # mass causing the isothermal sphere potential
#ISO_R200=0.05                # radius of the isothermal sphere potential
#ISO_Eps=0.1                   # softening of isothermal sphere potential
#ISO_FRACTION=0.9              # fraction in dark matter in isothermal sphere potential

SELFGRAVITY
COOLING
GRAVITY_NOT_PERIODIC
REFINEMENT_SPLIT_CELLS
REFINEMENT_MERGE_CELLS

USE_SFR
SFR_KEEP_CELLS
#--------------------------------------- Mesh motion and regularization
REGULARIZE_MESH_CM_DRIFT                 # Mesh regularization; Move mesh generating point towards center of mass to make cells rounder.
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED  # Limit mesh regularization speed by local sound speed
REGULARIZE_MESH_FACE_ANGLE               # Use maximum face angle as roundness criterion in mesh regularization

#--------------------------------------- Time integration options
TREE_BASED_TIMESTEPS                     # non-local timestep criterion (take 'signal speed' into account)

#---------------------------------------- Single/Double Precision
DOUBLEPRECISION=1                        # Mode of double precision: not defined: single; 1: full double precision 2: mixed, 3: mixed, fewer single precisions; unless short of memory, use 1.
INPUT_IN_DOUBLEPRECISION                 # initial conditions are in double precision
OUTPUT_CENTER_OF_MASS                    # output centers of cells

#--------------------------------------- Output/Input options
HAVE_HDF5                                # needed when HDF5 I/O support is desired; should this be standard?

#--------------------------------------- Testing and Debugging options
DEBUG                                    # enables core-dumps, should this be standard?
