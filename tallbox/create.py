""" @package ./examples/Noh_3d/create.py
Code that creates 3d Noh test problem initial conditions

created by Rainer Weinberger, last modified 24.02.2019
"""

""" load libraries """
import sys    # system specific calls
import numpy as np    # scientific computing package
import h5py    # hdf5 format

simulation_directory = str(sys.path[0])
print("tallbox/create.py: creating ICs in directory " + simulation_directory)

""" initial condition parameters """
FilePath = simulation_directory + '/IC.hdf5'

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32

Boxsize = FloatType(1.0)
CellsPerDimension = IntType(30)
NumberOfCells = CellsPerDimension * CellsPerDimension * CellsPerDimension*6

## initial state
code_units_dens = 2.468e7

gamma = 5./3.  ## note: this has to be consistent with the parameter settings for Arepo!


""" set up grid: cartesian 3d grid """
## spacing
dx = Boxsize / FloatType(CellsPerDimension)
## position of first and last cell
pos_first, pos_last = 0.5 * dx, Boxsize - 0.5 * dx
pos_firstz, pos_lastz = 0.5*dx, 6*Boxsize - 0.5*dx

## set up grid
Grid1d = np.linspace(pos_first, pos_last, CellsPerDimension, dtype=FloatType)
Grid1dz = np.linspace(pos_firstz, pos_lastz, 6*CellsPerDimension, dtype=FloatType)
xx, yy, zz = np.meshgrid(Grid1d, Grid1d, Grid1dz)
Pos = np.zeros([NumberOfCells, 3], dtype=FloatType)
Pos[:,0] = xx.reshape(NumberOfCells)
Pos[:,1] = yy.reshape(NumberOfCells)
Pos[:,2] = zz.reshape(NumberOfCells)
## calculate distance from center
#xPosFromCenter = (Pos[:,0] - 0.5 * Boxsize)
#Radius = np.sqrt( xPosFromCenter**2 + yPosFromCenter**2 + zPosFromCenter**2 )
zPosFromCenter = (Pos[:,2] - 0.5 * 6*Boxsize)
#yPosFromCenter = (Pos[:,1] - 0.5 * Boxsize)

def dl(z):
    term1 = 0.47*np.exp(-0.5*(z/0.09)**2)
    term2 = 0.13*np.exp(-0.5*(z/0.225)**2)
    term3 = 0.077*np.exp(-1*(np.abs(z)/0.403))
    term4 = 0.025*np.exp(-np.abs(z)/1.0)
    return term1+term2+term3+term4

""" set up hydrodynamical quantitites """
## mass insetad of density
Masses = np.zeros(NumberOfCells, dtype=FloatType)
Masses = code_units_dens*dl(zPosFromCenter)*dx*dx*dx
## velocity
Velocity = np.zeros([NumberOfCells,3], dtype=FloatType)
#Velocity[:,0] = velocity_radial_0 * xPosFromCenter / Radius
#Velocity[:,1] = velocity_radial_0 * yPosFromCenter / Radius
#Velocity[:,2] = velocity_radial_0 * zPosFromCenter / Radius
## specific internal energy
#Uthermal = np.full(NumberOfCells, utherm_0, dtype=FloatType)

""" write *.hdf5 file; minimum number of fields required by Arepo """
IC = h5py.File(simulation_directory+'/IC.hdf5', 'w')

## create hdf5 groups
header = IC.create_group("Header")
part0 = IC.create_group("PartType0")

## header entries
NumPart = np.array([NumberOfCells, 0, 0, 0, 0, 0], dtype=IntType)
header.attrs.create("NumPart_ThisFile", NumPart)
header.attrs.create("NumPart_Total", NumPart)
header.attrs.create("NumPart_Total_HighWord", np.zeros(6, dtype=IntType) )
header.attrs.create("MassTable", np.zeros(6, dtype=IntType) )
header.attrs.create("Time", 0.0)
header.attrs.create("Redshift", 0.0)
header.attrs.create("BoxSize", Boxsize)
header.attrs.create("NumFilesPerSnapshot", 1)
header.attrs.create("Omega0", 0.0)
header.attrs.create("OmegaB", 0.0)
header.attrs.create("OmegaLambda", 0.0)
header.attrs.create("HubbleParam", 1.0)
header.attrs.create("Flag_Sfr", 0)
header.attrs.create("Flag_Cooling", 0)
header.attrs.create("Flag_StellarAge", 0)
header.attrs.create("Flag_Metals", 0)
header.attrs.create("Flag_Feedback", 0)
if Pos.dtype == np.float64:
    header.attrs.create("Flag_DoublePrecision", 1)
else:
    header.attrs.create("Flag_DoublePrecision", 0)

## copy datasets
part0.create_dataset("ParticleIDs", data=np.arange(1, NumberOfCells+1) )
part0.create_dataset("Coordinates", data=Pos)
part0.create_dataset("Masses", data=Masses)
part0.create_dataset("Velocities", data=Velocity)
#part0.create_dataset("InternalEnergy", data=Uthermal)

## close file
IC.close()