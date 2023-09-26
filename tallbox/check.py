""" @package ./examples/Noh_3d/check.py
Code that checks results of 3d Noh problem

created by Rainer Weinberger, last modified 24.02.2019
"""

""" load libraries """
import sys    ## load sys; needed for exit codes
import numpy as np    ## load numpy
import h5py    ## load h5py; needed to read snapshots
import os      # file specific calls
import matplotlib.pyplot as plt    # needs to be active for plotting!


makeplots = True
if len(sys.argv) > 2:
  if sys.argv[2] == "True":
    makeplots = True
  else:
    makeplots = False

simulation_directory = str(sys.path[0])
print("tallbox/check.py: checking simulation output in directory " + simulation_directory) 

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32

## open initial conditiions to get parameters
try:
    data = h5py.File(simulation_directory + "/IC.hdf5", "r")
except:
    print("could not open initial  conditions!")
    exit(1)
Boxsize = FloatType(data["Header"].attrs["BoxSize"])
NumberOfCells = IntType(data["Header"].attrs["NumPart_Total"][0]) 
CellsPerDimension = (NumberOfCells/6)**(1./3.) ## 3d sim

## parameters for initial state
#density_0 = 1.0
#velocity_radial_0 = -1.0    ## radial inflow velocity
#pressure_0 = 1.0e-4
gamma = 5./3.  ## note: this has to be consistent with the parameter settings for Arepo!
#utherm_0 = pressure_0 / ( gamma - 1.0 ) / density_0

## maximum L1 error after one propagation; empirically based
DeltaMaxAllowed = 0.25

""" loop over all output files """
i_file = 0
while True:
    """ try to read in snapshot """
    directory = simulation_directory+"/snaps/"
    filename = "snap_%03d.hdf5" % (i_file)
    try:
        data = h5py.File(directory+filename, "r")
    except:
        break
    
    """ get simulation data """
    
    time = FloatType(data["Header"].attrs["Time"])
    ## simulation data
    Pos = np.array(data["PartType0"]["Coordinates"], dtype = FloatType) # CenterOfMass
    VoronoiPos = np.array(data["PartType0"]["Coordinates"], dtype = FloatType)
    Density = np.array(data["PartType0"]["Density"], dtype = FloatType)
    Mass = np.array(data["PartType0"]["Masses"], dtype = FloatType)
    Velocity = np.array(data["PartType0"]["Velocities"], dtype = FloatType)
    Uthermal = np.array(data["PartType0"]["InternalEnergy"], dtype = FloatType)
    CellVolume = Mass / Density
    
    xPosFromCenter = (Pos[:,0] - 0.5 * Boxsize)
    yPosFromCenter = (Pos[:,1] - 0.5 * Boxsize)
    zPosFromCenter = (Pos[:,2] - 0.5 * 6*Boxsize)
    Radius = np.sqrt( xPosFromCenter**2 + yPosFromCenter**2 + zPosFromCenter**2 )
    
    vRad = Velocity[:,0] * xPosFromCenter / Radius \
      + Velocity[:,1] * yPosFromCenter / Radius \
      + Velocity[:,2] * zPosFromCenter / Radius
    
    
    

    

    
    #### plot profiles
    if makeplots:
        fig = plt.figure( figsize=np.array([7.5,6.0]), dpi=300 )

        
        Nplot = 256
        Nplotz = 4*Nplot
        from scipy import spatial # needed for KDTree that we use for nearest neighbour search and Voronoi mesh
        Edges1d = np.linspace(0, Boxsize, Nplot+1, endpoint=True, dtype=FloatType)
        Edgesz = np.linspace(1,5,Nplotz+1, endpoint = True, dtype = FloatType)
        Grid1d = 0.5 * (Edges1d[1:] + Edges1d[:-1])
        Grid1dz = 0.5* (Edgesz[1:] + Edgesz[:-1])
        xx, yy = np.meshgrid(Grid1d, Grid1d)
        xx2,zz2 = np.meshgrid(Grid1d,Grid1dz)
        Grid2D = np.array( [xx.reshape(Nplot**2), yy.reshape(Nplot**2), np.ones(Nplot**2)*0.5*6*Boxsize] ).T
        dist, cells = spatial.KDTree( VoronoiPos[:] ).query( Grid2D, k=1 )

        Grid2D2 = np.array( [xx2.reshape(Nplot*Nplotz), np.ones(Nplot*Nplotz)*0.5*Boxsize, zz2.reshape(Nplot*Nplotz)] ).T
        dist2, cells2 = spatial.KDTree( VoronoiPos[:] ).query( Grid2D2, k=1 )

      
        ax  = plt.axes( [0.65,0.70,0.20,0.25] )
        pc  = ax.pcolormesh( Edges1d, Edges1d, Density[cells].reshape((Nplot,Nplot)), rasterized=True, cmap=plt.get_cmap('viridis') )
        cax = plt.axes( [0.88,0.70,0.02,0.25] )
        plt.colorbar( pc, cax=cax )
        ax.set_xlim( 0., Boxsize)
        ax.set_ylim( 0., Boxsize)

        ax  = plt.axes( [0.1,0.1,0.2,0.8] )
        pc  = ax.pcolormesh( Edges1d, Edgesz, Density[cells2].reshape((Nplotz,Nplot)), rasterized=True, cmap=plt.get_cmap('viridis') )
        cax = plt.axes( [0.3,0.1,0.02,0.8] )
        plt.colorbar( pc, cax=cax )
        #ax.set_xlim( 0., Boxsize)
        #ax.set_ylim( 0., 4*Boxsize)
       
        ax  = plt.axes( [0.65,0.38,0.20,0.25] )
        pc  = ax.pcolormesh( Edges1d, Edges1d, vRad[cells].reshape((Nplot,Nplot)), rasterized=True, cmap=plt.get_cmap('plasma'))
        cax = plt.axes( [0.88,0.38,0.02,0.25] )
        plt.colorbar( pc, cax=cax )
        ax.set_xlim( 0., Boxsize)
        ax.set_ylim( 0., Boxsize)

        #ax  = plt.axes( [0.1,0.38,0.2,0.25] )
        #pc  = ax.pcolormesh( Edgesz, Edges1d, vRad[cells2].reshape((Nplot,Nplotz)), rasterized=True, cmap=plt.get_cmap('viridis') )
        #cax = plt.axes( [0.3,0.38,0.02,0.25] )
        #plt.colorbar( pc, cax=cax )
        #ax.set_xlim( 0., Boxsize)
        #ax.set_ylim( 0., 4*Boxsize)
        
        
        ax  = plt.axes( [0.65,0.07,0.20,0.25] )
        pc  = ax.pcolormesh( Edges1d, Edges1d, Uthermal[cells].reshape((Nplot,Nplot)), rasterized=True, cmap=plt.get_cmap('magma') )
        cax = plt.axes( [0.88,0.07,0.02,0.25] )
        plt.colorbar( pc, cax=cax )
        ax.set_xlim( 0., Boxsize)
        ax.set_ylim( 0., Boxsize)

       # ax  = plt.axes( [0.1,0.07,0.2,0.25] )
       # pc  = ax.pcolormesh( Edgesz, Edges1d, Uthermal[cells2].reshape((Nplot,Nplotz)), rasterized=True, cmap=plt.get_cmap('viridis') )
       # cax = plt.axes( [0.3,0.07,0.02,0.25] )
       # plt.colorbar( pc, cax=cax )
       # ax.set_xlim( 0., Boxsize)
       # ax.set_ylim( 0., 4*Boxsize)
        
        if not os.path.exists( simulation_directory+"/plots" ):
          os.mkdir( simulation_directory+"/plots" )

        print(simulation_directory+"/plots/figure_%03d.pdf" % (i_file) )
        fig.savefig(simulation_directory+"/plots/figure_%03d.pdf" % (i_file), dpi=300)
        plt.close(fig)

    #### check against analytic solution
    
    print("tallbox/check snapshot %d!" % (i_file) )
    
    
    
    i_file += 1
    

## if everything is ok
sys.exit(0)
