.. _Parameterfile:

Parameterfile
***************** 

The parameter-file, often named param.txt, is a file containing run-time 
options for AREPO. These include things like the input and output directory,
maximum runtimes, memory limits and all kind of freely choosable simulaiton 
and model parameters. In general, the code will output an error message if 
there are either missing paramters for the given configuration AREPO was 
compiled with, or if there are obsolete paramters. The latter can be 
deactivated by setting the compile-time flag ``ALLOWEXTRAPARAMS``.
Unlike changing Config.sh, changing the parameters does not require 
re-compilation of the code.


Initial conditions 
==================

**InitCondFile**

  The filename of the initial conditions file. Can be a relative or absolute 
  path. The ICs can be distributed in one or more files, as with snapshots. 
  If using ICs with multiple files, only the basename without the trailing 
  ".n" should be specified here. Similarly, the ".hdf5" 
  extension can be omitted, and the code will automatically append it for 
  ``ICFormat=3``. If a restart from a snapshot with the "2" option is desired, 
  one needs to specify the snapshot file here.

-----

**ICFormat**

  The file format of the initial conditions. Currently, three different formats
  are supported, selected by one of the choices "1", "2", or "3". Format "1" 
  is the traditional fortran-style unformatted format familiar from GADGET. 
  Format "2" is a variant of this format, where each block of data is 
  preceeded by a 4-character block-identifier. Finally, format "3" selects the 
  HDF-5 format.

-----

**InitGasTemp**

  This sets the initial gas temperature (assuming either a mean molecular 
  weight corresponding to full ionization or full neutrality, depending on 
  whether the temperature is above or below 10^4 K) in Kelvin when initial 
  conditions are read. However, the gas temperature is only set to a certain
  temperature if ``InitGasTemp>0``, and if the temperature of the gas particles
  in the initial conditions file is zero, otherwise the initial gas 
  temperature is left at the value stored in the IC file.

-----

**MinimumDensityOnStartUp**

  This sets a lower limit to the density of gas cells after reading in. All 
  cells that have a lower density are set to this value.

-----

**MHDSeedDir**

``MHD_SEEDFIELD``

  Direction of the uniform B field that is set before starting the simulation.
  The direction is encoded by sum(2^k), where k is the index of direciton
  (0, 1, 2 for x, y, z, respectively). E.g. 3 is a diagonal field in 
  the xy plane, parallel to the z axis. This allows only orientations along
  coordinate axis, perpendicular to it or diagonal. Note that the equations of 
  ideal MHD do not change if the field is reversed.

-----

**MHDSeedValue**

``MHD_SEEDFIELD``

  Value of the uniform initial magnetic field in comoving Gauss.

-----

**TileICsFactor**

``TILE_ICS``

  Factor by which the ICs are dublicated in each dimension. Should be an 
  integer.

-----

**GridSize**

``ADDBACKGROUNDGRID``

  Initial guess of grid size for ``ADDBACKGROUNDGRID``. The input will be
  rounded to the nearest power of two.

-----



Output file names and formats 
=============================

**OutputDir**

  Pathname of the output directory of the code. Can be a relative or absolute 
  path. Path must exist.

-----

**SnapshotFileBase**

  Basename of snapshot files produced by the code. e.g. ``snap`` for output files
  ``snap_000.hdf5`` etc.

-----

**NumFilesPerSnapshot**

  The number of separate files requested for each snapshot dump. Each file of 
  the snapshot will hold the data of one or several processors, up to all of 
  them. ``NumFilesPerSnapshot`` must hence lie between 1 and the number of 
  processors used. Distributing a snapshot onto several files can be done 
  in parallel and may lead to much better I/O performance, depending on the 
  hardware configuration. It can also help to avoid problems due to big 
  files for large simulations. Note that initial conditions may also 
  be distributed into several files, the number of which is automatically 
  recognised by the code and does not have to be equal to 
  ``NumFilesPerSnapshot`` (it may also be larger than the number of 
  processors).

-----

**OutputListOn**

  If set to "1", the code tries to read a list of desired output times from the
  file given in ``OutputListFilename``. Otherwise, output times are generated 
  equally spaced from the values assigned from ``TimeOfFirstSnapshot`` onwards
  with spacing ``TimeBetSnapshot``.
  
-----
  
**OutputListFilename**

  File with a list of the desired output times. Can be specified with a 
  relative or absolute path.

-----

**SnapFormat**

  Similar to ``ICFormat``, this parameter selects the file-format of snapshot 
  dumps produced by the code. 1 and 2 are two binary formats, identical to 
  the ones in GADGET. 3 is an HDF5 output, which is recommended unless 
  there are good reasons not to use it.

-----

**NumFilesWrittenInParallel**

  The number of files the code may read or write simultaneously when writing 
  or reading snapshot/restart files. If the value of this parameter is larger 
  than the number of processors, it is capped by that. This parameter is only
  important for very large runs, where the file-system can be significantly 
  affected by too many tasks writing (restart files) at the same time.

-----

**AlternativeOutputDir**

``TOLERATE_WRITE_ERROR``

  Path name of an alternative output directory which is used in case output
  to OutputDir fails.

-----

Output frequency 
================

**CpuTimeBetRestartFile**

  The value specfied here gives the time in seconds the code will run before it
  writes regularly produced restart files. This can be useful to protect 
  against unexpected interruptions (for example due to a hardware problem) of 
  a simulation, particularly if it is run for a long time. It is then possible 
  to resume a simulation from the last restart file, reducing the potential 
  loss to the elapsed CPU-time since this was produced.

-----

**TimeBetSnapshot**

  The time interval in code units between two subsequent snapshot files in 
  case a file with output times is not specified. For cosmological 
  simulations, this is a multiplicative factor applied to the time of the 
  last snapshot, such that the snapshots will have a constant logarithmic 
  spacing in the scale factor. Otherwise, the parameter is an additive 
  constant that gives the linear spacing between snapshot times.

-----

**TimeBetStatistics**

  The code can be asked to measure the total kinetic, thermal, and potential 
  energy in regular intervals, and to write the results to the file given in 
  ``energy.txt``. The time interval between two such measurements
  is given by this parameter, in an analogous way as with ``TimeBetSnapshot``. 
  Note that the compile time option ``EVALPOTENTIAL`` needs to be activated to 
  obtain a measurement of the gravitational potential energy.

-----

**TimeOfFirstSnapshot**

  The time of the first desired snapshot file in code units in case a file with
  output times is not specified. For cosmological simulations, the value given 
  here is the scale factor of the first desired output.
  
-----
  
**FlushCpuTimeDiff**

  ``REDUCE_FLUSH``

  Time interval (in seconds) for flush calls on all log-files (i.e. save 
  write buffer to disk). In case ``Reduce_Flush`` is not set, this is done 
  during every sync-point.

-----
  

CPU-time limit and restarts 
===========================

**TimeLimitCPU**

  CPU-time limit for the present submission of the code. If 85 percent of this 
  time have been reached at the end of a timestep, the code terminates itself 
  and produces restart files. The extra 15% is used to guarantee that there is 
  enough time to safely finish the current time step and write the restart 
  files. This CPU time refers to the wall-lock time on a single 
  processor only.

-----

**ResubmitOn**

  If set to "1", the code will try to resubmit itself to the queuing system 
  when an interruption of the run due to the CPU-time limit occurs. The 
  resubmission itself is done by executing the program/script given with 
  ``ResubmitCommand``.

-----

**ResubmitCommand**

  The name of a script file or program that is executed for automatic 
  resubmission of the job to the queuing system. Note that the file given here 
  needs to be executable.

-----

Memory allocation 
=================

**MaxMemSize**

  The memory allocate per MPI task, in megabytes. A contiguous memory arena of 
  this total size is allocated at startup, and then partitioned internally 
  within AREPO for memory allocation and deallocation requests. Can generally 
  be set to ~95% of the total available, e.g. (memory per node / number of MPI 
  tasks per node), to leave room for operating system tasks and MPI buffers. 
  This value can be changed on a restart to increase the amount of memory 
  available to each task.

-----

Simulated time and spatial extent 
=================================

**BoxSize**

  The boxsize for the simulation, in internal code units. 
  All particles and gas cells in the ICs must have Coordinates 
  within the range ``[0,BoxSize]`` in each dimension. The only exception from 
  this is for collisionless particles in a tree-only gravity mode (no 
  ``PMGRID``) and ``GRAVITY_NONPERIODIC``.

-----

**PeriodicBoundariesOn**

  If set to "1", periodic boundary conditions are assumed, with a cubical 
  box-size of side-length ``BoxSize``. Particle coordinates are expected to 
  be in the range ``[0,BoxSize[``. Can only be set to zero if 
  ``GRAVITY_NOT_PERIODIC`` is set. Note: refers to gravity only! Hydrodynamic
  boundary conditions are handled by ``REFLECTIVE_X``, 
  ``REFLECTIVE_Y`` and ``REFLECTIVE_Z`` in ``Config.sh``

-----

**TimeBegin**

  This sets the starting time of a simulation when the code is started from 
  initial conditions in internal code units. For cosmological integrations, 
  the value specified here is taken as the initial scale factor.

-----

**TimeMax**

  This sets the final time for the simulation. The code normally tries to run 
  until this time is reached. For cosmological integrations, the value given 
  here is the final scale factor.

-----

Cosmological parameters 
=======================

**ComovingIntegrationOn**

  If set to "0", the code assumes plain Newtonian physics, with time, 
  positions, velocities, and masses measured in the internal system of units.
  If set to "1", the code assumes that a cosmological integration in comoving 
  coordinates should be carried out, assuming an expanding universe described 
  by the 'Cosmological parameters' below. In a cosmological integration, the 
  time variable is the scale factor.

-----

**Omega0**

  Gives the total matter density in units of the 
  critical density at z=0 for cosmological simulations. Relevant for comoving
  integration and halo/subhalo finder

-----

**OmegaBaryon**

  Gives the baryon density in units of the critical 
  densty at z=0 for cosmological simulations. Relevant for comoving 
  integration, halo/subhalo finder and star formation model.

-----

**OmegaLambda**

  Gives the vacuum energy density (cosmological 
  constant) at z=0 for cosmological simulations. Relevant for comoving 
  integration and halo/subhalo finder

-----

**HubbleParam**

  This gives the Hubble constant ("little h") at z=0 in units of 
  100 km/sec/Mpc.  Note that this parameter has been basically absorbed into 
  the definition of the internal code units, such that for gravitational 
  dynamics and adiabatic gas dynamics the actual value assigned for 
  ``HubbleParam`` is not used by the code. Only used when conversions to 
  physical cgs units are required (e.g. for radiative cooling physics). 
  In other cases, use 1.0.

-----

System of units 
===============

**UnitVelocity_in_cm_per_s**

  This sets the internal velocity unit in **cm/sec**. For example, the 
  choice of ``1e5`` sets the velocity unit to 1.0 *km/sec*. 
  Note that the specification of ``UnitLength_in_cm``, ``UnitMass_in_g``, and
  ``UnitVelocity_in_cm_per_s`` also determines the internal unit of time.

-----

**UnitLength_in_cm**

  This sets the internal length unit in **cm/h**, where H_0 = 100 h km/sec/Mpc. 
  For example, the choice of ``3.085678e21`` sets the length unit to 
  `1.0 kpc/h`.

-----

**UnitMass_in_g**

  This sets the internal mass unit in **g/h**, where H_0 = 100 h km/sec/Mpc.  
  For example, the choice of ``1.989e43`` sets the mass unit to 
  `10^10 Msun/h`.

-----

**GravityConstantInternal**

  The numerical value of the gravitational constant G in internal units depends
  on the system of units you choose. For example, for the choices above, 
  `G=43007.1` in internal units.  For ``GravityConstantInternal=0``, the code 
  calculates the value corresponding to the physical value of G automatically. 
  However, you might want to set G yourself.  For example, by specifying:
  ``GravityConstantInternal=1``, ``UnitLength_in_cm=1``, ``UnitMass_in_g=1``, 
  and ``UnitVelocity_in_cm_per_s=1``, one obtains a *natural* system of units. 
  Note that the code will nevertheless try to use the *correct* value of the 
  Hubble constant in this case, so you should not set 
  ``GravityConstantInternal`` in cosmological integrations.

-----

Gravitational force accuracy 
============================

**TypeOfOpeningCriterion**

  This selects the type of cell-opening criterion used in the tree walks. A 
  value of "1" selects the relative opening criterion. 
  *Required value: 1* (only implemented option).

-----

**ErrTolTheta**

  If the relative opening criterion is used, a first force estimate is 
  computed using the Barnes-Hut tree algorithm, which is then recomputed 
  with the relative opening criterion.

-----

**ErrTolForceAcc**

  The accuracy parameter for the relative opening criterion for the tree walk.
  Only used if ``ErrTolTheta 0``.

-----

Time integration accuracy 
=========================

**TypeOfTimestepCriterion**

  This parameter can in principle be used to select different kinds of 
  timestep criteria for gravitational dynamics. However, AREPO presently only 
  supports the criterion "0".

-----

**ErrTolIntAccuracy**

  This dimensionless parameter controls the accuracy of the timestep criterion 
  selected by ``TypeOfTimestepCriterion``. It is the variable ``eta``, 
  where the cosmological timestep for collisionless particles scales as 
  `dt ~ eta^0.5`.

-----

**MaxSizeTimestep**

  This gives the maximum timestep a particle may take. This should be set to a 
  sensible value in order to protect against too large timesteps for particles 
  with very small acceleration. For cosmological simulations, the parameter 
  given here is the maximum allowed step in the logarithm of the expansion 
  factor. For comoving runs, this is in units of ``dln(a)``.

-----

**MinSizeTimestep**

  If a particle requests a timestep smaller than the value specified here, the 
  code will normally terminate with a warning message. If compiled with the 
  ``NOSTOP_WHEN_BELOW_MINTIMESTEP`` option, the code will instead force the 
  timesteps to be at least as large as ``MinSizeTimestep``.

-----

**CourantFac**

  This sets the value of the Courant parameter used in the determination of the
  hydrodynamical timestep of gas cells. The hydrodynamical timestep is this 
  value times the Courant–Friedrichs–Lewy (CFL) condition calculated for each 
  cell.

-----


Domain decomposition 
====================

**ActivePartFracForNewDomainDecomp**

  Fraction of particles that need to be at least active in order to trigger a 
  new domain decomposition at a sync-point. All sync-points where fewer
  particles are active will not perform a domain decomposition.

-----

**MultipleDomains**

  Number of domains per MPI task. Consequently, the domain decomposition will
  cut computational box in ``MultipleDomains`` times number of tasks chunks.
  Too few of them will lead to cpu load and memory inbalances, too many to more
  MPI communication as there are more domain boundaries.

-----

**TopNodeFactor**

  Determines how deep the top-level tree is extended (1.0: only as far as 
  necessary to split domain in the required number of chunks). The higher this 
  factor is, the more precise the Peano-Hilbert curve can be cut into equal 
  pieces of cpu and memory load. A higher value, however, increases the global 
  tree which is stored on every MPI task, thus the total memory requirements.

-----

Moving mesh 
===========

**DesNumNgb**

  This sets the desired number of nearest neighbors for an initial density/size
  estimate for gas cells during code startup.

-----

**MaxNumNgbDeviation**

  This sets the allowed variation of the number of neighbours around the target 
  value ``DesNumNgb``. A larger tolerance will reduce the number of iterations.
  to find the correct radius.

-----

**MaxVolumeDiff**

``REFINEMENT_VOLUME_LIMIT``

  Maximum difference in volume of two neighboring cells. This avoids large
  cell-size gradients in the mesh which might cause numerical inaccuracies. 
  In case where the volume of a cell exceeds ``MaxVolumeDiff`` times the volume
  of the smallest neighboring cell, the cell is refined, irrespective of other
  refinement criteria.

-----

**MinVolume**

``REFINEMENT_VOLUME_LIMIT``

  Global minimal volume a cell is allowed to have, irrespective of other 
  refinement and derefinement criteria.

-----

**MaxVolume**

``REFINEMENT_VOLUME_LIMIT``

  Global maximal volume a cell is allowed to have, irrespective of other 
  refinement and derefinement criteria.

-----

**MeanVolume**

``NODEREFINE_BACKGROUND_GRID``

  Mean volume of cells. In case ``NODEREFINE_BACKGROUND_GRID``, cells with 
  more than 10% or this volume will not be derefined. Used for cosmological 
  zoom simulations.

-----

**CellMaxAngleFactor**

not  ``VORONOI_STATIC_MESH``,
``REGULARIZE_MESH_FACE_ANGLE``

  Cell "roundness" criterion.
  The face angle of an interface of a cell is defined as the square root of the
  area divided by pi, divided by the distance to the face. The maximum face 
  angle is the maximum of this values for all faces in a given cell. This value
  is a measure for the "roundness" of a cell.
  If this value exceeds 1.5 times ``CellMaxAngleFactor``, the cell is not 
  allowed to be refined, i.e. highly elongated cells are not allowed to be 
  refined. In addition to this, if the cell exceeds 0.75 times 
  ``CellMaxAngleFactor``, the movement of the mesh-generating point will 
  deviate from the pure Lagrangian motion to make the cell rounder.

-----

**CellShapingFactor**

not  ``VORONOI_STATIC_MESH``,
not  ``REGULARIZE_MESH_FACE_ANGLE``

  Alternative "roundness" criterion. This criterion uses the distance between 
  center of mass and mesh-generating point as a measure for roundness. If this 
  distnce exceeds twice the cell radius times ``CellShapingFactor``, the cell 
  will not be refined.
  If this distance exceeds 0.75 times the cell radius times 
  ``CellShapingFactor`` the movement of the mesh-generating point will deviate 
  from pure Lagrangian motion to make the cell rounder.

-----

**CellShapingSpeed**

not  ``VORONOI_STATIC_MESH``

  Determines the speed of the regularization of the mesh.
  ``CellShapingSpeed`` (times a characteristic speed) is the speed by which 
  the mesh is regularized (i.e. speed by which the motion of a mesh-generating
  point can deviate from Lagrangian motion).
  Higher values will lead to round cells in fewer timesteps, but will also 
  introduce more numerical noise.

-----

Refinement and derefinement 
===========================

**ReferenceGasPartMass**

``REFINEMENT``

  For comoving runs, it can either be given a non-zero value, in which case 
  this value * ``MassFactor`` is used as the target mass, or it can be given 0,
  in which case the code calculates the mean cell mass. For non-comoving runs, 
  it must be given non-zero, otherwise the run will exit.
  If ``REFINEMENT_HIGH_RES_GAS`` is enabled, then: if 
  ReferenceGasPartMass==0 in the parameter file, then all gas present in the 
  ICs will be allowed to be (de-)refined (and the code calculates the reference
  mass as the mean mass of those cells for which (de-)refinement is allowed), 
  and if that is not desired, then ReferenceGasPartMass should be set to the 
  correct value, in which case only gas with initial 
  mass<1.2*ReferenceGasPartMass will be allowed to be (de-)refined.
  In case of ``GENERATE_GAS_IN_ICS``, only the gas cells split off from 
  particle type 1 (usually the high-res dark matter particles) are flagged for 
  (de-)refinement, i.e. only these gas cells will be considered for the 
  ``ReferenceGasPartMass`` calculation (only in case ``ReferenceGasPartMass=0``
  in parameter file).

-----

**TargetGasMassFactor**

``REFINEMENT``

  The target gas cell mass, where (de-)refinement is triggered if a given cell 
  deviates by more than a factor of 2.0 above or below this value. 
  Multiplicative factor with respect to the mean cell mass. 

-----

**RefinementCriterion**

``REFINEMENT``

  Selects the criterion for refinement; "0" no refinement, "1" target mass 
  refinement, "2" Jeans stability refinement cirtrion.

-----

**DerefinementCriterion**

``REFINEMENT``

  Selects the criterion for derefinement; "0" no derefinement, "1" target mass 
  derefinement, "2" Jeans stability derefinement cirtrion.

-----

Hydrodynamics 
=============

**LimitUBelowThisDensity**

  Density threshold for a specific thermal energy lower limit for low density 
  gas.

-----

**LimitUBelowCertainDensityToThisValue**

  Minimum specific thermal energy for low density gas.

-----

**MinGasTemp**

  A minimum temperature floor imposed by the code. This may be set to zero, 
  but it may be desirable to prevent the gas from becoming too cold, e.g. for 
  resolution reasons or because of lower limits in the implemented cooling 
  function. (This value is converted by the code to a minimum thermal 
  energy per unit mass assuming the mean molecular weight of neutral gas).

-----

**MinEgySpec**

  Minimum specific energy allowed in a gas cell. If specific energy is smaller
  than the value specified here, AREPO will add additional thermal energy in 
  this cell such that it reaches a specific thermal energy of ``MinEgySpec``.
  This is mainly as a protection against negative specific energies emerging 
  from numerical round-off errors in kinetically or magnetically dominated 
  cells (keep in mind that the thermal energy is recomputed from the total 
  energy in AREPO.
  In case this parameter is nonzero, it overrides ``MinGasTemp``.
  If this is zero, internally, ``MinEgySpec`` will be calculated via the value
  of ``MinGasTemp``. In case both ``MinEgySpec``and ``MinGasTemp`` are nonzero,
  ``MinGasTemp`` will only set a lower limit to the cooling.

-----

**IsoSoundSpeed**

``ISOTHERM_EQS``

Sound speed of gas in runs with isothermal hydrodynamics.

-----


Gravitational softening 
=======================

**GasSoftFactor**

  The gravitational softening length of a gas cell is this value times the 
  cellsize, which is calculated as the radius of the volume-equilvalent-sphere.

-----

**SofteningComovingTypeX**

  A Plummer equivalent gravitational softening length, to be referenced by 
  one or more specific particle types. For cosmological simulations in 
  comoving coordinates, this is interpreted as a comoving softening length in
  code length units.

-----

**SofteningMaxPhysTypeX**

  When comoving integration is used, this parameter gives the maximum physical 
  gravitational softening length corresponding to ``SofteningComovingTypeX`` 
  (referenced by one or more specific 
  particle types depending on the entries of ``SofteningTypeOfPartTypeN``). 
  Depening on the relative settings of the *Comoving* and *MaxPhys* softenings,
  the code will hence switch from a softening constant in comoving units to 
  one constant in physical units. For example, if the *MaxPhys* value is 
  exactly half the *Comoving* value, then particles using this softening type 
  will have comoving softening until`z=1` and fixed physical softenings 
  after that point in time. Code length units.

-----

**SofteningTypeOfPartTypeX**

  For each particle type in the simulation which is involved gravitational 
  calculations, it must be assigned to a "softening type", a 0-based integer 
  index corresponding to one of the above 
  ``SofteningComovingTypeX``/``SofteningMaxPhysTypeX`` entry pairs.
  
-----

**MinimumComovingHydroSoftening**

``ADAPTIVE_HYDRO_SOFTENING``

  If this treatment for gas softenings is based used, a discrete spectrum of 
  possible softening lengths for gas cells is created at 
  startup. It contains ``NSOFTTYPES_HYDRO`` entries, controlled by 
  this 'minimum' parameter and the following 'spacing' parameter (as a 
  multiplicative factor). Code length units.

-----

**AdaptiveHydroSofteningSpacing**

``ADAPTIVE_HYDRO_SOFTENING``

  The logarithmic spacing for the adaptive gas softenings table, as described 
  above. Must be larger than one.

-----

Subfind parameters 
==================

**DesLinkNgb**

``SUBFIND``

  The (integer) minimum number of particles/cells, of all types, for Subfind 
  groups. If a Subfind group is identified with fewer than this number of 
  total particles/cells, it is discarded. Note that this means many small 
  friends-of-friends groups (with a nomimal minimum number of 32 member 
  particles) may frequently have no sufficiently large Subfind groups, and 
  so will have ``GroupFirstSub==-1`` indicating that that FoF has no central 
  subhalo in addition to no satellite subhalos.

-----

**ErrTolThetaSubfind**

``SUBFIND``

  This has the same meaning as the ``ErrTolTheta`` parameter, i.e. the tree 
  opening angle used to control the accuracy of the gravity calculation, for 
  uses within the Subfind algorithm.

-----

Cooling and star formation  
==========================

**CoolingOn**

  If set to "1", gas looses energy through a (optically-thin) radiative 
  cooling model at each timestep. Can only be set to zero if ``COOLING`` is 
  not set.

-----

**StarformationOn**

  If set to "1", gas can (stochastically) convert into collisionless star 
  particles based on a star formation model.
  Can only be set to zero if ``USE_SFR`` is not set.

-----

**TreecoolFile**

``COOLING``

  File path to cooling file. Possible files are available under ``./data/``.

-----

**CritOverDensity**

``USE_SFR``

  The critical (over-)density above which star formation may take place, where 
  the threshold density is then  
  
  rho_th = ``CritOverDensity`` 3 ``Omega_b`` H^2 / (8 pi G) 
  
  (redshift independent). Used in place of a critical physical density for 
  comoving integrations.

-----

**TemperatureThresh**

``USE_SFR``

  Star formation is prevented for cells which are hotter than the eEOS and 
  hotter than the TemperatureThresh parameter (in Kelvin). If this parameter is
  very large (e.g. 1e20), then nothing is changed compared to the base model. 
  If this parameter is small (e.g. 0, 1e4, or 1e5) then star-formation will be 
  prevented in hot halo gas.

-----

**CritPhysDensity**

``USE_SFR``

  The critical physical density above which star formation may take place 
  (in `cm^-3`). Used instead of ``CritOverDensity`` for non-comoving runs.

-----

**FactorSN**

``USE_SFR``

  The variable giving the mass fraction of massive stars 
  (``> 8 Msun``) formed for each initial population of stars. This is 
  thus determined by the stellar IMF.

-----

**FactorEVP**

``USE_SFR``

  The variable ``A``, giving the efficiency of the cloud evaporation 
  process.

-----

**TempSupernova**

``USE_SFR``

  The "supernova temperature" ``T_SN`` of the hot intercloud medium in Kelvin.

-----

**TempClouds**

``USE_SFR``

  The "cold cloud temperature" ``T_c``, in Kelvin.

-----

**MaxSfrTimescale**

``USE_SFR``

  This is the star-formation timescale ``t_0`` at the threshold density, such that 
  the local star-formation timescale is then calculated as 
  
  t_star(rho) = t_0 (rho / rho_th)^-0.5

-----

Sphericaly symmetric simulations 
================================

**CoreRadius**

``ONEDIMS_SPHERICAL``

  Inner radius (position of boundary conditions) of 1d spherical simulation.

-----

**CoreMass**

``ONEDIMS_SPHERICAL``

  Mass enclosed within the inner boundary radius in 1d spherical simulations.
  Required for gravity calculation.