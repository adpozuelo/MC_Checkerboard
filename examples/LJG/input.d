! ==============================================================================
! MONTE CARLO SIMULATION CONTROL PARAMETERS
! ==============================================================================
&Control_Params
  istep_ini   = 0,          ! Initial step number (0 for new simulation, >0 for restart)
  istep_fin   = 1000,       ! Total number of MC sweeps (steps) to execute
  Neq         = 0,          ! Number of equilibration sweeps before accumulating averages
  Nmove       = 20,         ! Trial moves per particle attempted in each MC sweep (10 trans + 10 rot)
  Nsave       = 100,        ! Frequency (in sweeps) for writing progress stdout & block averages
  Ndump       = 100,        ! Frequency (in sweeps) for writing trajectory frames
  Nrestart    = 100,        ! Frequency (in sweeps) for saving restart checkpoints
  imovie      = .true.,     ! Enable trajectory output (.true. / .false.)
  traj_format = 'both',     ! Trajectory output format: 'plain' (text), 'netcdf' (binary), or 'both'
  Npart_types = 1,          ! Number of particle species in system
  lammps      = .false.     ! Master switch: Enable hybrid LAMMPS HMC moves (disabled for LJG)
/

! ==============================================================================
! MONTE CARLO SAMPLING & THERMODYNAMIC PARAMETERS
! ==============================================================================
&MC_Params
  hmax         = 0.0125,    ! Maximum particle translation step size (in units of sigma)
  omax         = 0.5,       ! Maximum particle rotation step size (in radians)
  vmax         = 0.1, 0.1, 0.1, ! Maximum box edge length displacement (Lx, Ly, Lz) for NpT moves
  displ_update = .true.,    ! Enable automatic adjustment of hmax/omax/vmax (~40% acceptance)
  temp0        = 0.154,     ! Initial reduced temperature T0* = kB*T/epsilon
  temp1        = 0.154,     ! Final reduced temperature T1* (set T1* != T0* for linear annealing)
  npt          = .false.,   ! Isobaric-isothermal ensemble (.true. for NpT, .false. for NVT)
  Nvolf        = 1,         ! Frequency (in sweeps) of NpT volume change attempts
  iscale       = 0,         ! Isotropic box scaling flag (0 = independent Lx,Ly,Lz, 1 = isotropic)
  pres         = 0.1,       ! Reduced pressure P* = P*sigma^3/epsilon (for NpT ensemble)
  seed         = 345        ! Random number generator seed
/

! ==============================================================================
! POTENTIAL MODEL SELECTION & PATCHY INTERACTION PARAMETERS
! ==============================================================================
&Potential_Params
  model         = 'LJG',    ! Potential model: 'LJG' (Lennard-Jones-Gauss with patchy modulations)
  sigma_jon_aux = -0.3,     ! Angular patch width scaling factor sigma_ang
  sigma_tor_jon = 0.6,      ! Torsional patch alignment width parameter
  rangepp       = 2.5,      ! Radial potential cutoff distance in units of sigma (Rcut = rangepp*sigma)
  Bool_tor      = 0         ! Enable torsional patch alignment (1 = active, 0 = disabled)
/

! ==============================================================================
! TOPOLOGY AND PATCH GEOMETRY DATA
! ==============================================================================
--- TOPOLOGY AND PATCH GEOMETRY DATA ---
          7                 ! Number of patches per particle for species 1 (7-patch particle)
! Line format for each patch: px  py  pz  vop_x  vop_y  vop_z  bool_tor  ref_phi_1  ref_phi_2  ref_phi_3
 -0.80901699  0.30901699  0.50000000  0.85065807  0.52573111  0.00000000     0    0.000000    3.141592   -3.141592
 -0.30901699 -0.50000000  0.80901699  0.85065807  0.52573111  0.00000000     0    0.000000    3.141592   -3.141592
 -0.80901706  0.30901733 -0.50000000  0.85065807  0.52573111  0.00000000     0    0.000000    3.141592   -3.141592
  0.00000000 -1.00000000  0.00000000  0.85065807  0.52573111  0.00000000     0    0.000000    3.141592   -3.141592
 -0.30901699 -0.50000000 -0.80901699  0.85065807  0.52573111  0.00000000     0    0.000000    3.141592   -3.141592
  0.50000000  0.80901699  0.30901700  0.85065807  0.52573111  0.00000000     0    0.000000    3.141592   -3.141592
  1.00000000  0.00000000  0.00000000  0.85065807  0.52573111  0.00000000     0    0.000000    3.141592   -3.141592

! Patch-patch pairwise well-depth interaction matrix Vpot(alpha, beta) (7 x 7)
-- vpot matrix --
 1.00 1.00 1.00 1.00 1.00 1.10 1.10 
 1.00 1.00 1.00 1.00 1.00 1.10 1.10 
 1.00 1.00 1.00 1.00 1.00 1.10 1.10 
 1.00 1.00 1.00 1.00 1.00 1.10 1.10 
 1.00 1.00 1.00 1.00 1.00 1.10 1.10 
 1.10 1.10 1.10 1.10 1.10 1.20 1.20 
 1.10 1.10 1.10 1.10 1.10 1.20 1.20

! Patch angular interaction width parameters sigmas_ang for each patch (1 x 7)
-- sigmas_ang --
 0.300 0.300 0.300 0.300 0.300 0.300 0.300

! Colloid core LJ diameter parameters sigma_LJ (1 x 1)
-- sigma_LJ --
 1.000 

 
