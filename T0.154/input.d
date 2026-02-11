0.1.dat             
 Step_ini Step_fin Neq  Nmove  Nsave Nsave2 Nrestart
        0  500000   0    20    100  5000    5000
 Hmax Max_or vmax   adjust_displ (Lx, Ly,Lz)   temp_initial temp_final npt iscale pres seed
  0.0125  0.5000   0.1 0.1 0.1    T        0.15400     0.15400 T  0 0.1         345
 parametros para energias libres (not implemented)
 XXXXXXXXXXXXXX
 imovie
 T
 Trajectory format (plain, netcdf, both)
 netcdf
 Model  (K-F or LJG)
 LJG 
 Parametros modelo LJG       
 cosmax ctorsiona_max deltaa  bool_torsion
  -0.3000  0.6000  2.5000    0      
 parametros para AVB (not implemented)
 XXXXXXXXXXXXXXXXXXXXX
  Number  of particle types     
           1                
 Number of patches in each particle
           7         
 -0.80901699  0.30901699  0.50000000  0.85065807  0.52573111  0.00000000     0    0.000000    3.141592   -3.141592
 -0.30901699 -0.50000000  0.80901699  0.85065807  0.52573111  0.00000000     0    0.000000    3.141592   -3.141592
 -0.80901706  0.30901733 -0.50000000  0.85065807  0.52573111  0.00000000     0    0.000000    3.141592   -3.141592
  0.00000000 -1.00000000  0.00000000  0.85065807  0.52573111  0.00000000     0    0.000000    3.141592   -3.141592
 -0.30901699 -0.50000000 -0.80901699  0.85065807  0.52573111  0.00000000     0    0.000000    3.141592   -3.141592
  0.50000000  0.80901699  0.30901700  0.85065807  0.52573111  0.00000000     0    0.000000    3.141592   -3.141592
  1.00000000  0.00000000  0.00000000  0.85065807  0.52573111  0.00000000     0    0.000000    3.141592   -3.141592
 -vpot matrix --
 1.00 1.00 1.00 1.00 1.00 1.10 1.10 
 1.00 1.00 1.00 1.00 1.00 1.10 1.10 
 1.00 1.00 1.00 1.00 1.00 1.10 1.10 
 1.00 1.00 1.00 1.00 1.00 1.10 1.10 
 1.00 1.00 1.00 1.00 1.00 1.10 1.10 
 1.10 1.10 1.10 1.10 1.10 1.20 1.20 
 1.10 1.10 1.10 1.10 1.10 1.20 1.20 
sigmas_ang
 0.300 0.300 0.300 0.300 0.300 0.300 0.300 0.300
sigma_LJ
 1.000 


#  
# iscale=0 cubico
# isclae =1 tetragonal (Lx=ly != Lz)
# icale 2  Lx != Ly = !Lz (angulos simpre rectos) 
