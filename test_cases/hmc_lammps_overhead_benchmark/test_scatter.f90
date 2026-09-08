program test_scatter
    use, intrinsic :: iso_c_binding
    use liblammps
    implicit none

    type(lammps) :: lmp
    character(len=32) :: args(14)
    character(len=256) :: cmd
    real(c_double), allocatable :: pos(:)
    real(c_double) :: t1, t2, dt_scatter, dt_run, dt_gather
    integer :: n_atoms, i

    args(1) = "liblammps"
    args(2) = "-screen"; args(3) = "none"
    args(4) = "-log";    args(5) = "none"
    args(6) = "-pk";    args(7) = "gpu"; args(8) = "1"
    args(9) = "gpuID";  args(10) = "0"
    args(11) = "neigh"; args(12) = "yes"
    args(13) = "-sf";   args(14) = "gpu"

    lmp = lammps(args)

    call lmp%command("units lj")
    call lmp%command("atom_style full")
    call lmp%command("boundary p p p")
    call lmp%command("read_data lammps_hmc.data")
    call lmp%command("include potential_ssp_table.lmp")
    call lmp%command("pair_modify shift yes")
    call lmp%command("fix hmc_rigid all rigid/nve molecule")
    call lmp%command("compute hmc_ke all ke")
    call lmp%command("compute hmc_pe all pe")
    call lmp%command("velocity all create 0.15 12345 loop geom")
    call lmp%command("run 0")

    call lmp%gather_atoms("x", 3, pos)
    n_atoms = size(pos) / 3
    print *, "Successfully gathered atoms from LAMMPS: N_atoms = ", n_atoms

    ! Modify positions slightly (shift x by 0.01)
    pos(1:3*n_atoms:3) = pos(1:3*n_atoms:3) + 0.01_c_double

    ! Time scatter_atoms
    call cpu_time(t1)
    call lmp%scatter_atoms("x", pos)
    call cpu_time(t2)
    dt_scatter = t2 - t1
    print *, "scatter_atoms time = ", dt_scatter * 1000.0, " ms"

    ! Time run 50
    call cpu_time(t1)
    call lmp%command("run 50")
    call cpu_time(t2)
    dt_run = t2 - t1
    print *, "run 50 time = ", dt_run, " s"

    ! Time gather_atoms
    call cpu_time(t1)
    call lmp%gather_atoms("x", 3, pos)
    call cpu_time(t2)
    dt_gather = t2 - t1
    print *, "gather_atoms time = ", dt_gather * 1000.0, " ms"

    call lmp%close()
    print *, "TEST PASSED SUCCESSFULLY!"
end program test_scatter
