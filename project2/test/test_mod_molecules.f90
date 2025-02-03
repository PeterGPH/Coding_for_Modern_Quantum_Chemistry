program test_mod_molecules
    use mod_vibrations
    use iso_fortran_env, only: dp => real64
    implicit none

    type(vibrations) :: mol
    integer :: i, j, k, l, io_status
    character(len=100) :: filename1, filename2
    real(dp), parameter :: tolerance = 1.0e-4_dp
    ! Write test file
    filename1 = "molecule_test.txt"
    filename2 = 'molecule_hessian_test.txt'
    call write_test_file(filename1, filename2)

    ! Read geometry from file
    call mol%read_geometry(filename1)
    call mol%print_geometry()
    call mol%read_cartesian_hessian(filename2)
    call mol%print_cartesian_hessian()
    print *, "--------------------------------------------------"
    call mol%mass_weighted_hessian()
    call mol%print_cartesian_hessian()
    call mol%diagolization_ordered_mass_weighted_hessian()
    call mol%harmoinc_vib_freqs()
    ! Print the geometry
    ! print *, "Geometry read from file:"
contains

    subroutine write_test_file(filename1, filename2)
        character(len=*), intent(in) :: filename1, filename2
        integer :: unit

        open(newunit=unit, file=filename1, status='replace', action='write')
        write(unit, '(I1)') 3
        write(unit, '(I1, 3F20.12)') 8,  0.000000000000_dp,  0.000000000000_dp, -0.134503695264_dp
        write(unit, '(I1, 3F20.12)') 1,  0.000000000000_dp, -1.684916670000_dp,  1.067335684736_dp
        write(unit, '(I1, 3F20.12)') 1,  0.000000000000_dp,  1.684916670000_dp,  1.067335684736_dp
        close(unit)

        open(newunit=unit, file=filename2, status='replace', action='write')
        write(unit, '(I1)') 3
        write(unit, '(3F20.10)') 0.0927643390,  0.0000000000,  0.0000000000
        write(unit, '(3F20.10)') -0.0463821695,  0.0000000000,  0.0000000000
        write(unit, '(3F20.10)') -0.0463821695,  0.0000000000,  0.0000000000
        write(unit, '(3F20.10)') 0.0000000000,  0.3171327134,  0.0000000000
        write(unit, '(3F20.10)') 0.0000000000, -0.1585663567,  0.0800202030
        write(unit, '(3F20.10)') 0.0000000000, -0.1585663567, -0.0800202030
        write(unit, '(3F20.10)') 0.0000000000,  0.0000000000,  0.2800907293
        write(unit, '(3F20.10)') 0.0000000000,  0.0347765865, -0.1400453646
        write(unit, '(3F20.10)') 0.0000000000, -0.0347765865, -0.1400453646
        write(unit, '(3F20.10)') -0.0463821695,  0.0000000000,  0.0000000000
        write(unit, '(3F20.10)') 0.0514668232,  0.0000000000,  0.0000000000
        write(unit, '(3F20.10)') -0.0050846537,  0.0000000000,  0.0000000000
        write(unit, '(3F20.10)') 0.0000000000, -0.1585663567,  0.0347765865
        write(unit, '(3F20.10)') 0.0000000000,  0.1730075524, -0.0573983947
        write(unit, '(3F20.10)') 0.0000000000, -0.0144411957,  0.0226218083
        write(unit, '(3F20.10)') 0.0000000000,  0.0800202030, -0.1400453646
        write(unit, '(3F20.10)') 0.0000000000, -0.0573983947,  0.1268373488
        write(unit, '(3F20.10)') 0.0000000000, -0.0226218083,  0.0132080159
        write(unit, '(3F20.10)') -0.0463821695,  0.0000000000,  0.0000000000
        write(unit, '(3F20.10)') -0.0050846537,  0.0000000000,  0.0000000000
        write(unit, '(3F20.10)') 0.0514668232,  0.0000000000,  0.0000000000
        write(unit, '(3F20.10)') 0.0000000000, -0.1585663567, -0.0347765865
        write(unit, '(3F20.10)') 0.0000000000, -0.0144411957, -0.0226218083
        write(unit, '(3F20.10)') 0.0000000000,  0.1730075524,  0.0573983947
        write(unit, '(3F20.10)') 0.0000000000, -0.0800202030, -0.1400453646
        write(unit, '(3F20.10)') 0.0000000000,  0.0226218083,  0.0132080159
        write(unit, '(3F20.10)') 0.0000000000,  0.0573983947,  0.1268373488
        close(unit)
    end subroutine write_test_file
end program test_mod_molecules