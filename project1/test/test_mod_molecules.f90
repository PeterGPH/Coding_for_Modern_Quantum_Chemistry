program test_mod_molecules
    use mod_molecules
    use iso_fortran_env, only: dp => real64
    implicit none

    type(molecules) :: mol
    integer :: i, j, k, l, io_status
    character(len=100) :: filename
    real(dp), parameter :: tolerance = 1.0e-4_dp
    logical :: all_passed
    ! Write test file
    filename = "molecule_test.xyz"
    call write_test_file(filename)

    ! Read geometry from file
    call mol%read_geometry(filename)

    ! Print the geometry
    print *, "Geometry read from file:"
    call mol%print_geometry()


contains

    subroutine write_test_file(filename)
        character(len=*), intent(in) :: filename
        integer :: unit

        open(newunit=unit, file=filename, status='replace', action='write')
        write(unit, '(I1)') 7
        write(unit, '(I1, 3F20.12)') 6,  0.000000000000_dp,  0.000000000000_dp,  0.000000000000_dp
        write(unit, '(I1, 3F20.12)') 6,  0.000000000000_dp,  0.000000000000_dp,  2.845112131228_dp
        write(unit, '(I1, 3F20.12)') 8,  1.899115961744_dp,  0.000000000000_dp,  4.139062527233_dp
        write(unit, '(I1, 3F20.12)') 1, -1.894048308506_dp,  0.000000000000_dp,  3.747688672216_dp
        write(unit, '(I1, 3F20.12)') 1,  1.942500819960_dp,  0.000000000000_dp, -0.701145981971_dp
        write(unit, '(I1, 3F20.12)') 1, -1.007295466862_dp, -1.669971842687_dp, -0.705916966833_dp
        write(unit, '(I1, 3F20.12)') 1, -1.007295466862_dp,  1.669971842687_dp, -0.705916966833_dp
        close(unit)
    end subroutine write_test_file
end program test_mod_molecules