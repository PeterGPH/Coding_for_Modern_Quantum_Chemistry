program test_molecules
  use mod_molecules
  implicit none

  type(molecules) :: mol
  character(len=100) :: filename

  ! Prompt user for input file name
  write(*, '(A)', advance='no') "Enter the name of the geometry file: "
  read(*, '(A)') filename

  ! Read the geometry from the file
  call mol%read_geometry(filename)

  ! Print the geometry to verify the information
  ! call mol%print_geometry()

end program test_molecules