module mod_vibrations 
  use iso_fortran_env, only: dp => real64
  use mod_molecules
  use stdlib_linalg, only: eigh
  use stdlib_sorting, only: ord_sort

  implicit none

  real(dp), parameter :: an2masses(0:19) = [0.0_dp, 1.00784_dp, 4.00260_dp, 6.941_dp, 9.01218_dp, &
                                      10.81_dp, 12.01_dp, 14.01_dp, 16.00_dp, 19.00_dp, 20.18_dp, &
                                      22.99_dp, 24.31_dp, 26.98_dp, 28.09_dp, 30.97_dp, 32.06_dp, &
                                      35.45_dp, 39.95_dp, 39.10_dp]
  real(dp), parameter ::  pi = 3.14159265358979323846_dp

  private

  public :: vibrations

  type, extends(molecules) :: vibrations
    real(dp), dimension(:,:), allocatable :: cartesian_hessian
    real(dp), dimension(:), allocatable :: hessian_eigvals
  contains
    procedure :: read_cartesian_hessian
    procedure :: mass_weighted_hessian
    procedure :: print_cartesian_hessian
    procedure :: diagolization_ordered_mass_weighted_hessian
    procedure :: harmoinc_vib_freqs
  end type vibrations

contains
subroutine read_cartesian_hessian(self, filename)
  class(vibrations), intent(inout) :: self
  character(len=*), intent(in) :: filename
  integer :: io, i, j, fileunit, natom  

  open(newunit=fileunit, file=filename, status='old', action='read', iostat=io)
  if (io /= 0) then
      print *, "Error opening file: ", filename
      return
  end if
  
  read(fileunit, *, iostat=io) natom
  if (io /= 0) then
      print *, "Error reading number of atoms"
      close(fileunit)
      return
  end if

  if (natom /= self%natom) then
      print *, "Number of atoms in geometry and hessian do not match"
      close(fileunit)
      return
  end if

  if (allocated(self%cartesian_hessian)) deallocate(self%cartesian_hessian)
  allocate(self%cartesian_hessian(1:3*natom, 1:3*natom))
  self%cartesian_hessian = 0.0_dp
  do i = 1, self%natom*3
    do j = 1, self%natom
      read(fileunit, *, iostat=io) self%cartesian_hessian(i, 3*j-2), self%cartesian_hessian(i,3*j-1), self%cartesian_hessian(i,3*j)
      if (io /= 0) then
          print *, "Error reading atom data for atom ", i
          close(fileunit)
          return
      end if
    end do
  end do

  close(fileunit)
  print *, "Geometry read successfully from ", filename
end subroutine read_cartesian_hessian

subroutine mass_weighted_hessian(self)
  class(vibrations), intent(inout) :: self
  integer :: i, j, real_i
  real(dp) :: massi, massj

  do i = 1, 3*self%natom
    massi = an2masses(self%zvals(int((i-1)/3)+1))
    do j = 1, self%natom
      massj = an2masses(self%zvals(j))
      self%cartesian_hessian(i, 3*j-2:3*j) = self%cartesian_hessian(i, 3*j-2:3*j) / (sqrt(massi)*sqrt(massj))
    end do
  end do
end subroutine mass_weighted_hessian

subroutine diagolization_ordered_mass_weighted_hessian(self)
  class(vibrations), intent(inout) :: self
  integer :: i, j, k, l
  real(dp), dimension(:), allocatable :: lambda, work
  real(dp), dimension(:,:), allocatable :: vectors

  if (.not. allocated(self%cartesian_hessian)) then
      print *, "Hessian not allocated"
      return
  end if
  allocate(lambda(1:3*self%natom), work(1:3*self%natom), vectors(1:3*self%natom, 1:3*self%natom))
  call eigh(self%cartesian_hessian, lambda, vectors)
  call ord_sort(lambda, work)
  
  do i = 1, size(lambda)
    print *, lambda(i)
  end do
  allocate(self%hessian_eigvals(1:3*self%natom))
  self%hessian_eigvals = lambda
end subroutine diagolization_ordered_mass_weighted_hessian

subroutine harmoinc_vib_freqs(self)
  !How to determine the constant?
  class(vibrations), intent(in) :: self
  integer :: i
  real(dp) :: freq

  do i = 1, 3*self%natom
    freq = 5140.484157*sqrt(abs(self%hessian_eigvals(i)))
    print *, freq
  end do
end subroutine harmoinc_vib_freqs

subroutine print_cartesian_hessian(self)
  class(vibrations), intent(in) :: self
  integer :: i

  do i = 1, 3*self%natom
    print *, self%cartesian_hessian(i, :)
  end do
end subroutine print_cartesian_hessian
end module mod_vibrations