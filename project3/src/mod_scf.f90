module mod_scf
  use iso_fortran_env, only: dp => real64
  use mod_molecules
  use stdlib_linalg, only: eigh, operator(.inv.)
  implicit none
  integer, parameter :: npgs = 7 ! 7 number of primitive gaussians for STO-3G
  integer, parameter :: nelec = 10 ! 10 for H2O
  real(dp), parameter :: conv_threshold = 1.0e-8_dp ! Convergence threshold for SCF
  integer, parameter :: max_iter = 100 ! Maximum number of SCF iterations
  private

  public :: scf

  type, extends(molecules) :: scf
    real(dp)  ::  enuc, E_elec, E_total
    real(dp), allocatable :: S_C(:,:), T_C(:,:), V_C(:,:), H_C(:,:), TEI(:,:,:,:), S_C_inv(:, :), F_C(:,:), D_C(:,:)
  contains
    procedure :: read_enuc
    procedure :: read_overlap
    procedure :: read_kinetic
    procedure :: read_nuclear_attraction
    procedure :: cons_core_hamiltonian
    procedure :: read_two_electron
    procedure :: cons_orthogonalization_matrix
    procedure :: cons_initial_guess
    procedure :: cons_initial_scf_energy
    procedure :: update_fock_matrix
    procedure :: update_density_matrix
    procedure :: run_scf
  end type scf
contains

subroutine read_enuc(self, filename)
  class(scf), intent(inout) :: self
  character(len=*), intent(in) :: filename
  integer :: unit, io

  open(newunit=unit, file=filename, status='old', action='read')
  read(unit, '(F20.12)', iostat=io) self%enuc
  if (io /= 0) then
    print *, "Error reading number of atoms"
    close(unit)
    return
  end if
  close(unit)
end subroutine read_enuc

subroutine read_overlap(self, filename)
  class(scf), intent(inout) :: self
  character(len=*), intent(in) :: filename
  real(dp) :: value
  integer :: unit, io, i, j
  if (allocated(self%S_C)) deallocate(self%S_C)
  allocate(self%S_C(npgs, npgs))
  self%S_C(npgs, npgs) = 0.0_dp
  open(newunit=unit, file=filename, status='old', action='read')
  do 
    read(unit, *, iostat=io) i, j, value
    if (io /= 0) exit
    self%S_C(i, j) = value
    if (self%S_C(j, i) == 0.0_dp) then
      self%S_C(j, i) = value
    end if
  end do
  close(unit)
end subroutine read_overlap

subroutine read_kinetic(self, filename)
  class(scf), intent(inout) :: self
  character(len=*), intent(in) :: filename
  real(dp) :: value
  integer :: unit, io, i, j
  if (allocated(self%T_C)) deallocate(self%T_C)
  allocate(self%T_C(npgs, npgs))
  self%T_C(npgs, npgs) = 0.0_dp
  open(newunit=unit, file=filename, status='old', action='read')
  do 
    read(unit, *, iostat=io) i, j, value
    if (io /= 0) exit
    self%T_C(i, j) = value
  end do
  close(unit)
end subroutine read_kinetic

subroutine read_nuclear_attraction(self, filename)
  class(scf), intent(inout) :: self
  character(len=*), intent(in) :: filename
  real(dp) :: value
  integer :: unit, io, i, j
  if (allocated(self%V_C)) deallocate(self%V_C)
  allocate(self%V_C(npgs, npgs))
  self%V_C(npgs, npgs) = 0.0_dp
  open(newunit=unit, file=filename, status='old', action='read')
  do 
    read(unit, *, iostat=io) i, j, value
    if (io /= 0) exit
    self%V_C(i, j) = value
  end do
  close(unit)
end subroutine read_nuclear_attraction

subroutine read_two_electron(self, filename)
  class(scf), intent(inout) :: self
  character(len=*), intent(in) :: filename
  real(dp) :: value
  integer :: unit, io, mu, nu, lambda, sigma
  if (allocated(self%TEI)) deallocate(self%TEI)
  allocate(self%TEI(npgs, npgs, npgs, npgs))
  self%TEI(npgs, npgs, npgs, npgs) = 0.0_dp
  open(newunit=unit, file=filename, status='old', action='read')
  do 
    read(unit, *, iostat=io) mu, nu, lambda, sigma, value
    if (io /= 0) exit
    self%TEI(mu, nu, lambda, sigma) = value
    self%TEI(nu, mu, lambda, sigma) = value
    self%TEI(mu, nu, sigma, lambda) = value
    self%TEI(nu, mu, sigma, lambda) = value
    self%TEI(lambda, sigma, mu, nu) = value
    self%TEI(sigma, lambda, mu, nu) = value
    self%TEI(lambda, sigma, nu, mu) = value
    self%TEI(sigma, lambda, nu, mu) = value
  end do
  close(unit)
end subroutine read_two_electron

subroutine cons_core_hamiltonian(self)
  class(scf), intent(inout) :: self
  integer :: i, j
  if (allocated(self%H_C)) deallocate(self%H_C)
  allocate(self%H_C(npgs, npgs))
  self%H_C = self%T_C + self%V_C
  do i=1, npgs
    do j=1, npgs
      if (self%H_C(j, i) == 0.0_dp) then 
        self%H_C(j, i) = self%H_C(i, j)
      end if 
    end do
  end do
end subroutine cons_core_hamiltonian

subroutine cons_orthogonalization_matrix(self)
  class(scf), intent(inout) :: self
  real(dp), dimension(:), allocatable :: lambda, work, inverse_lambda
  real(dp), dimension(:,:), allocatable :: vectors 
  integer :: i
  if (.not. allocated(self%S_C)) then
    print *, "Overlap matrix not read"
    return
  end if
  allocate(lambda(npgs), work(npgs), vectors(npgs, npgs), inverse_lambda(npgs))
  if (allocated(self%S_C_inv)) deallocate(self%S_C_inv)
  allocate(self%S_C_inv(npgs, npgs))
  self%S_C_inv = 0.0_dp
  call eigh(self%S_C, lambda, vectors)
  inverse_lambda = sqrt(1.0_dp/lambda)
  do i = 1, size(lambda)
    self%S_C_inv(:, i) = vectors(:, i) * inverse_lambda(i)
  end do
  self%S_C_inv = matmul(self%S_C_inv, transpose(vectors))
  end subroutine cons_orthogonalization_matrix

  subroutine cons_initial_guess(self)
    class(scf), intent(inout) :: self
    real(dp), dimension(:), allocatable :: lambda, work
    real(dp), dimension(:,:), allocatable :: vectors 
    logical :: debug = .true.
    integer :: i, mu, nu
    if (.not. allocated(self%H_C)) then
      print *, "Core Hamiltonian not constructed"
      return
    end if
    if (allocated(self%F_C)) deallocate(self%F_C)
    allocate(self%F_C(npgs, npgs))
    allocate(lambda(npgs), work(npgs), vectors(npgs, npgs))
    self%F_C = matmul(transpose(self%S_C_inv), self%H_C)
    self%F_C = matmul(self%F_C, self%S_C_inv)
    call eigh(self%F_C, lambda, vectors)
    vectors = matmul(self%S_C_inv, vectors)
    if (debug) then
      print *, "-----------------------------------------------------------------"
      print *, "The Initial MO Coefficients is:"
      call dump2(vectors, npgs)
    end if
    if (allocated(self%D_C)) deallocate(self%D_C)
    allocate(self%D_C(npgs, npgs))
    self%D_C = 0.0_dp
    do i = 1, nelec/2
      do mu = 1, size(self%D_C, 1)
        do nu = 1, size(self%D_C, 2)
          self%D_C(mu, nu) = self%D_C(mu, nu) + vectors(mu, i) * vectors(nu, i)
        end do
      end do
    end do 
    if (debug) then
      print *, "-----------------------------------------------------------------"
      print *, "The Initial Density Matrix is:"
      call dump2(self%D_C, npgs)
    end if 

  end subroutine cons_initial_guess

  subroutine cons_initial_scf_energy(self)
    class(scf), intent(inout) :: self
    integer :: mu, nu, n 

    ! Check matrix dimensions
    n = size(self%D_C, 1)
    if (size(self%D_C, 2) /= n .or. size(self%H_C, 1) /= n .or. size(self%H_C, 2) /= n .or. &
        size(self%F_C, 1) /= n .or. size(self%F_C, 2) /= n) then
        print *, "Error: Matrix dimensions do not match"
        return
    end if
    ! print *, "D_C(1,1) = ", self%D_C(1, 1)
    ! print *, "H_C(1,1) = ", self%H_C(1, 1)
    ! print *, "F_C(1,1) = ", self%F_C(1, 1)
    if (.not. allocated(self%H_C)) then
      print *, "Core Hamiltonian is not constructed"
      return
    end if
    if (.not. allocated(self%F_C)) then
      print *, "Fock matrix is not constructed"
      return
    end if
    if (.not. allocated(self%D_C)) then
      print *, "Fock matrix is not constructed"
      return
    end if
    self%E_elec = 0.0_dp
    do mu = 1, npgs
      do nu = 1, npgs
        ! print *, "mu = ", mu, ", nu = ", nu
        ! print *, "D_C(mu, nu) = ", self%D_C(mu, nu)
        ! print *, "H_C(mu, nu) = ", self%H_C(mu, nu)
        ! print *, "F_C(mu, nu) = ", self%F_C(mu, nu)
        ! print *, "Term = ", self%D_C(mu, nu) * (self%H_C(mu, nu) + self%F_C(mu, nu))
        ! print *, "E_elec = ", self%E_elec
        self%E_elec = self%E_elec + (self%D_C(mu, nu)*(self%H_C(mu, nu) + self%F_C(mu, nu)))
      end do
    end do
    self%E_total = self%E_elec + self%enuc
  end subroutine cons_initial_scf_energy

  subroutine update_fock_matrix(self)
    class(scf), intent(inout) :: self
    integer :: mu, nu, lambda, sigma
    if (.not. allocated(self%D_C)) then
      print *, "Density matrix is not constructed"
      return
    end if
    do mu = 1, size(self%F_C, 1)
      do nu = 1, size(self%F_C, 2)
        self%F_C(mu, nu) = self%H_C(mu, nu)
        do lambda = 1, size(self%F_C, 1)
          do sigma = 1, size(self%F_C, 2)
            self%F_C(mu, nu) = self%F_C(mu, nu) + self%D_C(lambda, sigma) * (2.0_dp * self%TEI(mu, nu, lambda, sigma) - self%TEI(mu, lambda, nu, sigma))
          end do
        end do
      end do
    end do
  end subroutine update_fock_matrix

  subroutine update_density_matrix(self)
    class(scf), intent(inout) :: self
    real(dp), dimension(:), allocatable :: lambda, work
    real(dp), dimension(:,:), allocatable :: vectors 
    integer :: i, mu, nu
    if (.not. allocated(self%F_C)) then
      print *, "Fock matrix is not constructed"
      return
    end if
    allocate(lambda(npgs), work(npgs), vectors(npgs, npgs))
    self%F_C = matmul(transpose(self%S_C_inv), self%F_C)
    self%F_C = matmul(self%F_C, self%S_C_inv)
    call eigh(self%F_C, lambda, vectors)
    vectors = matmul(self%S_C_inv, vectors)
    self%D_C = 0.0_dp
    do i = 1, nelec/2
      do mu = 1, size(self%D_C, 1)
        do nu = 1, size(self%D_C, 2)
          self%D_C(mu, nu) = self%D_C(mu, nu) + vectors(mu, i) * vectors(nu, i)
        end do
      end do
    end do 
  end subroutine update_density_matrix

  subroutine run_scf(self)
    class(scf), intent(inout) :: self
    integer :: iter
    real(dp):: E_old, delta_E
    
    ! Initialize SCF
    call self%read_enuc('enuc.dat')
    call self%read_overlap('s.dat')
    call self%read_kinetic('t.dat')
    call self%read_nuclear_attraction('v.dat')
    call self%read_two_electron('eri.dat')
    call self%cons_core_hamiltonian()
    call self%cons_orthogonalization_matrix()
    call self%cons_initial_guess()
    call self%cons_initial_scf_energy()

    E_old = self%E_total
    print *, "Initial SCF Energy = ", E_old
    
    ! SCF Loop
    print *, "-----------------------------------------------------------------"
    print '(A5, A15, A15, A15)', "Iter", "E(elec)", "E(tot)", "Delta(E)"
    print *, "-----------------------------------------------------------------"
    print '(I5, F15.10, F15.10)', 0, self%E_elec, self%E_total
    do iter = 1, max_iter
      call self%update_fock_matrix()
      call self%update_density_matrix()
      call self%cons_initial_scf_energy()

      delta_E = abs(self%E_total - E_old)
      print '(I5, F15.10, F15.10, F15.10)', iter, self%E_elec, self%E_total, delta_E

      if (delta_E < conv_threshold) then
        print *, "SCF converged in ", iter, " iterations."
        exit
      end if

      E_old = self%E_total
    end do

    if (iter >= max_iter) then
      print *, "SCF did not converge in ", max_iter, " iterations."
    end if

    print *, "Final SCF Energy: ", self%E_total
  end subroutine run_scf
  
  subroutine dump2(W, N)
    integer, intent(in) :: N
    real(dp), dimension(N, N), intent(in) :: W
    integer :: i
  
    do i = 1, N
      print *, W(i, :)
    end do
  end subroutine dump2
end module mod_scf
