module mod_molecules
    use iso_fortran_env, only: dp => real64
    use stdlib_linalg, only: eigh
    use stdlib_sorting, only: ord_sort
    implicit none
    private
    real(dp), parameter :: an2masses(0:19) = [0.0_dp, 1.00784_dp, 4.00260_dp, 6.941_dp, 9.01218_dp, &
                                            10.81_dp, 12.01_dp, 14.01_dp, 16.00_dp, 19.00_dp, 20.18_dp, &
                                            22.99_dp, 24.31_dp, 26.98_dp, 28.09_dp, 30.97_dp, 32.06_dp, &
                                            35.45_dp, 39.95_dp, 39.10_dp]
    real(dp), parameter ::  pi = 3.14159265358979323846_dp

    type, public :: molecules
        integer :: natom
        integer :: charge
        integer, public, allocatable :: zvals(:)
        real(dp), allocatable :: geom(:,:)
        real(dp) :: com(3), inertia_eigvals(3)
        character(len=50) :: point_group

    contains
        procedure :: init_molecule
        procedure :: unit
        procedure :: read_geometry
        procedure :: bond_length
        procedure :: bond_angles
        procedure :: out_of_plane_angle
        procedure :: torsion
        procedure :: center_of_mass
        procedure :: principle_moments_of_inertia
        procedure :: rotational_constants
        procedure :: print_geometry
    end type molecules

contains
    subroutine init_molecule(self, natom, charge) 
        class(molecules), intent(inout) :: self
        integer, intent(in) :: natom, charge

        self%natom = natom
        self%charge = charge
        if (allocated(self%zvals)) deallocate(self%zvals)
        if (allocated(self%geom)) deallocate(self%geom)
        allocate(self%zvals(natom))
        self%zvals = 0
        allocate(self%geom(3, natom))
        self%geom = 0.0_dp
        self%point_group = ''
    end subroutine init_molecule

    subroutine read_geometry(self, filename)
        class(molecules), intent(inout) :: self
        character(len=*), intent(in) :: filename
        integer :: io, i, fileunit

        open(newunit=fileunit, file=filename, status='old', action='read', iostat=io)
        if (io /= 0) then
            print *, "Error opening file: ", filename
            return
        end if
        
        read(fileunit, *, iostat=io) self%natom
        if (io /= 0) then
            print *, "Error reading number of atoms"
            close(fileunit)
            return
        end if

        call self%init_molecule(self%natom, 0)
    
        do i = 1, self%natom
            read(fileunit, *, iostat=io) self%zvals(i), self%geom(:,i)
            if (io /= 0) then
                print *, "Error reading atom data for atom ", i
                close(fileunit)
                return
            end if
        end do
    
        close(fileunit)
        print *, "Geometry read successfully from ", filename
    end subroutine read_geometry

    function bond_length(self, i, j) result(length)
        class(molecules), intent(in) :: self
        integer, intent(in) :: i, j
        real(dp) :: length

        length = sqrt((self%geom(1,i) - self%geom(1,j))**2 + &
                      (self%geom(2,i) - self%geom(2,j))**2 + &
                      (self%geom(3,i) - self%geom(3,j))**2)
    end function bond_length


    function unit(self, i, j) result(unit_vector)
        class(molecules), intent(in) :: self
        integer, intent(in) :: i, j
        real(dp) :: unit_vector(3)

        unit_vector = (self%geom(:,j) - self%geom(:,i)) / self%bond_length(i, j)
    end function 

    function bond_angles(self, i, j, k) result(angle)
        class(molecules), intent(in) :: self
        integer, intent(in) :: i, j, k
        real(dp) :: angle
        real(dp) :: e_ji(3), e_jk(3)

        e_ji = self%unit(i, j)
        e_jk = self%unit(k, j)
        angle = acos(dot_product(e_ji, e_jk)) * 180.0_dp / pi
    end function bond_angles

    function out_of_plane_angle(self, i, j, k, l) result(angle)
        ! No idea why e_ki = self%geom(:,i) - self%geom(:,k) / self%bond_length(i, k) don't work, resulting in nan
        class(molecules), intent(in) :: self
        integer, intent(in) :: i, j, k, l
        real(dp) :: angle
        real(dp) :: e_kjl(3), e_kj(3), e_kl(3), e_ki(3)
        real(dp) :: exx, eyy, ezz

        e_kj = self%unit(j, k)
        e_kl = self%unit(l, k)
        e_ki = self%unit(i, k)
        
        e_kjl(1) = e_kj(2)*e_kl(3) - e_kj(3)*e_kl(2)
        e_kjl(2) = e_kj(3)*e_kl(1) - e_kj(1)*e_kl(3)
        e_kjl(3) = e_kj(1)*e_kl(2) - e_kj(2)*e_kl(1)

        exx = e_kjl(1) * e_ki(1)
        eyy = e_kjl(2) * e_ki(2)
        ezz = e_kjl(3) * e_ki(3)
        angle = (exx + eyy + ezz) / sin(self%bond_angles(j, k, l) * pi / 180.0_dp )

        if (angle > 1.0_dp) then
            angle = asin(1.0_dp) * 180.0_dp / pi
        else if (angle < -1.0_dp) then
            angle = asin(-1.0_dp) * 180.0_dp / pi
        else
            angle = asin(angle) * 180.0_dp / pi
        end if 
        !Not sure why?
        angle = -angle
    end function out_of_plane_angle

    function torsion(self, i, j, k, l) result(angle)
        ! acos(x) in fortran result in [pi, 0] with [-1, 1]
        ! However the acos(x) in C++ result in [0, pi] wiht [-1, 1]
        class(molecules), intent(in) :: self
        integer, intent(in) :: i, j, k, l
        real(dp) :: angle
        real(dp) :: e_ji(3), e_jk(3), e_kl(3)
        real(dp) :: e_jkjl(3), e_kjkl(3)
        real(dp) :: exx, eyy, ezz, cross_x, cross_y, cross_z, norm, sign

        e_ji = self%unit(i, j)
        e_jk = self%unit(k, j)
        e_kl = self%unit(l, k)
        e_jkjl(1) = e_ji(2)*e_jk(3) - e_ji(3)*e_jk(2)
        e_jkjl(2) = e_ji(3)*e_jk(1) - e_ji(1)*e_jk(3)
        e_jkjl(3) = e_ji(1)*e_jk(2) - e_ji(2)*e_jk(1)
        e_kjkl(1) = e_jk(2)*e_kl(3) - e_jk(3)*e_kl(2)
        e_kjkl(2) = e_jk(3)*e_kl(1) - e_jk(1)*e_kl(3)
        e_kjkl(3) = e_jk(1)*e_kl(2) - e_jk(2)*e_kl(1)

        exx = e_jkjl(1) * e_kjkl(1)
        eyy = e_jkjl(2) * e_kjkl(2)
        ezz = e_jkjl(3) * e_kjkl(3)
        angle = (exx + eyy + ezz) / (sin(self%bond_angles(i, j, k)* pi / 180.0_dp ) * sin(self%bond_angles(j, k, l)* pi / 180.0_dp ))

        if (angle > 1.0_dp) then
            angle = 180.0_dp - acos(1.0_dp) * 180.0_dp / pi
        else if (angle < -1.0_dp) then
            angle = 180.0_dp - acos(-1.0_dp) * 180.0_dp / pi
        else
            angle = 180.0_dp - acos(angle) * 180.0_dp / pi
        end if

        cross_x = e_jkjl(2) * e_kjkl(3) - e_jkjl(3) * e_kjkl(2)
        cross_y = e_jkjl(3) * e_kjkl(1) - e_jkjl(1) * e_kjkl(3)
        cross_z = e_jkjl(1) * e_kjkl(2) - e_jkjl(2) * e_kjkl(1)
        norm = sqrt(cross_x**2 + cross_y**2 + cross_z**2)
        cross_x = cross_x / norm
        cross_y = cross_y / norm
        cross_z = cross_z / norm
        sign = 1.0_dp
        if (cross_x*e_ji(1) + cross_y*e_ji(2) + cross_z*e_ji(3) < 0.0_dp) sign = -1.0_dp
        angle = angle * sign
    end function torsion

    subroutine center_of_mass(self)
        class(molecules), intent(inout) :: self
        real(dp) :: mass, total_mass
        integer :: i

        total_mass = 0.0_dp
        self%com = 0.0_dp

        do i = 1, self%natom
            mass = an2masses(self%zvals(i))
            ! write(*,*) mass
            total_mass = total_mass + mass
            self%com = self%com + mass * self%geom(:,i)
        end do
        self%com = self%com / total_mass
    end subroutine center_of_mass

    subroutine principle_moments_of_inertia(self)
        class(molecules), intent(inout) :: self
        real(dp) :: inertia(3,3), lambda(3), vectors(3,3)
        real(dp), allocatable :: relative_geom(:,:)
        real(dp) :: mass
        integer :: i

        inertia = 0.0_dp
        call self%center_of_mass()
        allocate(relative_geom(3, self%natom))
        do i = 1, self%natom
            relative_geom(:,i) = self%geom(:,i) - self%com
        end do
        do i = 1, self%natom
            mass = an2masses(self%zvals(i))
            inertia(1,1) = inertia(1,1) + mass * (relative_geom(2,i)**2 + relative_geom(3,i)**2)
            inertia(2,2) = inertia(2,2) + mass * (relative_geom(1,i)**2 + relative_geom(3,i)**2)
            inertia(3,3) = inertia(3,3) + mass * (relative_geom(1,i)**2 + relative_geom(2,i)**2)
            inertia(1,2) = inertia(1,2) - mass * relative_geom(1,i) * relative_geom(2,i)
            inertia(1,3) = inertia(1,3) - mass * relative_geom(1,i) * relative_geom(3,i)
            inertia(2,3) = inertia(2,3) - mass * relative_geom(2,i) * relative_geom(3,i)
        end do
        inertia(2,1) = inertia(1,2)
        inertia(3,1) = inertia(1,3)
        inertia(3,2) = inertia(2,3)
        print *, "Momemnt of inertia tensor (amu bohr^2)"
        do i=1,3
            print *, inertia(:, i)
        end do

        call eigh(inertia, lambda, vectors)

        ! print *, 'Real matrix'
        ! do i=1,3
        !    print *, 'eigenvalue  ',i,': ',lambda(i)
        !    print *, 'eigenvector ',i,': ',vectors(:,i)
        ! end do

        print *, "Principle moments of inertia (amu bohr^2)"
        print *, lambda(1), lambda(2), lambda(3)
        self%inertia_eigvals = lambda
        deallocate(relative_geom)
        if (self%natom .EQ. 2) then
            print *, "Molecule is diatomic."
        else if ((ABS(lambda(1) - lambda(2)) < 1.0e-4_dp) .AND. (ABS(lambda(2) - lambda(3)) < 1.0e-4_dp)) then
            print *, "Molecule is spherical top."
        else if ((ABS(lambda(1) - lambda(2)) < 1.0e-4_dp) .AND. (ABS(lambda(2) - lambda(3)) > 1.0e-4_dp)) then
            print *, "Molecule is oblate top."
        else if ((ABS(lambda(1) - lambda(2)) > 1.0e-4_dp) .AND. (ABS(lambda(2) - lambda(3)) < 1.0e-4_dp)) then
            print *, "Molecule is prolate top."
        else
            print *, "Molecule is asymmetric top."
        end if
    end subroutine principle_moments_of_inertia

    function rotational_constants(self) result(constants)
        class(molecules), intent(inout) :: self
        real(dp) :: inertia(3,3), lambda(3), work(3), constants(3)
        real(dp) :: mass
        integer :: i

        inertia = 0.0_dp
        lambda = self%inertia_eigvals
        call ord_sort(lambda, work)
        constants = 6.6260755e-34_dp / (8.0_dp * pi * pi * lambda)
        !Unit conversion to MHz
        constants = constants / (1.6605402e-27_dp * 0.529177249e-10_dp * 0.529177249e-10_dp) * 1.0e-6_dp
    end function rotational_constants

    subroutine print_geometry(self) 
        class(molecules), intent(inout) :: self
        integer :: i, j, k, l

        print *, "Number of atoms: ", self%natom
        print *, "Charge: ", self%charge
        print *, "Point group: ", self%point_group
        print *, "Atomic numbers and coordinates:"
        do i = 1, self%natom
            print *, self%zvals(i), self%geom(:,i)
        end do
        
        print *, "Interatomic distances (bohr):        "
        do i = 1, self%natom
            do j = 1, i-1
                print *, i,  j, self%bond_length(i, j)
            end do
        end do

        print *, "Bond angles (degrees):"
        do i = 1, self%natom
            do j = 1, i-1
                do k = 1, j-1
                    if ((self%bond_length(i, j) < 4.0_dp) .and. (self%bond_length(j, k) < 4.0_dp)) then
                        print *, i, j, k, self%bond_angles(i, j, k)
                    end if 
                end do
            end do
        end do

        print *, "Out-of-plane angles (degrees):"
        do i = 1, self%natom
            do k = 1, self%natom
                do j = 1, self%natom
                    do l = 1, j-1
                        if ((i .NE. j) .and. (i .NE. k) .and. (i .NE. l) .and. (j .NE. k) .and. (k .NE. l) .and. &
                            (self%bond_length(i, k) < 4.0_dp) .and. (self%bond_length(j, k) < 4.0_dp) .and. (self%bond_length(k, l) < 4.0_dp)) then
                            print *, i-1, j-1, k-1, l-1, self%out_of_plane_angle(i, j, k, l)
                        end if
                    end do
                end do
            end do
        end do

        print *, "Torsion angles (degrees):"
        do i = 1, self%natom
            do j = 1, i-1
                do k = 1, j-1
                    do l = 1, k-1
                        if ((self%bond_length(i, j) < 4.0_dp) .and. (self%bond_length(j, k) < 4.0_dp) .and. (self%bond_length(k, l) < 4.0_dp)) then
                            print *, i-1, j-1, k-1, l-1, self%torsion(i, j, k, l)
                        end if
                    end do
                end do
            end do
        end do
        
        call self%center_of_mass()
        print *, "Center of mass ", self%com(1), self%com(2), self%com(3)

        call self%principle_moments_of_inertia()
        print *, "rotational constants (MHz): ", self%rotational_constants()
    end subroutine print_geometry
end module mod_molecules