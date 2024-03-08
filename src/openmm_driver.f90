subroutine openmm(i3n, v)
    implicit real(kind=8) (a-h, o-z)
    
    integer, intent(in) :: i3n
    real(kind=8), intent(out) :: v

    integer, parameter :: nd1 = 100
    integer, parameter :: ndp = 10

    common /qpdot/ q(3 * nd1), pdot(3 * nd1)
    common /constn/ c1, c2, c3, c4, c5, c6, c7, pi, halfpi, twopi

    integer :: i

    !> write positions of atoms to file
    open(10, file="positions.txt")

    do i = 1, int(i3n / 3)
        write(10, "(3f15.9)") q(3 * i - 2), q(3 * i - 1), q(3 * i)
    end do

    close(10)

    !> calculate potential energy by OpenMM
    call execute_command_line("python3 openmm_driver.py -t energy -f positions.txt")

    !> read energy from file, in kcal/mol
    open(20, file="energy.txt", status="old")

    read(20, *) v

    close(20)

    !> convert to VENUS96 integration unit
    v = v * c1

end subroutine openmm

subroutine dopenmm(i3n)
    implicit real(kind=8) (a-h, o-z)

    integer, intent(in) :: i3n

    integer, parameter :: nd1 = 100
    integer, parameter :: ndp = 10

    common /qpdot/ q(3 * nd1), pdot(3 * nd1)
    common /constn/ c1, c2, c3, c4, c5, c6, c7, pi, halfpi, twopi

    integer :: i

    !> write positions of atoms to file
    open(10, file="positions.txt")

    do i = 1, int(i3n / 3)
        write(10, "(3f15.9)") q(3 * i - 2), q(3 * i - 1), q(3 * i)
    end do

    close(10)

    !> calculate force by OpenMM
    call execute_command_line("python3 openmm_driver.py -t forces -f positions.txt")

    !> read force from file, in kcal/mol/Angstrom
    open(20, file="forces.txt", status="old")

    do i = 1, int(i3n / 3)
        read(20, *) pdot(3 * i - 2), pdot(3 * i - 1), pdot(3 * i)
    end do

    close(20)

    !> convert to integration unit
    pdot(:) = -1.0d0 * pdot(:) * c1

end subroutine dopenmm
