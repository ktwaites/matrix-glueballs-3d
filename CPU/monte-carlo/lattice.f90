module lattice
    use matrix_glueball_3d_parameters
    use elementary_functions

    implicit none

    ! -------------------------- Conventions for lattice moves and plane indexing -------------------------- !
    ! Direction vectors: x == 1
    !                    y == 2
    !                    t == 3
    !
    ! Planes:            yt-plane == 1
    !                    xt-plane == 2
    !                    xy-plane == 3
    ! ------------------------------------------------------------------------------------------------------ !

    ! Array containing all neighbours of all lattice sites
    integer(int32), allocatable :: move(:,:,:)

    contains

    ! Populate move with lattice site neighbours
    subroutine lattice_setup()
        implicit none
        integer :: site, mu, step, next

        ! Allocate move appropriately
        if (LS > LT) then
            allocate(move(LATTICE_VOLUME, 3, -LS:LS))
        else
            allocate(move(LATTICE_VOLUME, 3, -LT:LT))
        endif

        ! Calculate neighbours of all sites
        do site = 1, LATTICE_VOLUME
            mu = 1
            do step = -LS, LS
                next = site + step
                if (ceiling(real(next)/real(LS)) == ceiling(real(site)/real(LS))) then
                    move(site,1,step) = next
                else
                    move(site,1,step) = next - sgn(step)*LS
                endif
            enddo

            mu = 2
            do step = -LS, LS
                next = site + step*LS
                if (ceiling(real(next)/real(SLICE_VOLUME)) == ceiling(real(site)/real(SLICE_VOLUME))) then
                    move(site,2,step) = next
                else
                    move(site,2,step) = next - sgn(step)*SLICE_VOLUME
                endif
            enddo

            mu = 3
            do step = -LT, LT
                next = site + step*SLICE_VOLUME
                if ((next > LATTICE_VOLUME).or.(next < 1)) then
                    move(site,3,step) = next - sgn(step)*LATTICE_VOLUME
                else
                    move(site,3,step) = next
                endif
            enddo
        enddo
    end subroutine

    ! Translate a given number of time steps up the lattice
    integer function translate_time_slice(time_slice, step)
        implicit none
        integer, intent(in) :: time_slice, step

        translate_time_slice = 1 + mod(time_slice + step - 1, LT)
    end function

    ! Return order of sites on which the Metropolis algorithm is carried out
    ! on for a single sweep of the lattice using the checkerboard prescription
    function checkerboard()
        implicit none

        integer, dimension(LATTICE_VOLUME) :: checkerboard
        integer, dimension(LATTICE_VOLUME/2) :: odd, even

        integer :: n, X, Y, T, next_odd, next_even

        next_odd = 1
        next_even = 1

        do n = 1, LATTICE_VOLUME
            X = mod(n, LS)
            if (X==0) X = LS
            Y = mod(ceiling(real(n)/real(LS)),LS)
            if (Y==0) Y = LS
            T = ceiling(real(n)/real(SLICE_VOLUME))
            if (mod(X+Y+T, 2) == 0) then
                even(next_even) = n
                next_even = next_even + 1
            else
                odd(next_odd) = n
                next_odd = next_odd + 1
            endif
        enddo
        checkerboard = [odd, even]
    end function

    ! Implement periodic boundary conditions in coordinates
    elemental function periodic(n)
        implicit none

        integer, intent(in) :: n
        integer :: periodic

        periodic = 1 + modulo(n-1, LS)
    end function

    ! Return site index from coordinates
    function site_index(coordinates)
        implicit none

        integer, intent(in) :: coordinates(2)
        integer :: site_index

        site_index = coordinates(1) + LS*(coordinates(2)-1)
    end function

    ! Return site coordinates from index
    function site_coordinates(index)
        implicit none

        integer, intent(in) :: index
        integer :: site_coordinates(2)

        site_coordinates(1) = periodic(modulo(index, LS))
        site_coordinates(2) = 1 + (index - site_coordinates(1))/LS
    end function
end module