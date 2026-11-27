module D4_loops
    use matrix_glueball_3d_parameters
    use elementary_functions, only: exp_phase
    use linalg, only: trace
    use lattice, only: site_coordinates, site_index, periodic
    
    implicit none

    complex(real64), allocatable :: momentum_phase(:,:,:)

    contains

    ! Pre-calculate all momentum phases to be used in momentum sums
    subroutine setup_momentum()
        implicit none

        integer :: q, x

        ! Allocate momentum phases array
        allocate(momentum_phase(0:LS-1, 0:LS-1, -(LS/2 - 1) : LS/2))

        ! Calculate phase factor for each site for each momentum index q
        do q = -(LS/2 - 1), LS/2
            do x = 0, LS-1
                ! Calculate phase factor
                momentum_phase(x, :, q) &
                = exp_phase(2.0_real64 * PI * real(q, kind=real64) * real(x, kind=real64) &
                / real(LS, kind=real64))
            enddo
        enddo
    end subroutine

    ! ------------------------------------- Calculate matrix operators ------------------------------------- !
    !
    ! An operator with a given set of symmetry transformation coefficients {c_g}, with a given winding w and
    ! a given lattice x-momentum q is given by:
    ! O({c_g}, w, q) = \sum_{g \in G[w]) c_g \sum_{x} e^{2 PI i q x / LS} Re[\sum_{n} e^{-2 PI i w n / NCOL}
    !     phi(n, x, g)]
    !
    ! This function therefore calculates phi(n, x, g) for all x, g on a given time slice from a given twisted
    ! configuration, labelled by n.
    !
    ! Inputs:
    !  * vertex_coordinates(2, <number of sites>): the coordinates of the vertices of the calculated operator
    !    wrt the origin
    !  * propagator_matrix(NCOL, NCOL, SLICE_VOLUME, SLICE_VOLUME): matrix containing propagators from any
    !    site to any other site on the time slice, measured on a single twisted configuration
    !  * calculate_loop(8): boolean array containing whether the symmetry transformation should be calculated 
    !    or not
    !
    ! Outputs:
    !  * phi(LS, LS, 8): translations and symmetry transformations of matrix operator constructed from
    !    vertex_coordinates
    !
    ! -----------------------------------------------------------------------------------------------------  !

    ! Helper function - calculate matrix operator between given coordinates
    pure function single_loop(coordinates, num_sites, propagator_matrix)
        implicit none
        integer(int32), intent(in) :: coordinates(2, num_sites), num_sites
        complex(real32), intent(in) :: propagator_matrix(NCOL, NCOL, SLICE_VOLUME, SLICE_VOLUME)
        complex(real64) :: single_loop

        integer(int32) :: coords(2, num_sites), sites(num_sites), i
        complex(real64) :: loop_matrix(NCOL, NCOL)

        ! Ensure coordinates lie in the range [1, LS]
        coords = periodic(coordinates)

        ! Convert vertex coordinates to site indices
        do i = 1, num_sites
            sites(i) = site_index(coords(:, i))
        enddo

        ! Initialise matrix as propagator from last site to first site
        loop_matrix = cmplx(propagator_matrix(:, :, sites(num_sites), sites(1)), kind=real64)

        ! Multiply matrices in path-order around the loop
        do i = 1, num_sites-1
            loop_matrix = matmul(loop_matrix, propagator_matrix(:, :, sites(i), sites(i+1)))
        enddo

        ! Trace over total propagator matrix to define operator phi[n](x, g)
        single_loop = trace(loop_matrix)
    end function

    ! Calculate all translations and symmetry transformations of a given matrix operator defined through vertex_coordinates
    function loop_propagator(vertex_coordinates, propagator_matrix, calculate_loop) result(phi)
        implicit none
        integer(int32), intent(in) :: vertex_coordinates(:,:)
        complex(real32), intent(in) :: propagator_matrix(NCOL, NCOL, SLICE_VOLUME, SLICE_VOLUME)
        logical, intent(in) :: calculate_loop(8)
        complex(real64) :: phi(LS, LS, 8)

        complex(real64) :: loop_matrix(NCOL, NCOL)
        integer(int32), allocatable :: translated_coords(:,:), coords(:,:), sites(:)
        integer(int32) :: n, g, x, y, i

        ! Get number of sites on path
        n = size(vertex_coordinates, dim=2)

        ! Allocate site coordiante and index arrays
        allocate(translated_coords(2, n), coords(2, n), sites(n))

        ! Calculate all symmetry transformations and translations of the loop defined by vertex_coordinates
        do concurrent (x = 0:LS-1, y = 0:LS-1, g = 1:8)
            ! Translate operator by given amount
            translated_coords(1, :) = periodic(vertex_coordinates(1, :) + x)
            translated_coords(2, :) = periodic(vertex_coordinates(2, :) + y)

            ! Select correct symmetry transformation
            select case (g)
            case(1)
                ! g=1 -> Id map
                coords = translated_coords
            case(2)
                ! g=2 -> Rotation, 90 degrees anticlockwise = R
                coords(1, :) = -translated_coords(2, :)
                coords(2, :) = translated_coords(1, :)
            case(3)
                ! g=3 -> Rotation, 180 degrees anticlockwise = R2
                coords = -translated_coords
            case(4)
                ! g=4 -> Rotation, 270 degrees anitclockwise = R3
                coords(1, :) = translated_coords(2, :)
                coords(2, :) = -translated_coords(1, :)
            case(5)
                ! g=5 -> Parity inversion (x -> - x) = P
                coords(1, :) = -translated_coords(1, :)
                coords(2, :) = translated_coords(2, :)
            case(6)
                ! g=6 -> PR
                coords(1, :) = translated_coords(2, :)
                coords(2, :) = translated_coords(1, :)
            case(7)
                ! g=7 -> PR2
                coords(1, :) = translated_coords(1, :)
                coords(2, :) = -translated_coords(2, :)
            case(8)
                ! g=8 -> PR3
                coords(1, :) = -translated_coords(2, :)
                coords(2, :) = -translated_coords(1, :)
            end select

            ! Calculate propagator
            if (calculate_loop(g)) then
                ! Ensure coordinates lie in the range [1, LS]
                coords = periodic(coords)

                ! Convert vertex coordinates to site indices
                do i = 1, n
                    sites(i) = site_index(coords(:, i))
                enddo

                ! Initialise matrix as propagator from last site to first site
                loop_matrix = cmplx(propagator_matrix(:, :, sites(n), sites(1)), kind=real64)

                ! Multiply matrices in path-order around the loop
                do i = 1, n-1
                    loop_matrix = matmul(loop_matrix, propagator_matrix(:, :, sites(i), sites(i+1)))
                enddo

                ! Trace over total propagator matrix to define operator phi[n](x, g)
                phi(x, y, g) = trace(loop_matrix)
            else
                ! Set loop to zero as it will not be used
                phi(x, y, g) = 0.0_real64
            endif
        enddo
    end function
end module D4_loops