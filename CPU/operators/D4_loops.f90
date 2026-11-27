module D4_loops
    use matrix_glueball_3d_parameters
    use elementary_functions, only: exp_phase
    use linalg, only: trace
    use lattice, only: site_coordinates, site_index, periodic
    
    implicit none

    complex(real64), allocatable :: momentum_phase(:,:,:)
    integer(int32) :: multiply_D4_elements(8,8)

    contains

    ! Setup function for loop calculations
    ! Pre-calculate all momentum phases to be used in momentum sums
    ! Assign D4 group multiplication table elements
    ! table(i, j) = g(i) * g(j)
    ! Group element indices:
    !  * g=1 -> id
    !  * g=2 -> R
    !  * g=3 -> R2
    !  * g=4 -> R3
    !  * g=5 -> P
    !  * g=6 -> PR
    !  * g=7 -> PR2
    !  * g=8 -> PR3
    subroutine setup_loops()
        implicit none

        integer(int32) :: q, x, Id(2,2), R(2,2), P(2,2), result(2,2), group_elements(2,2,8), i, j, k

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

        !! Define D4 group multiplication table by explicitely multiplying matrices
        ! Define group elements
        Id(1, :) = [1, 0]
        Id(2, :) = [0, 1]
        R(1, :) = [0, -1]
        R(2, :) = [1, 0]
        P(1, :) = [-1, 0]
        P(2, :) = [0, 1]
        group_elements(:,:,1) = Id
        group_elements(:,:,2) = R
        group_elements(:,:,3) = matmul(R, R)
        group_elements(:,:,4) = matmul(R, matmul(R, R))
        group_elements(:,:,5) = P
        group_elements(:,:,6) = matmul(P, R)
        group_elements(:,:,7) = matmul(P, matmul(R, R))
        group_elements(:,:,8) = matmul(P, matmul(R, matmul(R, R)))

        ! Multiply group elements and compare with look-up table to find index
        do j = 1, 8
            do i = 1, 8
                ! result = g(i) * g(j)
                result = matmul(group_elements(:,:,i), group_elements(:,:,j))

                ! Compare with group elements to find what group element this is
                multiply_D4_elements(i, j) = 0
                do k = 1, 8
                    if (all(result == group_elements(:,:,k))) then
                        multiply_D4_elements(i, j) = k
                        exit
                    endif
                enddo
                
                ! Check a group element was found
                if (multiply_D4_elements(i, j) == 0) then
                    write(*, '(a, 4x, a)') "[Error][setup_loops]", &
                    "Product of group elements not found in list of group elements"
                    error stop
                endif
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


    ! ------------------------------------------ Loop degeneracies ------------------------------------------ !
    !
    ! There are two types of loop degeneracy we must find:
    !  * Intra-loop: degeneracies between transformations of the same loop. These affect which states a given
    !    loop will project onto, and the weighting each symmetry transformation should be given when
    !    constructing states from that loop
    !  * Inter-loop: degeneracies between different loops. We say two loops are degenerate with each other if
    !    one is a symmetry transformation of the other. These affect which loops should be included in the
    !    operator basis for a given state, as degenerate loops will lead to a singular correlation matrix.
    !
    ! ------------------------------------------------------------------------------------------------------- !

    ! Check two loops are close to each other
    pure function is_close(x, y)
        implicit none
        complex(real64), intent(in) :: x, y
        logical :: is_close

        is_close = (abs(y - x) < 1e-12)
    end function

    ! Determine which sets of loops are degenerate with each other
    ! phi(g) = \sum_{x,y} loop_propagator(x,y,g) 
    function interloop_degeneracy(phi) result(distinct_loops)
        implicit none
        complex(real64), intent(in) :: phi(:,:)
        logical :: distinct_loops(size(phi, dim=2))

        integer(int32) :: num_sets_loops, i, j, g, g_prime
        logical :: loop_equivilance_matrix(8,8)

        ! Get number of sets of loops
        num_sets_loops = size(phi, dim=2)

        ! Default all loops to be distinct from one another
        distinct_loops = .true.

        ! Compare loop j with all previous loops. If it is equivilant to any previous loop, distinct_loops(j) is set to false
        do j = 2, num_sets_loops
            do i = 1, j-1
                ! Compute the logical matrix with elements g, g' as true if phi(g, i) = phi(gg', j) and false otherwise
                do concurrent (g = 1:8, g_prime = 1:8)
                    loop_equivilance_matrix(g, g_prime) = is_close(&
                        phi(g, i), &
                        phi(multiply_D4_elements(g, g_prime), j)&
                    )
                enddo

                ! Loops phi(i) and phi(j) are defined as equivilant if there exists a g' such that phi(g, i) = phi(gg', j) for all g. This means phi(i) and phi(j) are equivilant if the logical matrix contains a column where all elements are true.
                if (any(all(loop_equivilance_matrix, dim=1))) then
                    distinct_loops(j) = .false.
                    exit
                endif
            enddo
        enddo
    end function
end module D4_loops