module staggered_fermion_matrix
    use matrix_glueball_3d_parameters
    use linalg, only: hermc, solve_linear_system, trace
    use lattice

    implicit none

    ! Arrays containing neighbours
    integer(int32), allocatable :: odd(:), even(:), inverse_sites(:), &
    neighbours(:,:), next_neighbours(:,:), sub_neighbours(:,:,:)

    contains

    ! ----------------------------------------- Neighbour functions ----------------------------------------- !

    ! Return site one away from input site according to index given
    integer elemental function nearest_neighbour(site, index)
        implicit none
        integer(int32), intent(in) :: site, index
        select case(index)
        case(1)
            nearest_neighbour = move(site,1,1)
        case(2)
            nearest_neighbour = move(site,2,1)
        case(3)
            nearest_neighbour = move(site,1,-1)
        case(4)
            nearest_neighbour = move(site,2,-1)
        end select
    end function

    ! Return site two away from input site according to index given
    integer elemental function next_nearest_neighbour(site, index)
        implicit none
        integer(int32), intent(in) :: site, index
        select case(index)
        case(1)
            next_nearest_neighbour = move(site,1,2)
        case(2)
            next_nearest_neighbour = move(move(site,2,1), 1, 1)
        case(3)
            next_nearest_neighbour = move(site,2,2)
        case(4)
            next_nearest_neighbour = move(move(site,1,-1), 2, 1)
        case(5)
            next_nearest_neighbour = move(site,1,-2)
        case(6)
            next_nearest_neighbour = move(move(site,2,-1), 1, -1)
        case(7)
            next_nearest_neighbour = move(site,2,-2)
        case(8)
            next_nearest_neighbour = move(move(site,1,1), 2, -1)
        end select
    end function


    ! ------------------------------------------- Setup functions ------------------------------------------- !

    ! Complete all necessary setup for matrix module
    subroutine setup_matrix()
        implicit none
        integer(int32) :: site, i, sublattice, comp, k, n, X, Y, T, next_odd, next_even

        ! Allocate arrays
        allocate(neighbours(4,SLICE_VOLUME), next_neighbours(8,SLICE_VOLUME))

        ! Find nearest and next nearest neighbours for all sites
        do site = 1, SLICE_VOLUME
            do i = 1, 4
                neighbours(i,site) = nearest_neighbour(site,i)
            enddo
            do i = 1, 8
                next_neighbours(i,site) = next_nearest_neighbour(site,i)
            enddo
        enddo

        ! Setup odd/even lattice basis functions
        SLICE_SUBSET_VOLUME = SLICE_VOLUME/2

        next_odd = 1
        next_even = 1

        allocate(odd(SLICE_VOLUME/2), even(SLICE_VOLUME/2), inverse_sites(SLICE_VOLUME))

        do n = 1, SLICE_VOLUME
            X = mod(n, LS)
            if (X==0) X = LS
            Y = mod(ceiling(real(n)/real(LS)),LS)
            if (Y==0) Y = LS
            T = 1
            if (mod(X+Y+T, 2) == 0) then
                even(next_even) = n
                next_even = next_even + 1
            else
                odd(next_odd) = n
                next_odd = next_odd + 1
            endif
        enddo

        do k = 1, SLICE_VOLUME/2
            inverse_sites(even(k)) = k
            inverse_sites(odd(k)) = k
        enddo

        ! Allocate arrays
        allocate(sub_neighbours(8,SLICE_SUBSET_VOLUME,2))

        ! Find next nearest neighbours in odd/even basis
        do sublattice = 1, 2
            do i = 1, SLICE_SUBSET_VOLUME
                ! Find i in full lattice basis
                if (sublattice == 1) then
                    site = odd(i)
                else
                    site = even(i)
                endif

                do comp = 1, 8
                    sub_neighbours(comp,i,sublattice) = inverse_sites(next_neighbours(comp, site))
                enddo
            enddo
        enddo
    end subroutine


    ! ------------------------------------- Fill sparse representations ------------------------------------- !

    !! Fill sparse representation of M^dagger
    ! This is done by storing only the non-zero elements of M^dagger = -M_0, in the basis of the full lattice,
    ! in the following order:
    ! M^dagger(colour index, colour index, column index, site index)
    function define_M_dagger(links) result(M_dagger)
        implicit none
        complex(real32), dimension(NCOL,NCOL,SLICE_VOLUME,2), intent(in) :: links
        complex(real32), dimension(NCOL,NCOL,4,SLICE_VOLUME) :: M_dagger

        integer(int32) :: site
        logical :: oscillate

        oscillate = .true. ! used to apply osicllating factors appropriately
        ! NOTE: this method depends on LS being even!
        do site = 1, SLICE_VOLUME
            ! M(n,n+x)
            M_dagger(:,:,1,site) = -links(:,:,site,1)

            ! M(n,n+y)
            M_dagger(:,:,2,site) = links(:,:,site,2)
            if (.not.oscillate) then
                M_dagger(:,:,2,site) = -M_dagger(:,:,2,site)
            endif

            ! M(n,n-x)
            M_dagger(:,:,3,site) = hermc(links(:,:,neighbours(3,site),1))

            ! M(n,n-y)
            M_dagger(:,:,4,site) = hermc(links(:,:,neighbours(4,site),2))
            if (oscillate) then
                M_dagger(:,:,4,site) = -M_dagger(:,:,4,site)
            endif

            oscillate = .not.oscillate
        enddo
    end function

    !! Fill sparse representation of K = M_0^dagger M_0
    ! This is done by only storing the non-zero elements of K, seperated by odd/even sublattice, in the
    ! following order:
    ! K(colour index, colour index, column index, odd/even site index, odd/even lattice subset)
    function define_K(links) result(K)
        implicit none
        complex(real32), dimension(NCOL, NCOL, SLICE_VOLUME, 2), intent(in) :: links
        complex(real32), dimension(NCOL, NCOL, 8, SLICE_SUBSET_VOLUME, 2) :: K

        complex(real32), dimension(NCOL,NCOL) :: right_link, up_link, left_link, down_link
        integer(int32) :: r, u, l, d, ul, ll, ld, dd, dr
        logical :: oscillate

        integer(int32) :: sublattice, slice, site_index, site, site_oe

        ! Calculate odd and even sublattices individually
        do sublattice = 1, 2 ! 1 = odd, 2 = even
            ! Initialise oscillate for given sublattice (oscillate = eta_s(0))
            ! osicllate = false -> eta_s(0) = +1 -> odd
            ! oscillate = true -> eta_s(0) = -1 -> even
            oscillate = (sublattice == 2)

            ! Iterate over horizontal slices
            do slice = 1, LS
                ! Update osicllate for this slice
                oscillate = .not. oscillate

                ! Iterate over sublattice sites on given slice
                do site_index = 1, LS/2
                    ! Calculate site number in odd/even basis
                    site_oe = (slice-1)*LS/2 + site_index

                    ! Find site number in whole lattice basis
                    if (sublattice == 1) then
                        site = odd(site_oe)
                    else
                        site = even(site_oe)
                    endif

                    !! Calculate K(:,:,:,site,sublattice) using oscillate
                    ! compute nearest neighbour locations
                    r = neighbours(1,site)
                    u = neighbours(2,site)
                    l = neighbours(3,site)
                    d = neighbours(4,site)

                    ! compute second nearest neighbour locations
                    ul = next_neighbours(4,site)
                    ll = next_neighbours(5,site)
                    ld = next_neighbours(6,site)
                    dd = next_neighbours(7,site)
                    dr = next_neighbours(8,site)

                    ! find inner link variables (the ones used more than once)
                    right_link = links(:,:,site,1)
                    up_link = links(:,:,site,2)
                    left_link = hermc(links(:,:,l,1))
                    down_link = hermc(links(:,:,d,2))

                    ! K(n,n+2x)
                    K(:,:,1,site_oe,sublattice) = -matmul(right_link, links(:,:,r,1))

                    ! K(n,n+x+y)
                    K(:,:,2,site_oe,sublattice) = &
                    matmul(right_link, links(:,:,r,2)) - matmul(up_link, links(:,:,u,1))
                    if (oscillate) then
                        K(:,:,2,site_oe,sublattice) = -K(:,:,2,site_oe,sublattice)
                    endif

                    ! K(n,n+2y)
                    K(:,:,3,site_oe,sublattice) = -matmul(up_link, links(:,:,u,2))

                    ! K(n,n-x+y)
                    K(:,:,4,site_oe,sublattice) = &
                    matmul(up_link, hermc(links(:,:,ul,1))) - matmul(left_link, links(:,:,l,2))
                    if (oscillate) then
                        K(:,:,4,site_oe,sublattice) = -K(:,:,4,site_oe,sublattice)
                    endif

                    ! K(n,n-2x)
                    K(:,:,5,site_oe,sublattice) = -matmul(left_link, hermc(links(:,:,ll,1)))

                    ! K(n,n-x-y)
                    K(:,:,6,site_oe,sublattice) = &
                    matmul(left_link, hermc(links(:,:,ld,2))) - matmul(down_link, hermc(links(:,:,ld,1)))
                    if (oscillate) then
                        K(:,:,6,site_oe,sublattice) = -K(:,:,6,site_oe,sublattice)
                    endif

                    ! Z(n,n-2y)
                    K(:,:,7,site_oe,sublattice) = -matmul(down_link, hermc(links(:,:,dd,2)))

                    ! K(n,n+x-y)
                    K(:,:,8,site_oe,sublattice) = &
                    matmul(down_link, links(:,:,d,1)) - matmul(right_link, hermc(links(:,:,dr,2)))
                    if (oscillate) then
                        K(:,:,8,site_oe,sublattice) = -K(:,:,8,site_oe,sublattice)
                    endif
                enddo
            enddo
        enddo
    end function


    ! ----------------------------------- Sparse matrix - vector product ----------------------------------- !

    ! Multiply K by a vector
    function K_vector_product(K, v, mass, sublattice) result(Kv)
        implicit none
        complex(real32), dimension(NCOL, NCOL, 8, SLICE_SUBSET_VOLUME, 2), intent(in) :: K
        complex(real32), dimension(NCOL, NCOL, SLICE_SUBSET_VOLUME), intent(in) :: v
        integer(int32), intent(in) :: sublattice
        real(real32), intent(in) :: mass

        complex(real32), dimension(NCOL, NCOL, SLICE_SUBSET_VOLUME) :: Kv

        complex(real32), dimension(NCOL, NCOL, 8) :: Kv_site
        integer(int32) :: i, j, comp

        ! Iterate through components of resulting vector Kv
        do i = 1, SLICE_SUBSET_VOLUME
            ! Iterate through 8 non-zero components in row of K
            do comp = 1, 8
                ! Calculate component of v to be added
                j = sub_neighbours(comp,i,sublattice)

                ! Add relevant component to Kv_site
                Kv_site(:,:,comp) = matmul(K(:,:,comp,i,sublattice), v(:,:,j))
            enddo

            ! Add to Kv
            Kv(:,:,i) = sum(Kv_site, dim=3)
        enddo

        ! Add diagonal term
        Kv = Kv + (4.0 + mass**2)*v
    end function


    ! --------------------------------------- Block matrix functions --------------------------------------- !

    ! Dot product of block vectors
    function block_dot_product(v, w) result(vTw)
        implicit none
        complex(real32), dimension(NCOL, NCOL, SLICE_SUBSET_VOLUME), intent(in) :: v, w
        complex(real32), dimension(NCOL, NCOL) :: vTw

        integer(int32) :: i

        vTw = cmplx(0.0,0.0)
        do i = 1, SLICE_SUBSET_VOLUME
            vTw = vTw + matmul(hermc(v(:,:,i)), w(:,:,i))
        enddo
    end function


    ! ------------------------------ Apply Block Conjugate Gradient algorithm ------------------------------ !

    ! Use block Conjugate Gradient to solve for one column of the inverse of K_o or K_e, starting from an
    ! initial guess x_0
    subroutine solve_BCG(K, x, column, mass, epsilon, sublattice, iterations_completed)
        implicit none
        complex(real32), dimension(NCOL, NCOL, 8, SLICE_SUBSET_VOLUME, 2), intent(in) :: K
        complex(real32), dimension(NCOL, NCOL, SLICE_SUBSET_VOLUME), intent(out) :: x
        integer(int32), intent(in) :: column, sublattice
        integer(int32), intent(out), optional :: iterations_completed
        real(real32), intent(in) :: mass, epsilon

        complex(real32), dimension(NCOL,NCOL,SLICE_SUBSET_VOLUME) :: direction, residual, Kd
        complex(real32), dimension(NCOL,NCOL) :: alpha_BCG, beta_BCG, delta_old, delta_new
        real(real32) :: matrix_norm, old_matrix_norm, delta_0

        integer(int32) :: iteration, colour, i

        ! Initialize vectors to d = r = b-K_s*x, where b is the col'th unit vector
        residual = -K_vector_product(K, x, mass, sublattice)
        
        ! Set r = b-Ax
        do colour = 1, NCOL
            residual(colour,colour,column) = residual(colour,colour,column) + cmplx(1.0,0.0)
        enddo
        direction = residual
        matrix_norm = 100.0

        ! Find magnitude of initial residue
        delta_new = block_dot_product(residual, residual)

        ! Save initial residual
        delta_0 = sqrt(real(trace(delta_new), kind=real32))

        ! Conduct iterations
        do iteration = 1, SLICE_VOLUME
            ! First check if the residual is within a ball of radius epsilon around the solution
            ! We do this by checking if all elements of residual_squared have fallen below epsilon

            ! Calculate Frobenius norm of residual
            old_matrix_norm = matrix_norm
            matrix_norm = sqrt(real(trace(delta_new), kind=real32))

            if (matrix_norm < epsilon*delta_0) then
                ! The residual is small enough, so the problem is approximately solved
                ! print *, "Number of iterations completed: ", iteration-1
                if (present(iterations_completed)) iterations_completed = iteration-1
                return
            endif

            if ((matrix_norm > 1.1*old_matrix_norm)) then
                if (present(iterations_completed)) iterations_completed = SLICE_VOLUME
                return
            endif

            ! Check if the solution has blown up
            if (matrix_norm > 1e16*delta_0) then
                ! The residual is too large, and CG has not worked
                ! print *, "Column ", col, " diverges after ", iteration-1, " iterations"
                if (present(iterations_completed)) iterations_completed = SLICE_VOLUME
                return
            endif

            ! Compute K*d
            Kd = K_vector_product(K, direction, mass, sublattice)

            ! Find alpha (factor by which to update x along diretion d)
            alpha_BCG = solve_linear_system(block_dot_product(direction, Kd), delta_new)

            ! Update x
            do i = 1, SLICE_SUBSET_VOLUME
                x(:,:,i) = x(:,:,i) + matmul(direction(:,:,i), alpha_BCG)
            enddo

            ! Check if this is the last iteration
            if (iteration == SLICE_SUBSET_VOLUME) then
                print *, "Column ", column, " required the full number of iterations"
                if (present(iterations_completed)) iterations_completed = SLICE_VOLUME
                Return
            else
                ! Update r. Every LS iterations, recalculate r to avoid round-off error
                if (mod(iteration, LS) == 0) then
                    residual = -K_vector_product(K, x, mass, sublattice)
                    ! Set r = b-Ax
                    do colour = 1, NCOL
                        residual(colour,colour,column) = residual(colour,colour,column)&
                        + cmplx(1.0,0.0)
                    enddo
                else
                    do i = 1, SLICE_SUBSET_VOLUME
                        residual(:,:,i) = residual(:,:,i) - matmul(Kd(:,:,i), alpha_BCG)
                    enddo
                endif

                ! Find new magnitude of residual
                delta_old = delta_new
                delta_new = block_dot_product(residual, residual)

                ! Find beta (factor by which to update direction vector)
                beta_BCG = solve_linear_system(delta_old, delta_new)

                ! Update direcion vector
                do i = 1, SLICE_SUBSET_VOLUME
                    direction(:,:,i) = residual(:,:,i) + matmul(direction(:,:,i), beta_BCG)
                enddo
            endif
        enddo
    end subroutine

    ! Solve for the inverse of K_o or K_e
    subroutine invert_BCG(K, K_s_inv, mass, epsilon, sublattice)
        implicit none
        complex(real32), dimension(NCOL, NCOL, 8, SLICE_SUBSET_VOLUME, 2), intent(in) :: K
        complex(real32), dimension(NCOL, NCOL, SLICE_SUBSET_VOLUME, SLICE_SUBSET_VOLUME), intent(out) :: K_s_inv
        integer(int32), intent(in) :: sublattice
        real(real32), intent(in) :: mass, epsilon

        integer(int32) :: column

        ! Solve for each column of the inverse sequentially
        do column = 1, SLICE_SUBSET_VOLUME            
            ! Use conjugate gradient to solve for the next column of A_inv
            call solve_BCG(K, K_s_inv(:,:,:,column), column, mass, epsilon, sublattice)
        enddo
    end subroutine


    ! ------------------------ Solve for the inverse of M using M^-1 = M^dagger K^-1 ------------------------ !

    ! Construct M_inv using M_dagger and K^-1
    function construct_M_inv(M_dagger, K_o_inv, K_e_inv, mass) result(M_inv)
        implicit none
        complex(real32), dimension(NCOL, NCOL, 4, SLICE_VOLUME), intent(in) :: M_dagger
        complex(real32), dimension(NCOL, NCOL, SLICE_SUBSET_VOLUME, SLICE_SUBSET_VOLUME), intent(in) :: K_o_inv, K_e_inv
        complex(real32), dimension(NCOL, NCOL, SLICE_VOLUME, SLICE_VOLUME) :: M_inv
        real(real32), intent(in) :: mass

        complex(real32), dimension(NCOL, NCOL, 4) :: M_inv_site_site
        integer(int32) :: n, m, i, j, k, comp

        !! M_inv splits into block form in the odd/even basis.
        ! The diagonal blocks (elements (n,m) of M_inv for which type(n) = type(m)) and found by scaling
        ! K_inv by the mass: M_inv(n,m) = mass*K_inv(n,m)
        ! The off-diagonal blocks are given by M_inv = M_0^dagger*K_inv

        ! Iterate through odd columns of M_inv, then through even columns
        do j = 1, SLICE_SUBSET_VOLUME
            !! Odd columns
            m = odd(j)

            ! Iterate through odd rows of M_inv (diagonal odd-odd block), then through even rows
            ! (off-diagonal even-odd block)
            do i = 1, SLICE_SUBSET_VOLUME
                !! Odd rows
                n = odd(i)

                ! For these elements, M_inv = mass*K_inv
                M_inv(:,:,n,m) = mass * K_o_inv(:,:,i,j)

                !! Even rows
                n = even(i)

                ! Iterate through neighbours of n
                do comp = 1, 4
                    ! Find index of K_o_inv
                    k = inverse_sites(neighbours(comp, n))

                    ! Add to inverse matrix element
                    M_inv_site_site(:,:,comp) = matmul(M_dagger(:,:,comp,n), K_o_inv(:,:,k,j))
                enddo

                ! Sum over neighbours to find inverse matrix element
                M_inv(:,:,n,m) = sum(M_inv_site_site, dim=3)
            enddo

            !! Even columns
            m = even(j)

            ! Iterate through odd rows of M_inv (off-diagonal odd-even block), then through even rows
            ! (diagonal even-even block)
            do i = 1, SLICE_SUBSET_VOLUME
                !! Odd rows
                n = odd(i)

                ! Iterate through neighbours of n
                do comp = 1, 4
                    ! Find index of K_e_inv
                    k = inverse_sites(neighbours(comp, n))

                    ! Add to inverse matrix element
                    M_inv_site_site(:,:,comp) = matmul(M_dagger(:,:,comp,n), K_e_inv(:,:,k,j))
                enddo

                ! Sum over neighbours to find inverse matrix element
                M_inv(:,:,n,m) = sum(M_inv_site_site, dim=3)

                !! Even rows
                n = even(i)

                ! For these elements, M_inv = mass*K_inv
                M_inv(:,:,n,m) = mass * K_e_inv(:,:,i,j)
            enddo
        enddo
    end function

    ! Invert matrix M at single mass using vanilla method
    function invert_matrix_vanilla(links, mass, epsilon) result(M_inv)
        implicit none
        complex(real32), intent(in) :: links(NCOL, NCOL, SLICE_VOLUME, 2)
        real(real32), intent(in) :: mass, epsilon
        complex(real32) :: M_inv(NCOL, NCOL, SLICE_VOLUME, SLICE_VOLUME)

        complex(real32), dimension(NCOL, NCOL, 4, SLICE_VOLUME) :: M_dagger
        complex(real32), dimension(NCOL, NCOL, 8, SLICE_SUBSET_VOLUME, 2) :: K
        complex(real32), dimension(NCOL, NCOL, SLICE_SUBSET_VOLUME, SLICE_SUBSET_VOLUME) :: K_o_inv, K_e_inv
        integer(int32) :: column

        ! Define sparse representation (valid for all masses and epsilons)
        M_dagger = define_M_dagger(links)
        K = define_K(links)

        ! Initialise all vectors to 0
        K_o_inv = cmplx(0.0,0.0)
        K_e_inv = cmplx(0.0,0.0)

        ! Invert odd and even blocks of K seperately
        do column = 1, SLICE_SUBSET_VOLUME
            call solve_BCG(K, K_o_inv(:,:,:,column), column, mass, epsilon, 1)
            call solve_BCG(K, K_e_inv(:,:,:,column), column, mass, epsilon, 2)
        enddo

        ! Construct M_inv
        M_inv = construct_M_inv(M_dagger, K_o_inv, K_e_inv, mass)
    end function
end module staggered_fermion_matrix