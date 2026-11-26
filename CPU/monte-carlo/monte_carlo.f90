module monte_carlo
    use matrix_glueball_3d_parameters
    use elementary_functions
    use linalg, only: hermc, det_2x2, determinant
    use lattice
    use luxury

    implicit none

    contains

    ! -------------------------------------- Initialisation functions -------------------------------------- !

    ! Initialise ranlux generator using system clock
    subroutine initialise_ranlux(offset)
        implicit none
        integer(int32), optional, intent(in) :: offset
        integer(int32) :: initialize

        ! Initialize ranlux generator
        call system_clock(initialize)

        if (present(offset)) then
            call rluxgo(3, initialize+offset, 0, 0)
        else
            call rluxgo(3, initialize, 0, 0)
        endif
    end subroutine

    ! Load gauge field into memory from file
    subroutine load_gauge_field(gauge_field)
        complex(real32), intent(out) :: gauge_field(NCOL, NCOL, LATTICE_VOLUME, 3)

        character(len=127) :: gauge_field_filename, msg
        character(len=15) :: file_beta, file_LS, file_LT, file_ncol, file_index
        integer(int32) :: stat
        logical :: file_exists

        ! Write gauge field parameters to strings
        write(file_beta, '(f0.6)') BETA
        write(file_LS, '(i0)') LS
        write(file_LT, '(i0)') LT
        write(file_ncol, '(i0)') NCOL
        write(file_index, '(i0)') CONFIG_START

        ! Construct file name
        gauge_field_filename = &
            trim(FILENAME) // &
            trim(file_LS) // "x" // trim(file_LS) // "x" // trim(file_LT) // &
            "nc" // trim(file_ncol) // &
            "b" // trim(file_beta) // &
            "n" // trim(file_index) // &
            ".dat"

        ! Check file exists at specified directory
        inquire(file = trim(FILEPATH) // trim(gauge_field_filename), exist=file_exists)
        if (.not. file_exists) then
            write(*, '(a, 4x, a)') "[Error][load_gauge_field]", &
            "Gauge field file could not be found"
            write(*, '(a, a)') "File: ", trim(FILEPATH) // trim(gauge_field_filename)
            error stop
        endif
        
        ! Open gauge file from file in specified directory
        open(&
            10, &
            file = trim(FILEPATH) // trim(gauge_field_filename), &
            form='unformatted', &
            access='stream', &
            status='old', &
            action='read', &
            iostat=stat, &
            iomsg=msg &
        )
        
        ! Check gauge field was opened correctly
        if (stat /= 0) then
            write(*, '(a, 4x, a)') "[Error][load_gauge_field]", trim(msg)
            error stop
        endif

        ! Read gauge field into container
        read(10) gauge_field
        close(10)
    end subroutine

    ! Save gauge field to file
    subroutine save_gauge_field(gauge_field, gauge_field_index)
        complex(real32), intent(in) :: gauge_field(NCOL, NCOL, LATTICE_VOLUME, 3)
        integer(int32), intent(in), optional :: gauge_field_index

        character(len=127) :: gauge_field_filename, msg
        character(len=15) :: file_beta, file_LS, file_LT, file_ncol, file_index
        integer(int32) :: stat, config_index
        logical :: file_exists

        ! Write gauge field parameters to strings
        write(file_beta, '(f0.6)') BETA
        write(file_LS, '(i0)') LS
        write(file_LT, '(i0)') LT
        write(file_ncol, '(i0)') NCOL

        ! Construct file name
        gauge_field_filename = &
            trim(FILENAME) // &
            trim(file_LS) // "x" // trim(file_LS) // "x" // trim(file_LT) // &
            "nc" // trim(file_ncol) // &
            "b" // trim(file_beta)
        
        ! Determine configuration index to assign to this gauge field by looking for the first missing configuration number, or use value passed
        if (present(gauge_field_index)) then
            ! Set configuration index
            config_index = gauge_field_index

            ! Write trial configurations index to string
            write(file_index, '(i0)') config_index
        else
            config_index = 1
            do
                ! Write trial configurations index to string
                write(file_index, '(i0)') config_index

                ! Check if file exists
                inquire(&
                    file = trim(FILEPATH) // trim(gauge_field_filename) // "n" // trim(file_index) // ".dat", &
                    exist=file_exists &
                )

                ! Update trial configuration index if file already exists
                if (file_exists) then
                    config_index = config_index + 1
                else
                    exit
                endif

                ! Throw error if configuration index becomes very large
                if (config_index >= 1e6) then
                    write(*, '(a, 4x, a)') "[Error][save_gauge_field]", &
                    "Error in configration index determination"
                    error stop
                endif
            end do
        endif

        ! Update gauge field file name
        gauge_field_filename = &
            trim(gauge_field_filename) // &
            "n" // trim(file_index) // &
            ".dat"

        ! Check file exists at specified directory
        inquire(file = trim(FILEPATH) // trim(gauge_field_filename), exist=file_exists)
        if (.not. file_exists) then
            write(*, '(a, 4x, a)') "[Error][save_gauge_field]", &
            "Gauge field file could not be found"
            write(*, '(a, a)') "File: ", trim(FILEPATH) // trim(gauge_field_filename)
            error stop
        endif
        
        ! Open gauge field file
        open(&
            10, &
            file = trim(FILEPATH) // trim(gauge_field_filename), &
            form='unformatted', &
            access='stream', &
            status='new', &
            action='write', &
            iostat=stat, &
            iomsg=msg &
        )

        ! Check gauge field was opened correctly
        if (stat /= 0) then
            write(*, '(a, 4x, a)') "[Error][save_gauge_field]", trim(msg)
            error stop
        endif

        ! Save gauge field to file
        write(10) gauge_field
        close(10)
    end subroutine

    ! Setup monte carlo data
    subroutine setup_monte_carlo(gauge_field, sweep_sequence, random_start, offset)
        implicit none

        complex(real32), intent(out) :: gauge_field(NCOL, NCOL, LATTICE_VOLUME, 3)
        integer(int32), intent(out) :: sweep_sequence(LATTICE_VOLUME)
        logical, intent(in), optional :: random_start
        integer(int32), intent(in), optional :: offset

        logical :: random_start_

        ! Default ranlux generator to a random seed using system clock
        if (present(random_start)) then
            random_start_ = random_start
        else
            random_start_ = .true.
        endif

        ! Initialise ranlux generator
        if (random_start_) then
            if (present(offset)) then
                call initialise_ranlux(offset)
            else
                call initialise_ranlux()
            endif
        else
            if (present(offset)) then
                call rluxgo(3, 314159265+offset, 0, 0)
            else
                call rluxgo(3, 314159265, 0, 0)
            endif
        endif

        ! Define sweep sequence used in Monte Carlo sweeps
        sweep_sequence = checkerboard()

        ! Load gauge field
        call load_gauge_field(gauge_field)
    end subroutine


    ! ------------------------------------- Monte Carlo step functions ------------------------------------- !

    ! Find staples to a given link (n, mu)
    function find_staples(n, mu, links)
        implicit none

        integer(int32), intent(in) :: n, mu
        complex(real32), dimension(NCOL, NCOL, LATTICE_VOLUME, 3), intent(in) :: links

        complex(real32), dimension(NCOL, NCOL, 4) :: find_staples

        integer(int32) :: i, j, ni, nj, nk, n_i, n_j

        select case (mu)
        case(1)
            i = 2; j = 3
        case(2)
            i = 3; j = 1
        case(3)
            i = 1; j = 2
        case default
            error stop "Error: Invalid mu in action_change"
        end select

        ni = move(n,i,1)
        nj = move(n,j,1)
        nk = move(n,mu,1)
        n_i = move(n,i,-1)
        n_j = move(n,j,-1)

        ! Find staples to link (n,mu)
        find_staples(:,:,1) = matmul(links(:,:,nk,i), matmul(hermc(links(:,:,ni,mu)), hermc(links(:,:,n,i))))
        find_staples(:,:,2) = matmul(links(:,:,nk,j), matmul(hermc(links(:,:,nj,mu)), hermc(links(:,:,n,j))))
        find_staples(:,:,3) = matmul(hermc(links(:,:,move(n_i,mu,1),i)), matmul(hermc(links(:,:,n_i,mu)),&
        links(:,:,n_i,i)))
        find_staples(:,:,4) = matmul(hermc(links(:,:,move(n_j,mu,1),j)), matmul(hermc(links(:,:,n_j,mu)),&
        links(:,:,n_j,j)))
    end function

    ! Perform heatbath algorithm, given 2x2 input matrix w (W = U*Sigma)
    function heatbath_subgroup(w) result(a)
        implicit none

        complex(real32), dimension(2,2), intent(in) :: w
        complex(real32), dimension(2,2) :: a

        real(real32), dimension(0:3) :: chi_pauli
        complex(real32), dimension(2,2) :: w_hat, chi
        real(real32) :: w_hat_normalisation, rho, trial_sample, gauss_sample, exp_sample
        real(real32), dimension(6) :: random_numbers
        real(real32), dimension(3) :: sphere_vector

        ! Make an SU(2) variable out of w. Store this in w_hat.
        w_hat(1,1) = w(1,1) + conjg(w(2,2))
        w_hat(2,1) = w(2,1) - conjg(w(1,2))
        w_hat(1,2) = -conjg(w_hat(2,1))
        w_hat(2,2) = conjg(w_hat(1,1))
        w_hat_normalisation = sqrt((abs(w_hat(1,1)))**2 + (abs(w_hat(1,2)))**2)
        w_hat = w_hat/w_hat_normalisation
        w_hat_normalisation = w_hat_normalisation/2.0

        ! Define rho
        rho = 2.0*beta*w_hat_normalisation/NCOL

        ! Take 6 random numbers
        call ranlux(random_numbers, 6)

        !! Sample from a_0 distribution
        do
            ! Sample from Gaussian distribution
            gauss_sample = sqrt(-log(1.0-random_numbers(1)))*cos(real(PI, kind=real32)*random_numbers(2)/2.0)

            ! Sample from exponential ditribution
            exp_sample = -log(1.0-random_numbers(3))

            ! Construct trial sample
            trial_sample = gauss_sample**2 + exp_sample

            ! Test if trial sample is accepted
            if (2.0*random_numbers(4)**2 <= 2.0-trial_sample/rho) exit

            !! If the program has progressed here, the trial was rejected.
            ! Re-randomize the 3rd and 4th random numbers
            call ranlux(random_numbers(3:4), 2)

            ! Re-sample from Gaussian distribution
            gauss_sample = sqrt(-log(1.0-random_numbers(1)))*sin(real(PI, kind=real32)*random_numbers(2)/2.0)

            ! Re-sample from exponential ditribution
            exp_sample = -log(1.0-random_numbers(3))

            ! Reconstruct trial sample
            trial_sample = gauss_sample**2 + exp_sample

            ! Test if trial sample is accepted
            if (2.0*random_numbers(4)**2 <= 2.0-trial_sample/rho) exit

            !! If the program has progessed here, both samples have been rejected.
            ! Re-randomize all 4 random numbers used and start again
            call ranlux(random_numbers(1:4), 4)
        enddo

        !! Sample vector from the surface of a 2-sphere
        sphere_vector(1) = 1.0 - 2.0*random_numbers(5)
        sphere_vector(2) = sqrt(1.0-sphere_vector(1)**2)*cos(2.0*real(PI, kind=real32)*random_numbers(6))
        sphere_vector(3) = sqrt(1.0-sphere_vector(1)**2)*sin(2.0*real(PI, kind=real32)*random_numbers(6))

        ! Undo transformations to find coefficients for v
        chi_pauli(0) = 1.0 - trial_sample/rho
        chi_pauli(1:3) = sqrt(1.0-chi_pauli(0)**2)*sphere_vector
        ! print *, "Chi (pauli rep.): ", chi_pauli

        ! Put v into matrix form using Pauli matrices
        chi(1,1) = cmplx(chi_pauli(0),chi_pauli(3))
        chi(2,1) = cmplx(-chi_pauli(2),chi_pauli(1))
        chi(1,2) = cmplx(chi_pauli(2),chi_pauli(1))
        chi(2,2) = cmplx(chi_pauli(0),-chi_pauli(3))

        ! Update a using a = v*w_hat
        a = matmul(chi, conjg(transpose(w_hat)))
    end function

    ! Perform overrelaxation algorithm, given 2x2 input matrix w (W = U*Sigma)
    function overrelaxation_subgroup(w) result(a)
        implicit none

        complex(real32), dimension(2,2), intent(in) :: w
        complex(real32), dimension(2,2) :: a

        real(real32) :: w_pauli

        ! Extract real part of the Pauli coefficient 0 of w
        w_pauli = w(1,1)%re + w(2,2)%re

        ! Construct a as v in textbook
        a(1,1) = conjg(w(1,1)) + w(2,2)
        a(2,1) = conjg(w(1,2)) - w(2,1)
        a(1,2) = -conjg(a(2,1))
        a(2,2) = conjg(a(1,1))

        ! Correct to form SU(2) variable
        a = 2*w_pauli*a/det_2x2(a)
        a(1,1) = a(1,1) - cmplx(1.0,0.0)
        a(2,2) = a(2,2) - cmplx(1.0,0.0)
    end function

    ! Find staples and iterate SU(2) procedure over SU(2) subgroups of SU(N)
    subroutine monte_carlo_step(n, mu, links, algorithm)
        implicit none

        interface
            function algorithm(w) result(a)
                use iso_fortran_env, only: real32
                use luxury
                implicit none
                complex(real32), dimension(2,2), intent(in) :: w
                complex(real32), dimension(2,2) :: a
            end function
        end interface

        integer(int32), intent(in) :: n, mu
        complex(real32), dimension(NCOL, NCOL, LATTICE_VOLUME, 3), intent(inout) :: links

        complex(real32), dimension(NCOL, NCOL, 4) :: staples
        complex(real32), dimension(NCOL, NCOL) :: Sigma, U
        complex(real32), dimension(2,2) :: w_SU2, a_SU2
        complex(real32), dimension(2, NCOL) :: U_prime
        integer(int32) :: i, j, k

        ! Find staples to link (n, mu)
        staples = find_staples(n, mu, links)

        ! Sum staples and multiply by U(n,mu) to find W
        U = links(:,:,n,mu)
        Sigma = cmplx(0.0,0.0)
        do k = 1, 4
            Sigma = Sigma + staples(:,:,k)
        enddo

        !! Iterate through SU(2) subgroups
        do j = 1, NCOL
            do i = 1, j-1
                ! (i,j) give indices of SU(2) subgroups in SU(N) matrices
                ! elements are located at (i,i), (i,j), (j,i) and (j,j)

                ! Update W and real(PI, kind=real32)ck out SU(2) subgroup
                w_SU2 = cmplx(0.0,0.0)
                do k = 1, NCOL
                    w_SU2(1,1) = w_SU2(1,1) + U(i,k)*Sigma(k,i)
                    w_SU2(2,1) = w_SU2(2,1) + U(j,k)*Sigma(k,i)
                    w_SU2(1,2) = w_SU2(1,2) + U(i,k)*Sigma(k,j)
                    w_SU2(2,2) = w_SU2(2,2) + U(j,k)*Sigma(k,j)
                enddo

                ! Carry out Monte Carlo procedure on this SU(2) subgroup
                a_SU2 = algorithm(w_SU2)

                ! Find new link matrix using U' = AU
                U_prime(1,:) = a_SU2(1,1)*U(i,:) + a_SU2(1,2)*U(j,:)
                U_prime(2,:) = a_SU2(2,1)*U(i,:) + a_SU2(2,2)*U(j,:)
                U(i,:) = U_prime(1,:)
                U(j,:) = U_prime(2,:)
            enddo
        enddo

        links(:,:,n,mu) = U
    end subroutine

    ! Conduct a full sweep of the lattice (paralelliseable)
    subroutine sweep(links, sequence, algorithm)
        implicit none

        interface
            function algorithm(w) result(a)
                use iso_fortran_env, only: real32
                use luxury
                implicit none
                complex(real32), dimension(2,2), intent(in) :: w
                complex(real32), dimension(2,2) :: a
            end function
        end interface

        integer(int32), dimension(LATTICE_VOLUME), intent(in) :: sequence
        complex(real32), dimension(NCOL, NCOL, LATTICE_VOLUME, 3), intent(inout) :: links

        integer(int32) :: i, direction

        ! Sweep over lattice according to sequence defined
        do direction = 1, 3
            do i = 1, LATTICE_VOLUME
                call monte_carlo_step(sequence(i), direction, links, algorithm)
            enddo
        enddo
    end subroutine

    ! Conduct 1 heatbath sweep followed by several overrelaxation sweeps
    subroutine compound_sweep(links, sequence)
        implicit none
        integer(int32), dimension(LATTICE_VOLUME), intent(in) :: sequence
        complex(real32), dimension(NCOL, NCOL, LATTICE_VOLUME, 3), intent(inout) :: links

        integer(int32) :: i

        ! Conduct heatbath sweep
        call sweep(links, sequence, heatbath_subgroup)

        ! Conduct L/2a overrelaxation sweeps
        do i = 1, 3
            call sweep(links, sequence, overrelaxation_subgroup)
        enddo
    end subroutine

    ! Renormalise 2x2 matrix to lie in SU(2)
    subroutine renormalise_SU2(w)
        implicit none
        complex(real32), dimension(2,2), intent(out) :: w

        complex(real32), dimension(2,2) :: w_hat
        real(real32) :: w_hat_normalisation

        ! Make an SU(2) variable out of w
        w_hat(1,1) = w(1,1) + conjg(w(2,2))
        w_hat(2,1) = w(2,1) - conjg(w(1,2))
        w_hat(1,2) = -conjg(w_hat(2,1))
        w_hat(2,2) = conjg(w_hat(1,1))
        w_hat_normalisation = sqrt((abs(w_hat(1,1)))**2 + (abs(w_hat(1,2)))**2)
        w = w_hat/w_hat_normalisation
    end subroutine

    ! Renormalise links to remain in SU(N) using Gram-Schmidt
    subroutine renormalise_links_Gram_Schmidt(links)
        complex(real32), dimension(NCOL, NCOL, LATTICE_VOLUME, 3), intent(out) :: links

        complex(real32), dimension(NCOL, NCOL) :: U
        complex(real32), dimension(NCOL) :: orthogonal_vector, input_vector
        complex(real32) :: det
        integer(int32) :: n, mu, i, j

        do mu = 1, 3
            do n = 1, LATTICE_VOLUME
                !! Use Gram-Schmidt to orthonormalise columns
                ! Make columns of U orthogonal
                do j = 1, NCOL
                    input_vector = links(:,j,n,mu)
                    orthogonal_vector = input_vector

                    ! Subtract projections of other orthogonal vectors (columns of U) from vector
                    do i = 1, j-1
                        ! orthogonal_vector = orthogonal_vector - &
                        ! (dot_product(input_vector,U(:,i))/dot_product(U(:,i),U(:,i)))*U(:,i)
                        orthogonal_vector = orthogonal_vector - dot_product(U(:,i), orthogonal_vector)*U(:,i)
                    enddo

                    ! Set column of U as new orthogonal vector
                    U(:,j) = orthogonal_vector/sqrt(real(dot_product(orthogonal_vector, orthogonal_vector)))
                enddo

                ! Orthonormalise columns by multiplying U by phase
                det = determinant(U)
                det = det/abs(det) ! ensure det(U) is a phase as expected
                ! Scale U to force det(U) = 1 whilst maintaining unitarity
                U = U/(det**(1.0/real(NCOL)))

                ! Replace U in links matrix
                links(:,:,n,mu) = U
            enddo
        enddo
    end subroutine
end module