program thermalise
    use matrix_glueball_3d_parameters
    use lattice, only: setup_lattice
    use monte_carlo, only: &
        setup_monte_carlo, &
        save_gauge_field, &
        compound_sweep, &
        renormalise_links_Gram_Schmidt, &
        average_plaquette

    implicit none

    complex(real32), allocatable :: gauge_field(:,:,:,:)
    integer(int32), allocatable :: sequence(:)
    integer(int32) :: config, renorm
    complex(real64) :: plaquette

    ! --------------------------------------- Thermalise gauge field --------------------------------------- !
    !
    ! Read gauge field with identitfier CONFIG_START (n=0 corresponds to the frozen configuration)
    ! Thermalise gauge field for NMEASURE * CONFIG_STEP Monte Carlo sweeps
    ! Renormalise links every CONFIG_STEP configurations
    ! Save gauge field to file with next available identifier
    !
    ! ------------------------------------------------------------------------------------------------------ !

    ! Load all parameters
    call initialise_parameters("parameter_file.txt")

    ! Define number of configurations to thermalise over
    NCONFIG = NMEASURE * CONFIG_STEP

    ! Setup monte carlo data
    allocate(gauge_field(NCOL, NCOL, LATTICE_VOLUME, 3), sequence(LATTICE_VOLUME))
    call setup_lattice()
    call setup_monte_carlo(gauge_field, sequence, random_start=.true.)

    ! Thermalise gauge field
    do renorm = 1, NMEASURE
        ! Evolve CONFIG_STEP configurations
        do config = 1, CONFIG_STEP
            call compound_sweep(gauge_field, sequence)
        enddo

        ! Renormalise links
        call renormalise_links_Gram_Schmidt(gauge_field)
    enddo

    ! Save gauge field
    call save_gauge_field(gauge_field)

    ! Output average plaquette
    plaquette = average_plaquette(gauge_field)
    if (plaquette%im >= 0) then
        write(*, '(a, f8.6, a, f8.6, a)') "Average plaquette: ", plaquette%re, "+", plaquette%im, "i"
    else
        write(*, '(a, f8.6, a, f8.6, a)') "Average plaquette: ", plaquette%re, "-", abs(plaquette%im), "i"
    endif
end program thermalise