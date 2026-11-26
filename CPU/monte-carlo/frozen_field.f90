program identity_field
    use matrix_glueball_3d_parameters
    use lattice, only: setup_lattice
    use monte_carlo, only: setup_monte_carlo, save_gauge_field, compound_sweep, renormalise_links_Gram_Schmidt

    implicit none

    complex(real32), allocatable :: gauge_field(:,:,:,:)
    integer(int32) :: n, mu, i

    ! ----------------------------------------- Save identity field ----------------------------------------- !
    !
    ! Save the frozen gauge configuration with parameters given by the input file to a file with n=0
    !
    ! ------------------------------------------------------------------------------------------------------- !

    ! Import parameters
    call initialise_parameters("parameter_file.txt")

    ! Set frozen gauge field configuration
    allocate(gauge_field(NCOL, NCOL, LATTICE_VOLUME, 3))
    do concurrent (n = 1:LATTICE_VOLUME, mu = 1:3, i = 1:NCOL)
        gauge_field(i, i, n, mu) = cmplx(1.0_real32, 0.0_real32, kind=real32)
    enddo

    ! Save gauge field to file with identitfier n=0
    call save_gauge_field(gauge_field, gauge_field_index=0)
end program identity_field