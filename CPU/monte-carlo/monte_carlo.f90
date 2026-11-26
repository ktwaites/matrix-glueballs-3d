module monte_carlo
    use matrix_glueball_3d_parameters
    use elementary_functions
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
end module