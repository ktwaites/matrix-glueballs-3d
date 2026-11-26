program convert_csv
    use matrix_glueball_3d_parameters
    use monte_carlo, only: save_gauge_field

    implicit none

    complex(real32), allocatable :: gauge_field(:,:,:,:)
    integer(int32) :: run_id
    logical :: run_passed

    ! Define output filepath and filename
    FILEPATH = "/mnt/users/twaitesk/gauge_fields/"
    FILENAME = "gauge_field_"

    ! Define where to start enumerating field configurations from
    CONFIG_START = 1

    ! Get parameters from command line
    call get_command_line_parameters(run_id, run_passed)

    ! Allocate gauge field
    allocate(gauge_field(NCOL, NCOL, LATTICE_VOLUME, 3))

    ! Read gauge field from csv
    if (run_passed) then
        gauge_field = load_field_config(run_id)
    else
        gauge_field = load_field_config()
    endif

    ! Save gauge field to file
    if (run_passed) then
        call save_gauge_field(gauge_field, gauge_field_index=run_id)
    else
        call save_gauge_field(gauge_field)
    endif

    contains

    ! Get NCOL, BETA, LS & LT from the command line
    subroutine get_command_line_parameters(run_id, run_passed)
        implicit none
        integer(int32), intent(out) :: run_id
        logical, intent(out) :: run_passed

        character(len=128) :: arg
        logical :: ncol_passed, beta_passed, ls_passed, lt_passed
        integer(int32) :: i, nargs

        ! Get number of arguments
        nargs = command_argument_count()

        ! Initialise logical parameters
        ncol_passed = .false.
        beta_passed = .false.
        ls_passed = .false.
        lt_passed = .false.
        run_passed = .false.
        i = 1

        ! Get options
        do while (i <= nargs)
            call get_command_argument(i, arg)

            select case (trim(arg))
            case("--ncol")
                i = i + 1
                call get_command_argument(i, arg)
                read(arg, *) NCOL
                ncol_passed = .true.
            case("--beta")
                i = i + 1
                call get_command_argument(i, arg)
                read(arg, *) BETA
                beta_passed = .true.
            case("--ls")
                i = i + 1
                call get_command_argument(i, arg)
                read(arg, *) LS
                ls_passed = .true.
            case("--lt")
                i = i + 1
                call get_command_argument(i, arg)
                read(arg, *) LT
                lt_passed = .true.
            case("--run")
                i = i + 1
                call get_command_argument(i, arg)
                read(arg, *) run_id
                run_passed = .true.
            case default
                write(*, '(a, 4x, a)') "[Error]", "Unknown option"
            end select

            i = i + 1
        enddo

        ! Check all options were passed
        if (.not. run_passed) run_id = 0
        if (all([ncol_passed, beta_passed, ls_passed, lt_passed])) then
            SLICE_VOLUME = LS ** 2
            LATTICE_VOLUME = SLICE_VOLUME * LT
            return
        else
            write(*, '(a, 4x, a)') "[Error]", &
            "Insufficient options passed"
            write(*, '(a, a)') "Run with ", &
            "./convert_csv --ncol <NCOL> --beta <BETA> --ls <LS> --lt <LT> --run <run_id>"
            error stop
        endif
    end subroutine

    ! Load field configuration from csv file
    function load_field_config(run_id) result(U)
        implicit none

        integer, optional, intent(in) :: run_id
        complex, dimension(NCOL,NCOL,LATTICE_VOLUME,3) :: U

        character(len=128) :: file_name, file_path
        character(len=10) :: file_beta, file_N_s, file_N_t, file_N_col, file_id
        real, dimension(2*NCOL**2, LATTICE_VOLUME*3) :: input
        integer :: pos_x, pos_y, i, j, direction, site

        ! Construct file name
        write(file_beta, '(f4.1)') BETA
        write(file_N_s, '(i0)') LS
        write(file_N_t, '(i0)') LT
        write(file_N_col, '(i0)') NCOL
        if (present(run_id)) then
            write(file_id, '(i0)') run_id
            file_name = "field_config_SU" // trim(adjustl(file_N_col)) // "_" // trim(adjustl(file_beta)) // "_"&
            // trim(adjustl(file_N_s)) // "x" // trim(adjustl(file_N_s)) // "x" // trim(adjustl(file_N_t)) //&
            "_" // trim(adjustl(file_id)) // ".csv"
        else
            file_name = "field_config_SU" // trim(adjustl(file_N_col)) // "_" // trim(adjustl(file_beta)) // "_"&
            // trim(adjustl(file_N_s)) // "x" // trim(adjustl(file_N_s)) // "x" // trim(adjustl(file_N_t)) //&
            ".csv"
        endif

        ! Construct file path
        file_path = "/mnt/users/twaitesk/SUN/field_configs/"

        ! Read in field configuration from file
        ! open(10, file=trim(adjustl(file_name)), status='old', action='read')
        open(10, file=trim(adjustl(file_path)) // trim(adjustl(file_name)), status='old', action='read')
        do i = 1, LATTICE_VOLUME*3
            read(10, *) input(:,i)
        enddo
        close(10)

        ! Reshape input into U
        pos_y = 1
        do direction = 1, 3
            do site = 1, LATTICE_VOLUME
                pos_x = 1
                ! Insert matrix U(n,mu) in real form into output
                do i = 1, NCOL
                    do j = 1, NCOL
                        U(i,j,site,direction)%re = input(pos_x, pos_y)
                        U(i,j,site,direction)%im = input(pos_x+1, pos_y)
                        pos_x = pos_x + 2
                    enddo
                enddo
                pos_y = pos_y + 1
            enddo
        enddo
    end function
end program
