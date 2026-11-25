module matrix_glueball_3d_parameters
    use iso_fortran_env, only : int8, int32, int64, real32, real64
    implicit none

    !! Define pi
    real(real64), parameter :: PI = 4.0_real64 * atan(1.0_real64)

    !! Declare input parameters
    ! Lattice parameters
    integer(int32) :: NCOL, LS, LT
    real(real32) :: BETA

    ! Monte Carlo parameters
    integer(int32) :: CONFIG_STEP ! Minimum number of configurations between measurements
    character(len=256) :: FILEPATH, FILENAME ! Configuration file path and file name
    integer :: CONFIG_START ! Index of first configuration to be read

    ! Matrix method parameters
    real(real32) :: BARE_MASS ! Bare mass of matter field
    real(real32) :: TOL_CG ! Tolerance on Conjugate Gradient algorithm. Set between 1e-2 - 1e-6.
    real(real32) :: TOL_SVD ! Tolerance on SVD-unitarisation algorithm. Set to between 1e-6 - 1e-10.

    ! Measurement parameters
    integer(int32) :: MEASURE_PER_BIN ! Number of measurements per jackknife bin
    integer :: NUM_BINS ! Number of jackknife bins for this calculation

    !! Declare dependant parameters
    ! Lattice parameters
    integer(int32) :: SLICE_VOLUME ! Number of sites in one slice, renamed from LSIZEB
    integer(int32) :: LATTICE_VOLUME ! Total number of sites, renamed from LSIZE
    integer(int32) :: NCOL2 ! Square number of colours
    integer(int32) :: MAX_DELTA_T ! Maximum time extent of correlation functions, renamed from MAXDTLS

    ! Monte Carlo parameters
    integer(int32) :: NCONFIG ! Total number of configurations generated during simulation

    ! Measurement parameters
    integer(int32) :: NMEASURE ! Total number of measurements to make during calculation

    contains

    ! Read all parameters from file and set dependant parameters
    subroutine initialise_parameters(parameter_filename)
        implicit none
        character(len=*), intent(in) :: parameter_filename

        ! Define variables to read in from parameter file
        namelist /params/ FILEPATH, FILENAME, NCOL, LS, LT, BETA, CONFIG_STEP, &
        CONFIG_START, BARE_MASS, TOL_CG, TOL_SVD, MEASURE_PER_BIN, NUM_BINS

        ! Read parameters from parameter file
        open(10, file=trim(parameter_filename))
        read(10, params)
        close(10)

        ! Set dependant parameters
        SLICE_VOLUME = LS ** 2
        LATTICE_VOLUME = SLICE_VOLUME * LT
        NCOL2 = NCOL * NCOL
        MAX_DELTA_T = LT / 2
        NMEASURE = MEASURE_PER_BIN * NUM_BINS
    end subroutine initialise_parameters
end module matrix_glueball_3d_parameters