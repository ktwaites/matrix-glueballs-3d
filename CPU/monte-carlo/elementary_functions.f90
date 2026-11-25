module elementary_functions
    use matrix_glueball_3d_parameters

    implicit none

    interface sgn
        module procedure sgn_int32
        module procedure sgn_real32
        module procedure sgn_real64
    end interface

    interface exp_phase
        module procedure exp_phase_real32
        module procedure exp_phase_real64
    end interface

    contains

    ! Calculate error on quantity from jackknife bins
    function jackknife_error(data)
        implicit none
        real(real64), intent(in) :: data(:)
        real(real64) :: jackknife_error

        real(real64) :: jackknife_average
        integer :: nbins

        ! Get number of bins
        nbins = size(data)

        ! Get mean of data
        jackknife_average = sum(data) / nbins

        ! Get jackknife error from data
        jackknife_error = sqrt(real(nbins-1, kind=real64) / real(nbins, kind=real64) &
        * sum((data - jackknife_average) ** 2))
    end function

    ! Compute sign of input
    pure elemental function sgn_int32(x) result(sgn)
        implicit none
        integer(int32), intent(in) :: x
        integer(int32) :: sgn

        if (x > 0) then
            sgn = 1
        elseif (x < 0) then
            sgn = -1
        else
            sgn = 0
        endif
    end function

    pure elemental function sgn_real32(x) result(sgn)
        implicit none
        real(real32), intent(in) :: x
        integer(int32) :: sgn

        if (x > 0) then
            sgn = 1
        elseif (x < 0) then
            sgn = -1
        else
            sgn = 0
        endif
    end function

    pure elemental function sgn_real64(x) result(sgn)
        implicit none
        real(real64), intent(in) :: x
        integer(int32) :: sgn

        if (x > 0) then
            sgn = 1
        elseif (x < 0) then
            sgn = -1
        else
            sgn = 0
        endif
    end function

    ! Display progress bar
    subroutine print_progress_bar(progress, runtime)
        implicit none
        real(real32), intent(in) :: progress
        real(real32), intent(in), optional :: runtime

        ! Variables to print progress bar
        integer(int32), parameter :: width = 50
        character(len=width) :: str = repeat("=", width-1) // ">"
        character(len=width) :: display_str
        integer(int32) :: int_progress, lpad, rpad

        ! Variables to compute time remaining
        real(real32) :: remaining_runtime
        integer(int32) :: hours, minutes, seconds

        ! Compute quantities needed to print progress bar
        int_progress = int(progress*100, kind=int32)
        lpad = int(progress * (width+1), kind=int32)
        rpad = width - lpad
        display_str = str(width-lpad+1:width)
        display_str(lpad+1:width) = " "

        ! Compute time remaining
        if (present(runtime)) then
            remaining_runtime = (1.0/progress - 1)*runtime
            hours = int(remaining_runtime, kind=int32)/3600
            remaining_runtime = remaining_runtime - 3600*hours
            minutes = int(remaining_runtime, kind=int32)/60
            seconds = int(remaining_runtime, kind=int32) - 60*minutes
        endif

        ! Print progress bar and time remaining
        write(*, '(a, a, f5.1, a2, a1, a1)', advance='no') char(13), "Progress: ", progress*100, "%"
        write(*, '(a)', advance='no') "["
        write(*, '(a)', advance='no') display_str
        write(*, '(a)', advance='no') "]"
        if (present(runtime)) then
            write(*, '(a, i0, a, i0, a, i0, a)', advance='no') " Time remaining: ", hours, " hours, ",&
            minutes, " minutes, ", seconds, " seconds"
        endif
        flush(6)
    end subroutine

    ! Compute exponential of an angle
    pure elemental function exp_phase_real32(theta) result(exp_phase)
        implicit none
        real(real32), intent(in) :: theta
        complex(real32) :: exp_phase

        exp_phase = cmplx(cos(theta), sin(theta), kind=real32)
    end function

    pure elemental function exp_phase_real64(theta) result(exp_phase)
        implicit none
        real(real64), intent(in) :: theta
        complex(real64) :: exp_phase

        exp_phase = cmplx(cos(theta), sin(theta), kind=real64)
    end function
end module