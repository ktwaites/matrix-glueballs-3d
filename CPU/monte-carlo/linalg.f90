module linalg
    use iso_fortran_env, only : real32, real64
    implicit none

    interface eigenpairs
        module procedure eigenpairs_real_symmetric
    end interface

    interface generalised_eigenpairs
        module procedure generalised_eigenpairs_real_symmetric
        module procedure generalised_eigenpairs_complex_Hermitian
    end interface

    interface solve_linear_system
        module procedure solve_linear_system_real32
        module procedure solve_linear_system_complex32
        module procedure solve_linear_system_real64
        module procedure solve_linear_system_complex64
    end interface

    interface determinant
        module procedure determinant_real32
        module procedure determinant_complex32
        module procedure determinant_real64
        module procedure determinant_complex64
    end interface

    interface hermc
        module procedure hermc_complex32
        module procedure hermc_complex64
    end interface

    contains

    ! Output hermitian conjugate of matrix
    function hermc_complex32(matrix) result(hermc)
        implicit none
        complex(real32), intent(in) :: matrix(:,:)
        complex(real32) :: hermc(size(matrix, dim=1), size(matrix, dim=2))

        hermc = conjg(transpose(matrix))
    end function

    function hermc_complex64(matrix) result(hermc)
        implicit none
        complex(real64), intent(in) :: matrix(:,:)
        complex(real64) :: hermc(size(matrix, dim=1), size(matrix, dim=2))

        hermc = conjg(transpose(matrix))
    end function

    ! Solve a general real system Ax=b using LAPACK routine SGESV.
    function solve_linear_system_real32(A, b) result(x)
        implicit none

        real(real32), intent(in) :: A(:,:), b(:)
        real(real32) :: x(size(b))

        real(real32), allocatable :: Acopy(:,:)
        integer, allocatable :: ipiv(:)
        integer :: n, info

        n = size(A, 1)

        if (n < 1 .or. size(A, 2) /= n) &
        error stop "[Error][solve_linear_system]    A must be a non-empty square matrix"

        if (size(b) /= n) &
        error stop "[Error][solve_linear_system]    b must have size n"

        allocate(Acopy(n,n), ipiv(n))
        Acopy = A
        x = b

        ! SGESV overwrites A and b, so solve using local copies.
        call SGESV(n, 1, Acopy, n, ipiv, x, n, info)

        if (info /= 0) then
            if (info < 0) then
                write(*,'(a,i0)') '[Error][solve_linear_system]    SGESV: illegal argument ', -info
            else
                write(*,'(a,i0,a,i0,a)') '[Error][solve_linear_system]    SGESV: U(', info, ',', info, ') is zero'
            endif
            error stop
        endif
    end function solve_linear_system_real32

    ! Solve a general complex system Ax=b using LAPACK routine CGESV.
    function solve_linear_system_complex32(A, b) result(x)
        implicit none

        complex(real32), intent(in) :: A(:,:), b(:)
        complex(real32) :: x(size(b))

        complex(real32), allocatable :: Acopy(:,:)
        integer, allocatable :: ipiv(:)
        integer :: n, info

        n = size(A, 1)

        if (n < 1 .or. size(A, 2) /= n) &
        error stop "[Error][solve_linear_system]    A must be a non-empty square matrix"

        if (size(b) /= n) &
        error stop "[Error][solve_linear_system]    b must have size n"

        allocate(Acopy(n,n), ipiv(n))
        Acopy = A
        x = b

        ! CGESV overwrites A and b, so solve using local copies.
        call CGESV(n, 1, Acopy, n, ipiv, x, n, info)

        if (info /= 0) then
            if (info < 0) then
                write(*,'(a,i0)') '[Error][solve_linear_system]    CGESV: illegal argument ', -info
            else
                write(*,'(a,i0,a,i0,a)') '[Error][solve_linear_system]    CGESV: U(', info, ',', info, ') is zero'
            endif
            error stop
        endif
    end function solve_linear_system_complex32

    ! Solve a general real system Ax=b using LAPACK routine DGESV.
    function solve_linear_system_real64(A, b) result(x)
        implicit none

        real(real64), intent(in) :: A(:,:), b(:)
        real(real64) :: x(size(b))

        real(real64), allocatable :: Acopy(:,:)
        integer, allocatable :: ipiv(:)
        integer :: n, info

        n = size(A, 1)

        if (n < 1 .or. size(A, 2) /= n) &
        error stop "[Error][solve_linear_system]    A must be a non-empty square matrix"

        if (size(b) /= n) &
        error stop "[Error][solve_linear_system]    b must have size n"

        allocate(Acopy(n,n), ipiv(n))
        Acopy = A
        x = b

        ! DGESV overwrites A and b, so solve using local copies.
        call DGESV(n, 1, Acopy, n, ipiv, x, n, info)

        if (info /= 0) then
            if (info < 0) then
                write(*,'(a,i0)') '[Error][solve_linear_system]    DGESV: illegal argument ', -info
            else
                write(*,'(a,i0,a,i0,a)') '[Error][solve_linear_system]    DGESV: U(', info, ',', info, ') is zero'
            endif
            error stop
        endif
    end function solve_linear_system_real64

    ! Solve a general complex system Ax=b using LAPACK routine ZGESV.
    function solve_linear_system_complex64(A, b) result(x)
        implicit none

        complex(real64), intent(in) :: A(:,:), b(:)
        complex(real64) :: x(size(b))

        complex(real64), allocatable :: Acopy(:,:)
        integer, allocatable :: ipiv(:)
        integer :: n, info

        n = size(A, 1)

        if (n < 1 .or. size(A, 2) /= n) &
        error stop "[Error][solve_linear_system]    A must be a non-empty square matrix"

        if (size(b) /= n) &
        error stop "[Error][solve_linear_system]    b must have size n"

        allocate(Acopy(n,n), ipiv(n))
        Acopy = A
        x = b

        ! ZGESV overwrites A and b, so solve using local copies.
        call ZGESV(n, 1, Acopy, n, ipiv, x, n, info)

        if (info /= 0) then
            if (info < 0) then
                write(*,'(a,i0)') '[Error][solve_linear_system]    ZGESV: illegal argument ', -info
            else
                write(*,'(a,i0,a,i0,a)') '[Error][solve_linear_system]    ZGESV: U(', info, ',', info, ') is zero'
            endif
            error stop
        endif
    end function solve_linear_system_complex64

    ! Calculate a real32 determinant using LAPACK routine SGETRF.
    function determinant_real32(A) result(det)
        implicit none

        real(real32), intent(in) :: A(:,:)
        real(real32) :: det

        real(real32), allocatable :: Acopy(:,:)
        integer, allocatable :: ipiv(:)
        integer :: n, i, info

        n = size(A, 1)
        if (n < 1 .or. size(A, 2) /= n) &
        error stop "[Error][determinant]    A must be a non-empty square matrix"

        allocate(Acopy(n,n), ipiv(n))
        Acopy = A
        call SGETRF(n, n, Acopy, n, ipiv, info)

        if (info < 0) then
            write(*,'(a,i0)') '[Error][determinant]    SGETRF: illegal argument ', -info
            error stop
        endif

        det = 1.0_real32
        do i = 1, n
            det = det * Acopy(i,i)
            if (ipiv(i) /= i) det = -det
        enddo
    end function determinant_real32

    ! Calculate a complex32 determinant using LAPACK routine CGETRF.
    function determinant_complex32(A) result(det)
        implicit none

        complex(real32), intent(in) :: A(:,:)
        complex(real32) :: det

        complex(real32), allocatable :: Acopy(:,:)
        integer, allocatable :: ipiv(:)
        integer :: n, i, info

        n = size(A, 1)
        if (n < 1 .or. size(A, 2) /= n) &
        error stop "[Error][determinant]    A must be a non-empty square matrix"

        allocate(Acopy(n,n), ipiv(n))
        Acopy = A
        call CGETRF(n, n, Acopy, n, ipiv, info)

        if (info < 0) then
            write(*,'(a,i0)') '[Error][determinant]    CGETRF: illegal argument ', -info
            error stop
        endif

        det = cmplx(1.0_real32, 0.0_real32, kind=real32)
        do i = 1, n
            det = det * Acopy(i,i)
            if (ipiv(i) /= i) det = -det
        enddo
    end function determinant_complex32

    ! Calculate a real64 determinant using LAPACK routine DGETRF.
    function determinant_real64(A) result(det)
        implicit none

        real(real64), intent(in) :: A(:,:)
        real(real64) :: det

        real(real64), allocatable :: Acopy(:,:)
        integer, allocatable :: ipiv(:)
        integer :: n, i, info

        n = size(A, 1)
        if (n < 1 .or. size(A, 2) /= n) &
        error stop "[Error][determinant]    A must be a non-empty square matrix"

        allocate(Acopy(n,n), ipiv(n))
        Acopy = A
        call DGETRF(n, n, Acopy, n, ipiv, info)

        if (info < 0) then
            write(*,'(a,i0)') '[Error][determinant]    DGETRF: illegal argument ', -info
            error stop
        endif

        det = 1.0_real64
        do i = 1, n
            det = det * Acopy(i,i)
            if (ipiv(i) /= i) det = -det
        enddo
    end function determinant_real64

    ! Calculate a complex64 determinant using LAPACK routine ZGETRF.
    function determinant_complex64(A) result(det)
        implicit none

        complex(real64), intent(in) :: A(:,:)
        complex(real64) :: det

        complex(real64), allocatable :: Acopy(:,:)
        integer, allocatable :: ipiv(:)
        integer :: n, i, info

        n = size(A, 1)
        if (n < 1 .or. size(A, 2) /= n) &
        error stop "[Error][determinant]    A must be a non-empty square matrix"

        allocate(Acopy(n,n), ipiv(n))
        Acopy = A
        call ZGETRF(n, n, Acopy, n, ipiv, info)

        if (info < 0) then
            write(*,'(a,i0)') '[Error][determinant]    ZGETRF: illegal argument ', -info
            error stop
        endif

        det = cmplx(1.0_real64, 0.0_real64, kind=real64)
        do i = 1, n
            det = det * Acopy(i,i)
            if (ipiv(i) /= i) det = -det
        enddo
    end function determinant_complex64

    ! Return a list of evenly spaced real numbers in the interval [start, stop]
    function linspace(start, stop, num)
        implicit none

        real(8), intent(in) :: start, stop
        integer, intent(in) :: num

        real(8) :: linspace(num)

        real(8) :: step
        integer :: i

        step = (stop - start) / real(num - 1, kind=8)
        do i = 1, num
            linspace(i) = start + real(i - 1, kind=8) * step
        enddo
    end function

    ! Solve real symmetric eigenvalue problem using LAPACK routine DSYEV.
    ! All eigenpairs are returned in ascending order of eigenvalue.
    ! The columns of eigenvectors contain the corresponding eigenvectors.
    subroutine eigenpairs_real_symmetric(A, eigenvalues, eigenvectors)
        implicit none

        real(8), intent(in)  :: A(:,:)
        real(8), intent(out) :: eigenvalues(:), eigenvectors(:,:)

        real(8), allocatable :: Acopy(:,:), work(:)
        integer :: lwork, info, n

        ! Get size of matrix
        n = size(A, 1)

        ! Check dimensions
        if (n < 1 .or. size(A, 2) /= n) &
        error stop "[Error][eigenpairs_real_symmetric]    A must be a non-empty square matrix"

        if (size(eigenvalues) /= n) &
        error stop "[Error][eigenpairs_real_symmetric]    eigenvalues must have shape (n)"

        if (size(eigenvectors, 1) /= n .or. size(eigenvectors, 2) /= n) &
        error stop "[Error][eigenpairs_real_symmetric]    eigenvectors must have shape (n,n)"

        ! Copy matrix, since DSYEV overwrites it
        allocate(Acopy(n,n))
        Acopy = A

        ! Workspace query
        allocate(work(1))
        lwork = -1

        call DSYEV('V', 'U', n, Acopy, n, eigenvalues, work, lwork, info)

        if (info /= 0) &
        error stop "[Error][eigenpairs_real_symmetric]    DSYEV workspace query failed"

        ! Allocate optimal workspace
        lwork = max(1, 3*n-1, int(work(1)))
        deallocate(work)
        allocate(work(lwork))

        ! Restore matrix before solving
        Acopy = A

        ! Solve symmetric eigenvalue problem
        call DSYEV('V', 'U', n, Acopy, n, eigenvalues, work, lwork, info)

        if (info /= 0) then
            if (info < 0) then
                write(*,'(a,i0)') &
                    '[Error][eigenpairs_real_symmetric]    DSYEV: illegal argument ', -info
            else
                write(*,'(a,i0)') &
                    '[Error][eigenpairs_real_symmetric]    DSYEV failed to converge, INFO = ', info
            endif
            error stop
        endif

        ! DSYEV returns eigenvectors as columns of Acopy
        eigenvectors = Acopy

    end subroutine eigenpairs_real_symmetric


    ! Solve real symmetric generalised eigenvalue problem using LAPACK routine DSYGVX. Eigenpairs are returned in ascending order (so eigenvalues(k) refers to the largest eigenvalue)
    subroutine generalised_eigenpairs_real_symmetric(A, B, eigenvalues, eigenvectors)
        implicit none

        real(8), intent(in) :: A(:,:), B(:,:)
        real(8), intent(out) :: eigenvalues(:), eigenvectors(:,:)

        real(8), allocatable :: Acopy(:,:), Bcopy(:,:), work(:)
        integer, allocatable :: iwork(:), ifail(:)
        integer :: lwork, info, n, m, il, iu, k
        real(8) :: abstol

        ! Get size of matrices
        n = size(A, 1)

        ! Get requested number of eigenpairs
        k = size(eigenvalues)

        ! Check dimensions
        if (size(A, 2) /= n .or. size(B, 1) /= n .or. size(B, 2) /= n) &
        error stop "[Error][generalised_eigenpairs]    A and B must be square and of equal size"

        if (k < 1 .or. k > n) &
        error stop "[Error][generalised_eigenpairs]    k must be between 1 and n"

        if (size(eigenvectors, 1) /= n .or. size(eigenvectors, 2) /= k) &
        error stop "[Error][generalised_eigenpairs]    eigenvectors must have shape (n,k)"

        ! Allocate arrays
        allocate(Acopy(n,n), Bcopy(n,n))
        allocate(iwork(5*n), ifail(n))

        ! Copy matrices, since DSYGVX overwrites them
        Acopy = A
        Bcopy = B

        ! Select largest k eigenvalues
        il = n - k + 1
        iu = n

        ! Default absolute tolerance
        abstol = 0.0d0

        ! Workspace query
        allocate(work(1))
        lwork = -1

        call DSYGVX(1, 'V', 'I', 'U', n, Acopy, n, Bcopy, n, &
                    0.0d0, 0.0d0, il, iu, abstol, m, &
                    eigenvalues, eigenvectors, n, &
                    work, lwork, iwork, ifail, info)

        if (info /= 0) &
        error stop "[Error][generalised_eigenpairs]    DSYGVX workspace query failed"

        ! Allocate optimal workspace
        lwork = max(8*n, int(work(1)))
        deallocate(work)
        allocate(work(lwork))

        ! Solve generalized eigenvalue problem
        call DSYGVX(1, 'V', 'I', 'U', n, Acopy, n, Bcopy, n, &
                    0.0d0, 0.0d0, il, iu, abstol, m, &
                    eigenvalues, eigenvectors, n, &
                    work, lwork, iwork, ifail, info)

        if (info /= 0) then
            if (info < 0) then
                write(*,'(a,i0)') '[Error][generalised_eigenpairs]    DSYGVX: illegal argument ', -info
            else
                write(*,'(a,i0)') '[Error][generalised_eigenpairs]    DSYGVX failed, INFO = ', info
            endif
            error stop
        endif

        if (m /= k) &
        error stop "[Error][generalised_eigenpairs]    DSYGVX returned an unexpected number of eigenpairs"
    end subroutine

    subroutine generalised_eigenpairs_complex_Hermitian(A, B, eigenvalues, eigenvectors)
        implicit none

        complex(8), intent(in) :: A(:,:), B(:,:)
        real(8), intent(out) :: eigenvalues(:)
        complex(8), intent(out) :: eigenvectors(:,:)

        complex(8), allocatable :: Acopy(:,:), Bcopy(:,:), work(:)
        real(8), allocatable :: rwork(:)
        integer, allocatable :: iwork(:), ifail(:)

        integer :: lwork, info, n, m, il, iu, k
        real(8) :: abstol

        ! Get size of matrices
        n = size(A, 1)

        ! Get requested number of eigenpairs
        k = size(eigenvalues)

        ! Check dimensions
        if (size(A, 2) /= n .or. size(B, 1) /= n .or. size(B, 2) /= n) &
        error stop "[Error][generalised_eigenpairs]    A and B must be square and of equal size"

        if (k < 1 .or. k > n) &
        error stop "[Error][generalised_eigenpairs]    k must be between 1 and n"

        if (size(eigenvectors, 1) /= n .or. size(eigenvectors, 2) /= k) &
        error stop "[Error][generalised_eigenpairs]    eigenvectors must have shape (n,k)"

        ! Allocate arrays
        allocate(Acopy(n,n), Bcopy(n,n))
        allocate(rwork(7*n), iwork(5*n), ifail(n))

        ! Copy matrices, since ZHEGVX overwrites them
        Acopy = A
        Bcopy = B

        ! Select largest k eigenvalues
        il = n - k + 1
        iu = n

        ! Default absolute tolerance
        abstol = 0.0d0

        ! Workspace query
        allocate(work(1))
        lwork = -1

        call ZHEGVX(1, 'V', 'I', 'U', n, Acopy, n, Bcopy, n, &
                    0.0d0, 0.0d0, il, iu, abstol, m, &
                    eigenvalues, eigenvectors, n, &
                    work, lwork, rwork, iwork, ifail, info)

        if (info /= 0) &
        error stop "ZHEGVX workspace query failed"

        ! Allocate optimal workspace
        lwork = max(2*n, int(real(work(1))))
        deallocate(work)
        allocate(work(lwork))

        ! Solve generalized eigenvalue problem
        call ZHEGVX(1, 'V', 'I', 'U', n, Acopy, n, Bcopy, n, &
                    0.0d0, 0.0d0, il, iu, abstol, m, &
                    eigenvalues, eigenvectors, n, &
                    work, lwork, rwork, iwork, ifail, info)

        if (info /= 0) then
            if (info < 0) then
                write(*,'(a,i0)') '[Error][generalised_eigenpairs]    ZHEGVX: illegal argument ', -info
            else if (info <= n) then
                write(*,'(a,i0)') &
                '[Error][generalised_eigenpairs]    ZHEGVX: eigenvectors failed to converge, INFO = ', info
            else
                write(*,'(a,i0)') &
                '[Error][generalised_eigenpairs]    ZHEGVX: B is not positive definite, INFO = ', info
            endif
            error stop
        endif

        if (m /= k) &
        error stop "[Error][generalised_eigenpairs]    ZHEGVX returned an unexpected number of eigenpairs"
    end subroutine
end module
