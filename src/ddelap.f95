!----------------------------------------------------------------------------------------
!
! MODULE: Delaporte
!
! AUTHOR: Avraham Adler <Avraham.Adler@gmail.com>
!
! DESCRIPTION: Probability mass, cumulative distribution, and quantile
!              functions for the Delaporte distribution. Random number
!              generation and method of moments functions as well.
!
! HISTORY:
!          Version 0.1: 2016-11-11
!                       Initial porting from C++ code in Delaporte package
!                       for R.
!----------------------------------------------------------------------------------------
module delaporte
    use, intrinsic :: iso_c_binding
    use utils
    use lgam

    implicit none
    private
    public :: ddelap_f                      !Vectorized functions called from C

contains

!----------------------------------------------------------------------------------------
! FUNCTION: ddelap_f_s
!
! DESCRIPTION: Calculate the Delaporte probability mass function for a single observation
!              and return the value or its log. Calculated through explicit summation.
!              Will coerce real observations to integer by calling ceiling.
!----------------------------------------------------------------------------------------

    function ddelap_f_s (x, alpha, beta, lambda) result (pmf)

    external set_nan                                              ! C-based Nan
    real(kind = c_double)               :: pmf                    ! Result
    real(kind = c_double), intent(in)   :: x, alpha, beta, lambda ! Observation & Parms
!    logical(kind = c_bool), intent(in)  :: lg                     ! Log flag
    integer                             :: i, k                   ! Integers

    if (alpha < EPS .or. beta < EPS .or. lambda < EPS .or. x < 0) then
        call set_nan(pmf)
    else
        k = ceiling(x)
        pmf = 0_c_double
        do i = 0, k, 1
            pmf = pmf + exp(gamln(alpha + i) + i * log(beta) &
                      + (k - i) * log(lambda) - lambda &
                      - gamln(alpha) - gamln(i + ONE) &
                      - (alpha + i) * log1p(beta) &
                      - gamln(k - i + ONE))
        end do
!        if (lg) then
!          pmf = log(pmf)
!        end if
    end if

end function ddelap_f_s

    !------------------------------------------------------------------------------------
    ! ROUTINE: ddelap_f
    !
    ! DESCRIPTION: Vector-based PMF allowing parameter vector recycling and called from C
    !------------------------------------------------------------------------------------

    subroutine ddelap_f (x, nx, a, na, b, nb, l, nl, lg, pmfv) &
                       bind(C, name="ddelap_f")

    integer(kind = c_int), intent(in), value         :: nx, na, nb, nl     ! Sizes
    real(kind = c_double), intent(in), dimension(nx) :: x                  ! Observations
    real(kind = c_double), intent(out), dimension(nx):: pmfv               ! Result
    real(kind = c_double), intent(in)                :: a(na), b(nb), l(nl)! Parameters
    logical(kind = c_bool), intent(in)               :: lg                 ! Log flag
    integer                                          :: i                  ! Integer

        do i = 1, nx
            pmfv(i) = ddelap_f_s(x(i), a(mod(i - 1, na) + 1), b(mod(i - 1, nb) + 1), &
                                 l(mod(i - 1, nl) + 1))
        end do
        if (lg) then
            pmfv = log(pmfv)
        end if

    end subroutine ddelap_f


    !------------------------------------------------------------------------------------
    ! FUNCTION: pdelap_f_s
    !
    ! DESCRIPTION: Calculate the Delaporte probability mass function for a single
    !              observation and return the value or its log. Calculated through
    !              explicit summation. Will coerce real observations to integer
    !              by calling ceiling.
    !------------------------------------------------------------------------------------

    function pdelap_f_s (x, alpha, beta, lambda) result (cdf)
    external set_nan                                              ! C-based Nan
    real(kind = c_double)               :: cdf                    ! Result
    real(kind = c_double), intent(in)   :: x, alpha, beta, lambda ! Observation & Parms
    integer                             :: i, k                   ! Integers

    end function pdelap_f_s

end module delaporte
