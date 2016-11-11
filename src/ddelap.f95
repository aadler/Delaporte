module delaporte_f
  use, intrinsic :: iso_c_binding
  use utils
  use lgam
  implicit none

  real(kind = c_double), parameter :: EPS = 2.2204460492503131e-16_c_double
  real(kind = c_double), parameter :: ONE = 1_c_double
  real(kind = c_double), parameter :: NAN = -ONE

  contains

  function ddelap_f_s (x, alpha, beta, lambda) result (ddp) &
                     bind(C, name="ddelap_f_s")

    external set_nan
    real(kind = c_double)                      :: ddp
    real(kind = c_double), intent(in)          :: x, alpha, beta, lambda
    integer(kind = c_int)                      :: i, k


    if (alpha < EPS .or. beta < EPS .or. lambda < EPS .or. x < 0) then
        call set_nan(ddp)
    else
        k = ceiling(x)
        ddp = 0_c_double
        do i = 0, k, 1
            ddp = ddp + exp(gamln(alpha + i) + i * log(beta) &
                            + (k - i) * log(lambda) - lambda &
                            - gamln(alpha) - gamln(i + ONE) &
                           - (alpha + i) * log1p(beta) &
                            - gamln(k - i + ONE))
        end do
    end if
  end function ddelap_f_s

end module delaporte_f
