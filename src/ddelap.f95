module delaporte_f
use, intrinsic :: iso_c_binding
use utils
use lgam
implicit none

contains

function ddelap_f_s (x, alpha, beta, lambda, lg) result (ddp)

external set_nan
real(kind = c_double)                      :: ddp
real(kind = c_double), intent(in)          :: x, alpha, beta, lambda
logical(kind = c_bool), intent(in)         :: lg
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
    if (lg) then
      ddp = log(ddp)
    end if
  end if

end function ddelap_f_s

subroutine ddelap_f (x, nx, a, na, b, nb, l, nl, lg, ddpv) &
                   bind(C, name="ddelap_f")

integer(kind = c_int), intent(in), value          :: nx, na, nb, nl
real(kind = c_double), intent(in), dimension(nx)  :: x
real(kind = c_double), intent(out), dimension(nx) :: ddpv
real(kind = c_double), intent(in), dimension(na)  :: a
real(kind = c_double), intent(in), dimension(nb)  :: b
real(kind = c_double), intent(in), dimension(nl)  :: l
logical(kind = c_bool), intent(in)                :: lg
real(kind = c_double), dimension(nx)              :: ax, bx, lx
integer                                           :: i

  call extend_v(a, nx, na, ax)
  call extend_v(b, nx, nb, bx)
  call extend_v(l, nx, nl, lx)
  do i = 1, nx
    ddpv(i) = ddelap_f_s(x(i), ax(i), bx(i), lx(i), lg)
  end do

end subroutine ddelap_f

end module delaporte_f
