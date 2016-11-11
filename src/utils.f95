module utils
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  implicit none

  real(kind = c_double), parameter :: ONE = 1_c_double
  real(kind = c_double), parameter :: EPS = 2.2204460492503131e-16_c_double

  contains

  elemental function log1p(x) result (y)
    real(kind = c_double), intent(in) :: x
    real(kind = c_double) :: y, z


    z = x + ONE
    y = log(z) - ((z - ONE) - x) / z

  end function log1p

  subroutine extend_v(y, nx, ny, y_out)
    integer(kind = c_int), intent(in), value          :: nx, ny
    real(kind = c_double), dimension(ny), intent(in)  :: y
    real(kind = c_double), dimension(nx), intent(out) :: y_out

    if (ny >= nx) then
      y_out = y(1:nx)
    else
      y_out = reshape(y, [nx], y)
    end if
  end subroutine extend_v
end module utils
