module utils
  use, intrinsic :: iso_c_binding
  implicit none

  real(kind = c_double), parameter :: ONE = 1_c_double
  real(kind = c_double), parameter :: ZERO = 0_c_double
  real(kind = c_double), parameter :: EPS = 2.2204460492503131e-16_c_double

  contains

  elemental function log1p(x) result (y)
    real(kind = c_double), intent(in) :: x
    real(kind = c_double) :: y, z


    z = x + ONE
    y = log(z) - ((z - ONE) - x) / z

  end function log1p

  pure function position(x, a) result (k)
      real(kind = c_double), intent(in)                :: x
      real(kind = c_double), intent(in), dimension(:)  :: a
      integer(kind = c_int)                            :: k

      k = 1
      do while (a(k) < x)
          k = k + 1
      end do
  end function position

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
