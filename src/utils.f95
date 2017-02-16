!----------------------------------------------------------------------------------------
!
! MODULE: Utils
!
! AUTHOR: Avraham Adler <Avraham.Adler@gmail.com>
!
! DESCRIPTION: Utility functions and definitions for Delaporte package
!
! HISTORY:
!          Version 1.0: 2016-11-20
!
! LICENSE:
!   Copyright (c) 2016, Avraham Adler
!   All rights reserved.
!
!   Redistribution and use in source and binary forms, with or without modification, are
!   permitted provided that the following conditions are met:
!       1. Redistributions of source code must retain the above copyright notice, this
!          list of conditions and the following disclaimer.
!       2. Redistributions in binary form must reproduce the above copyright notice,
!          this list of conditions and the following disclaimer in the documentation
!          and/or other materials provided with the distribution.
!
!   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
!   EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
!   OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
!   SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
!   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
!   TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
!   BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
!   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
!   DAMAGE.
!----------------------------------------------------------------------------------------

module utils
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none

    real(kind = c_double), parameter :: ZERO = 0._c_double
    real(kind = c_double), parameter :: HALF = 0.5_c_double
    real(kind = c_double), parameter :: ONE = 1._c_double
    real(kind = c_double), parameter :: THREEHALFS = 1.5_c_double
    real(kind = c_double), parameter :: TWO = 2._c_double
    real(kind = c_double), parameter :: THREE = 3._c_double
    real(kind = c_double), parameter :: EPS = 2.2204460492503131e-16_c_double
    real(kind = c_double), parameter :: NAN = TRANSFER(z'7FF0000000000001', ONE)
    real(kind = c_double), parameter :: INFTY = TRANSFER(z'7FF0000000000000', ONE)
  
contains

!----------------------------------------------------------------------------------------
! FUNCTION: log1p
!
! DESCRIPTION: Fortran (2003/8) does not have 1og1p as an intrinsic. This serves as such.
!----------------------------------------------------------------------------------------


    elemental function log1p(x) result(y)

        real(kind = c_double), intent(in) :: x
        real(kind = c_double) :: y, z


        z = x + ONE
        y = log(z) - ((z - ONE) - x) / z   !Serves to eliminate catastrophic subtraction

    end function log1p
  
!----------------------------------------------------------------------------------------
! FUNCTION: Position
!
! DESCRIPTION: Returns the position in vector a of greatest number less than x. Used in
!              singleton versions of pdelap and qdelap.
!----------------------------------------------------------------------------------------  

    pure function position(x, a) result(k)
    
        real(kind = c_double), intent(in)                :: x
        real(kind = c_double), intent(in), dimension(:)  :: a
        integer                                          :: k

        k = 1
        do while (a(k) < x)
            k = k + 1
        end do
        
    end function position

end module utils
