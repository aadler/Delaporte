!-------------------------------------------------------------------------------
!
! MODULE: Utils
!
! AUTHOR: Avraham Adler <Avraham.Adler@gmail.com>
!
! DESCRIPTION: Utility functions and definitions for Delaporte package
!
! HISTORY:
!          Version 1.0: 2016-11-20
!          Version 1.1: 2017-03-01
!          Version 1.2: 2017-11-20
!                       Reformatted to 80 columns
!          Version 1.3: 2018-11-20
!                       Added cleanzeros function to handle EPS issues for right
!                       tail. See Issue #1
!          Version 1.4: 2018-12-10
!                       Cleanzeros removed in favor of floor & ceiling of 0 & 1
!          Version 1.5: 2021-01-03
!                       Added parameter for pdelap max vector size
!          Version 2.0: 2023-01-29
!                       Updated to rely on Fortran 2008 intrinsics
!          Version 2.1: 2023-08-08
!                       Added OpenMP control functions
!          Version 3.0: 2023-09-28
!                       Converted log1p to one based on its Taylor expansion.
!                       This is also the degree 2 polynomial minimax
!                       approximation.
!          Version 3.1: 2024-05-21
!                       Added pure header to log1p.
!          Version 4.0: 2024-06-17
!                       Added imk helper function. A smidgen faster—I'm not sure
!                       why, perhaps due to pre-compilation in module—and easier
!                       to read. Turn FP error cleaning into a function.
!          Version 5.0: 2026-07-07
!                       Removed sOMPT_f as that was a global setting and could
!                       interfere with other OMP packages. Threads are now
!                       handled in a package-specific environment. See zzz.R and
!                       omp.R for more.
!                       Changed binding names for header/source refactor.
!                       Use specific "only" lists to prevent scope infractions.
!                       Added interface to C unifrnd and drop need for
!                       external and F77_SUB calls.
!
! LICENSE:
!   Copyright (c) 2016, Avraham Adler
!   All rights reserved.
!
!   Redistribution and use in source and binary forms, with or without
!   modification, are permitted provided that the following conditions are met:
!       1. Redistributions of source code must retain the above copyright
!          notice, this list of conditions and the following disclaimer.
!       2. Redistributions in binary form must reproduce the above copyright
!          notice, this list of conditions and the following disclaimer in the
!          documentation and/or other materials provided with the distribution.
!
!   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
!   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
!   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!   POSSIBILITY OF SUCH DAMAGE.
!-------------------------------------------------------------------------------

module utils
    use, intrinsic :: iso_c_binding,   only: c_int, c_double
    use, intrinsic :: iso_fortran_env, only: INT64 
    !$ use omp_lib
    implicit none

    real(kind = c_double), parameter :: ZERO = 0._c_double
    real(kind = c_double), parameter :: HALF = 0.5_c_double
    real(kind = c_double), parameter :: ONE = 1._c_double
    real(kind = c_double), parameter :: THREEHALFS = 1.5_c_double
    real(kind = c_double), parameter :: TWO = 2._c_double
    real(kind = c_double), parameter :: THREE = 3._c_double
    real(kind = c_double), parameter :: EPS = 2.2204460492503131e-16_c_double
    real(kind = c_double), parameter :: MAXD = REAL(HUGE(1_INT64), c_double)
    integer, parameter               :: MAXVECSIZE = 16384
    
! ------------------------------------------------------------------------------
! Interface to C-side RNG bridge (defined in utils_and_wrappers.c). Declared
! here at module scope so any procedure that uses this module gets a
! compiler-checked, explicit interface -- replaces the implicit "external
! unifrnd" declaration formerly local to rdelap_f. bind(C) pins the linker
! symbol to literally "unifrnd", so the C definition needs no F77_SUB/mangling
! macro.
! NOTE: interface bodies are their own scoping unit and do NOT inherit the
! module's use statements, so iso_c_binding must be re-imported inside the body
! (or brought in via an IMPORT statement).
! ---------------------------------------------------------------------
    interface
        subroutine unifrnd(n, x) bind(C, name = "unifrnd")
            use, intrinsic :: iso_c_binding, only: c_int, c_double
            
            integer(kind = c_int), intent(in), value :: n
            real(kind = c_double), intent(out)       :: x(n)
        end subroutine unifrnd
    end interface

contains

!-------------------------------------------------------------------------------
! FUNCTION: log1p
!
! DESCRIPTION: Fortran 2008 does not have log1p as an intrinsic. This uses the
!              Taylor expansion for small x to reduce relative error.
!-------------------------------------------------------------------------------

    pure elemental function log1p(x) result(y)

        real(kind = c_double), intent(in) :: x
        real(kind = c_double)             :: y

        if (abs(x) <= 1.e-4_c_double) then
            y = (-x * HALF + ONE) * x
        else
            y = log(x + ONE)
        end if
            
    end function log1p
    
!-------------------------------------------------------------------------------
! FUNCTION: imk (i mod k)
!
! DESCRIPTION: Calculates mod(i - 1, k) + 1 for vector recyling.
!-------------------------------------------------------------------------------

    pure elemental function imk(i, k) result(j)

    integer(kind = c_int), intent(in) :: i, k
    integer(kind = c_int)             :: j
    
        j = mod(i - 1, k) + 1
    
    end function imk
    
!-------------------------------------------------------------------------------
! FUNCTION: cFPe (clearFPerrors)
!
! DESCRIPTION: Restricts solutions to [0, 1] and eliminates spurious FP errors.
!-------------------------------------------------------------------------------

    pure elemental function cFPe(x) result(y)

    real(kind = c_double), intent(in) :: x
    real(kind = c_double)             :: y
    
        y = max(min(x, ONE), ZERO)
    
    end function cFPe
    
!-------------------------------------------------------------------------------
! FUNCTION: gOMPT
!
! DESCRIPTION: Gets the OpenMP runtime's maximum thread count. Called once at
!              package load to seed the package-local thread setting; the
!              per-call thread count is passed explicitly to each parallel
!              region via its num_threads clause. The former companion
!              sOMPT_f (omp_set_num_threads) was removed in 9.0.0 because it
!              mutated process-global OpenMP state shared with other packages.
!-------------------------------------------------------------------------------

    subroutine gOMPT_f(n) bind(C, name="gOMPT_f")
    
    integer(kind = c_int), intent(out) :: n
    
        n = 1_c_int
        !$ n = omp_get_max_threads()
    
    end subroutine gOMPT_f    
    
end module utils ! # nocov covr doesn't always pick up the end module
