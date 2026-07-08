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
!                       Added "lower_bound" which replaces minloc in qdelap.
!                       This is a binary search, O(log n), and not a linear
!                       scan, O(n).
!                       Changed MAXVECSIZE to 2^24 given ddelap enhancements.
!                       Fixed "bug" in log1p. Cutoff too high. See note.
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
    
    ! alpha*beta above which pdelap_f routes around ddelap_table to the legacy
    ! per-element summation build. See routing comment in pdelap_f.
    real(kind = c_double), parameter :: TBLMAXCOEF = 1.e30_c_double
    
    ! Floor on the guaranteed per-step PMF ratio below which ddelap_f routes
    ! log-space requests around ddelap_table to the elemental path. The
    ! table's downward rescale can rescue a decaying scaled mass only within
    ! the band between the rescale floor 2**(-900) and hard underflow at
    ! 2**(-1074); a single step whose ratio is below that 2**(-174)-wide band
    ! jumps over it, lands on exact zero, and the log-PMF would wrongly
    ! return -Inf (e.g. ddelap(1e4, 1e-50, 1e-100, 1e-50, log = TRUE), whose
    ! true value is about -1.23e6). Because the Delaporte is a convolution of
    ! two log-concave PMFs:
    ! (f⊛g)(n+1) = Σ f(i)·g(n+1−i) ≥ min_ratio(g) · Σ f(i)·g(n−i)
    ! its ratio p(n+1)/p(n) is bounded below by the larger of the components'
    ! minimum ratios: λ/(k+1) for the Poisson factor and
    ! (α+n)/(n+1) · β/(1+β) ≥ min(α,1)·β/(1+β) or min(alpha, 1)*beta/(1+beta)
    ! for the NB factor. Requiring that bound >= 2**(-120) keeps every 
    ! intermediate at or above 2**(-1020)---nonzero AND full-precision
    ! normal---with a 2**54 safety margin over the band. Any statistically
    ! meaningful distribution passes by hundreds of orders of magnitude; only
    ! point-mass-plus-dust degeneracies route.
    ! (AA & Claude: 2026-07-07)
    real(kind = c_double), parameter :: TBLMINRATIO = 2._c_double ** (-120)
    
    ! Maximum allowable q for pdelap's singleton fast path. Raised from 2**14
    ! = 16384 to 2**24 now that the CDF table is built by an O(K) three-term
    ! recurrence; the binding constraint is now the two K+1-length work
    ! vectors (16 bytes per support point, ~270 MB at 2**24), not compute
    ! time.
    ! (AA & Claude: 2026-07-07)
    integer, parameter               :: MAXVECSIZE = 16777216
    
! ------------------------------------------------------------------------------
! Interface to C-side RNG bridge (defined in delaporteC.c). Declared here at
! module scope so any procedure that uses this module gets a compiler-checked,
! explicit interface. It replaces the implicit "external unifrnd" declaration
! formerly local to rdelap_f. bind(C) pins the linker symbol to literally
! "unifrnd", so the C definition needs no F77_SUB/mangling macro.
! NOTE: interface bodies are their own scoping unit and do NOT inherit the
! module's use statements, so iso_c_binding must be re-imported inside the body
! (or brought in via an IMPORT statement).
! ------------------------------------------------------------------------------

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
! DESCRIPTION: Fortran 2008 does not have log1p as an intrinsic. This uses a
!              three-term Taylor expansion for small x to reduce relative error.
!
! GENERAL NOTE: The degree-3 Taylor expansion for log(x) around 1 is 
!               x - x^2/2 + x^3/3. The old degree-2 form, x - x^2/2, truncated
!               at O(x^3), so its absolute error was ~x^3/3 ~ 3.3e-13 at the
!               switch x = 1e-4 -- ~3e4x worse there than plain log(1+x)
!              (~1e-16). log1p enters only as log1p(beta) in
!              c = -lambda - alpha*log1p(beta), the seed of ddelap_table, so
!              the error became a uniform ~alpha*3.3e-13 relative error on every
!              d/p/q value (e.g. ddelap(300, 3e6, 1e-4, 1) was off by 1e-6 vs
!              the NB (x) Poisson oracle). The cubic term cuts truncation to
!              ~x^4/4 ~ 2.5e-17 < EPS at the switch (verified: abs err 2.5e-17
!              at 1e-4, 0 for x <= 1e-6), so the polynomial is now at least as
!              accurate as log(1+x) across all of [0, 1e-4] while keeping its
!              small-x edge. ONE/THREE folds to a compile-time constant (both
!              are parameters), so no runtime division is added.
!-------------------------------------------------------------------------------

    pure elemental function log1p(x) result(y)

        real(kind = c_double), intent(in) :: x
        real(kind = c_double)             :: y

        if (abs(x) <= 1.e-4_c_double) then
           
            ! (AA & Claude: 2026-07-07)
            y = ((x * (ONE / THREE) - HALF) * x + ONE) * x
        else
            y = log(x + ONE)
        end if

    end function log1p
    
!-------------------------------------------------------------------------------
! FUNCTION: imk (i mod k)
!
! DESCRIPTION: Calculates mod(i - 1, k) + 1 for vector recyling.
! 
! GENERAL NOTE: mod(i - 1, k) has NO internal k == 0 guard: k == 0 is integer
!               division by zero -> SIGFPE -> the R process dies (the original
!               zero-length-parameter crash). Safe ONLY because every caller
!               passes a validated parameter-vector length (na/nb/nl >= 1),
!               enforced in the R wrappers AND in the C entry points
!               (delaporteC.c) exported via R_RegisterCCallable. Any new call
!               site MUST keep k fenced >= 1. 
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
! FUNCTION: lower_bound
!
! DESCRIPTION: Index of the first element of the non-decreasing vector v that
!              is >= p, or 0 if no element is >= p.
!
! GENERAL NOTE: Drop-in replacement for minloc(v, dim = 1, mask = v >= p) on
!               sorted data - it returns the identical index - but in O(log n)
!               instead of a full O(n) masked scan. With the O(K) table build in
!               place, the masked minloc scan had become the dominant cost of
!               quantile lookups: qdelap/rdelap on np variates cost O(np * K)
!               there, e.g. ~11s for rdelap(1e6, 50, 100, 100) versus ~0.06s
!               with the binary search. Requires v non-decreasing, which the
!               cFPe-clamped cumulative sums built in qdelap_f always are.
!-------------------------------------------------------------------------------

    pure function lower_bound(v, p) result(lo)

    real(kind = c_double), intent(in)   :: v(:), p
    integer                             :: lo, hi, mid

        if (v(size(v)) < p) then
            ! Unreachable from qdelap_f in ordinary operation: every p sent
            ! here satisfies p <= x, and both table builds only exit once
            ! svec(size) >= x, so v(size(v)) >= p by transitivity on the
            ! identical stored doubles. The sole breach is the MAXTBL cap
            ! exit (an ~8GB table), the same contingency as the j == 0
            ! defence at the call site, which carries the same marker. The
            ! guard stays because without it this input violates the loop
            ! invariant below and the search would silently return
            ! size(v)---a plausible wrong answer---rather than a value the
            ! caller detects and clamps.
            ! (AA & Claude: 2026-07-07)
            lo = 0                          ! No element reaches p.  ! # nocov
        else
            lo = 1
            hi = size(v)
            do while (lo < hi)              ! Invariant: v(hi) >= p, and
                mid = (lo + hi) / 2         ! every index below lo is < p.
                if (v(mid) >= p) then
                    hi = mid
                else
                    lo = mid + 1
                end if
            end do
        end if

    end function lower_bound
    
!-------------------------------------------------------------------------------
! FUNCTION: gOMPT
!
! DESCRIPTION: Gets the OpenMP runtime's maximum thread count.
!
! GENERAL NOTE: Called once at package load to seed the package-local thread
!              setting; the per-call thread count is passed explicitly to each
!              parallel region via its num_threads clause. The former companion
!              sOMPT_f (omp_set_num_threads) was removed in 9.0.0 because it
!              mutated process-global OpenMP state shared with other packages.
!-------------------------------------------------------------------------------

    subroutine gOMPT_f(n) bind(C, name="gOMPT_f")
    
    integer(kind = c_int), intent(out) :: n
    
        n = 1_c_int
        !$ n = omp_get_max_threads()
    
    end subroutine gOMPT_f    
    
end module utils ! # nocov covr doesn't always pick up the end module
