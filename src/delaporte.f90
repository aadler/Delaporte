!-------------------------------------------------------------------------------
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
!          Version 1.0: 2016-11-20
!                       Porting from C++ code in Delaporte package for R.
!          Version 1.1: 2017-03-01
!                       Various tweaks.
!          Version 1.2: 2017-08-13
!                       Corrected MoMdelap code.
!          Version 1.3: 2017-11-20
!                       Updates.
!          Version 1.4: 2018-06-18
!                       Added skew bias correction option to MoMdelap.
!          Version 1.5: 2018-11-20
!                       Zapping absolute values <= EPS to 0.
!          Version 1.6: 2018-12-10
!                       Replaced zapping with setting min to 0 and max to 1
!                       as appropriate. Less monkeying with values this way.
!          Version 2.0: 2021-01-03
!                       Setting limits as < 0 to be more consistent with R
!                       defaults for d/p/q/r functions. Returning NaN for NaN
!                       inputs per R defaults. Trapping for INFTY more
!                       consistently with base R. Should use ieee_arithmetic
!                       once current oldrelease gets deprecated and min GCC
!                       version is > 5. Using iso_fortran_env for INT64 to allow
!                       wider domain for d/pdelap.
!          Version 3.0: 2023-01-29
!                       Updated to rely on Fortran 2008 intrinsics and use
!                       ieee_arithmetic.
!          Version 4.0: 2023-08-08
!                       Added OpenMP thread control functionality.
!          Version 4.1: 2024-04-04
!                       Use "source" when allocating arrays in qdelap_f.
!          Version 5.0: 2024-06-17
!                       OpenMP is still significantly faster than extending
!                       parameter vectors and applying the elemental singleton
!                       function. Use imk helper function to calculate position
!                       for vector recycling. Turn floating point error cleaner
!                       into a function. Move lt and lg for p/delap into
!                       existing loops. While elemental functions can be even
!                       faster than OpenMP loops, they are still another loop.
!                       Saving the overhead by calling the conditionals inside
!                       the necessary loops is still faster. This is not useful
!                       for qdelap as we need to know every "real" percentile
!                       first. So run those before anything else. Other minor
!                       tweaks to code for efficiency. SIMD instructions do not
!                       save that much time and prevent compilation with current
!                       versions of flang (through 19) as it does not have full
!                       OpenMP 4.5 implementation yet.
!          Version 5.1: 2025-07-17
!                       Made checks of the lt and lg variables passed from C to
!                       be against _c_int variables, which they should be.
!          Version 5.2: 2025-12-31
!                       Declared intent of threads variable in rdelap_f.
!          Version 6.0  2026-07-08
!                       Changed binding names for header/source refactor.
!                       Use specific "only" lists to prevent scope infractions.
!                       Change unifrnd to interface and drop "external".
!                       Squashed a whole bunch of bugs/errors.
!                       Function Specific:
!                         ddelap:
!                               1) Added ddelap_f_s_log. ddelap with log = TRUE
!                                  now accumulates the PMF in log space via a
!                                  streaming log-sum-exp, so deep-tail
!                                  log-probabilities no longer return -Inf when
!                                  the linear-space PMF underflows.
!                               2) Added ddelap_table which uses a recurrence
!                                  relation based on the probability generating
!                                  function to calculate PMF values instead of
!                                  the nested loops. This is now the new fast
!                                  hot loop. Existing machinery retained to
!                                  handle exception cases.
!                               3) Bring into compliance with base R in that
!                                  negative values have 0 probability.
!                               4) Caught catatstrophic cancellation bug in
!                                  ddelap_f_s & ddelap_f_s_log.
!                         pdelap:
!                               1) When log.p = TRUE, new functions exist which
!                                  accumuate all values in log space. These end
!                                  in _log, like their ddelap analogue.
!                               2) When lower.tail = FALSE, new functions exists
!                                  that will calculate the upper tail from the
!                                  top to prevent catastrophic cancelation of
!                                  1 - CDF when CDF is very small (near or below
!                                  machine precision). These begin with "s" for
!                                  "survival".
!                               3) Prevented from overwriting passed p vector.
!                               4) Uses ddelap_table where possible for speed.
!                               5) Bring into compliance with base R in that
!                                  negative values have 0 probability so have
!                                  0 CDF as well.
!                         qdelap:
!                               1) Rewrite of the lower.tail = FALSE and
!                                  log.p = TRUE search paths (the default path
!                                  is retained verbatim and verified
!                                  bitwise-identical) to achieve full deep-tail
!                                  accuracy and establish the round-trip
!                                  identity qdelap(pdelap(x)) = x across all
!                                  modes. It is exact in log space and it is
!                                  accurate, to documented representation limits
!                                  in linear space. Extreme parameter sets now
!                                  refuse with NaN and a warning instead of
!                                  hanging. Squashed several small bugs along
!                                  the way, most of which would not surface in
!                                  normal use but were found when stress-testing
!                                  edge cases.
!                         rdelap:
!                               1) Trap singleton NaN error which resulted in
!                                  function hanging until memory was exhausted.
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

module delaporte
    use, intrinsic :: iso_c_binding,   only: c_int, c_double
    use, intrinsic :: iso_fortran_env, only: INT64
    use, intrinsic :: ieee_arithmetic, only: ieee_positive_inf, ieee_value, &
                                             ieee_quiet_nan, ieee_is_nan, &
                                             ieee_is_finite, ieee_negative_inf
    !$ use omp_lib
    use utils, only: imk, cFPe, log1p, unifrnd, lower_bound, logaddexp, ZERO, &
                     HALF, ONE, THREEHALFS, TWO, THREE, EPS, MAXD, MAXVECSIZE, &
                     TBLMAXCOEF, TBLMINRATIO, TAILSWITCH

    implicit none
    private
    public :: ddelap_f, pdelap_f, qdelap_f, rdelap_f, momdelap_f

contains

!-------------------------------------------------------------------------------
! FUNCTION:       ddelap_f_s
!
! DESCRIPTION:    Calculate the Delaporte probability mass function for a single
!                 observation and return the value or its log.
!
! GENERAL NOTE:   This function uses explicit summation. Follows R convention
!                 real observations, non-finite observations, and observations
!                 outside the support have 0 probability. Calls cFPe to
!                 implement a hard floor of 0 and a hard ceiling of 1 to prevent
!                 spurious floating point errors.
!
! ADDITIONAL:     The four loop-invariant transcendental calls---log(beta),
!                 log(lambda), log_gamma(alpha), log1p(beta)---are manually
!                 hoisted out of the summation loop . gfortran under R's default
!                 IEEE-strict flags will not hoist libm calls on its own, as
!                 they may set errno and floating-point exception flags. The
!                 measured speedup from manual hoisting is roughly 10%. The
!                 hoisted scalars replace the calls in-place with the operand
!                 order and association of the original expression untouched, so
!                 results are bitwise identical to the pre-hoist code. The
!                 summand must stay in sync with the identical expression in
!                 ddelap_f_s_log; it is deliberately duplicated there rather
!                 than shared through a helper function, since hiding it behind
!                 a call boundary blocks all invariant reuse and was measured to
!                 slow down this hot loop by roughly 24%.
!
! PARAMETER GATE: Parameters must be strictly positive AND finite. is_finite
!                 fails on both NaN and Inf, and a Delaporte with any infinite
!                 parameter has infinite mean, so no finite point carries 
!                 positive mass thus there is no PMF to report. Without the
!                 screen, Inf-driven NaNs from the log-space summands are
!                 laundered into hard 0s and 1s by the cFPe clamps. For example,
!                 ddelap(1, Inf, Inf, Inf) returned exactly 1.
!
! LOOP SUMMAND:   (log_gamma(alpha + ii) - lga) is parenthesized to force it to
!                 resolve BEFORE lambda and the other terms are added. For small
!                 alpha, log_gamma(alpha) can be large (e.g. ~115 for alpha =
!                 1e-50), and at ii = 0 the two terms are identical and must
!                 cancel exactly. Adding lambda to the running sum before that
!                 cancellation resolves means lambda gets rounded away whenever
!                 it is smaller than the ULP of that still-large intermediate
!                 value. This silently loses a real contribution before the
!                 cancellation ever happens. Isolating the cancellation first
!                 means the running sum is already at its true, potentially 0,
!                 magnitude before lambda is added, so nothing is lost
!                 regardless of how large log_gamma(alpha) is. This was
!                 invisible here because exp() of the resulting log-space error
!                 is itself unrepresentably close to exp(0) = 1 in linear space
!                 for the affected magnitudes, but the identical formula in
!                 ddelap_f_s_log makes the error directly visible. For example,
!                 ddelap_f_s_log(0, 1e-50, 1e-100, 1e-50) returned -1e-150
!                 instead of the true -1e-50 before this fix).
!-------------------------------------------------------------------------------

    pure elemental function ddelap_f_s(x, alpha, beta, lambda) result(pmf)

    real(kind = c_double), intent(in)   :: x, alpha, beta, lambda
    real(kind = c_double)               :: pmf, ii, kk
    real(kind = c_double)               :: lb, ll, lga, l1pb
    integer(INT64)                      :: i, k

        if (alpha <= ZERO .or. beta <= ZERO .or. lambda <= ZERO &
            .or. ieee_is_nan(x) &
            .or. .not. ieee_is_finite(alpha + beta + lambda)) then
            pmf = ieee_value(x, ieee_quiet_nan)
        else if (x < ZERO) then
            pmf = ZERO     ! x is below (outside) the support
        else
            pmf = ZERO     ! x >= MAXD means x is functionally +Inf
            if (x < MAXD) then
                k = floor(x, INT64)
                kk = real(k, c_double)
                if (x == kk) then ! x != kk means x is non-integral
                    ! Hoist terms that do not depend on the summation index.
                    lb = log(beta)
                    ll = log(lambda)
                    lga = log_gamma(alpha)
                    l1pb = log1p(beta)
                    do i = 0_INT64, k
                        ii = real(i, c_double)
                        pmf = pmf + &
                              exp((log_gamma(alpha + ii) - lga) + ii * lb &
                             + (kk - ii) * ll - lambda - log_gamma(ii + ONE) &
                             - (alpha + ii) * l1pb - log_gamma(kk - ii + ONE))
                    end do
                    pmf = cFPe(pmf)
                end if
            end if
        end if

    end function ddelap_f_s

!-------------------------------------------------------------------------------
! FUNCTION:       ddelap_f_s_log
!
! DESCRIPTION:    Calculate the LOG of the Delaporte probability mass function
!                 for a single observation directly in log space.
!
! GENERAL NOTE:   The linear-space function underflows to 0 once the PMF drops
!                 below the smallest double (~1e-308), so log(ddelap_f_s(...))
!                 returns -Inf for deep-tail log-probabilities that are actually
!                 representable, e.g. log P(X = 2000 | 1, 1, 1) is about -1386.
!                 Since every summand is already assembled in log space, this
!                 accumulates them with a streaming log-sum-exp: it tracks the
!                 running maximum log-term mx and the sum s of exp(term - mx),
!                 rescaling s whenever a new maximum appears. The result
!                 mx + log(s) never underflows while the true log-PMF is finite,
!                 and s lies in [1, k + 1] so it can neither underflow nor
!                 overflow. Guards and structure deliberately mirror ddelap_f_s:
!                 NaN for invalid parameters, negative, or NaN x; -Inf (the log
!                 of 0) for non-integer or over-large x. The hoisted invariants
!                 and the summand are the identical expression, token for token,
!                 as in ddelap_f_s to keep the two in sync.
!
! PARAMETER GATE: Parameters must be strictly positive AND finite. is_finite
!                 fails on both NaN and Inf, and a Delaporte with any infinite
!                 parameter has infinite mean, so no finite point carries
!                 positive mass thus there is no PMF to report. Without the
!                 screen, Inf-driven NaNs from the log-space summands are
!                 laundered into hard 0s and 1s by the cFPe clamps. For example,
!                 ddelap(1, Inf, Inf, Inf) returned exactly 1). Also, base R
!                 convention is that x below the support, including x = -Inf,
!                 has probability 0 and thus a log-probability pf -Inf. This has
!                 to be handled before the call to floor(x, INT64), which would
!                 overflow the integer kind on -Inf.
!
! LOOP SUMMAND:   (log_gamma(alpha+ii) - lga) is parenthesized to force it to
!                 resolve BEFORE lambda and the other terms are added. For small
!                 alpha, log_gamma(alpha) can be large (e.g. ~115 for alpha =
!                 1e-50), and at ii = 0 the two terms are identical and must
!                 cancel exactly. Adding lambda to the running sum before that
!                 cancellation resolves means lambda gets rounded away whenever
!                 it is smaller than the ULP of that still-large intermediate
!                 value. This silently loses a real contribution before the
!                 cancellation ever happens. This was confirmed by execution.
!                 With alpha = 1e-50, beta = 1e-100, and lambda = 1e-50, the old
!                 ordering returned log-PMF(0) = -1e-150 instead of the true
!                 -1e-50 (off by 60 orders of magnitude), which then corrupted
!                 pdelap's log-space CDF accumulation, as logaddexp of the
!                 corrupted PMF(0) with PMF(1) yielded a positive "log-CDF", an
!                 impossibility that surfaced the bug). Isolating the
!                 cancellation first means the running sum is already at its
!                 true,potentially 0, magnitude before lambda is added, so
!                 nothing is lost regardless of how large log_gamma(alpha) is.!
!-------------------------------------------------------------------------------

    pure elemental function ddelap_f_s_log(x, alpha, beta, lambda) result(lpmf)

    real(kind = c_double), intent(in)   :: x, alpha, beta, lambda
    real(kind = c_double)               :: lpmf, ii, kk, lt, mx, s
    real(kind = c_double)               :: lb, ll, lga, l1pb
    integer(INT64)                      :: i, k

        if (alpha <= ZERO .or. beta <= ZERO .or. lambda <= ZERO &
            .or. ieee_is_nan(x) &
            .or. .not. ieee_is_finite(alpha + beta + lambda)) then
            lpmf = ieee_value(x, ieee_quiet_nan)
        else if (x < ZERO) then
            lpmf = ieee_value(x, ieee_negative_inf) ! log(0)
        else
            lpmf = ieee_value(x, ieee_negative_inf)
            if (x < MAXD) then
                k = floor(x, INT64)
                kk = real(k, c_double)
                if (x == kk) then
                    ! Hoist terms that do not depend on the summation index.
                    lb = log(beta)
                    ll = log(lambda)
                    lga = log_gamma(alpha)
                    l1pb = log1p(beta)
                    mx = ieee_value(x, ieee_negative_inf)
                    s = ZERO
                    do i = 0_INT64, k
                        ii = real(i, c_double)
                        lt = (log_gamma(alpha + ii) - lga) + ii * lb &
                             + (kk - ii) * ll - lambda - log_gamma(ii + ONE) &
                             - (alpha + ii) * l1pb - log_gamma(kk - ii + ONE)
                        if (lt > mx) then
                            ! New running maximum: rescale the accumulated sum
                            ! to the new base (equivalent of multiplying by
                            ! mx(old) / lt(new)) and add ONE representing the
                            ! new "largest value" that enters the sum. On the
                            ! very first term mx is -Inf, so exp(mx - lt) is 0
                            ! and thus s is 1.
                            s = s * exp(mx - lt) + ONE
                            mx = lt
                        else
                            s = s + exp(lt - mx)
                        end if
                    end do
                    ! Log-space analogue of cFPe's ceiling of 1: a log-PMF
                    ! cannot exceed log(1) = 0. The floor of 0 needs no analogue
                    ! as its log-space image is -Inf itself.
                    lpmf = min(mx + log(s), ZERO)
                end if
            end if
        end if

    end function ddelap_f_s_log

!-------------------------------------------------------------------------------
! ROUTINE:        ddelap_table
!
! DESCRIPTION:    Fill pv(1:k+1) with the Delaporte PMF at 0, 1, ..., k in O(k)
!                 total operations using the distribution's own three-term
!                 recurrence.
!
! DERIVATION:     Differentiating the Delaporte probability generating function:
!                 P(z) = exp(lambda(z - 1)) * (1 + beta - beta * z) ** (-alpha)
!                 gives (1 + beta - beta * z) P'(z) =
!                 (lambda(1 + beta - beta * z) + alpha * beta) P(z);
!                 matching coefficients of z ** n yields:
!
!                 (1 + b)(n + 1) p(n + 1) =
!                 (b * n + lambda(1 + b) + alpha * b) p(n) - lambda * b p(n-1),
!                 with p(-1) = 0.
!
!                 Each support point therefore costs a handful of flops instead
!                 of an O(n) summation of log_gamma/exp terms, turning the CDF
!                 table build in pdelap_f from O(K**2) transcendental calls into
!                 O(K) arithmetic. Forward recursion is numerically stable here:
!                 the wanted PMF is the *dominant* solution of the recurrence
!                 (geometric tail ~ (b / (1 + b)) ** n, versus a recessive
!                 solution decaying factorially like a Poisson tail), so
!                 rounding errors are damped rather than amplified.
!
! SCALING:        p(0) = exp(-lambda) * (1+b) ** (-alpha) underflows to 0 for
!                 lambda + alpha * log(1 + b) > ~745 even though
!                 mid-distribution masses are perfectly representable, which
!                 would zero the whole forward pass. The recurrence therefore
!                 runs on scaled values ps(n) = p(n) * exp(-c), seeded with
!                 ps(0) = 1 and c = -lambda - alpha * log1p(b). Whenever the
!                 scaled value climbs past CAP = 2**900 it and its predecessor
!                 are divided by CAP and c increases by log(CAP); rescaling can
!                 only occur while the true masses are still far below the
!                 smallest double, so the masses stored as pv = ps * exp(c) in
!                 that region are correctly 0 and no CDF accuracy is lost. Once
!                 no rescale has happened (every lambda/alpha/beta is of
!                 ordinary size), c never changes and results carry no scaling
!                 error at all; with rescales the log bookkeeping costs about
!                 (lambda + alpha * log1p(b)) * EPS relative error, i.e. ~2e-13
!                 even at lambda = 1000.
!
!                 A spurious negative from the single subtraction is clamped to
!                 zero; the subtraction cannot cancel catastrophically because
!                 the positive term dominates by construction once n exceeds the
!                 mode, and below the mode both terms are of the same modest
!                 magnitude as the result.
!
!                 When lg == 1, pv receives the LOG of the PMF instead, computed
!                 as c + log(ps) so that deep-tail log-masses whose linear
!                 values underflow to 0 still come back finite - preserving the
!                 log-space guarantee ddelap_f_s_log provides on the elemental
!                 path (e.g. log P(X = 2000 | 1, 1, 1) ~ -1387, not -Inf). The
!                 log variant additionally rescales DOWNWARD whenever the scaled
!                 mass decays below CAPINV = 2**(-900): without it, ps itself
!                 underflows in a long decaying tail (the scale c only ever
!                 climbed) and the log would hit -Inf exactly like the linear
!                 path. Downward rescaling is deliberately restricted to
!                 lg == 1 so the lg == 0 path remains byte-identical to the
!                 version validated against the convolution oracle and relied
!                 on, bitwise, by the pdelap/qdelap round trip.
!
! PRECONDITION:   Callers must have validated alpha, beta, lambda as strictly
!                 positive, finite, non-NaN, and small enough that the bracketed
!                 coefficient below cannot overflow when multiplied by CAP. 
!                 pdelap_f checks coefmax < TBLMAXCOEF before selecting this
!                 path.
!-------------------------------------------------------------------------------

    pure subroutine ddelap_table(k, alpha, beta, lambda, lg, pv)

    integer, intent(in)                 :: k
    real(kind = c_double), intent(in)   :: alpha, beta, lambda
    integer(kind = c_int), intent(in)   :: lg
    real(kind = c_double), intent(out)  :: pv(k + 1)
    real(kind = c_double)               :: ob, cf0, lb, c, sc
    real(kind = c_double)               :: psm1, ps, psp1
    integer                             :: n

    real(kind = c_double), parameter    :: CAP = 2._c_double ** 900
    real(kind = c_double), parameter    :: CAPINV = 2._c_double ** (-900)
    real(kind = c_double), parameter    :: LCAP = log(CAP)
    real(kind = c_double), parameter    :: SCMIN = -690._c_double

        ob = ONE + beta                    ! Hoisted loop invariants:
        cf0 = lambda * ob + alpha * beta   ! n-independent coefficient piece
        lb = lambda * beta                 ! weight of the trailing term
        c = -lambda - alpha * log1p(beta)  ! log of the true p(0)
        sc = exp(c)                        ! current scale; may be 0 (see above)
        psm1 = ZERO                        ! scaled p(-1)
        ps = ONE                           ! scaled p(0)
        if (lg == ONEi) then
            pv(1) = c                      ! log(p(0)); log(ps) = log(1) = 0
        else
            pv(1) = ps * sc
        end if
        do n = 0, k - 1
            psp1 = ((beta * real(n, c_double) + cf0) * ps - lb * psm1) &
                   / (ob * real(n + 1, c_double))
            psp1 = max(psp1, ZERO)         ! Clamp rounding-level negatives.
            if (psp1 > CAP) then           ! Climbing out of the underflowed
                psp1 = psp1 / CAP          ! region: rescale the two live
                ps = ps / CAP              ! values and grow the scale.
                c = c + LCAP
                sc = exp(c)
            else if (lg == ONEi .and. psp1 < CAPINV .and. &
                     psp1 > ZERO) then
                
                ! Decaying tail, log variant only: rescale downward so ps
                ! never underflows and c + log(ps) stays finite as long as
                ! the true log-mass is. The linear variant skips this on
                ! purpose - see header - since its values are correctly 0
                ! down there anyway.
                psp1 = psp1 * CAP
                ps = ps * CAP
                c = c - LCAP
                sc = exp(c)
            end if
            psm1 = ps
            ps = psp1
            if (lg == ONEi) then
                if (ps > ZERO) then
                    pv(n + 2) = min(c + log(ps), ZERO) ! log version of cFPe
                else
                    ! Defence only: callers admit the log variant solely when
                    ! the guaranteed step ratio is >= TBLMINRATIO, which keeps
                    ! every ps at or above 2**(-1020) (see the constant's
                    ! comment), and the clamp above cannot produce exact zero
                    ! because forward recursion tracks the dominant solution,
                    ! keeping the subtrahend a bounded fraction of the
                    ! minuend. Retained because it is the only correct value
                    ! should a future caller relax that admission rule.
                    pv(n + 2) = ieee_value(ps, ieee_negative_inf)   ! # nocov
                end if
            
            ! While c is so negative that exp(c) is zero or subnormal, the
            ! product ps * sc would return 0 (or lose significand bits) even
            ! when the true mass exp(c + log(ps)) is a perfectly good normal
            ! double, so assemble the value in log space instead. SCMIN is
            ! -690 > log(2.2e-308) ~ -708.4, guaranteeing sc in the multiply
            ! branch is a full-precision normal number. The log branch can
            ! only run in the pre-rescale climb (c only ever increases when
            ! lg == 0), so ordinary parameter values never pay for the
            ! transcendentals. ps == 0 (hard underflow deep in a decaying
            ! tail, where the true mass is below the smallest double and 0 is
            ! the correct stored value) is folded into the multiply branch:
            ! 0 * sc is exactly 0, since sc = exp(c) is finite and in [0, 1)
            ! - c < 0 is an invariant, as c starts negative and every upward
            ! rescale leaves exp(c) <= p/ps < 1 with ps > 1 - so no 0 * Inf
            ! is possible. Folding keeps the assignment total (no branch of
            ! this chain can fall through leaving pv(n + 2) undefined) and
            ! keeps log(ps) unevaluated at zero, without a dead-code arm
            ! whose reachability would need a proof.
            else if (c > SCMIN .or. ps == ZERO) then
                pv(n + 2) = ps * sc
            else
                pv(n + 2) = exp(c + log(ps))
            end if
        end do

    end subroutine ddelap_table

!-------------------------------------------------------------------------------
! ROUTINE:        ddelap_f
!
! DESCRIPTION:    Vector-based PMF allowing parameter vector recycling.
!
! GENERAL NOTE:   As Fortran starts its indices at 1, for the mod function to
!                 properly recycle the vectors, the index needs to be reduced by
!                 one, mod applied, and then increased by one again. This is
!                 handled by the imk function found in the utils module. Follows
!                 R convention that real observations are errors and have 0
!                 probability, so returns 0 for non-integer without calling
!                 summation loop.
!
!                 When every parameter is a singleton and the observations are
!                 well-behaved, the PMF (or log-PMF) at 0..max(x) is built once
!                 by the O(K) three-term recurrence in ddelap_table and every
!                 element is answered by an O(1) lookup, replacing an O(x_i)
!                 log_gamma summation per element. This is the same
!                 O(K**2) -> O(K) restructure which pdelap_f received. One test
!                 bed showed the timing of ddelap(0:5000, 4, 6, 10) dropping
!                 from ~0.8s to under a millisecond. The routing test mirrors
!                 pdelap_f in that any vector parameter, NaN, negative, or
!                 over-large observation; invalid parameters; or parameters past
!                 TBLMAXCOEF fall through to the per-element path, which is
!                 unchanged. Non-integer x inside the fast path needs no table
!                  -Inf (lg = 1) by convention, and its floor is still <= k so
!                 sizing remains unaffected.
!-------------------------------------------------------------------------------

    subroutine ddelap_f(x, nx, a, na, b, nb, l, nl, lg, threads, pmfv) &
               bind(C, name="ddelap_f")
                        
    integer(kind = c_int), intent(in), value     :: nx, na, nb, nl
    real(kind = c_double), intent(in)            :: x(nx), a(na), b(nb), l(nl)
    integer(kind = c_int), intent(in)            :: lg, threads
    real(kind = c_double), intent(out)           :: pmfv(nx)
    real(kind = c_double), allocatable           :: pv(:)
    real(kind = c_double)                        :: xi
    integer                                      :: i, k
    
        if (na > 1 .or. nb > 1 .or. nl > 1 .or. minval(x) < ZERO .or. &
            maxval(x) > REAL(MAXVECSIZE, c_double) .or. &
            any(ieee_is_nan(x))) then
            !$omp parallel do num_threads(threads) default(shared) private(i) &
            !$omp schedule(static)
            do i = 1, nx
                if (lg == ONEi) then
                    pmfv(i) = ddelap_f_s_log(x(i), a(imk(i, na)), &
                    b(imk(i, nb)), l(imk(i, nl)))
                else
                    pmfv(i) = ddelap_f_s(x(i), a(imk(i, na)), b(imk(i, nb)), &
                    l(imk(i, nl)))
                end if
            end do
            !$omp end parallel do
        else if (a(1) <= ZERO .or. b(1) <= ZERO .or. l(1) <= ZERO .or. &
                 .not. ieee_is_finite(a(1) + b(1) + l(1))) then
            pmfv = ieee_value(x, ieee_quiet_nan)
        else if (b(1) * (maxval(x) + ONE) + l(1) * (ONE + b(1)) &
                 + a(1) * b(1) >= TBLMAXCOEF) then
            do i = 1, nx
                if (lg == ONEi) then
                    pmfv(i) = ddelap_f_s_log(x(i), a(1), b(1), l(1))
                else
                    pmfv(i) = ddelap_f_s(x(i), a(1), b(1), l(1))
                end if
            end do
        else if (lg == ONEi .and. &
                 max(l(1) / (real(floor(maxval(x)), c_double) + ONE), &
                     min(a(1), ONE) * b(1) / (ONE + b(1))) < TBLMINRATIO) then
            
            ! Degenerate point-mass-plus-dust parameters whose PMF can fall
            ! by more than the rescue band in a single recurrence step (see
            ! TBLMINRATIO); the log-space table could wrongly return -Inf
            ! where the true log-mass is finite, so compute elementally. Only
            ! lg == 1 routes: in linear space (lg == 0) a hard-zero deep in the
            ! tail is the correct stored double, so pdelap_f/qdelap_f/rdelap_f
            ! tables need no change. 
            do i = 1, nx
                pmfv(i) = ddelap_f_s_log(x(i), a(1), b(1), l(1))
            end do
        else
            k = floor(maxval(x))
            allocate(pv(k + 1))
            call ddelap_table(k, a(1), b(1), l(1), lg, pv)
            if (lg /= ONEi) pv = cFPe(pv)
            !$omp parallel do num_threads(threads) default(shared) &
            !$omp private(i, xi) schedule(static)
            do i = 1, nx
                xi = x(i)
                if (xi == real(floor(xi), c_double)) then
                    pmfv(i) = pv(floor(xi) + 1)
                else if (lg == ONEi) then
                    pmfv(i) = ieee_value(xi, ieee_negative_inf)
                else
                    pmfv(i) = ZERO
                end if
            end do
            !$omp end parallel do
        
            deallocate(pv)
        end if
        
        if (any(ieee_is_nan(pmfv))) call rwarn("NaNs produced")

    end subroutine ddelap_f
    
!-------------------------------------------------------------------------------
! FUNCTION:       pdelap_f_s
!
! DESCRIPTION:    Calculate the Delaporte cumulative distribution function for a
!                 single observation and return the value or its log.
!
! GENERAL NOTE:   Calculated through explicit summation. Follows R convention
!                 that real observations, non-finite observations, and
!                 observations outside the support have 0 probability. So
!                 non-integers are "floored" and values below the support have
!                 CDF = 0. However, +Inf has the entire support beneath it, so
!                 it has CDF = 1. The function calls cFPe to implement a hard
!                 floor of 0 and a hard ceiling of 1 to prevent spurious
!                 floating point errors.
!-------------------------------------------------------------------------------

    pure elemental function pdelap_f_s(q, alpha, beta, lambda) result(cdf)

    real(kind = c_double), intent(in)   :: q, alpha, beta, lambda
    real(kind = c_double)               :: cdf
    integer(INT64)                      :: i, k

        if (alpha <= ZERO .or. beta <= ZERO .or. lambda <= ZERO &
            .or. ieee_is_nan(q) &
            .or. .not. ieee_is_finite(alpha + beta + lambda)) then
            cdf = ieee_value(q, ieee_quiet_nan)
        else if (q < ZERO) then
            cdf = ZERO
        else if (.not. ieee_is_finite(q)) then
            cdf = ONE
        else
            k = floor(q, INT64)
            cdf = exp(-lambda) / ((beta + ONE) ** alpha)
            do i = 1_INT64, k
                cdf = cdf + ddelap_f_s(real(i, c_double), alpha, beta, lambda)
            end do
            cdf = cFPe(cdf)
        end if

    end function pdelap_f_s
    
!-------------------------------------------------------------------------------
! FUNCTION:       pdelap_f_s_log
!
! DESCRIPTION:    Calculate the LOG of the Delaporte cumulative distribution
!                 function for a single observation directly in log space.
!
! GENERAL NOTE:   Log-space image of pdelap_f_s: log(pdelap_f_s(...)) returns
!                 -Inf once the linear-space CDF underflows below ~1e-308 (e.g.
!                 deep in the lower tail with a very large lambda), even though
!                 the true log-CDF is finite. Every summand is obtained from
!                 ddelap_f_s_log (already in log space) and accumulated with the
!                 same streaming log-sum-exp ddelap_f_s_log itself uses: running
!                 maximum mx, running sum s of exp(term - mx). Used only by
!                 pdelap_f's vector-recycling/irregular-input fallback path;
!                 the singleton fast path builds an equivalent log-space table
!                 once for the whole vector instead of calling this per element.
!-------------------------------------------------------------------------------

    pure elemental function pdelap_f_s_log(q, alpha, beta, lambda) result(lcdf)

    real(kind = c_double), intent(in)   :: q, alpha, beta, lambda
    real(kind = c_double)               :: lcdf, lterm, mx, s
    integer(INT64)                      :: i, k

        if (alpha <= ZERO .or. beta <= ZERO .or. lambda <= ZERO &
            .or. ieee_is_nan(q) &
            .or. .not. ieee_is_finite(alpha + beta + lambda)) then
            lcdf = ieee_value(q, ieee_quiet_nan)
        else if (q < ZERO) then
            lcdf = ieee_value(q, ieee_negative_inf)    ! log(0): below support
        else if (.not. ieee_is_finite(q)) then
            lcdf = ZERO                                ! log(1): q = +Inf
        else
            k = floor(q, INT64)
            mx = ieee_value(q, ieee_negative_inf)
            s = ZERO
            do i = 0_INT64, k
                lterm = ddelap_f_s_log(real(i, c_double), alpha, beta, lambda)
                if (lterm > mx) then
                    s = s * exp(mx - lterm) + ONE
                    mx = lterm
                else
                    s = s + exp(lterm - mx)
                end if
            end do
            lcdf = min(mx + log(s), ZERO)         ! log-space ceiling; see cFPe
        end if

    end function pdelap_f_s_log    

!-------------------------------------------------------------------------------
! FUNCTION:       sdelap_f_s
!
! DESCRIPTION:    Calculate the Delaporte survival function P(X > q) for a
!                 single observation by summing the PMF upwards from
!                 floor(q) + 1, instead of computing 1 - CDF.
!
! GENERAL NOTE:   When the survival probability is below about 1e-16, computing
!                 it as 1 - CDF loses all significant digits to catastrophic
!                 cancellation (the CDF rounds to 1). Direct summation of the
!                 tail keeps full relative precision, mirroring base R's
!                 practice of computing the smaller tail directly.
!
!                 The summation stops when a provable bound on the remaining
!                 tail is negligible relative to the accumulated sum. The
!                 Delaporte PMF is unimodal---it is the convolution of a
!                 negative binomial with a Poisson, and the Poisson is
!                 log-concave---so once terms decrease they keep decreasing, and
!                 successive term ratios approach beta / (1 + beta) - the
!                 geometric decay rate of the dominant negative binomial tail.
!                 Bounding all future ratios by rb = max(observed ratio,
!                 beta / (1 + beta)) < 1 bounds the uncomputed remainder by
!                 term * rb / (1 - rb).
!
! PRECONDITION:   Callers must invoke this only when the upper tail is the
!                 smaller tail, i.e. when CDF(q) > 0.5. In practice every caller
!                 (the pdelap_f and qdelap_f table seeds, the saturated-band
!                 repair, and the elemental survival branches) applies the far
!                 stricter threshold of CDF > TAILSWITCH, or fires only once the
!                 forward CDF has plateaued within rounding of 1. This
!                 guarantees floor(q) is at or past the distribution's mode, so
!                 the first term cannot have underflowed while real mass remains
!                 further out, and the loop is guaranteed to terminate as terms
!                 decay to zero past the mode.
!-------------------------------------------------------------------------------

    pure elemental function sdelap_f_s(q, alpha, beta, lambda) result(sf)
    
    real(kind = c_double), intent(in)   :: q, alpha, beta, lambda
    real(kind = c_double)               :: sf, term, ptrm, rb
    integer(INT64)                      :: i

        if (alpha <= ZERO .or. beta <= ZERO .or. lambda <= ZERO .or. q < ZERO &
            .or. ieee_is_nan(q) &
            .or. .not. ieee_is_finite(alpha + beta + lambda)) then
            
            ! Defence in depth, unreachable through R: pdelap_f only calls this
            ! function when the CDF it just computed for the very same arguments
            ! exceeds one half, which invalid parameters and NaNs can never
            ! satisfy. For in that case, the CDF is NaN, and NaN comparisons are
            ! always false. Retained nevertheless, because without it a future
            ! caller passing invalid parameters would hang rather return NaN as
            ! every term from ddelap_f_s would be NaN, and NaN fails both exit
            ! comparisons in the summation loop below, so it would never end!
            sf = ieee_value(q, ieee_quiet_nan) ! # nocov
        else if (.not. ieee_is_finite(q)) then
            sf = ZERO
        else
            sf = ZERO
            ptrm = ieee_value(q, ieee_positive_inf)
            i = floor(q, INT64) + 1_INT64
            do
                term = ddelap_f_s(real(i, c_double), alpha, beta, lambda)
                sf = sf + term
               
               ! If tail fully underflowed; per the precondition we are past the
               ! mode, so every subsequent term is smaller still and contributes
               ! nothing.
                if (term <= ZERO) exit
               
               ! Only test convergence once terms are strictly decreasing (past
               ! the mode); rb < 1 is then guaranteed and the geometric
               ! remainder bound is valid.
                if (term < ptrm) then
                    rb = max(term / ptrm, beta / (beta + ONE))
                    if (term * (rb / (ONE - rb)) <= sf * EPS) exit
                end if
                ptrm = term
                i = i + 1_INT64
            end do
            sf = cFPe(sf)
        end if

    end function sdelap_f_s
    
!-------------------------------------------------------------------------------
! FUNCTION:       sdelap_f_s_log
!
! DESCRIPTION:    Calculate the LOG of the Delaporte survival function
!                 log P(X > q) for a single observation directly in log space,
!                 by log-sum-exp accumulating the tail PMF upward from
!                 floor(q) + 1.
!
! GENERAL NOTE:   The linear-space sdelap_f_s underflows to exactly 0 once the
!                 survival probability drops below ~1e-308, so log(sdelap_f_s
!                 (...)) returns -Inf for deep-tail log-survival values that are
!                 actually finite and representable (e.g. log P(X > 2000 | 1, 1,
!                 1) is about -1386, not -Inf). This mirrors sdelap_f_s exactly
!                 - same tail walk, same geometric-remainder termination rule -
!                 but each term is ddelap_f_s_log(i, ...) (already in log space,
!                 so it never underflows before the true mass does) and the
!                 running total is kept as a streaming log-sum-exp (running
!                 max mx, running sum s of exp(term - mx)), the identical
!                 pattern ddelap_f_s_log uses to accumulate its own summation.
!                 The remainder bound is evaluated as a log-space inequality
!                 (comparing logs directly) rather than by exponentiating the
!                 bound itself, so no intermediate can overflow or underflow
!                 even when accumulated terms are astronomically small.
!
! PRECONDITION:   Identical to sdelap_f_s: callers must invoke this only when
!                 the survival probability is small enough to need direct
!                 summation, i.e. when the lower-tail log-CDF at q exceeds
!                 log(TAILSWITCH) (see TAILSWITCH's definition in utils.f90).
!-------------------------------------------------------------------------------

    pure elemental function sdelap_f_s_log(q, alpha, beta, lambda) result(lsf)

    real(kind = c_double), intent(in)   :: q, alpha, beta, lambda
    real(kind = c_double)               :: lsf, lterm, plterm, rb, mx, s
    integer(INT64)                      :: i

        if (alpha <= ZERO .or. beta <= ZERO .or. lambda <= ZERO .or. q < ZERO &
            .or. ieee_is_nan(q) &
            .or. .not. ieee_is_finite(alpha + beta + lambda)) then
            
            ! Defence in depth; see the identical note in sdelap_f_s - pdelap_f
            ! only reaches here after a log-CDF computed from the same
            ! arguments exceeded log(0.5), which invalid parameters or NaN q
            ! can never satisfy (NaN comparisons are always false).
            lsf = ieee_value(q, ieee_quiet_nan) ! # nocov
        else if (.not. ieee_is_finite(q)) then
            lsf = ieee_value(q, ieee_negative_inf)   ! log(0): no mass above Inf
        else
            mx = ieee_value(q, ieee_negative_inf)
            s = ZERO
            plterm = ieee_value(q, ieee_positive_inf)
            i = floor(q, INT64) + 1_INT64
            do
                lterm = ddelap_f_s_log(real(i, c_double), alpha, beta, lambda)
                if (lterm > mx) then
                    s = s * exp(mx - lterm) + ONE
                    mx = lterm
                else
                    s = s + exp(lterm - mx)
                end if

                ! Tail term hit -Inf (hard underflow); precondition guarantees
                ! we are past the mode, so every subsequent term is smaller
                ! still and the remainder is exactly zero.
                if (.not. ieee_is_finite(lterm)) exit

                ! Only test convergence once terms are strictly decreasing
                ! (past the mode); rb < 1 is then guaranteed and the geometric
                ! remainder bound is valid. Comparison is done entirely in log
                ! space: lterm + log(rb / (1 - rb)) <= (running log-sum) +
                ! log(EPS), the log-space image of
                ! term * (rb / (1 - rb)) <= sf * EPS in sdelap_f_s.
                if (lterm < plterm) then
                    rb = max(exp(lterm - plterm), beta / (beta + ONE))
                    if (lterm + log(rb / (ONE - rb)) <= mx + log(s) + log(EPS)) &
                        exit
                end if
                plterm = lterm
                i = i + 1_INT64
            end do
            ! Log-space analogue of cFPe's ceiling of 1 (see ddelap_f_s_log).
            lsf = min(mx + log(s), ZERO)
        end if

    end function sdelap_f_s_log    

!-------------------------------------------------------------------------------
! ROUTINE:        pdelap_f
!
! DESCRIPTION:    Vector-based CDF allowing parameter vector recycling.
!
! GENERAL NOTE:   If parameters are all singletons, not vectors, the idea is to
!                 find the largest value in the vector and build the PDF up to
!                 that point. Building the vector has each succesive value
!                 piggyback off of the prior instead of calling the singleton
!                 function each time which increases the speed dramatically.
!                 Once created, remaining values are simple lookups off of the
!                 svec vector. For vector parameters, or if the R wrapper senses
!                 NaNs or other weirdness, each entry will be built by calling
!                 the appropriate singleton function. Calls cFPe to implement a
!                 hard floor of 0 and hard ceiling of 1 to prevent spurious
!                 floating point errors.
!
! LOG BRANCH:     The vector-recycling/irregular-input branches had the
!                 identical deep-tail log.p underflow as the table-built branch:
!                 the linear CDF/survival was computed first and log() taken
!                 afterward, so any value that underflowed to exactly 0 came
!                 back -Inf even though the true log value is finite. lg == 1
!                 now calls pdelap_f_s_log/sdelap_f_s_log directly, which never
!                 form linear-space values at all; lg == 0 is unchanged.
!
! UPPER TAIL:     For the upper tail, the function uses the direct tail
!                 summation (sdelap_f_s / sdelap_f_s_log) only once the survival
!                 probability is small enough that 1 - CDF would lose precision
!                 to catastrophic cancellation: when CDF > TAILSWITCH, the
!                 constant set at 1 - sqrt(EPS). Below that, the complement is
!                 accurate to ~1e-13 relative, already at the precision of the
!                 direct sum, but whose time is O(1) instead of the direct sum's
!                 O(beta**2) worst case (see TAILSWITCH's definition in
!                 utils.f90).  NaN CDFs fail both threshold tests and propagate
!                 through the complement unchanged. The lg == 1 complement
!                 recovers the linear CDF via exp() of the already-accurate
!                 log-CDF. This is safe because this branch only runs when CDF
!                 <= TAILSWITCH, nowhere near 1 or 0. So it takes the identical
!                 linear complement, and re-logs it. This avoids needing a
!                 general log(1 - exp(.)) routine, which would require a
!                 specialized EXPM1 to stay accurate as log-CDF approaches 0.
!                 Fortran 2008 does not have EXPM1 intrinsically.
!
!                 Specifically for the log version, the backwards accumulation
!                 calls new function logaddexp, the log-space image of the
!                 linear backward accumulation.
!
! TABLE MAX:      When parameters are so extreme that the ddelap_table bracketed 
!                 coefficient, at magnitude up to coefmax and multiplied by a
!                 scaled mass as large as CAP = 2**900 (~8.5e270), could
!                 overflow a double (~1.8e308), the function cannot use the fast
!                 table. TBLMAXCOEF of 1e30 leaves seven orders of magnitude of
!                 headroom. Parameters beyond it (including infinities) take the
!                 legacy O(K**2) singleton-function-based  summation build,
!                 which reproduces the pre-recurrence behavior exactly.
!
!                 Specifically for the log version, there is a TBLMINRATIO
!                 degeneracy check: a single recurrence step whose ratio is too
!                 small to survive the table's rescue band would wrongly hand
!                 back -Inf even though the true log-mass is finite. A failure
!                 here also falls back to the elemental ddelap_f_s_log, which
!                 has no scaled-recurrence overflow/underflow exposure at all.
!
! SATURATED BAND: Saturated-band repair: forward logaddexp accumulation carries
!                 an absolute offset of O(K * EPS), so log-CDF values above
!                 log(TAILSWITCH)---a CDF within ~sqrt(EPS) of 1---are pure
!                 noise on the log scale (a true -1e-25 can surface as ~-5e-10)
!                 even though they are excellent as probabilities. In the band
!                 the complement is tiny and a survival accumulation computes it
!                 to relative accuracy, so each band value is replaced by
!                 log1p(-exp(log S)); with S < sqrt(EPS) the module log1p is
!                 always inside its Taylor branch. The band is a suffix of a
!                 nondecreasing table, so a single scalar accumulator seeded
!                 with the direct tail sum walks down from the top and stops at
!                 the first value at or below the switch. This is
!                 pdelap_f_s_log's endpoint TAILSWITCH rule extended to the
!                 whole band, and makes qdelap_f's Taylor inversion of these
!                 values exact rather than noise-limited.
!-------------------------------------------------------------------------------

    subroutine pdelap_f(q, nq, a, na, b, nb, l, nl, lt, lg, threads, pmfv) &
               bind(C, name="pdelap_f")
                        
    integer(kind = c_int), intent(in), value    :: nq, na, nb, nl
    real(kind = c_double), intent(in)           :: q(nq), a(na), b(nb), l(nl)
    integer(kind = c_int), intent(in)           :: lg, lt, threads
    real(kind = c_double), intent(out)          :: pmfv(nq)
    real(kind = c_double), allocatable          :: svec(:), pv(:)
    real(kind = c_double), allocatable          :: logpv(:), logsvec(:)
    real(kind = c_double)                       :: srun
    integer                                     :: i, k

! If there are any complications at all, don't use the fast version. pdelap_f_s
! and ddelap_f_s are more robust to improper entries.

        if (na > 1 .or. nb > 1 .or. nl > 1 .or. minval(q) < ZERO .or. &
            maxval(q) > REAL(MAXVECSIZE, c_double) &
            .or. any(ieee_is_nan(q))) then
            !$omp parallel do num_threads(threads) default(shared) private(i) &
            !$omp schedule(static)
                do i = 1, nq
                    if (lg == ONEi) then
                        pmfv(i) = pdelap_f_s_log(q(i), a(imk(i, na)), &
                        b(imk(i, nb)), l(imk(i, nl)))
                    else
                        pmfv(i) = pdelap_f_s(q(i), a(imk(i, na)), &
                        b(imk(i, nb)), l(imk(i, nl)))
                    end if
                    if (lt == ONEi .and. lg == ONEi .and. &
                        pmfv(i) > log(TAILSWITCH)) then !saturated-band repair
                        pmfv(i) = log1p(-exp(sdelap_f_s_log(q(i), &
                                  a(imk(i, na)), b(imk(i, nb)), &
                                  l(imk(i, nl)))))
                    end if
                    if (lt == ZEROi) then
                        if (lg == ONEi) then
                            if (pmfv(i) > log(TAILSWITCH)) then
                                pmfv(i) = sdelap_f_s_log(q(i), a(imk(i, na)), &
                                b(imk(i, nb)), l(imk(i, nl)))
                            else if (pmfv(i) < HALF * log(EPS)) then
                                ! Mirror band: the CDF is below sqrt(EPS),
                                ! so the linear flip would round the
                                ! survival to 1 and its log to 0; the exact
                                ! complement of the log-CDF keeps the
                                ! sub-EPS resolution, as in the table path.
                                pmfv(i) = log1p(-exp(pmfv(i)))
                            else
                                pmfv(i) = log(HALF - exp(pmfv(i)) + HALF)
                            end if
                        else
                            if (pmfv(i) > TAILSWITCH) then
                                pmfv(i) = sdelap_f_s(q(i), a(imk(i, na)), &
                                b(imk(i, nb)), l(imk(i, nl)))
                            else
                                pmfv(i) = HALF - pmfv(i) + HALF   ! See dpq.h
                            end if
                        end if
                    end if
                end do
            !$omp end parallel do
        else if (a(1) <= ZERO .or. b(1) <= ZERO .or. l(1) <= ZERO .or. &
                 .not. ieee_is_finite(a(1) + b(1) + l(1))) then
            pmfv = ieee_value(q, ieee_quiet_nan)
        else
            k = floor(maxval(q))

            if (lg == ONEi) then
                allocate (logpv(k + 1))
                allocate (logsvec(k + 1))
                
                if (b(1) * (real(k, c_double) + ONE) + l(1) * (ONE + b(1)) &
                    + a(1) * b(1) < TBLMAXCOEF .and. &
                    max(l(1) / (real(k, c_double) + ONE), &
                        min(a(1), ONE) * b(1) / (ONE + b(1))) &
                        >= TBLMINRATIO) then
                    call ddelap_table(k, a(1), b(1), l(1), ONEi, logpv)
                else
                    do i = 1, k + 1
                        logpv(i) = ddelap_f_s_log(real(i - 1, c_double), &
                                                   a(1), b(1), l(1))
                    end do
                end if

                ! min(.., ZERO) is the log-space image of cFPe's ceiling of 1
                ! (as in pdelap_f_s_log): near saturation the logaddexp
                ! cumulative sum can drift a few ulp above log(1) = 0, which
                ! is not a representable log-probability and which qdelap_f
                ! correctly rejects as NaN if fed back in.
                logsvec(1) = min(logpv(1), ZERO)
                do i = 2, k + 1
                    logsvec(i) = min(logaddexp(logsvec(i - 1), logpv(i)), &
                                     ZERO)
                end do

                if (lt == ONEi) then
                    if (logsvec(k + 1) > log(TAILSWITCH)) then
                        srun = sdelap_f_s_log(real(k, c_double), a(1), &
                                              b(1), l(1))
                        logsvec(k + 1) = log1p(-exp(srun))
                        do i = k, 1, -1
                            if (logsvec(i) <= log(TAILSWITCH)) exit
                            srun = logaddexp(srun, logpv(i + 1))
                            logsvec(i) = log1p(-exp(srun))
                        end do
                    end if
                else
                    if (logsvec(k + 1) > log(TAILSWITCH)) then
                        logsvec(k + 1) = sdelap_f_s_log(real(k, c_double), &
                                                         a(1), b(1), l(1))
                    else if (logsvec(k + 1) < HALF * log(EPS)) then
                        ! Same mirror-band rule as the loop below: with the
                        ! CDF under sqrt(EPS) the linear flip rounds the
                        ! survival to 1 (log to 0); the exact complement of
                        ! the log-CDF keeps the sub-EPS resolution.
                        logsvec(k + 1) = log1p(-exp(logsvec(k + 1)))
                    else
                        logsvec(k + 1) = log(HALF - exp(logsvec(k + 1)) + HALF)
                    end if
                    do i = k, 1, -1
                        ! Mirror band repair for the survival direction:
                        ! where the forward log-CDF is below log(sqrt(EPS))
                        ! the log-survival is within ~sqrt(EPS) of 0 and the
                        ! backward accumulation can no longer resolve it,
                        ! while the still-intact forward value at this index
                        ! gives it exactly as log1p(-exp(log F)). The band
                        ! is a prefix, so once entered every remaining index
                        ! substitutes; the one accumulation step that
                        ! consumes a substituted neighbour is the same
                        ! pattern as the seed rule above and stays capped.
                        if (logsvec(i) < HALF * log(EPS)) then
                            logsvec(i) = log1p(-exp(logsvec(i)))
                        else
                            logsvec(i) = min(logaddexp(logsvec(i + 1), &
                                             logpv(i + 1)), ZERO)
                        end if
                    end do
                end if

                do i = 1, nq
                    pmfv(i) = logsvec(floor(q(i)) + 1)
                end do
                deallocate(logsvec)
                deallocate(logpv)
            else
            ! Retain the individual PMF values while building the CDF; the
            ! upper-tail branch reuses them to build the survival vector without
            ! recomputation.
            allocate (pv(k + 1))
            allocate (svec(k + 1))

            if (b(1) * (real(k, c_double) + ONE) + l(1) * (ONE + b(1)) &
                + a(1) * b(1) < TBLMAXCOEF) then
                call ddelap_table(k, a(1), b(1), l(1), ZEROi, pv)
                pv(1) = cFPe(pv(1))
            else
                pv(1) = cFPe(exp(-l(1)) / ((b(1) + ONE) ** a(1)))
                do i = 2, k + 1
                    pv(i) = ddelap_f_s(real(i - 1, c_double), a(1), &
                                       b(1), l(1))
                end do
            end if
            svec(1) = pv(1)
            do i = 2, k + 1
                svec(i) = cFPe(svec(i - 1) + pv(i))
            end do
            if (lt == ZEROi) then
                if (svec(k + 1) > TAILSWITCH) then
                    svec(k + 1) = sdelap_f_s(real(k, c_double), a(1), b(1), &
                                             l(1))
                else
                    svec(k + 1) = HALF - svec(k + 1) + HALF     ! See dpq.h
                end if
                do i = k, 1, -1
                    svec(i) = cFPe(svec(i + 1) + pv(i + 1))
                end do
            end if
            do i = 1, nq
                pmfv(i) = svec(floor(q(i)) + 1)
            end do
            deallocate(svec)
            deallocate(pv)
            end if
        end if
        
        if (any(ieee_is_nan(pmfv))) call rwarn("NaNs produced")
        
    end subroutine pdelap_f

!-------------------------------------------------------------------------------
! FUNCTION:       qdelap_f_s
!
! DESCRIPTION:    Calculate the Delaporte quantile function for a single 
!                 observation and return the value.
!
! GENERAL NOTE:   Calculated through explicit summation. Returns NaN and Inf
!                 where appropriate.
!
! QUADRATICMAX:   Explicit summation costs O(value ** 2): each step's ddelap_f_s
!                 carries ! an O(value) inner sum (measured: 0.03/0.23/0.88 s
!                 at answers of 1e3/3e3/6e3---cleanly quadratic). A truncated
!                 sum holds no information about a quantile beyond it (the
!                 accumulated CDF is still ~ 0), so past the cap the honest
!                 answer is NaN, not a guess. The parameter sets that reach the
!                 cap generally put the true quantile beyond 2 ** 53, where a
!                 count is no longer exactly representable in a double in any
!                 case. The primary refusal (with its specific warning) lives in
!                 qdelap_f; this cap is the airtight backstop for the
!                 R_RegisterCCallable entry points that bypass the driver, and
!                 must stay in sync with QMAXLEGACY there.
!-------------------------------------------------------------------------------
 
    pure elemental function qdelap_f_s(p, alpha, beta, lambda) result(value)
 
    real(kind = c_double), intent(in)   :: p, alpha, beta, lambda
    real(kind = c_double)               :: testcdf, value
    real(kind = c_double), parameter    :: QSMAX = 2._c_double ** 15
 
        if (alpha <= ZERO .or. beta <= ZERO .or. lambda <= ZERO .or. p < ZERO &
          .or. p > ONE .or. ieee_is_nan(p) &
          .or. .not. ieee_is_finite(alpha + beta + lambda)) then
            value = ieee_value(p, ieee_quiet_nan)
        else if (p == ONE) then
            value = ieee_value(p, ieee_positive_inf)
        else
            value = ZERO
            testcdf = exp(-lambda) / ((beta + ONE) ** alpha)
            do while (p > testcdf)
                value = value + ONE
                if (value > QSMAX) then
                    value = ieee_value(p, ieee_quiet_nan)
                    exit
                end if
                testcdf = testcdf + ddelap_f_s(value, alpha, beta, lambda)
            end do
        end if
 
    end function qdelap_f_s

!-------------------------------------------------------------------------------
! ROUTINE:        qdelap_f
!
! DESCRIPTION:    Vector-based quantile function w/ parameter vector recycling.
!
! GENERAL NOTE:   If parameters are all singletons (not vectors) then the idea
!                 is to find the largest value in the vector and build the PDF
!                 up to that point. Building the vector has each successive
!                 value piggyback off of the prior instead of calling pdelap_f_s
!                 each time which increases the speed dramatically. Remaining
!                 values are lookups off of the table vectors. Otherwise, each
!                 entry will need to build its own pmf value by calling
!                 qdelap_f_s on each entry.
!
! NATIVE-SPACE:   Each of the four (lower.tail, log.p) modes searches in the
!                 space its targets arrive in. Transforming first (exp of a
!                 log-probability, 1 - p of a survival probability) and then
!                 searching the linear CDF collapses every log-probability below
!                 ~ -745 to 0 and every survival probability below EPS to
!                 exactly 1 - so qdelap returned Inf where base R's qpois
!                 resolves a finite quantile - and inflates near-1 targets by an
!                 ulp, stepping round-trip quantiles up by one. Instead:
!                 (T,F) searches the linear CDF (unchanged, bitwise-stable);
!                 (F,F) searches a linear survival table built exactly as
!                 pdelap_f's upper-tail branch builds one; (T,T) splits interior
!                 targets at log(TAILSWITCH) between the forward log-CDF and,
!                 via ts = log(-expm1(p)) (2-term Taylor, |p| < 1.5e-8 there),
!                 the log-survival table; (F,T) searches the log-survival table
!                 directly. Base R >= 4.1 searches discrete quantiles in native
!                 log space for the same reason.
!
! CRITICAL:       The tables must be built with the same routines and the same
!                 accumulation order as pdelap_f (existing test suite caught it
!                 immediately). The shared ddelap_table + identical cumsum
!                 guarantees bitwise-identical CDFs on the lower.tail = TRUE,
!                 log.p = FALSE path, and survival/log tables consistent with
!                 pdelap_f's to within table-length context (see pdelap_f).
!
! TABLE GROWTH:   The table length is unknown in advance, so it starts from a
!                 moment-based estimate of mean + 10 standard deviations, which
!                 reaches any practical percentile directly. If the vector
!                 exhausts, it doubles in size so long as the accumulated table
!                 cannot yet resolve the most extreme requested target. Every
!                 sufficiency test is made against the very array the lookup
!                 loop searches, so coverage is exact, never estimated.
!                 Geometric doubling with full rebuilds costs at most twice the
!                 final build, O(K), total work and O(1) allocations. This
!                 replaces the old grow-by-one allocate/copy/move_alloc dance
!                 whose copying alone was O(K**2).
!
! ARCHITECTURE:   It is preferable to place the copy on the heap to prevent
!                 blowing out the stack. So the accumulation vectors are created
!                 as allocatable and not as fixed size.
!
! TOLERANCES:     Interior search targets carry one-sided widenings sized to
!                 their actual noise sources, in the spirit of base R's 64 * EPS
!                 fuzz in qpois.c: relative 2 * K * EPS for survival-table
!                 chain-length noise, plus an absolute 2 * K * EPS component
!                 only for targets that can have come through pdelap_f's
!                 1 - CDF flip (linear targets >= sqrt(EPS); log images
!                 likewise). The (T,F) and (T,T) forward-CDF searches are exact
!                 (bitwise prefix property). The elemental vector-parameter
!                 path retains base R's verbatim 64 * EPS multiplicative shrink
!                 on transformed interior targets.
!
! EXTREMES:       The TBLMAXCOEF legacy route costs O(answer ** 2); past
!                 QMAXLEGACY steps the remaining targets are refused as NaN
!                 with a specific warning rather than chased (a truncated CDF
!                 holds no information about a quantile beyond it, and such
!                 quantiles generally exceed 2 ** 53, beyond exact integer
!                 representation in a double). qdelap_f_s carries the matching
!                 QSMAX backstop for the C-callable entry points. MAXTBL bounds
!                 table memory at 2 ** 27 entries (~1GB per vector).
!-------------------------------------------------------------------------------
 
    subroutine qdelap_f(pp, np, a, na, b, nb, l, nl, lt, lg, threads, obsv) &
               bind(C, name="qdelap_f")
 
    integer(kind = c_int), intent(in), value        :: np, na, nb, nl
    real(kind = c_double), intent(in)               :: a(na), b(nb), l(nl)
    integer(kind = c_int), intent(in)               :: lg, lt, threads
    real(kind = c_double), intent(in)               :: pp(np)
    real(kind = c_double), intent(out)              :: obsv(np)
    real(kind = c_double), allocatable              :: p(:), svec(:), pv(:)
    real(kind = c_double), allocatable              :: surv(:)
    real(kind = c_double), allocatable              :: logpv(:), logsvec(:)
    real(kind = c_double), allocatable              :: logsurv(:)
    real(kind = c_double)                           :: x, mu, prevCeil, ts
    real(kind = c_double)                           :: xc, xs
    logical                                         :: needC, needS, needBigger
    logical                                         :: legacyCapped
    integer                                         :: i, j, k, qlo, qhi, qmid
 
    ! Hard ceiling on the lookup table (2**27 support points ~ 1GB per
    ! vector); unreachable for any parameters of practical size, present so
    ! integer arithmetic in the doubling below can never overflow and so a
    ! pathological deep-log target cannot allocate multi-GB tables before
    ! its best-effort return.
    integer, parameter                              :: MAXTBL = 2 ** 27
    ! Ceiling on the legacy grow-by-one build in the extreme-parameter
    ! (TBLMAXCOEF) route, whose cost is O(answer ** 2); must stay in sync
    ! with QSMAX in qdelap_f_s (see the cap note there). At the measured
    ! quadratic rate the worst time spent before refusing is ~half a minute.
    integer, parameter                              :: QMAXLEGACY = 2 ** 15
 
        allocate(p, source = pp)
        legacyCapped = .false.
 
        if(na == 1 .and. nb == na .and. nl == nb) then
            if (a(1) <= ZERO .or. b(1) <= ZERO .or. l(1) <= ZERO .or. &
                .not. ieee_is_finite(a(1) + b(1) + l(1))) then
                obsv = ieee_value(p, ieee_quiet_nan)
            else
                mu = a(1) * b(1) + l(1)
                k = int(min(mu + 10._c_double * &
                    sqrt(a(1) * b(1) * (ONE + b(1)) + l(1)) + 9._c_double, &
                    real(MAXTBL, c_double)))
                if (lt == ONEi .and. lg == ZEROi) then
                x = maxval(p, 1, p < ONE)
                ! Sentinel: no completed build to compare against yet, so the
                ! stagnation check below cannot fire on the very first pass.
                prevCeil = -ONE
                do
                    ! EXIT 1: extreme parameters. svec is NOT allocated here -
                    ! either this is the first iteration and it never was, or a
                    ! previous iteration deallocated it at the bottom before
                    ! doubling k (and the doubled k is what pushed the
                    ! coefficient past TBLMAXCOEF).
                    if (b(1) * (real(k, c_double) + ONE) + l(1) * &
                        (ONE + b(1)) + a(1) * b(1) >= TBLMAXCOEF) exit
                    allocate(pv(k + 1))
                    allocate(svec(k + 1))
                    call ddelap_table(k, a(1), b(1), l(1), ZEROi, pv)
                    pv(1) = cFPe(pv(1))
                    svec(1) = pv(1)
                    do i = 2, k + 1
                        svec(i) = cFPe(svec(i - 1) + pv(i))
                    end do
                    deallocate(pv)
                    ! EXIT 2: the table reached the largest target x (the
                    ! normal success), or -- defensively -- hit the MAXTBL
                    ! ceiling. Either way svec is complete for every target
                    ! it can resolve; targets beyond its saturated ceiling
                    ! are handled by the survival search in the lookup loop.
                    if (svec(k + 1) >= x .or. k >= MAXTBL) exit
                    ! EXIT 3: bitwise stagnation. x sits closer to 1 than this
                    ! linear-space cumulative sum can ever resolve: cFPe caps
                    ! every partial sum at 1, and once the tail mass being
                    ! added on each doubling falls below the ULP of the
                    ! running total, svec(k + 1) stops changing at all, no
                    ! matter how large k grows. This is PURELY a performance
                    ! optimization: it stops the doubling as soon as the CDF
                    ! plateaus instead of grinding all the way to MAXTBL (a
                    ! ~1e9 element, multi-GB table -- ~20s and ~16GB in
                    ! testing). It does not change any returned value: the
                    ! lookup loop below resolves every saturated-band target
                    ! from the accurate upper-tail survival, not from svec,
                    ! so the answer is identical whether the doubling stopped
                    ! here or was allowed to run on. The comparison is exact
                    ! equality on purpose: ddelap_table's recurrence is
                    ! deterministic, so a doubled build reproduces the shorter
                    ! build's partial sums bit for bit through the shorter
                    ! table's own length; svec(k + 1) can only be unchanged
                    ! from the previous iteration's ceiling if every term
                    ! added since then was individually too small to move the
                    ! sum.
                    if (prevCeil >= ZERO .and. svec(k + 1) == prevCeil) exit
                    prevCeil = svec(k + 1)
                    ! Reached only when looping again: the table was too
                    ! short, so release it before rebuilding at double size.
                    k = min(2 * k, MAXTBL)
                    deallocate(svec)
                end do
                ! TRUE only via EXIT 1, i.e. the extreme-parameter route.
                if (.not. allocated(svec)) then
                    ! Legacy incremental build, retained for the
                    ! extreme-parameter route. The seed is assigned in a
                    ! separate statement rather than through source= on a
                    ! continued line: gcov attributes a continued statement's
                    ! execution to its final physical line and reports the
                    ! opening line as an uncovered phantom.
                    allocate(svec(1))
                    svec(1) = exp(-l(1)) / ((b(1) + ONE) ** a(1))
                    i = 1
                    do
                        if (svec(i) >= x) exit
                        ! O(i ** 2) growth; past QMAXLEGACY the deepest
                        ! target is declared unreachable. Targets the
                        ! partial CDF did resolve still answer normally;
                        ! the rest are NaN'd in the lookup below with a
                        ! specific warning. This exit can only fire with
                        ! svec(i) < x, i.e. with at least one target
                        ! genuinely unresolved.
                        if (i > QMAXLEGACY) then
                            legacyCapped = .true.
                            exit
                        end if
                        i = i + 1
                        allocate(pv(1:i), source = ZERO)
                        pv(1:i-1) = svec
                        call move_alloc(pv, svec)
                        svec(i) = svec(i - 1) + &
                            ddelap_f_s(real(i - 1, c_double), a(1), &
                                       b(1), l(1))
                    end do
                end if
                ! Any target above the forward CDF's saturated ceiling cannot
                ! be resolved from svec: within ~K * EPS of 1 the cumulative
                ! sum's accumulated rounding deficit exceeds 1 - p, so it
                ! plateaus below such p and lower_bound below returns 0. The
                ! upper tail is nonetheless finite and accurately summable from
                ! the top, so build a survival table by the same backward
                ! recurrence pdelap_f uses -- seed the top with the directly
                ! summed tail sdelap_f_s(k), then accumulate PMF mass moving
                ! down -- and resolve those targets against it. surv(i) holds
                ! P(X > i - 1) and is monotone decreasing. Built only when a
                ! saturated target actually exists; svec and the byte-identical
                ! CDF lookup used for every resolvable target are untouched.
                ! Never entered on the capped legacy route: ddelap_table is
                ! exactly what the TBLMAXCOEF gate protects against there,
                ! and the unresolved targets are NaN'd in the lookup instead.
                if (svec(size(svec)) < x .and. .not. legacyCapped) then
                    allocate(surv(size(svec)))
                    allocate(pv(size(svec)))
                    call ddelap_table(size(svec) - 1, a(1), b(1), l(1), &
                                      ZEROi, pv)
                    surv(size(svec)) = sdelap_f_s(real(size(svec) - 1, &
                                                  c_double), a(1), b(1), l(1))
                    do i = size(svec) - 1, 1, -1
                        surv(i) = cFPe(surv(i + 1) + pv(i + 1))
                    end do
                    deallocate(pv)
                end if
                do i = 1, np
                    if (p(i) < ZERO .or. p(i) > ONE .or. ieee_is_nan(p(i))) then
                        obsv(i) = ieee_value(p(i), ieee_quiet_nan)
                    else if (p(i) == ONE) then
                        obsv(i) = ieee_value(p(i), ieee_positive_inf)
                    else
                        j = lower_bound(svec, p(i))
                        if (j > 0) then
                            obsv(i) = real(j - 1, c_double)
                        else if (legacyCapped) then
                            ! The capped legacy table plateaued below this
                            ! target; a truncated CDF carries no information
                            ! about a quantile beyond it, so NaN rather than
                            ! a bound (see QMAXLEGACY above).
                            obsv(i) = ieee_value(p(i), ieee_quiet_nan)
                        else
                            ! Saturated target: the forward CDF plateaued below
                            ! p(i) (base R's qpois returns a finite quantile for
                            ! the analogous case, and so must this). Resolve it
                            ! from the accurate survival table: find the
                            ! smallest q = qlo - 1 with P(X > q) = surv(qlo)
                            ! <= 1 - p(i). surv is decreasing and surv(size) ~ 0,
                            ! so a satisfying index always exists. Computing
                            ! ts = 1 - p(i) can shed a little precision when
                            ! p(i) itself arrived as an upper-tail 1 - p input,
                            ! but by well under one integer of quantile (the
                            ! tail is steep here), so the returned integer is
                            ! unaffected.
                            ts = ONE - p(i)
                            qlo = 1
                            qhi = size(surv)
                            do while (qlo < qhi)
                                qmid = (qlo + qhi) / 2
                                if (surv(qmid) <= ts) then
                                    qhi = qmid
                                else
                                    qlo = qmid + 1
                                end if
                            end do
                            obsv(i) = real(qlo - 1, c_double)
                        end if
                    end if
                end do
                if (allocated(surv)) deallocate(surv)
                if (legacyCapped) then
                    call rwarn("quantile too large to compute exactly for &
                               &these parameter values; NaN returned")
                end if
                deallocate(svec)
                else if (lg == ZEROi) then
                    ! lower.tail = FALSE in linear space: the targets are
                    ! survival probabilities, searched directly against a
                    ! survival table. Flipping 1 - p first inflates targets
                    ! near 1 by an ulp (off-by-one round trips) and collapses
                    ! targets below EPS to exactly 1, returning Inf where
                    ! base R's qpois resolves a finite quantile. Targets at
                    ! or above TAILSWITCH are the one exception: there the
                    ! survival table's absolute O(K * EPS) noise exceeds the
                    ! distance to 1, while the complement 1 - p is exact
                    ! (Sterbenz, p >= 1/2) and the forward CDF resolves it
                    ! accurately, so those route to a CDF-side search -- the
                    ! mirror image of pdelap_f's TAILSWITCH seed rule.
                    needS = any(pp > ZERO .and. pp < TAILSWITCH)
                    needC = any(pp >= TAILSWITCH .and. pp < ONE)
                    if (needS) then
                        x = minval(pp, 1, pp > ZERO .and. pp < TAILSWITCH)
                    end if
                    if (needC) then
                        xc = HALF - minval(pp, 1, pp >= TAILSWITCH .and. &
                                           pp < ONE) + HALF
                    end if
                    do
                        ! Extreme parameters: fall through to the elemental
                        ! search below with surv never allocated (the exact
                        ! analogue of the lower-tail legacy route).
                        if (b(1) * (real(k, c_double) + ONE) + l(1) * &
                            (ONE + b(1)) + a(1) * b(1) >= TBLMAXCOEF) exit
                        allocate(pv(k + 1))
                        allocate(svec(k + 1))
                        allocate(surv(k + 1))
                        call ddelap_table(k, a(1), b(1), l(1), ZEROi, pv)
                        pv(1) = cFPe(pv(1))
                        svec(1) = pv(1)
                        do i = 2, k + 1
                            svec(i) = cFPe(svec(i - 1) + pv(i))
                        end do
                        ! Identical survival construction to pdelap_f's
                        ! upper-tail branch: an accurate direct tail sum
                        ! seeds the top when the forward CDF sits inside the
                        ! cancellation band, backward PMF accumulation fills
                        ! below. surv(i) = P(X > i - 1), monotone decreasing.
                        if (svec(k + 1) > TAILSWITCH) then
                            surv(k + 1) = sdelap_f_s(real(k, c_double), &
                                                     a(1), b(1), l(1))
                        else
                            surv(k + 1) = HALF - svec(k + 1) + HALF ! See dpq.h
                        end if
                        do i = k, 1, -1
                            surv(i) = cFPe(surv(i + 1) + pv(i + 1))
                        end do
                        deallocate(pv)
                        ! Sufficiency is tested on the very arrays the
                        ! lookups search, so coverage is exact. surv(k + 1)
                        ! shrinks as k doubles -- the seed is a direct tail
                        ! sum that decreases strictly until it underflows to
                        ! zero -- so any positive target is eventually
                        ! reached; MAXTBL is the defensive ceiling. No
                        ! short-circuit guarantee in Fortran, so each array
                        ! is referenced only under its own need flag.
                        needBigger = .false.
                        if (needS) then
                            if (surv(k + 1) > x) needBigger = .true.
                        end if
                        if (needC) then
                            if (svec(k + 1) < xc) needBigger = .true.
                        end if
                        if (.not. needBigger .or. k >= MAXTBL) exit
                        k = min(2 * k, MAXTBL)
                        deallocate(svec)
                        deallocate(surv)
                    end do
                    if (.not. allocated(surv)) then
                        ! TBLMAXCOEF corner: qdelap_f_s on the flipped target
                        ! performs the identical seed, identical ddelap_f_s
                        ! forward summation, and identical >= stopping
                        ! comparison the legacy incremental table performed,
                        ! so the historical extreme-parameter semantics
                        ! (including the 1 - p flip) are reproduced bit for
                        ! bit without a second copy of the legacy machinery.
                        do i = 1, np
                            obsv(i) = qdelap_f_s(HALF - pp(i) + HALF, a(1), &
                                                 b(1), l(1))
                        end do
                        ! NaN for a valid in-range input can only be
                        ! qdelap_f_s's QSMAX iteration cap (parameters were
                        ! screened above); surface it with the same specific
                        ! warning as the legacy route.
                        if (any(ieee_is_nan(obsv) .and. &
                                .not. ieee_is_nan(pp) .and. &
                                pp >= ZERO .and. pp <= ONE)) then
                            call rwarn("quantile too large to compute &
                                       &exactly for these parameter values; &
                                       &NaN returned")
                        end if
                    else
                        do i = 1, np
                            if (ieee_is_nan(pp(i)) .or. pp(i) < ZERO .or. &
                                pp(i) > ONE) then
                                obsv(i) = ieee_value(pp(i), ieee_quiet_nan)
                            else if (pp(i) == ZERO) then
                                obsv(i) = ieee_value(pp(i), &
                                                     ieee_positive_inf)
                            else if (pp(i) == ONE) then
                                obsv(i) = ZERO
                            else if (pp(i) >= TAILSWITCH) then
                                ! Near-1 survival target: search the forward
                                ! CDF for the exact complement instead (see
                                ! the routing note above). Same quantile
                                ! criterion, computed from the side that can
                                ! actually resolve it.
                                ts = HALF - pp(i) + HALF
                                j = lower_bound(svec, ts)
                                if (j > 0) then
                                    obsv(i) = real(j - 1, c_double)
                                else
                                    ! Only reachable if the doubling hit
                                    ! MAXTBL (2^30 entries, multi-GB) before
                                    ! covering the target; table end as a
                                    ! best effort. Untestable without that
                                    ! allocation, retained as the safety net
                                    ! for lower_bound's index-0 contract.
                                    obsv(i) = real(size(svec) - 1, &
                                                   c_double)   ! # nocov
                                end if
                            else
                                ! Smallest q with P(X > q) <= p; surv is
                                ! decreasing so the binary search descends.
                                ! The target is widened by a relative
                                ! 2 * K * EPS: the caller's survival
                                ! probability and this table approximate the
                                ! same quantity through accumulation chains
                                ! of different lengths, so they can disagree
                                ! by O(K * EPS) relative (measured: ~2.8e3
                                ! ulps at K ~ 3e3), and without the band a
                                ! round-tripped target lands one support
                                ! point high. Same rationale as base R's
                                ! 64 * EPS fuzz in qpois, scaled to the
                                ! actual chain length. If the doubling hit
                                ! MAXTBL before covering the target this
                                ! returns the table end as a best effort.
                                ts = pp(i) * (ONE + TWO * &
                                     real(size(surv), c_double) * EPS)
                                ! Targets at or above sqrt(EPS) can have come
                                ! through pdelap_f's 1 - CDF flip (taken when
                                ! CDF <= TAILSWITCH), whose error is absolute
                                ! O(K * EPS), not relative; widen by the same
                                ! absolute amount there. Below sqrt(EPS) the
                                ! flip is never the source (pdelap_f seeds
                                ! with the direct tail sum instead), and an
                                ! absolute band would wipe out deep-tail
                                ! resolution, so it is not applied.
                                if (pp(i) >= sqrt(EPS)) then
                                    ts = ts + TWO * &
                                         real(size(surv), c_double) * EPS
                                end if
                                qlo = 1
                                qhi = size(surv)
                                do while (qlo < qhi)
                                    qmid = (qlo + qhi) / 2
                                    if (surv(qmid) <= ts) then
                                        qhi = qmid
                                    else
                                        qlo = qmid + 1
                                    end if
                                end do
                                obsv(i) = real(qlo - 1, c_double)
                            end if
                        end do
                        deallocate(svec)
                        deallocate(surv)
                    end if
                else
                    ! log.p = TRUE (either tail), searched natively in log
                    ! space. exp() of the target first collapses every
                    ! log-probability below ~ -745 to 0 and every one above
                    ! ~ -EPS to 1, destroying exactly the deep-tail
                    ! resolution the log-space tables carry.
                    if (lt == ONEi) then
                        ! Interior targets split at log(TAILSWITCH), the log
                        ! image of pdelap_f's cancellation band: at or below
                        ! it the forward log-CDF resolves the target; above
                        ! it (within ~sqrt(EPS) of 0) the log-CDF is
                        ! saturated and the target converts to a
                        ! log-survival target log(1 - exp(p)) =
                        ! log(-expm1(p)), evaluated by the 2-term Taylor
                        ! -p * (1 + p / 2) since |p| < 1.5e-8 there
                        ! (relative error ~ p**2 / 6, far below EPS;
                        ! gfortran has no expm1 intrinsic). The Taylor
                        ! argument is monotone in p, so the sizing target xs
                        ! is the image of the largest survival-side p.
                        needC = any(ieee_is_finite(pp) .and. &
                                    pp <= log(TAILSWITCH))
                        needS = any(pp > log(TAILSWITCH) .and. pp < ZERO)
                        if (needC) then
                            xc = maxval(pp, 1, ieee_is_finite(pp) .and. &
                                        pp <= log(TAILSWITCH))
                        end if
                        if (needS) then
                            xs = maxval(pp, 1, pp > log(TAILSWITCH) .and. &
                                        pp < ZERO)
                            xs = log(-xs * (ONE + xs * HALF))
                        end if
                    else
                        ! lower.tail = FALSE: interior targets are
                        ! log-survival probabilities; the deepest one sizes
                        ! the table. Targets above log(TAILSWITCH) -- a
                        ! survival probability within ~sqrt(EPS) of 1 --
                        ! route to the CDF side through the exact complement
                        ! log(1 - exp(p)) = log(-expm1(p)), by the same
                        ! 2-term Taylor as the lower-tail split: there the
                        ! log-survival table's own absolute O(K * EPS) noise
                        ! exceeds the target's distance to 0, while the
                        ! forward log-CDF resolves the complement
                        ! accurately. Mirror image of pdelap_f's TAILSWITCH
                        ! rule, in the opposite direction to the lt == 1
                        ! split above.
                        needC = any(pp > log(TAILSWITCH) .and. pp < ZERO)
                        needS = any(ieee_is_finite(pp) .and. &
                                    pp <= log(TAILSWITCH))
                        if (needS) then
                            xs = minval(pp, 1, ieee_is_finite(pp) .and. &
                                        pp <= log(TAILSWITCH))
                        end if
                        if (needC) then
                            xc = minval(pp, 1, pp > log(TAILSWITCH) .and. &
                                        pp < ZERO)
                            xc = log(-xc * (ONE + xc * HALF))
                        end if
                    end if
                    do
                        allocate(logpv(k + 1))
                        allocate(logsvec(k + 1))
                        ! Same two gates as pdelap_f's log path: extreme
                        ! coefficients or a per-step PMF ratio below the
                        ! rescue band route to the elemental fill (see
                        ! TBLMINRATIO in utils.f90).
                        if (b(1) * (real(k, c_double) + ONE) + l(1) * &
                            (ONE + b(1)) + a(1) * b(1) < TBLMAXCOEF .and. &
                            max(l(1) / (real(k, c_double) + ONE), &
                                min(a(1), ONE) * b(1) / (ONE + b(1))) &
                                >= TBLMINRATIO) then
                            call ddelap_table(k, a(1), b(1), l(1), ONEi, &
                                              logpv)
                        else
                            do i = 1, k + 1
                                logpv(i) = ddelap_f_s_log(real(i - 1, &
                                           c_double), a(1), b(1), l(1))
                            end do
                        end if
                        ! Same log-space ceiling as pdelap_f's tables.
                        logsvec(1) = min(logpv(1), ZERO)
                        do i = 2, k + 1
                            logsvec(i) = min(logaddexp(logsvec(i - 1), &
                                             logpv(i)), ZERO)
                        end do
                        if (needS) then
                            ! Log image of pdelap_f's upper-tail seed rule.
                            allocate(logsurv(k + 1))
                            if (logsvec(k + 1) > log(TAILSWITCH)) then
                                logsurv(k + 1) = sdelap_f_s_log(real(k, &
                                                 c_double), a(1), b(1), l(1))
                            else
                                logsurv(k + 1) = log(HALF - &
                                                 exp(logsvec(k + 1)) + HALF)
                            end if
                            do i = k, 1, -1
                                logsurv(i) = min(logaddexp(logsurv(i + 1), &
                                                 logpv(i + 1)), ZERO)
                            end do
                        end if
                        deallocate(logpv)
                        ! Sufficiency is tested on the very arrays the
                        ! lookups search, so coverage is exact. Fortran does
                        ! not guarantee short-circuit evaluation, so each
                        ! array is only referenced under its own need flag.
                        needBigger = .false.
                        if (needC) then
                            if (logsvec(k + 1) < xc) needBigger = .true.
                        end if
                        if (needS) then
                            if (logsurv(k + 1) > xs) needBigger = .true.
                        end if
                        if (.not. needBigger .or. k >= MAXTBL) exit
                        k = min(2 * k, MAXTBL)
                        deallocate(logsvec)
                        if (allocated(logsurv)) deallocate(logsurv)
                    end do
                    if (lt == ONEi) then
                        do i = 1, np
                            if (ieee_is_nan(pp(i)) .or. pp(i) > ZERO) then
                                obsv(i) = ieee_value(pp(i), ieee_quiet_nan)
                            else if (pp(i) == ZERO) then
                                ! log(1): the entire mass, so q = Inf.
                                obsv(i) = ieee_value(pp(i), &
                                                     ieee_positive_inf)
                            else if (.not. ieee_is_finite(pp(i))) then
                                obsv(i) = ZERO         ! log(0): q = 0
                            else if (pp(i) <= log(TAILSWITCH)) then
                                j = lower_bound(logsvec, pp(i))
                                if (j > 0) then
                                    obsv(i) = real(j - 1, c_double)
                                else
                                    ! Only reachable if the doubling hit
                                    ! MAXTBL (2^30 entries, multi-GB) before
                                    ! covering the target; table end as a
                                    ! best effort. Untestable without that
                                    ! allocation, retained as the safety net
                                    ! for lower_bound's index-0 contract.
                                    obsv(i) = real(size(logsvec) - 1, &
                                                   c_double)   ! # nocov
                                end if
                            else
                                ! Survival-side target. With pdelap_f's
                                ! saturated-band values now derived from the
                                ! survival side, the 2-term Taylor inversion
                                ! recovers log S to relative accuracy, so
                                ! the tolerance is the same chain-length
                                ! band as the lower.tail = FALSE searches:
                                ! additive 2 * (K + |ts|) * EPS in log space
                                ! is a relative 2 * K * EPS band in linear
                                ! space plus the log value's representation
                                ! noise. (An absolute widening here would
                                ! flatten every target below it and destroy
                                ! the deep-band resolution the band repair
                                ! provides.)
                                ts = log(-pp(i) * (ONE + pp(i) * HALF))
                                ts = ts + TWO * EPS * &
                                     (real(size(logsurv), c_double) + &
                                      abs(ts))
                                qlo = 1
                                qhi = size(logsurv)
                                do while (qlo < qhi)
                                    qmid = (qlo + qhi) / 2
                                    if (logsurv(qmid) <= ts) then
                                        qhi = qmid
                                    else
                                        qlo = qmid + 1
                                    end if
                                end do
                                obsv(i) = real(qlo - 1, c_double)
                            end if
                        end do
                    else
                        do i = 1, np
                            if (ieee_is_nan(pp(i)) .or. pp(i) > ZERO) then
                                obsv(i) = ieee_value(pp(i), ieee_quiet_nan)
                            else if (pp(i) == ZERO) then
                                obsv(i) = ZERO         ! log(1): S = 1, q = 0
                            else if (.not. ieee_is_finite(pp(i))) then
                                ! log(0): S = 0, only reached in the limit.
                                obsv(i) = ieee_value(pp(i), &
                                                     ieee_positive_inf)
                            else if (pp(i) > log(TAILSWITCH)) then
                                ! Near-1 survival target: search the forward
                                ! log-CDF for the Taylor complement instead
                                ! (see the routing note above).
                                ts = log(-pp(i) * (ONE + pp(i) * HALF))
                                j = lower_bound(logsvec, ts)
                                if (j > 0) then
                                    obsv(i) = real(j - 1, c_double)
                                else
                                    ! Only reachable if the doubling hit
                                    ! MAXTBL (2^30 entries, multi-GB) before
                                    ! covering the target; table end as a
                                    ! best effort. Untestable without that
                                    ! allocation, retained as the safety net
                                    ! for lower_bound's index-0 contract.
                                    obsv(i) = real(size(logsvec) - 1, &
                                                   c_double)   ! # nocov
                                end if
                            else
                                ! Smallest q with log P(X > q) <= p, with
                                ! the same chain-length tolerance as the
                                ! linear survival search: additive
                                ! 2 * (K + |p|) * EPS in log space is the
                                ! image of a relative 2 * K * EPS band in
                                ! linear space plus the log value's own
                                ! representation noise.
                                ts = pp(i) + TWO * EPS * &
                                     (real(size(logsurv), c_double) + &
                                      abs(pp(i)))
                                ! Same flip-band absolute widening as the
                                ! linear search, expressed additively in log
                                ! space: an absolute band of 2 * K * EPS on a
                                ! survival value s is 2 * K * EPS / s in its
                                ! log. Bounded by 2 * K * sqrt(EPS) on this
                                ! branch, orders below any support step.
                                if (pp(i) >= HALF * log(EPS)) then
                                    ts = ts + TWO * &
                                         real(size(logsurv), c_double) * &
                                         EPS * exp(-pp(i))
                                end if
                                qlo = 1
                                qhi = size(logsurv)
                                do while (qlo < qhi)
                                    qmid = (qlo + qhi) / 2
                                    if (logsurv(qmid) <= ts) then
                                        qhi = qmid
                                    else
                                        qlo = qmid + 1
                                    end if
                                end do
                                obsv(i) = real(qlo - 1, c_double)
                            end if
                        end do
                    end if
                    deallocate(logsvec)
                    if (allocated(logsurv)) deallocate(logsurv)
                end if
            end if
        else
            if (lg == ONEi) p = exp(p)
            if (lt == ZEROi) p = HALF - p + HALF  ! See dpq.h in R source code
            if (lg == ONEi .or. lt == ZEROi) then
                ! Base R's fuzz (qpois.c and kin): a transformed target can
                ! arrive one ulp above the CDF value it round-trips from
                ! (1 - p inflation, exp/log noise), stepping the quantile up
                ! by one. Shrink strictly interior targets by 64 EPS exactly
                ! as base R does before its search; boundaries and the
                ! untransformed path are untouched.
                do i = 1, np
                    if (p(i) > ZERO .and. p(i) < ONE) then
                        p(i) = p(i) * (ONE - 64._c_double * EPS)
                    end if
                end do
            end if
            !$omp parallel do num_threads(threads) default(shared) private(i) &
            !$omp schedule(static)
            do i = 1, np
                obsv(i) = qdelap_f_s(p(i), a(imk(i, na)), b(imk(i, nb)), &
                          l(imk(i, nl)))
            end do
            !$omp end parallel do
        end if
        
        deallocate(p)
        
    end subroutine qdelap_f

!-------------------------------------------------------------------------------
! ROUTINE:        rdelap_f
!
! DESCRIPTION:    Vector-based random number generator with parameter vector
!                 recycling.
!
! GENERAL NOTE:   The routine calls a C procedure to generate uniform random
!                 variates which jibe with R's own internals and then calls
!                 qdelap_f on the uniforms. This allows qdelap's singleton mode
!                 to activate if appropriate.
!-------------------------------------------------------------------------------

    subroutine rdelap_f(n, a, na, b, nb, l, nl, threads, vars) &
               bind(C, name="rdelap_f")

    integer(kind = c_int), intent(in), value           :: n, na, nb, nl
    real(kind = c_double), intent(in)                  :: a(na), b(nb), l(nl)
    real(kind = c_double), intent(out)                 :: vars(n)
    integer(kind = c_int), intent(in)                  :: threads
    real(kind = c_double)                              :: p(n)

        call unifrnd(n, p)
        call qdelap_f(p, n, a, na, b, nb, l, nl, 1, 0, threads, vars)

    end subroutine rdelap_f

!-------------------------------------------------------------------------------
! ROUTINE:      momdelap_f
!
! DESCRIPTION:  Calculates method of moments estimates of parameters for a 
!               Delaporte distribution based on supplied vector. Based on
!               algorithms of Welford, Knuth, and Cook.
!               https://www.johndcook.com/blog/skewness_kurtosis/
!-------------------------------------------------------------------------------

    pure subroutine momdelap_f(obs, n, tp, params) bind(C, name="momdelap_f")

    integer(kind = c_int), intent(in), value :: n
    integer(kind = c_int), intent(in)        :: tp
    real(kind = c_double), intent(in)        :: obs(n)
    real(kind = c_double), intent(out)       :: params(3)
    real(kind = c_double)                    :: nnm1, P, Mu_D, M2, M3, T1
    real(kind = c_double)                    :: delta, delta_i, nn, Var_D
    real(kind = c_double)                    :: Skew_D, VmM_D, ii
    integer                                  :: i

        nn = real(n, c_double)
        nnm1 = nn - ONE
        select case (tp)
            case (ONEi)
                P = ONE
            case (2_c_int)
                P = sqrt(nn * nnm1) / (nn - TWO)
            case (3_c_int)
                P = (nnm1 / nn) ** THREEHALFS
            case default
                P = sqrt(nn * nnm1) / (nn - TWO)
        end select
        Mu_D = ZERO
        M2 = ZERO
        M3 = ZERO
        do i = 1, n
            ii = real(i, c_double)
            delta = obs(i) - Mu_D
            delta_i = delta / ii
            T1 = delta * delta_i * (ii - ONE)
            Mu_D = Mu_D + delta_i
            M3 = M3 + (T1 * delta_i * (ii - TWO) - THREE * delta_i * M2)
            M2 = M2 + T1
        end do
        Var_D = M2 / nnm1
        Skew_D = P * sqrt(nn) * M3 / (M2 ** THREEHALFS)
        VmM_D = Var_D - Mu_D
        params(2) = HALF * (Skew_D * (Var_D ** THREEHALFS) - Mu_D - THREE &
                            * VmM_D) / VmM_D
        params(1) = VmM_D / (params(2) ** 2)
        params(3) = Mu_D - params(1) * params(2)
 
    end subroutine momdelap_f

end module delaporte
