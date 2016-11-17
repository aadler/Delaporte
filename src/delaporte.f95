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
    public :: ddelap_f, pdelap_f, qdelap_f, rdelap_f !Vectorized functions called from C

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
    integer(kind = c_int)               :: i, k                   ! Integers

    if (alpha < EPS .or. beta < EPS .or. lambda < EPS .or. x < 0) then
        call set_nan(pmf)
    else
        k = ceiling(x)
        pmf = ZERO
        do i = 0, k, 1
            pmf = pmf + exp(gamln(alpha + i) + i * log(beta) &
                      + (k - i) * log(lambda) - lambda &
                      - gamln(alpha) - gamln(i + ONE) &
                      - (alpha + i) * log1p(beta) &
                      - gamln(k - i + ONE))
        end do
    end if

end function ddelap_f_s

!----------------------------------------------------------------------------------------
! ROUTINE: ddelap_f
!
! DESCRIPTION: Vector-based PMF allowing parameter vector recycling and called from C.
!              As Fortran starts its indices at 1, for the mod function to properly
!              recycle the vectors, the index needs to be reduced by one, mod applied,
!              and then increased by one again.
!----------------------------------------------------------------------------------------

    subroutine ddelap_f (x, nx, a, na, b, nb, l, nl, lg, pmfv) &
                       bind(C, name="ddelap_f")

    integer(kind = c_int), intent(in), value         :: nx, na, nb, nl     ! Sizes
    real(kind = c_double), intent(in), dimension(nx) :: x                  ! Observations
    real(kind = c_double), intent(out), dimension(nx):: pmfv               ! Result
    real(kind = c_double), intent(in)                :: a(na), b(nb), l(nl)! Parameters
    logical(kind = c_bool), intent(in)               :: lg                 ! Log flag
    integer(kind = c_int)                            :: i                  ! Integer

        do i = 1, nx
            pmfv(i) = ddelap_f_s(x(i), a(mod(i - 1, na) + 1), b(mod(i - 1, nb) + 1), &
                                 l(mod(i - 1, nl) + 1))
        end do

        if (lg) then
            pmfv = log(pmfv)
        end if

    end subroutine ddelap_f


!----------------------------------------------------------------------------------------
! FUNCTION: pdelap_f_s
!
! DESCRIPTION: Calculate the Delaporte probability mass function for a single observation
!              and return the value or its log. Calculated through explicit summation.
!              Will coerce real observations to integer by calling ceiling.
!----------------------------------------------------------------------------------------

    function pdelap_f_s (q, alpha, beta, lambda) result (cdf)

    external set_nan
    real(kind = c_double)               :: cdf                    ! Result
    real(kind = c_double), intent(in)   :: q, alpha, beta, lambda ! Observation & Parms
    integer(kind = c_int)               :: i, k                   ! Integers

        if (alpha < EPS .or. beta < EPS .or. lambda < EPS .or. q < 0) then
            call set_nan(cdf)
        else
            k = ceiling(q)
            cdf = exp(-lambda) / ((beta + ONE) ** alpha)
            do i = 1, k
                cdf = cdf + ddelap_f_s (real(i, c_double), alpha, beta, lambda)
            end do
        end if

    end function pdelap_f_s

!----------------------------------------------------------------------------------------
! ROUTINE: pdelap_f
!
! DESCRIPTION: Vector-based CDF allowing parameter vector recycling and called from C. If
!              parameters are all singletons (not vectors) then the idea is to find the
!              largest value in the vector and build the PDF up to that point. Building
!              the vector has each succesive value piggyback off of the prior instead of
!              calling p_delap_f_s each time which increases the speed dramatically. Once
!              created, remaining values are simple lookups off of the singlevec vector.
!              Otherwise, each entry will need to build its own pmf value by calling
!              p_delap_f_s on each entry.
!----------------------------------------------------------------------------------------

    subroutine pdelap_f (q, nq, a, na, b, nb, l, nl, lt, lg, pmfv) &
                       bind(C, name="pdelap_f")

    integer(kind = c_int), intent(in), value         :: nq, na, nb, nl     ! Sizes
    real(kind = c_double), intent(in), dimension(nq) :: q                  ! Observations
    real(kind = c_double), intent(out), dimension(nq):: pmfv               ! Result
    real(kind = c_double), intent(in)                :: a(na), b(nb), l(nl)! Parameters
    logical(kind = c_bool), intent(in)               :: lg, lt             ! Log/Tail flags
    integer(kind = c_int)                            :: i, k               ! Integers
    real(kind = c_double), allocatable, dimension(:) :: singlevec          ! holds pmf

        if(na == 1 .and. nb == na .and. nl == nb) then
            if (a(1) < EPS .or. b(1) < EPS .or. l(1) < EPS) then
                do i = 1, nq
                    call set_nan(pmfv(i))
                end do
            else
                k = ceiling(maxval(q))
                allocate (singlevec(k + 1))
                singlevec(1) = exp(-l(1)) / ((b(1) + ONE) ** a(1))
                do i = 2, k + 1
                    singlevec(i) = singlevec(i - 1) &
                                   + ddelap_f_s (real(i - 1, c_double), a(1), b(1), l(1))
                end do
                do i = 1, nq
                    k = ceiling(q(i))
                    pmfv(i) = singlevec(k + 1)
                end do
                deallocate(singlevec)
            end if
        else
            do i = 1, nq
                pmfv(i) = pdelap_f_s(q(i), a(mod(i - 1, na) + 1), b(mod(i - 1, nb) + 1), &
                                     l(mod(i - 1, nl) + 1))
            end do
        end if
        if (.not.(lt)) then
            pmfv = ONE - pmfv
        end if
        if (lg) then
            pmfv = log(pmfv)
        end if

    end subroutine pdelap_f


!----------------------------------------------------------------------------------------
! FUNCTION: qdelap_f_s
!
! DESCRIPTION: Calculate the Delaporte quantile function for a single observation and
!              return the value. Calculated through explicit summation. Will coerce real
!              observations to integer by calling ceiling but returns a real to handle
!              NaN and Inf.
!----------------------------------------------------------------------------------------

    function qdelap_f_s (p, alpha, beta, lambda) result (value)

    external set_nan
    external set_inf
    real(kind = c_double), intent(in)   :: p, alpha, beta, lambda ! Percentile & Parms
    real(kind = c_double)               :: testcdf, value         ! Result and testcdf


        if (alpha < EPS .or. beta < EPS .or. lambda < EPS .or. p < 0) then
            call set_nan(value)
        else if (p >= 1) then
            call set_inf(value)
        else
            value = 0
            testcdf = exp(-lambda) / ((beta + ONE) ** alpha)
            do while (p - testcdf > EPS)
                value = value + 1
                testcdf = testcdf + ddelap_f_s (real(value, c_double), alpha, beta, lambda)
            end do
        end if

    end function qdelap_f_s

!----------------------------------------------------------------------------------------
! ROUTINE: qdelap_f
!
! DESCRIPTION: Vector-based quantile function with parameter vector recycling. If
!              parameters are all singletons (not vectors) then the idea is to find the
!              largest value in the vector and build the PDF up to that point. Building
!              the vector has each succesive value piggyback off of the prior instead of
!              calling p_delap_f_s each time which increases the speed dramatically. Once
!              created, remaining values are lookups off of the singlevec vector.
!              Otherwise, each entry will need to build its own pmf value by calling
!              q_delap_f_s on each entry.
!----------------------------------------------------------------------------------------

    subroutine qdelap_f (p, np, a, na, b, nb, l, nl, lt, lg, obsv) &
                       bind(C, name="qdelap_f")
    external set_nan
    external set_inf
    integer(kind = c_int), intent(in), value           :: np, na, nb, nl     ! Sizes
    real(kind = c_double), intent(inout), dimension(np):: p                  ! %iles
    real(kind = c_double), intent(out), dimension(np)  :: obsv               ! Result
    real(kind = c_double), intent(in)                  :: a(na), b(nb), l(nl)! Parameters
    logical(kind = c_bool), intent(in)                 :: lg, lt             ! Flags
    integer(kind = c_int)                              :: i, j               ! Integers
    real(kind = c_double), allocatable, dimension(:)   :: svec, tvec         ! Results
    real(kind = c_double)                              :: x                  ! current %

        if (lg) then
            p = exp(p)
        end if

        if (.not. lt) then
            p = ONE - p
        end if

        if(na == 1 .and. nb == na .and. nl == nb) then
            if (a(1) < EPS .or. b(1) < EPS .or. l(1) < EPS) then
                do i = 1, np
                    call set_nan(obsv(i))
                end do
            else
                x = maxval(p, 1, p < 1)
                i = 1
                allocate (svec(i))
                svec(1) = exp(-l(1)) / ((b(1) + ONE) ** a(1))
                do
                    if (svec(i) >= x) then
                        exit
                    end if
                    i = i + 1
                    allocate(tvec(1:i))
                    tvec(1:i-1) = svec
                    call move_alloc(tvec, svec)
                    svec(i) = svec(i - 1) + ddelap_f_s(real(i - 1, c_double), a(1), b(1), l(1))
                end do
                do i = 1, np
                    if (p(i) < 0) then
                        call set_nan(obsv(i))
                    else if (p(i) >= 1) then
                        call set_inf(obsv(i))
                    else
                        obsv(i) = real(position(p(i), svec) - 1)
                    end if
                end do
                deallocate(svec)
            end if
        else
            do i = 1, np
                obsv(i) = qdelap_f_s(p(i), a(mod(i - 1, na) + 1), b(mod(i - 1, nb) + 1), &
                                     l(mod(i - 1, nl) + 1))
            end do
        end if

    end subroutine qdelap_f

!----------------------------------------------------------------------------------------
! ROUTINE: rdelap_f
!
! DESCRIPTION: Vector-based random number generator with parameter vector recycling. It
!              Calls a C procedure to generate uniform random variates which jibe with
!              R's own internals and then calls qdelap_f on the uniforms. This allows
!              qdelap's singleton mode to activate if appropriate. This is the one routine
!              slower in Fortran than C++, as the vector creation and pushback is more
!              efficient in C++ STL than the ballet between allocate and move_alloc in
!              Fortran. On vectors-valued parameters Fortran is faster than C++.
!----------------------------------------------------------------------------------------

    subroutine rdelap_f (n, a, na, b, nb, l, nl, vars) bind(C, name="rdelap_f")

    external unifrnd

    integer(kind = c_int), intent(in), value           :: n, na, nb, nl      ! Sizes
    real(kind = c_double), intent(out), dimension(n)   :: vars               ! Result
    real(kind = c_double), intent(in)                  :: a(na), b(nb), l(nl)! Parameters
    real(kind = c_double), dimension(n)                :: p                  ! %iles
    logical(kind = c_bool)                             :: lg, lt             ! Flags

    call unifrnd(n, p)
    lt = .TRUE.
    lg = .FALSE.
    call qdelap_f(p, n, a, na, b, nb, l, nl, lt, lg, vars)

    end subroutine rdelap_f

!----------------------------------------------------------------------------------------
! ROUTINE: momdelap_f
!
! DESCRIPTION: Calculates method of moments estimates of parameters for a Delaporte
!              distribution based on supplied vector.
!----------------------------------------------------------------------------------------

    subroutine momdelap_f (obs, n, params) bind(C, name="momdelap_f")

    integer(kind = c_int), intent(in), value           :: n                  ! Sizes
    real(kind = c_double), intent(in), dimension(n)    :: obs                ! Result
    real(kind = c_double), intent(out), dimension(3)   :: params
    integer(kind = c_int)                              :: i
    real(kind = c_double)                              :: nm1, P, Mu_D, M2, M3, T1, delta
    real(kind = c_double)                              :: delta_i, Var_D, Skew_D, VmM_D

    nm1 = n - ONE
    P = n * sqrt(nm1) / (n - 2_c_double)
    Mu_D = ZERO
    M2 = ZERO
    M3 = ZERO
    do i = 1, n
        delta = obs(i) - Mu_D
        delta_i = delta / i
        T1 = delta * delta_i * (i - ONE)
        Mu_D = Mu_D + delta_i
        M3 = M3 + (T1 * delta_i * (i - 2) - 3 * delta_i * M2)
        M2 = M2 + T1
    end do
    Var_D = M2 / nm1
    Skew_D = P * M3 / (M2 **1.5)
    VmM_D = Var_D - Mu_D
    params(2) = 0.5 * (Skew_D * (Var_D ** 1.5) / VmM_D - 3)
    params(1) = VmM_D / (params(2) ** 2)
    params(3) = Mu_D - params(1) * params(2)

    end subroutine momdelap_f

end module delaporte
