module lgam

! Taken from bratio.f90 found at http://jblevins.org/mirror/amiller/#nswc
! Both procedures written by ALFRED H. MORRIS as part of the
! Naval Surface Warfare Center Mathematical Library and released in public domain

  use, intrinsic :: iso_c_binding, only: c_double, c_int
  implicit none

  contains

  function gamln1 (a) result(fn_val)
  !-----------------------------------------------------------------------
  !     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 <= A <= 1.25
  !-----------------------------------------------------------------------
  ! Modified to interact with C using Fortran 2003 by Avraham Adler, 2016-11-07

    real(kind = c_double), parameter  :: ONE = 1_c_double
    real(kind = c_double), intent(in) :: a
    real(kind = c_double)             :: fn_val
    real(kind = c_double) :: w, x, &
               p0 =  .577215664901533_c_double,     &
               p1 =  .844203922187225_c_double,     &
               p2 = -.168860593646662_c_double,     &
               p3 = -.780427615533591_c_double,     &
               p4 = -.402055799310489_c_double,     &
               p5 = -.673562214325671e-01_c_double, &
               p6 = -.271935708322958e-02_c_double, &
               q1 =  .288743195473681e+01_c_double, &
               q2 =  .312755088914843e+01_c_double, &
               q3 =  .156875193295039e+01_c_double, &
               q4 =  .361951990101499_c_double,     &
               q5 =  .325038868253937e-01_c_double, &
               q6 =  .667465618796164e-03_c_double, &
               r0 = .422784335098467_c_double,      &
               r1 = .848044614534529_c_double,      &
               r2 = .565221050691933_c_double,      &
               r3 = .156513060486551_c_double,      &
               r4 = .170502484022650e-01_c_double,  &
               r5 = .497958207639485e-03_c_double,  &
               s1 = .124313399877507e+01_c_double,  &
               s2 = .548042109832463_c_double,      &
               s3 = .101552187439830_c_double,      &
               s4 = .713309612391000e-02_c_double,  &
               s5 = .116165475989616e-03_c_double

    if(a < 0.6_c_double) then
        w = ((((((p6 * a + p5) * a + p4) * a + p3) * a + p2) * a + p1) * a + p0) / &
            ((((((q6 * a + q5) * a + q4) * a + q3) * a + q2) * a + q1) * a + ONE)
        fn_val = -a * w
    else
        x = a - ONE
        w = (((((r5 * x + r4) * x + r3) * x + r2) * x + r1) * x + r0) / &
            (((((s5 * x + s4) * x + s3) * x + s2) * x + s1) * x + ONE)
        fn_val = x*w
    end if

  end function gamln1

  function gamln (a) result(fn_val) bind(C, name = "gamln")
  !-----------------------------------------------------------------------
  !            EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A
  !-----------------------------------------------------------------------
  !     WRITTEN BY ALFRED H. MORRIS
  !          NAVAL SURFACE WARFARE CENTER
  !          DAHLGREN, VIRGINIA
  !--------------------------
  !     D = 0.5*(LN(2*PI) - 1)
  !--------------------------
  ! Modified to interact with C using Fortran 2003 by Avraham Adler, 2016-11-07

    real(kind = c_double), parameter  :: ONE = 1_c_double
    real(kind = c_double), intent(in) :: a
    real(kind = c_double)             :: fn_val
    real(kind = c_double) :: c0 = 0.833333333333333e-01_c_double, &
                           c1 = -.277777777760991e-02_c_double, &
                           c2 = .793650666825390e-03_c_double,  &
                           c3 = -.595202931351870e-03_c_double, &
                           c4 = .837308034031215e-03_c_double,  &
                           c5 = -.165322962780713e-02_c_double, &
                            d = .418938533204673_c_double, t, w
    integer(kind = c_int) :: i, n

    if (a <= 0.8_c_double) then
        fn_val = gamln1(a) - log(a)
    else if(a <= 2.25_c_double) then
        fn_val = gamln1(a - ONE)
    else if (a < 10_c_double) then
        n = a - 1.25_c_double
        t = a
        w = ONE
        do  i = 1, n
            t = t - ONE
            w = t * w
        end do
        fn_val = gamln1(t - ONE) + log(w)
    else
        t = (ONE / a) ** 2
        w = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / a
        fn_val = (d + w) + (a - 0.5_c_double) * (LOG(a) - ONE)
    end if

  end function gamln

end module lgam
