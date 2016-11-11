module lgam

! Taken from bratio.f90 found at http://jblevins.org/mirror/amiller/#nswc
! Both procedures apparently written by ALFRED H. MORRIS as part of the
! Naval Surface Warfare Center Mathematical Library and released in public domain

use, intrinsic :: iso_c_binding, only: c_double, c_int
implicit none

contains

FUNCTION gamln1 (a) RESULT(fn_val)
!-----------------------------------------------------------------------
!     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 .LE. A .LE. 1.25
!-----------------------------------------------------------------------
! Modified to interact with C using Fortran 2003 by Avraham Adler, 2016-11-07

real(kind = c_double), intent(in) :: a
real(kind = c_double)             :: fn_val

REAL (kind = c_double) :: w, x, &
             p0 =  .577215664901533D+00, p1 =  .844203922187225D+00,  &
             p2 = -.168860593646662D+00, p3 = -.780427615533591D+00,  &
             p4 = -.402055799310489D+00, p5 = -.673562214325671D-01,  &
             p6 = -.271935708322958D-02,   &
             q1 =  .288743195473681D+01, q2 =  .312755088914843D+01,  &
             q3 =  .156875193295039D+01, q4 =  .361951990101499D+00,  &
             q5 =  .325038868253937D-01, q6 =  .667465618796164D-03,  &
             r0 = .422784335098467D+00,  r1 = .848044614534529D+00,  &
             r2 = .565221050691933D+00,  r3 = .156513060486551D+00,  &
             r4 = .170502484022650D-01,  r5 = .497958207639485D-03,  &
             s1 = .124313399877507D+01,  s2 = .548042109832463D+00,  &
             s3 = .101552187439830D+00,  s4 = .713309612391000D-02,  &
             s5 = .116165475989616D-03

IF (a >= 0.6D0) GO TO 10
w = ((((((p6*a + p5)*a + p4)*a + p3)*a + p2)*a + p1)*a + p0)/  &
    ((((((q6*a + q5)*a + q4)*a + q3)*a + q2)*a + q1)*a + 1.0D0)
fn_val = -a*w
RETURN

10 x = (a - 0.5D0) - 0.5D0
w = (((((r5*x + r4)*x + r3)*x + r2)*x + r1)*x + r0)/  &
    (((((s5*x + s4)*x + s3)*x + s2)*x + s1)*x + 1.0D0)
fn_val = x*w
RETURN
END FUNCTION gamln1

FUNCTION gamln (a) RESULT(fn_val) bind(C, name="gamln")
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

real(kind = c_double), intent(in) :: a
real(kind = c_double)             :: fn_val

real(kind = c_double) :: c0 = .833333333333333D-01, c1 = -.277777777760991D-02,  &
                         c2 = .793650666825390D-03, c3 = -.595202931351870D-03,  &
                         c4 = .837308034031215D-03, c5 = -.165322962780713D-02,  &
                          d = .418938533204673D0, t, w
integer(kind = c_int) :: i, n
!--------------------------
IF (a > 0.8D0) GO TO 10
fn_val = gamln1(a) - LOG(a)
RETURN
10 IF (a > 2.25D0) GO TO 20
t = (a - 0.5D0) - 0.5D0
fn_val = gamln1(t)
RETURN

20 IF (a >= 10.0D0) GO TO 30
n = a - 1.25D0
t = a
w = 1.0D0
DO i = 1, n
  t = t - 1.0D0
  w = t*w
END DO
fn_val = gamln1(t - 1.0D0) + LOG(w)
RETURN

30 t = (1.0D0/a)**2
w = (((((c5*t + c4)*t + c3)*t + c2)*t + c1)*t + c0)/a
fn_val = (d + w) + (a - 0.5D0)*(LOG(a) - 1.0D0)

RETURN
END FUNCTION gamln

end module lgam
