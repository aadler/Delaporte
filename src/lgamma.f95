!----------------------------------------------------------------------------------------
!
! MODULE: lgam
!
! AUTHOR: ALFRED H. MORRIS
!         Modified for Fortran 2003 and use with C by
!         Avraham Adler <Avraham.Adler@gmail.com>
!
! DESCRIPTION: Taken from bratio.f90 found at http://jblevins.org/mirror/amiller/#nswc
!              Both procedures written by ALFRED H. MORRIS as part of the Naval Surface
!              Warfare Center Mathematical Library and released in public domain. While
!              Fortran 2008 has an intrinsic GAMMA and LOG_GAMMA, Fortran 2003 does not.
!
! HISTORY:
!          Version 1.0: 2016-11-20
!                       Gently massaged into Fortran 2003
!          Version 1.1: 2017-02-04
!                       Adjusted code slightly to become elemental by explicitly
!                       declaring parameters as such.
! LICENSE: The original is a work of the US government and thus in the public domain.
!          The updated code below is released under the BSD-2 License below:
!
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
module lgam

    use, intrinsic :: iso_c_binding, only: c_double
    use utils
    implicit none

contains

!----------------------------------------------------------------------------------------
! FUNCTION: gamln1
!
! DESCRIPTION: EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 <= A <= 1.25
!----------------------------------------------------------------------------------------

    elemental function gamln1(a) result(fn_val)

    real(kind = c_double), intent(in) :: a
    real(kind = c_double)             :: fn_val
    real(kind = c_double)             :: w, x
    real(kind = c_double), parameter  :: &
        P0 =  .577215664901533_c_double,     &
        P1 =  .844203922187225_c_double,     &
        P2 = -.168860593646662_c_double,     &
        P3 = -.780427615533591_c_double,     &
        P4 = -.402055799310489_c_double,     &
        P5 = -.673562214325671e-01_c_double, &
        P6 = -.271935708322958e-02_c_double, &
        Q1 =  .288743195473681e+01_c_double, &
        Q2 =  .312755088914843e+01_c_double, &
        Q3 =  .156875193295039e+01_c_double, &
        Q4 =  .361951990101499_c_double,     &
        Q5 =  .325038868253937e-01_c_double, &
        Q6 =  .667465618796164e-03_c_double, &
        R0 = .422784335098467_c_double,      &
        R1 = .848044614534529_c_double,      &
        R2 = .565221050691933_c_double,      &
        R3 = .156513060486551_c_double,      &
        R4 = .170502484022650e-01_c_double,  &
        R5 = .497958207639485e-03_c_double,  &
        S1 = .124313399877507e+01_c_double,  &
        S2 = .548042109832463_c_double,      &
        S3 = .101552187439830_c_double,      &
        S4 = .713309612391000e-02_c_double,  &
        S5 = .116165475989616e-03_c_double

        if(a < 0.6_c_double) then
            w = ((((((P6 * a + P5) * a + P4) * a + P3) * a + P2) * a + P1) * a + P0) / &
                ((((((Q6 * a + Q5) * a + Q4) * a + Q3) * a + Q2) * a + Q1) * a + ONE)
            fn_val = -a * w
        else
            x = a - ONE
            w = (((((R5 * x + R4) * x + R3) * x + R2) * x + R1) * x + R0) / &
                (((((S5 * x + S4) * x + S3) * x + S2) * x + S1) * x + ONE)
            fn_val = x * w
        end if

    end function gamln1

!----------------------------------------------------------------------------------------
! FUNCTION: gamln
!
! DESCRIPTION: EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A
!              WRITTEN BY ALFRED H. MORRIS
!              NAVAL SURFACE WARFARE CENTER
!              DAHLGREN, VIRGINIA
!              Modified to interact with C using Fortran 2003
!              by Avraham Adler <Avraham.Adler@gmail.com>, 2016-11-07
!--------------------------
!     D = 0.5*(LN(2*PI) - 1)
!--------------------------

    elemental function gamln(a) result(fn_val)

    real(kind = c_double), intent(in) :: a
    real(kind = c_double)             :: fn_val, t, w
    real(kind = c_double), parameter  :: &
        C0 = 0.833333333333333e-01_c_double, &
        C1 = -.277777777760991e-02_c_double, &
        C2 = .793650666825390e-03_c_double,  &
        C3 = -.595202931351870e-03_c_double, &
        C4 = .837308034031215e-03_c_double,  &
        C5 = -.165322962780713e-02_c_double, &
        D = .418938533204673_c_double
    integer                           :: i, n

        if (a <= 0.8_c_double) then
            fn_val = gamln1(a) - log(a)
        else if(a <= 2.25_c_double) then
            fn_val = gamln1((a - HALF) - HALF)
        else if (a < 10._c_double) then
            n = int(a - 1.25_c_double)
            t = a
            w = ONE
            do  i = 1, n
                t = t - ONE
                w = t * w
            end do
            fn_val = gamln1(t - ONE) + log(w)
        else
            t = (ONE / a) ** 2
            w = (((((C5 * t + C4) * t + C3) * t + C2) * t + C1) * t + C0) / a
            fn_val = (D + w) + (a - HALF) * (LOG(a) - ONE)
        end if

    end function gamln

end module lgam
