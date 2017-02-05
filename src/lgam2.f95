!----------------------------------------------------------------------------------------
!
! MODULE: lgam2
!
! AUTHOR: W. J. Cody and L. Stoltz: Argonne National Laboratory
!         Modified for Fortran 2003 and use with C by
!         Avraham Adler <Avraham.Adler@gmail.com>
!
! DESCRIPTION: LogGamma function taken from http://www.netlib.org/specfun/algama
!              
! HISTORY:
!          Version 1.0: 2017-02-04
!                       Refactored Netlib code to be C-double precision and Fortran 2003.
!
! LICENSE: The original is was freely donated to the public domain.
!          The updated code below is released under the BSD-2 License below:
!
!   Copyright (c) 2017, Avraham Adler
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
!   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN if ADVISED OF THE POSSIBILITY OF SUCH
!   DAMAGE.
!----------------------------------------------------------------------------------------

!----------------------------------------------------------------------
! {ORIGINAL NETLIB NOTES - AA}
!
! This routine calculates the LOG(GAMMA) function for a positive real
!   argument X.  Computation is based on an algorithm outlined in
!   references 1 and 2.  The program uses rational functions that
!   theoretically approximate LOG(GAMMA) to at least 18 significant
!   decimal digits.  The approximation for X > 12 is from reference
!   3, while approximations for X < 12.0 are similar to those in
!   reference 1, but are unpublished.  The accuracy achieved depends
!   on the arithmetic system, the compiler, the intrinsic functions,
!   and proper selection of the machine-dependent constants.
!
!
!*********************************************************************
!*********************************************************************
!
! Explanation of machine-dependent constants
!
! beta   - radix for the floating-point representation
! maxexp - the smallest positive power of beta that overflows
! XBIG   - largest argument for which LN(GAMMA(X)) is representable
!          in the machine, i.e., the solution to the equation
!                  LN(GAMMA(XBIG)) = beta**maxexp
! XINF   - largest machine representable floating-point number;
!          approximately beta**maxexp.
! EPS    - The smallest positive floating-point number such that
!          1.0+EPS .GT. 1.0
! FRTBIG - Rough estimate of the fourth root of XBIG
!
!
!     Approximate values for some important machines are:
!
!                           beta      maxexp         XBIG
!
! CRAY-1        (S.P.)        2        8191       9.62E+2461
! Cyber 180/855
!   under NOS   (S.P.)        2        1070       1.72E+319
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)        2         128       4.08E+36
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)        2        1024       2.55D+305
! IBM 3033      (D.P.)       16          63       4.29D+73
! VAX D-Format  (D.P.)        2         127       2.05D+36
! VAX G-Format  (D.P.)        2        1023       1.28D+305
!
!
!                           XINF        EPS        FRTBIG
!
! CRAY-1        (S.P.)   5.45E+2465   7.11E-15    3.13E+615
! Cyber 180/855
!   under NOS   (S.P.)   1.26E+322    3.55E-15    6.44E+79
! IEEE (IBM/XT,
!   SUN, etc.)  (S.P.)   3.40E+38     1.19E-7     1.42E+9
! IEEE (IBM/XT,
!   SUN, etc.)  (D.P.)   1.79D+308    2.22D-16    2.25D+76
! IBM 3033      (D.P.)   7.23D+75     2.22D-16    2.56D+18
! VAX D-Format  (D.P.)   1.70D+38     1.39D-17    1.20D+9
! VAX G-Format  (D.P.)   8.98D+307    1.11D-16    1.89D+76
!
!**************************************************************
!**************************************************************
!
! Error returns
!
!  The program returns the value XINF for X .LE. 0.0 or when
!     overflow would occur.  The computation is believed to 
!     be free of underflow and overflow.
!
!
! Intrinsic functions required are:
!
!      LOG
!
!
! References:
!
!  1) W. J. Cody and K. E. Hillstrom, 'Chebyshev Approximations for
!     the Natural Logarithm of the Gamma Function,' Math. Comp. 21,
!     1967, pp. 198-203.
!
!  2) K. E. Hillstrom, ANL/AMD Program ANLC366S, DGAMMA/DLGAMA, May,
!     1969.
! 
!  3) Hart, Et. Al., Computer Approximations, Wiley and sons, New
!     York, 1968.
!
!
!  Authors: W. J. Cody and L. Stoltz
!           Argonne National Laboratory
!
!  Latest modification: June 16, 1988
!
!----------------------------------------------------------------------------------------

module lgam2

  use, intrinsic :: iso_c_binding, only: c_double
  use utils
  implicit none
  
  contains
  
  elemental function lgamma(x) result(res)

    integer :: i
    real(kind = c_double), intent(in) :: x
    real(kind = c_double)             :: res
    real(kind = c_double) :: corr, xden, xm1, xm2, xm4, xnum, y, ysq
    real(kind = c_double), parameter :: PNT68 = 6796875_c_double
    real(kind = c_double), parameter :: FRTBIG = 2.25e76_c_double

! Machine dependant parameters. Using IEEE
    real(kind = c_double), parameter :: XBIG = 2.55e305_c_double
    real(kind = c_double), parameter :: XINF = 1.797693e+308_c_double
    
! Numerator & denominator coefficients for rational minimax approximation over (0.5, 1.5)
!----------------------------------------------------------------------------------------
    real(kind = c_double), parameter :: D1 = -5.772156649015328605195174e-1_c_double
    real(kind = c_double), parameter, dimension(8) :: P1 = &
         [4.945235359296727046734888_c_double, &
          2.018112620856775083915565e2_c_double, &
          2.290838373831346393026739e3_c_double, &
          1.131967205903380828685045e4_c_double, &
          2.855724635671635335736389e4_c_double, &
          3.848496228443793359990269e4_c_double, &
          2.637748787624195437963534e4_c_double, &
          7.225813979700288197698961e3_c_double]
    real(kind = c_double), parameter, dimension(8) :: Q1 = &
         [6.748212550303777196073036e1_c_double, &
          1.113332393857199323513008e3_c_double, &
          7.738757056935398733233834e3_c_double, &
          2.763987074403340708898585e4_c_double, &
          5.499310206226157329794414e4_c_double, &
          6.161122180066002127833352e4_c_double, &
          3.635127591501940507276287e4_c_double, &
          8.785536302431013170870835e3_c_double]
  
! Numerator & denominator coefficients for rational minimax approximation over (1.5,4.0)
!----------------------------------------------------------------------------------------
    real(kind = c_double), parameter :: D2 = 4.227843350984671393993777e-1_c_double
    real(kind = c_double), parameter, dimension(8) :: P2 = &
         [4.974607845568932035012064e0_c_double, &
          5.424138599891070494101986e2_c_double, &
          1.550693864978364947665077e4_c_double, &
          1.847932904445632425417223e5_c_double, &
          1.088204769468828767498470e6_c_double, &
          3.338152967987029735917223e6_c_double, &
          5.106661678927352456275255e6_c_double, &
          3.074109054850539556250927e6_c_double]
    real(kind = c_double), parameter, dimension(8) :: Q2 = &
         [1.830328399370592604055942e2_c_double, &
          7.765049321445005871323047e3_c_double, &
          1.331903827966074194402448e5_c_double, &
          1.136705821321969608938755e6_c_double, &
          5.267964117437946917577538e6_c_double, &
          1.346701454311101692290052e7_c_double, &
          1.782736530353274213975932e7_c_double, &
          9.533095591844353613395747e6_c_double]

! Numerator & denominator coefficients for rational minimax approximation over (4.0,12.0)
!----------------------------------------------------------------------------------------
    real(kind = c_double), parameter :: D4 = 1.791759469228055000094023e0_c_double
    real(kind = c_double), parameter, dimension(8) :: P4 = &
         [1.474502166059939948905062e4_c_double, &
          2.426813369486704502836312e6_c_double, &
          1.214755574045093227939592e8_c_double, &
          2.663432449630976949898078e9_c_double, &
          2.940378956634553899906876e10_c_double, &
          1.702665737765398868392998e11_c_double, &
          4.926125793377430887588120e11_c_double, &
          5.606251856223951465078242e11_c_double]
    real(kind = c_double), parameter, dimension(8) :: Q4 = &
         [2.690530175870899333379843e3_c_double, &
          6.393885654300092398984238e5_c_double, &
          4.135599930241388052042842e7_c_double, &
          1.120872109616147941376570e9_c_double, &
          1.488613728678813811542398e10_c_double, &
          1.016803586272438228077304e11_c_double, &
          3.417476345507377132798597e11_c_double, &
          4.463158187419713286462081e11_c_double]
    
!  Coefficients for minimax approximation over (12, INF).
!----------------------------------------------------------------------------------------
    real(kind = c_double), parameter, dimension(7) :: C = &
             [-1.910444077728e-03_c_double, &
              8.4171387781295e-04_c_double, &
              -5.952379913043012e-04_c_double, &
              7.93650793500350248e-04_c_double, &
              -2.777777777777681622553e-03_c_double, &
              8.333333333333333331554247e-02_c_double, &
              5.7083835261e-03_c_double]

! Begin calculations
!----------------------------------------------------------------------------------------
      y = x
      if ((y > ZERO) .and. (y <= XBIG)) then
        if (y <= EPS) then
          res = -log(y)
        else if (y < THRHAL) then
!----------------------------------------------------------------------
!  EPS < X <= 1.5
!----------------------------------------------------------------------
          if (y < PNT68) then
            corr = -log(y)
            xm1 = y
          else
            corr = ZERO
            xm1 = (y - HALF) - HALF
          end if
          if ((y <= HALF) .or. (y >= PNT68)) then
            xden = ONE
            xnum = ZERO
            do i = 1, 8
              xnum = xnum * xm1 + P1(i)
              xden = xden * xm1 + Q1(i)
            end do
            res = corr + (xm1 * (D1 + xm1 * (xnum / xden)))
          else
            xm2 = (y - HALF) - HALF
            xden = ONE
            xnum = ZERO
            do i = 1, 8
              xnum = xnum * xm2 + P2(i)
              xden = xden * xm2 + Q2(i)
            end do
            res = corr + xm2 * (D2 + xm2 * (xnum / xden))
          end if
        else if (y <= FOUR) then
!----------------------------------------------------------------------
!  1.5 < X <= 4.0
!----------------------------------------------------------------------
          xm2 = y - TWO
          xden = ONE
          xnum = ZERO
          do i = 1, 8
            xnum = xnum * xm2 + P2(i)
            xden = xden * xm2 + Q2(i)
          end do
          res = xm2 * (D2 + xm2 * (xnum / xden))
        else if (y <= TWELVE) then
!----------------------------------------------------------------------
!  4.0 < X <= 12.0
!----------------------------------------------------------------------
          xm4 = y - FOUR
          xden = -ONE
          xnum = ZERO
          do i = 1, 8
            xnum = xnum * xm4 + P4(i)
            xden = xden * xm4 + Q4(i)
          end do
          res = D4 + xm4 * (xnum / xden)
        else 
!----------------------------------------------------------------------
!  Evaluate for argument >= 12.0
!----------------------------------------------------------------------
          res = ZERO
          if (y <= FRTBIG) then
            res = C(7)
            ysq = y * y
            do i = 1, 6
              res = res / ysq + C(i)
            end do
          end if
          res = res / y
          corr = log(y)
          res = res + SQRTPI - HALF * corr
          res = res + y * (corr - ONE)
        end if
      else
!----------------------------------------------------------------------
!  Return for bad arguments
!----------------------------------------------------------------------
        res = XINF
      end if
      
  end function lgamma

end module lgam2