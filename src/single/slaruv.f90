!*==slaruv.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLARUV returns a vector of n random real numbers from a uniform distribution.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLARUV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaruv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaruv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaruv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLARUV( ISEED, N, X )
!
!       .. Scalar Arguments ..
!       INTEGER            N
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       REAL               X( N )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLARUV returns a vector of n random real numbers from a uniform (0,1)
!> distribution (n <= 128).
!>
!> This is an auxiliary routine called by SLARNV and CLARNV.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          On entry, the seed of the random number generator; the array
!>          elements must be between 0 and 4095, and ISEED(4) must be
!>          odd.
!>          On exit, the seed is updated.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of random numbers to be generated. N <= 128.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is REAL array, dimension (N)
!>          The generated random numbers.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date December 2016
!
!> \ingroup OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  This routine uses a multiplicative congruential method with modulus
!>  2**48 and multiplier 33952834046453 (see G.S.Fishman,
!>  'Multiplicative congruential random number generators with modulus
!>  2**b: an exhaustive analysis for b = 32 and a partial analysis for
!>  b = 48', Math. Comp. 189, pp 331-344, 1990).
!>
!>  48-bit integers are stored in 4 integer array elements with 12 bits
!>  per element. Hence the routine is portable across machines with
!>  integers of 32 bits or more.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE SLARUV(Iseed,N,X)
      IMPLICIT NONE
!*--SLARUV99
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER N
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4)
      REAL X(N)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE
      PARAMETER (ONE=1.0E0)
      INTEGER LV , IPW2
      REAL R
      PARAMETER (LV=128,IPW2=4096,R=ONE/IPW2)
!     ..
!     .. Local Scalars ..
      INTEGER i , i1 , i2 , i3 , i4 , it1 , it2 , it3 , it4 , j
!     ..
!     .. Local Arrays ..
      INTEGER mm(LV,4)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MIN , MOD , REAL
!     ..
!     .. Data statements ..
      DATA (mm(1,j),j=1,4)/494 , 322 , 2508 , 2549/
      DATA (mm(2,j),j=1,4)/2637 , 789 , 3754 , 1145/
      DATA (mm(3,j),j=1,4)/255 , 1440 , 1766 , 2253/
      DATA (mm(4,j),j=1,4)/2008 , 752 , 3572 , 305/
      DATA (mm(5,j),j=1,4)/1253 , 2859 , 2893 , 3301/
      DATA (mm(6,j),j=1,4)/3344 , 123 , 307 , 1065/
      DATA (mm(7,j),j=1,4)/4084 , 1848 , 1297 , 3133/
      DATA (mm(8,j),j=1,4)/1739 , 643 , 3966 , 2913/
      DATA (mm(9,j),j=1,4)/3143 , 2405 , 758 , 3285/
      DATA (mm(10,j),j=1,4)/3468 , 2638 , 2598 , 1241/
      DATA (mm(11,j),j=1,4)/688 , 2344 , 3406 , 1197/
      DATA (mm(12,j),j=1,4)/1657 , 46 , 2922 , 3729/
      DATA (mm(13,j),j=1,4)/1238 , 3814 , 1038 , 2501/
      DATA (mm(14,j),j=1,4)/3166 , 913 , 2934 , 1673/
      DATA (mm(15,j),j=1,4)/1292 , 3649 , 2091 , 541/
      DATA (mm(16,j),j=1,4)/3422 , 339 , 2451 , 2753/
      DATA (mm(17,j),j=1,4)/1270 , 3808 , 1580 , 949/
      DATA (mm(18,j),j=1,4)/2016 , 822 , 1958 , 2361/
      DATA (mm(19,j),j=1,4)/154 , 2832 , 2055 , 1165/
      DATA (mm(20,j),j=1,4)/2862 , 3078 , 1507 , 4081/
      DATA (mm(21,j),j=1,4)/697 , 3633 , 1078 , 2725/
      DATA (mm(22,j),j=1,4)/1706 , 2970 , 3273 , 3305/
      DATA (mm(23,j),j=1,4)/491 , 637 , 17 , 3069/
      DATA (mm(24,j),j=1,4)/931 , 2249 , 854 , 3617/
      DATA (mm(25,j),j=1,4)/1444 , 2081 , 2916 , 3733/
      DATA (mm(26,j),j=1,4)/444 , 4019 , 3971 , 409/
      DATA (mm(27,j),j=1,4)/3577 , 1478 , 2889 , 2157/
      DATA (mm(28,j),j=1,4)/3944 , 242 , 3831 , 1361/
      DATA (mm(29,j),j=1,4)/2184 , 481 , 2621 , 3973/
      DATA (mm(30,j),j=1,4)/1661 , 2075 , 1541 , 1865/
      DATA (mm(31,j),j=1,4)/3482 , 4058 , 893 , 2525/
      DATA (mm(32,j),j=1,4)/657 , 622 , 736 , 1409/
      DATA (mm(33,j),j=1,4)/3023 , 3376 , 3992 , 3445/
      DATA (mm(34,j),j=1,4)/3618 , 812 , 787 , 3577/
      DATA (mm(35,j),j=1,4)/1267 , 234 , 2125 , 77/
      DATA (mm(36,j),j=1,4)/1828 , 641 , 2364 , 3761/
      DATA (mm(37,j),j=1,4)/164 , 4005 , 2460 , 2149/
      DATA (mm(38,j),j=1,4)/3798 , 1122 , 257 , 1449/
      DATA (mm(39,j),j=1,4)/3087 , 3135 , 1574 , 3005/
      DATA (mm(40,j),j=1,4)/2400 , 2640 , 3912 , 225/
      DATA (mm(41,j),j=1,4)/2870 , 2302 , 1216 , 85/
      DATA (mm(42,j),j=1,4)/3876 , 40 , 3248 , 3673/
      DATA (mm(43,j),j=1,4)/1905 , 1832 , 3401 , 3117/
      DATA (mm(44,j),j=1,4)/1593 , 2247 , 2124 , 3089/
      DATA (mm(45,j),j=1,4)/1797 , 2034 , 2762 , 1349/
      DATA (mm(46,j),j=1,4)/1234 , 2637 , 149 , 2057/
      DATA (mm(47,j),j=1,4)/3460 , 1287 , 2245 , 413/
      DATA (mm(48,j),j=1,4)/328 , 1691 , 166 , 65/
      DATA (mm(49,j),j=1,4)/2861 , 496 , 466 , 1845/
      DATA (mm(50,j),j=1,4)/1950 , 1597 , 4018 , 697/
      DATA (mm(51,j),j=1,4)/617 , 2394 , 1399 , 3085/
      DATA (mm(52,j),j=1,4)/2070 , 2584 , 190 , 3441/
      DATA (mm(53,j),j=1,4)/3331 , 1843 , 2879 , 1573/
      DATA (mm(54,j),j=1,4)/769 , 336 , 153 , 3689/
      DATA (mm(55,j),j=1,4)/1558 , 1472 , 2320 , 2941/
      DATA (mm(56,j),j=1,4)/2412 , 2407 , 18 , 929/
      DATA (mm(57,j),j=1,4)/2800 , 433 , 712 , 533/
      DATA (mm(58,j),j=1,4)/189 , 2096 , 2159 , 2841/
      DATA (mm(59,j),j=1,4)/287 , 1761 , 2318 , 4077/
      DATA (mm(60,j),j=1,4)/2045 , 2810 , 2091 , 721/
      DATA (mm(61,j),j=1,4)/1227 , 566 , 3443 , 2821/
      DATA (mm(62,j),j=1,4)/2838 , 442 , 1510 , 2249/
      DATA (mm(63,j),j=1,4)/209 , 41 , 449 , 2397/
      DATA (mm(64,j),j=1,4)/2770 , 1238 , 1956 , 2817/
      DATA (mm(65,j),j=1,4)/3654 , 1086 , 2201 , 245/
      DATA (mm(66,j),j=1,4)/3993 , 603 , 3137 , 1913/
      DATA (mm(67,j),j=1,4)/192 , 840 , 3399 , 1997/
      DATA (mm(68,j),j=1,4)/2253 , 3168 , 1321 , 3121/
      DATA (mm(69,j),j=1,4)/3491 , 1499 , 2271 , 997/
      DATA (mm(70,j),j=1,4)/2889 , 1084 , 3667 , 1833/
      DATA (mm(71,j),j=1,4)/2857 , 3438 , 2703 , 2877/
      DATA (mm(72,j),j=1,4)/2094 , 2408 , 629 , 1633/
      DATA (mm(73,j),j=1,4)/1818 , 1589 , 2365 , 981/
      DATA (mm(74,j),j=1,4)/688 , 2391 , 2431 , 2009/
      DATA (mm(75,j),j=1,4)/1407 , 288 , 1113 , 941/
      DATA (mm(76,j),j=1,4)/634 , 26 , 3922 , 2449/
      DATA (mm(77,j),j=1,4)/3231 , 512 , 2554 , 197/
      DATA (mm(78,j),j=1,4)/815 , 1456 , 184 , 2441/
      DATA (mm(79,j),j=1,4)/3524 , 171 , 2099 , 285/
      DATA (mm(80,j),j=1,4)/1914 , 1677 , 3228 , 1473/
      DATA (mm(81,j),j=1,4)/516 , 2657 , 4012 , 2741/
      DATA (mm(82,j),j=1,4)/164 , 2270 , 1921 , 3129/
      DATA (mm(83,j),j=1,4)/303 , 2587 , 3452 , 909/
      DATA (mm(84,j),j=1,4)/2144 , 2961 , 3901 , 2801/
      DATA (mm(85,j),j=1,4)/3480 , 1970 , 572 , 421/
      DATA (mm(86,j),j=1,4)/119 , 1817 , 3309 , 4073/
      DATA (mm(87,j),j=1,4)/3357 , 676 , 3171 , 2813/
      DATA (mm(88,j),j=1,4)/837 , 1410 , 817 , 2337/
      DATA (mm(89,j),j=1,4)/2826 , 3723 , 3039 , 1429/
      DATA (mm(90,j),j=1,4)/2332 , 2803 , 1696 , 1177/
      DATA (mm(91,j),j=1,4)/2089 , 3185 , 1256 , 1901/
      DATA (mm(92,j),j=1,4)/3780 , 184 , 3715 , 81/
      DATA (mm(93,j),j=1,4)/1700 , 663 , 2077 , 1669/
      DATA (mm(94,j),j=1,4)/3712 , 499 , 3019 , 2633/
      DATA (mm(95,j),j=1,4)/150 , 3784 , 1497 , 2269/
      DATA (mm(96,j),j=1,4)/2000 , 1631 , 1101 , 129/
      DATA (mm(97,j),j=1,4)/3375 , 1925 , 717 , 1141/
      DATA (mm(98,j),j=1,4)/1621 , 3912 , 51 , 249/
      DATA (mm(99,j),j=1,4)/3090 , 1398 , 981 , 3917/
      DATA (mm(100,j),j=1,4)/3765 , 1349 , 1978 , 2481/
      DATA (mm(101,j),j=1,4)/1149 , 1441 , 1813 , 3941/
      DATA (mm(102,j),j=1,4)/3146 , 2224 , 3881 , 2217/
      DATA (mm(103,j),j=1,4)/33 , 2411 , 76 , 2749/
      DATA (mm(104,j),j=1,4)/3082 , 1907 , 3846 , 3041/
      DATA (mm(105,j),j=1,4)/2741 , 3192 , 3694 , 1877/
      DATA (mm(106,j),j=1,4)/359 , 2786 , 1682 , 345/
      DATA (mm(107,j),j=1,4)/3316 , 382 , 124 , 2861/
      DATA (mm(108,j),j=1,4)/1749 , 37 , 1660 , 1809/
      DATA (mm(109,j),j=1,4)/185 , 759 , 3997 , 3141/
      DATA (mm(110,j),j=1,4)/2784 , 2948 , 479 , 2825/
      DATA (mm(111,j),j=1,4)/2202 , 1862 , 1141 , 157/
      DATA (mm(112,j),j=1,4)/2199 , 3802 , 886 , 2881/
      DATA (mm(113,j),j=1,4)/1364 , 2423 , 3514 , 3637/
      DATA (mm(114,j),j=1,4)/1244 , 2051 , 1301 , 1465/
      DATA (mm(115,j),j=1,4)/2020 , 2295 , 3604 , 2829/
      DATA (mm(116,j),j=1,4)/3160 , 1332 , 1888 , 2161/
      DATA (mm(117,j),j=1,4)/2785 , 1832 , 1836 , 3365/
      DATA (mm(118,j),j=1,4)/2772 , 2405 , 1990 , 361/
      DATA (mm(119,j),j=1,4)/1217 , 3638 , 2058 , 2685/
      DATA (mm(120,j),j=1,4)/1822 , 3661 , 692 , 3745/
      DATA (mm(121,j),j=1,4)/1245 , 327 , 1194 , 2325/
      DATA (mm(122,j),j=1,4)/2252 , 3660 , 20 , 3609/
      DATA (mm(123,j),j=1,4)/3904 , 716 , 3285 , 3821/
      DATA (mm(124,j),j=1,4)/2774 , 1842 , 2046 , 3537/
      DATA (mm(125,j),j=1,4)/997 , 3987 , 2107 , 517/
      DATA (mm(126,j),j=1,4)/2573 , 1368 , 3508 , 3017/
      DATA (mm(127,j),j=1,4)/1148 , 1848 , 3525 , 2141/
      DATA (mm(128,j),j=1,4)/545 , 2366 , 3801 , 1537/
!     ..
!     .. Executable Statements ..
!
      i1 = Iseed(1)
      i2 = Iseed(2)
      i3 = Iseed(3)
      i4 = Iseed(4)
!
      DO i = 1 , MIN(N,LV)
         DO
!
!
!        Multiply the seed by i-th power of the multiplier modulo 2**48
!
            it4 = i4*mm(i,4)
            it3 = it4/IPW2
            it4 = it4 - IPW2*it3
            it3 = it3 + i3*mm(i,4) + i4*mm(i,3)
            it2 = it3/IPW2
            it3 = it3 - IPW2*it2
            it2 = it2 + i2*mm(i,4) + i3*mm(i,3) + i4*mm(i,2)
            it1 = it2/IPW2
            it2 = it2 - IPW2*it1
            it1 = it1 + i1*mm(i,4) + i2*mm(i,3) + i3*mm(i,2)            &
     &            + i4*mm(i,1)
            it1 = MOD(it1,IPW2)
!
!        Convert 48-bit integer to a real number in the interval (0,1)
!
            X(i) = R*(REAL(it1)+R*(REAL(it2)+R*(REAL(it3)+R*REAL(it4))))
!
            IF ( X(i)==1.0 ) THEN
!           If a real number has n bits of precision, and the first
!           n bits of the 48-bit integer above happen to be all 1 (which
!           will occur about once every 2**n calls), then X( I ) will
!           be rounded to exactly 1.0. In IEEE single precision arithmetic,
!           this will happen relatively often since n = 24.
!           Since X( I ) is not supposed to return exactly 0.0 or 1.0,
!           the statistically correct thing to do in this situation is
!           simply to iterate again.
!           N.B. the case X( I ) = 0.0 should not be possible.
               i1 = i1 + 2
               i2 = i2 + 2
               i3 = i3 + 2
               i4 = i4 + 2
               CYCLE
            ENDIF
            EXIT
         ENDDO
!
      ENDDO
!
!     Return final value of seed
!
      Iseed(1) = it1
      Iseed(2) = it2
      Iseed(3) = it3
      Iseed(4) = it4
!
!     End of SLARUV
!
      END SUBROUTINE SLARUV
