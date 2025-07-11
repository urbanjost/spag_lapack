!*==slasq2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLASQ2 computes all the eigenvalues of the symmetric positive definite tridiagonal matrix associated with the qd Array Z to high relative accuracy. Used by sbdsqr and sstegr.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLASQ2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasq2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasq2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasq2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLASQ2( N, Z, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       REAL               Z( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLASQ2 computes all the eigenvalues of the symmetric positive
!> definite tridiagonal matrix associated with the qd array Z to high
!> relative accuracy are computed to high relative accuracy, in the
!> absence of denormalization, underflow and overflow.
!>
!> To see the relation of Z to the tridiagonal matrix, let L be a
!> unit lower bidiagonal matrix with subdiagonals Z(2,4,6,,..) and
!> let U be an upper bidiagonal matrix with 1's above and diagonal
!> Z(1,3,5,,..). The tridiagonal is L*U or, if you prefer, the
!> symmetric tridiagonal to which it is similar.
!>
!> Note : SLASQ2 defines a logical variable, IEEE, which is true
!> on machines which follow ieee-754 floating-point standard in their
!> handling of infinities and NaNs, and false otherwise. This variable
!> is passed to SLASQ3.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>        The number of rows and columns in the matrix. N >= 0.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is REAL array, dimension ( 4*N )
!>        On entry Z holds the qd array. On exit, entries 1 to N hold
!>        the eigenvalues in decreasing order, Z( 2*N+1 ) holds the
!>        trace, and Z( 2*N+2 ) holds the sum of the eigenvalues. If
!>        N > 2, then Z( 2*N+3 ) holds the iteration count, Z( 2*N+4 )
!>        holds NDIVS/NIN^2, and Z( 2*N+5 ) holds the percentage of
!>        shifts that failed.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>        = 0: successful exit
!>        < 0: if the i-th argument is a scalar and had an illegal
!>             value, then INFO = -i, if the i-th argument is an
!>             array and the j-entry had an illegal value, then
!>             INFO = -(i*100+j)
!>        > 0: the algorithm failed
!>              = 1, a split was marked by a positive value in E
!>              = 2, current block of Z not diagonalized after 100*N
!>                   iterations (in inner while loop).  On exit Z holds
!>                   a qd array with the same eigenvalues as the given Z.
!>              = 3, termination criterion of outer while loop not met
!>                   (program created more than N unreduced blocks)
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
!> \ingroup auxOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Local Variables: I0:N0 defines a current unreduced segment of Z.
!>  The shifts are accumulated in SIGMA. Iteration count is in ITER.
!>  Ping-pong is controlled by PP (alternates between 0 and 1).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE SLASQ2(N,Z,Info)
      IMPLICIT NONE
!*--SLASQ2116
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , N
!     ..
!     .. Array Arguments ..
      REAL Z(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL CBIAS
      PARAMETER (CBIAS=1.50E0)
      REAL ZERO , HALF , ONE , TWO , FOUR , HUNDRD
      PARAMETER (ZERO=0.0E0,HALF=0.5E0,ONE=1.0E0,TWO=2.0E0,FOUR=4.0E0,  &
     &           HUNDRD=100.0E0)
!     ..
!     .. Local Scalars ..
      LOGICAL ieee
      INTEGER i0 , i4 , iinfo , ipn4 , iter , iwhila , iwhilb , k ,     &
     &        kmin , n0 , nbig , ndiv , nfail , pp , splt , ttype , i1 ,&
     &        n1
      REAL d , dee , deemin , desig , dmin , dmin1 , dmin2 , dn , dn1 , &
     &     dn2 , e , emax , emin , eps , g , oldemn , qmax , qmin , s , &
     &     safmin , sigma , t , tau , temp , tol , tol2 , trace , zmax ,&
     &     tempe , tempq
!     ..
!     .. External Subroutines ..
      EXTERNAL SLASQ3 , SLASRT , XERBLA
!     ..
!     .. External Functions ..
      REAL SLAMCH
      EXTERNAL SLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , REAL , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!     (in case SLASQ2 is not called by SLASQ1)
!
      Info = 0
      eps = SLAMCH('Precision')
      safmin = SLAMCH('Safe minimum')
      tol = eps*HUNDRD
      tol2 = tol**2
!
      IF ( N<0 ) THEN
         Info = -1
         CALL XERBLA('SLASQ2',1)
         RETURN
      ELSEIF ( N==0 ) THEN
         RETURN
      ELSEIF ( N==1 ) THEN
!
!        1-by-1 case.
!
         IF ( Z(1)<ZERO ) THEN
            Info = -201
            CALL XERBLA('SLASQ2',2)
         ENDIF
         RETURN
      ELSEIF ( N==2 ) THEN
!
!        2-by-2 case.
!
         IF ( Z(1)<ZERO ) THEN
            Info = -201
            CALL XERBLA('DLASQ2',2)
            RETURN
         ELSEIF ( Z(2)<ZERO ) THEN
            Info = -202
            CALL XERBLA('SLASQ2',2)
            RETURN
         ELSEIF ( Z(3)<ZERO ) THEN
            Info = -203
            CALL XERBLA('SLASQ2',2)
            RETURN
         ELSEIF ( Z(3)>Z(1) ) THEN
            d = Z(3)
            Z(3) = Z(1)
            Z(1) = d
         ENDIF
         Z(5) = Z(1) + Z(2) + Z(3)
         IF ( Z(2)>Z(3)*tol2 ) THEN
            t = HALF*((Z(1)-Z(3))+Z(2))
            s = Z(3)*(Z(2)/t)
            IF ( s<=t ) THEN
               s = Z(3)*(Z(2)/(t*(ONE+SQRT(ONE+s/t))))
            ELSE
               s = Z(3)*(Z(2)/(t+SQRT(t)*SQRT(t+s)))
            ENDIF
            t = Z(1) + (s+Z(2))
            Z(3) = Z(3)*(Z(1)/t)
            Z(1) = t
         ENDIF
         Z(2) = Z(3)
         Z(6) = Z(2) + Z(1)
         RETURN
      ENDIF
!
!     Check for negative data and compute sums of q's and e's.
!
      Z(2*N) = ZERO
      emin = Z(2)
      qmax = ZERO
      zmax = ZERO
      d = ZERO
      e = ZERO
!
      DO k = 1 , 2*(N-1) , 2
         IF ( Z(k)<ZERO ) THEN
            Info = -(200+k)
            CALL XERBLA('SLASQ2',2)
            RETURN
         ELSEIF ( Z(k+1)<ZERO ) THEN
            Info = -(200+k+1)
            CALL XERBLA('SLASQ2',2)
            RETURN
         ENDIF
         d = d + Z(k)
         e = e + Z(k+1)
         qmax = MAX(qmax,Z(k))
         emin = MIN(emin,Z(k+1))
         zmax = MAX(qmax,zmax,Z(k+1))
      ENDDO
      IF ( Z(2*N-1)<ZERO ) THEN
         Info = -(200+2*N-1)
         CALL XERBLA('SLASQ2',2)
         RETURN
      ENDIF
      d = d + Z(2*N-1)
      qmax = MAX(qmax,Z(2*N-1))
      zmax = MAX(qmax,zmax)
!
!     Check for diagonality.
!
      IF ( e==ZERO ) THEN
         DO k = 2 , N
            Z(k) = Z(2*k-1)
         ENDDO
         CALL SLASRT('D',N,Z,iinfo)
         Z(2*N-1) = d
         RETURN
      ENDIF
!
      trace = d + e
!
!     Check for zero data.
!
      IF ( trace==ZERO ) THEN
         Z(2*N-1) = ZERO
         RETURN
      ENDIF
!
!     Check whether the machine is IEEE conformable.
!
!     IEEE = ( ILAENV( 10, 'SLASQ2', 'N', 1, 2, 3, 4 ).EQ.1 )
!
!     [11/15/2008] The case IEEE=.TRUE. has a problem in single precision with
!     some the test matrices of type 16. The double precision code is fine.
!
      ieee = .FALSE.
!
!     Rearrange data for locality: Z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...).
!
      DO k = 2*N , 2 , -2
         Z(2*k) = ZERO
         Z(2*k-1) = Z(k)
         Z(2*k-2) = ZERO
         Z(2*k-3) = Z(k-1)
      ENDDO
!
      i0 = 1
      n0 = N
!
!     Reverse the qd-array, if warranted.
!
      IF ( CBIAS*Z(4*i0-3)<Z(4*n0-3) ) THEN
         ipn4 = 4*(i0+n0)
         DO i4 = 4*i0 , 2*(i0+n0-1) , 4
            temp = Z(i4-3)
            Z(i4-3) = Z(ipn4-i4-3)
            Z(ipn4-i4-3) = temp
            temp = Z(i4-1)
            Z(i4-1) = Z(ipn4-i4-5)
            Z(ipn4-i4-5) = temp
         ENDDO
      ENDIF
!
!     Initial split checking via dqd and Li's test.
!
      pp = 0
!
      DO k = 1 , 2
!
         d = Z(4*n0+pp-3)
         DO i4 = 4*(n0-1) + pp , 4*i0 + pp , -4
            IF ( Z(i4-1)<=tol2*d ) THEN
               Z(i4-1) = -ZERO
               d = Z(i4-3)
            ELSE
               d = Z(i4-3)*(d/(d+Z(i4-1)))
            ENDIF
         ENDDO
!
!        dqd maps Z to ZZ plus Li's test.
!
         emin = Z(4*i0+pp+1)
         d = Z(4*i0+pp-3)
         DO i4 = 4*i0 + pp , 4*(n0-1) + pp , 4
            Z(i4-2*pp-2) = d + Z(i4-1)
            IF ( Z(i4-1)<=tol2*d ) THEN
               Z(i4-1) = -ZERO
               Z(i4-2*pp-2) = d
               Z(i4-2*pp) = ZERO
               d = Z(i4+1)
            ELSEIF ( safmin*Z(i4+1)<Z(i4-2*pp-2) .AND.                  &
     &               safmin*Z(i4-2*pp-2)<Z(i4+1) ) THEN
               temp = Z(i4+1)/Z(i4-2*pp-2)
               Z(i4-2*pp) = Z(i4-1)*temp
               d = d*temp
            ELSE
               Z(i4-2*pp) = Z(i4+1)*(Z(i4-1)/Z(i4-2*pp-2))
               d = Z(i4+1)*(d/Z(i4-2*pp-2))
            ENDIF
            emin = MIN(emin,Z(i4-2*pp))
         ENDDO
         Z(4*n0-pp-2) = d
!
!        Now find qmax.
!
         qmax = Z(4*i0-pp-2)
         DO i4 = 4*i0 - pp + 2 , 4*n0 - pp - 2 , 4
            qmax = MAX(qmax,Z(i4))
         ENDDO
!
!        Prepare for the next iteration on K.
!
         pp = 1 - pp
      ENDDO
!
!     Initialise variables to pass to SLASQ3.
!
      ttype = 0
      dmin1 = ZERO
      dmin2 = ZERO
      dn = ZERO
      dn1 = ZERO
      dn2 = ZERO
      g = ZERO
      tau = ZERO
!
      iter = 2
      nfail = 0
      ndiv = 2*(n0-i0)
!
      DO iwhila = 1 , N + 1
         IF ( n0<1 ) GOTO 200
!
!        While array unfinished do
!
!        E(N0) holds the value of SIGMA when submatrix in I0:N0
!        splits from the rest of the array, but is negated.
!
         desig = ZERO
         IF ( n0==N ) THEN
            sigma = ZERO
         ELSE
            sigma = -Z(4*n0-1)
         ENDIF
         IF ( sigma<ZERO ) THEN
            Info = 1
            RETURN
         ENDIF
!
!        Find last unreduced submatrix's top index I0, find QMAX and
!        EMIN. Find Gershgorin-type bound if Q's much greater than E's.
!
         emax = ZERO
         IF ( n0>i0 ) THEN
            emin = ABS(Z(4*n0-5))
         ELSE
            emin = ZERO
         ENDIF
         qmin = Z(4*n0-3)
         qmax = qmin
         DO i4 = 4*n0 , 8 , -4
            IF ( Z(i4-5)<=ZERO ) GOTO 50
            IF ( qmin>=FOUR*emax ) THEN
               qmin = MIN(qmin,Z(i4-3))
               emax = MAX(emax,Z(i4-5))
            ENDIF
            qmax = MAX(qmax,Z(i4-7)+Z(i4-5))
            emin = MIN(emin,Z(i4-5))
         ENDDO
         i4 = 4
!
 50      i0 = i4/4
         pp = 0
!
         IF ( n0-i0>1 ) THEN
            dee = Z(4*i0-3)
            deemin = dee
            kmin = i0
            DO i4 = 4*i0 + 1 , 4*n0 - 3 , 4
               dee = Z(i4)*(dee/(dee+Z(i4-2)))
               IF ( dee<=deemin ) THEN
                  deemin = dee
                  kmin = (i4+3)/4
               ENDIF
            ENDDO
            IF ( (kmin-i0)*2<n0-kmin .AND. deemin<=HALF*Z(4*n0-3) ) THEN
               ipn4 = 4*(i0+n0)
               pp = 2
               DO i4 = 4*i0 , 2*(i0+n0-1) , 4
                  temp = Z(i4-3)
                  Z(i4-3) = Z(ipn4-i4-3)
                  Z(ipn4-i4-3) = temp
                  temp = Z(i4-2)
                  Z(i4-2) = Z(ipn4-i4-2)
                  Z(ipn4-i4-2) = temp
                  temp = Z(i4-1)
                  Z(i4-1) = Z(ipn4-i4-5)
                  Z(ipn4-i4-5) = temp
                  temp = Z(i4)
                  Z(i4) = Z(ipn4-i4-4)
                  Z(ipn4-i4-4) = temp
               ENDDO
            ENDIF
         ENDIF
!
!        Put -(initial shift) into DMIN.
!
         dmin = -MAX(ZERO,qmin-TWO*SQRT(qmin)*SQRT(emax))
!
!        Now I0:N0 is unreduced.
!        PP = 0 for ping, PP = 1 for pong.
!        PP = 2 indicates that flipping was applied to the Z array and
!               and that the tests for deflation upon entry in SLASQ3
!               should not be performed.
!
         nbig = 100*(n0-i0+1)
         DO iwhilb = 1 , nbig
            IF ( i0>n0 ) GOTO 100
!
!           While submatrix unfinished take a good dqds step.
!
            CALL SLASQ3(i0,n0,Z,pp,dmin,sigma,desig,qmax,nfail,iter,    &
     &                  ndiv,ieee,ttype,dmin1,dmin2,dn,dn1,dn2,g,tau)
!
            pp = 1 - pp
!
!           When EMIN is very small check for splits.
!
            IF ( pp==0 .AND. n0-i0>=3 ) THEN
               IF ( Z(4*n0)<=tol2*qmax .OR. Z(4*n0-1)<=tol2*sigma ) THEN
                  splt = i0 - 1
                  qmax = Z(4*i0-3)
                  emin = Z(4*i0-1)
                  oldemn = Z(4*i0)
                  DO i4 = 4*i0 , 4*(n0-3) , 4
                     IF ( Z(i4)<=tol2*Z(i4-3) .OR. Z(i4-1)<=tol2*sigma )&
     &                    THEN
                        Z(i4-1) = -sigma
                        splt = i4/4
                        qmax = ZERO
                        emin = Z(i4+3)
                        oldemn = Z(i4+4)
                     ELSE
                        qmax = MAX(qmax,Z(i4+1))
                        emin = MIN(emin,Z(i4-1))
                        oldemn = MIN(oldemn,Z(i4))
                     ENDIF
                  ENDDO
                  Z(4*n0-1) = emin
                  Z(4*n0) = oldemn
                  i0 = splt + 1
               ENDIF
            ENDIF
!
         ENDDO
!
         Info = 2
!
!        Maximum number of iterations exceeded, restore the shift
!        SIGMA and place the new d's and e's in a qd array.
!        This might need to be done for several blocks
!
         i1 = i0
         n1 = n0
         DO
            tempq = Z(4*i0-3)
            Z(4*i0-3) = Z(4*i0-3) + sigma
            DO k = i0 + 1 , n0
               tempe = Z(4*k-5)
               Z(4*k-5) = Z(4*k-5)*(tempq/Z(4*k-7))
               tempq = Z(4*k-3)
               Z(4*k-3) = Z(4*k-3) + sigma + tempe - Z(4*k-5)
            ENDDO
!
!        Prepare to do this on the previous block if there is one
!
            IF ( i1>1 ) THEN
               n1 = i1 - 1
               DO WHILE ( (i1>=2) .AND. (Z(4*i1-5)>=ZERO) )
                  i1 = i1 - 1
               ENDDO
               IF ( i1>=1 ) THEN
                  sigma = -Z(4*n1-1)
                  CYCLE
               ENDIF
            ENDIF
 
            DO k = 1 , N
               Z(2*k-1) = Z(4*k-3)
!
!        Only the block 1..N0 is unfinished.  The rest of the e's
!        must be essentially zero, although sometimes other data
!        has been stored in them.
!
               IF ( k<n0 ) THEN
                  Z(2*k) = Z(4*k-1)
               ELSE
                  Z(2*k) = 0
               ENDIF
            ENDDO
            RETURN
         ENDDO
!
!        end IWHILB
!
!
 100  ENDDO
!
      Info = 3
      RETURN
!
!     end IWHILA
!
!
!     Move q's to the front.
!
 200  DO k = 2 , N
         Z(k) = Z(4*k-3)
      ENDDO
!
!     Sort and compute sum of eigenvalues.
!
      CALL SLASRT('D',N,Z,iinfo)
!
      e = ZERO
      DO k = N , 1 , -1
         e = e + Z(k)
      ENDDO
!
!     Store trace, sum(eigenvalues) and information on performance.
!
      Z(2*N+1) = trace
      Z(2*N+2) = e
      Z(2*N+3) = REAL(iter)
      Z(2*N+4) = REAL(ndiv)/REAL(N**2)
      Z(2*N+5) = HUNDRD*nfail/REAL(iter)
!
!     End of SLASQ2
!
      END SUBROUTINE SLASQ2
