!*==slasd3.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLASD3 finds all square roots of the roots of the secular equation, as defined by the values in D and Z, and then updates the singular vectors by matrix multiplication. Used by sbdsdc.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLASD3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLASD3( NL, NR, SQRE, K, D, Q, LDQ, DSIGMA, U, LDU, U2,
!                          LDU2, VT, LDVT, VT2, LDVT2, IDXC, CTOT, Z,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDQ, LDU, LDU2, LDVT, LDVT2, NL, NR,
!      $                   SQRE
!       ..
!       .. Array Arguments ..
!       INTEGER            CTOT( * ), IDXC( * )
!       REAL               D( * ), DSIGMA( * ), Q( LDQ, * ), U( LDU, * ),
!      $                   U2( LDU2, * ), VT( LDVT, * ), VT2( LDVT2, * ),
!      $                   Z( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLASD3 finds all the square roots of the roots of the secular
!> equation, as defined by the values in D and Z.  It makes the
!> appropriate calls to SLASD4 and then updates the singular
!> vectors by matrix multiplication.
!>
!> This code makes very mild assumptions about floating point
!> arithmetic. It will work on machines with a guard digit in
!> add/subtract, or on those binary machines without guard digits
!> which subtract like the Cray XMP, Cray YMP, Cray C 90, or Cray 2.
!> It could conceivably fail on hexadecimal or decimal machines
!> without guard digits, but we know of none.
!>
!> SLASD3 is called from SLASD1.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NL
!> \verbatim
!>          NL is INTEGER
!>         The row dimension of the upper block.  NL >= 1.
!> \endverbatim
!>
!> \param[in] NR
!> \verbatim
!>          NR is INTEGER
!>         The row dimension of the lower block.  NR >= 1.
!> \endverbatim
!>
!> \param[in] SQRE
!> \verbatim
!>          SQRE is INTEGER
!>         = 0: the lower block is an NR-by-NR square matrix.
!>         = 1: the lower block is an NR-by-(NR+1) rectangular matrix.
!>
!>         The bidiagonal matrix has N = NL + NR + 1 rows and
!>         M = N + SQRE >= N columns.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>         The size of the secular equation, 1 =< K = < N.
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is REAL array, dimension(K)
!>         On exit the square roots of the roots of the secular equation,
!>         in ascending order.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is REAL array, dimension (LDQ,K)
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>         The leading dimension of the array Q.  LDQ >= K.
!> \endverbatim
!>
!> \param[in,out] DSIGMA
!> \verbatim
!>          DSIGMA is REAL array, dimension(K)
!>         The first K elements of this array contain the old roots
!>         of the deflated updating problem.  These are the poles
!>         of the secular equation.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is REAL array, dimension (LDU, N)
!>         The last N - K columns of this matrix contain the deflated
!>         left singular vectors.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>         The leading dimension of the array U.  LDU >= N.
!> \endverbatim
!>
!> \param[in] U2
!> \verbatim
!>          U2 is REAL array, dimension (LDU2, N)
!>         The first K columns of this matrix contain the non-deflated
!>         left singular vectors for the split problem.
!> \endverbatim
!>
!> \param[in] LDU2
!> \verbatim
!>          LDU2 is INTEGER
!>         The leading dimension of the array U2.  LDU2 >= N.
!> \endverbatim
!>
!> \param[out] VT
!> \verbatim
!>          VT is REAL array, dimension (LDVT, M)
!>         The last M - K columns of VT**T contain the deflated
!>         right singular vectors.
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER
!>         The leading dimension of the array VT.  LDVT >= N.
!> \endverbatim
!>
!> \param[in,out] VT2
!> \verbatim
!>          VT2 is REAL array, dimension (LDVT2, N)
!>         The first K columns of VT2**T contain the non-deflated
!>         right singular vectors for the split problem.
!> \endverbatim
!>
!> \param[in] LDVT2
!> \verbatim
!>          LDVT2 is INTEGER
!>         The leading dimension of the array VT2.  LDVT2 >= N.
!> \endverbatim
!>
!> \param[in] IDXC
!> \verbatim
!>          IDXC is INTEGER array, dimension (N)
!>         The permutation used to arrange the columns of U (and rows of
!>         VT) into three groups:  the first group contains non-zero
!>         entries only at and above (or before) NL +1; the second
!>         contains non-zero entries only at and below (or after) NL+2;
!>         and the third is dense. The first column of U and the row of
!>         VT are treated separately, however.
!>
!>         The rows of the singular vectors found by SLASD4
!>         must be likewise permuted before the matrix multiplies can
!>         take place.
!> \endverbatim
!>
!> \param[in] CTOT
!> \verbatim
!>          CTOT is INTEGER array, dimension (4)
!>         A count of the total number of the various types of columns
!>         in U (or rows in VT), as described in IDXC. The fourth column
!>         type is any column which has been deflated.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is REAL array, dimension (K)
!>         The first K elements of this array contain the components
!>         of the deflation-adjusted updating row vector.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>         = 0:  successful exit.
!>         < 0:  if INFO = -i, the i-th argument had an illegal value.
!>         > 0:  if INFO = 1, a singular value did not converge
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
!> \date June 2017
!
!> \ingroup OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>     Ming Gu and Huan Ren, Computer Science Division, University of
!>     California at Berkeley, USA
!>
!  =====================================================================
      SUBROUTINE SLASD3(Nl,Nr,Sqre,K,D,Q,Ldq,Dsigma,U,Ldu,U2,Ldu2,Vt,   &
     &                  Ldvt,Vt2,Ldvt2,Idxc,Ctot,Z,Info)
      IMPLICIT NONE
!*--SLASD3227
!
!  -- LAPACK auxiliary routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      INTEGER Info , K , Ldq , Ldu , Ldu2 , Ldvt , Ldvt2 , Nl , Nr ,    &
     &        Sqre
!     ..
!     .. Array Arguments ..
      INTEGER Ctot(*) , Idxc(*)
      REAL D(*) , Dsigma(*) , Q(Ldq,*) , U(Ldu,*) , U2(Ldu2,*) ,        &
     &     Vt(Ldvt,*) , Vt2(Ldvt2,*) , Z(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO , NEGONE
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0,NEGONE=-1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER ctemp , i , j , jc , ktemp , m , n , nlp1 , nlp2 , nrp1
      REAL rho , temp
!     ..
!     .. External Functions ..
      REAL SLAMC3 , SNRM2
      EXTERNAL SLAMC3 , SNRM2
!     ..
!     .. External Subroutines ..
      EXTERNAL SCOPY , SGEMM , SLACPY , SLASCL , SLASD4 , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , SIGN , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
!
      IF ( Nl<1 ) THEN
         Info = -1
      ELSEIF ( Nr<1 ) THEN
         Info = -2
      ELSEIF ( (Sqre/=1) .AND. (Sqre/=0) ) THEN
         Info = -3
      ENDIF
!
      n = Nl + Nr + 1
      m = n + Sqre
      nlp1 = Nl + 1
      nlp2 = Nl + 2
!
      IF ( (K<1) .OR. (K>n) ) THEN
         Info = -4
      ELSEIF ( Ldq<K ) THEN
         Info = -7
      ELSEIF ( Ldu<n ) THEN
         Info = -10
      ELSEIF ( Ldu2<n ) THEN
         Info = -12
      ELSEIF ( Ldvt<m ) THEN
         Info = -14
      ELSEIF ( Ldvt2<m ) THEN
         Info = -16
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SLASD3',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( K==1 ) THEN
         D(1) = ABS(Z(1))
         CALL SCOPY(m,Vt2(1,1),Ldvt2,Vt(1,1),Ldvt)
         IF ( Z(1)>ZERO ) THEN
            CALL SCOPY(n,U2(1,1),1,U(1,1),1)
         ELSE
            DO i = 1 , n
               U(i,1) = -U2(i,1)
            ENDDO
         ENDIF
         RETURN
      ENDIF
!
!     Modify values DSIGMA(i) to make sure all DSIGMA(i)-DSIGMA(j) can
!     be computed with high relative accuracy (barring over/underflow).
!     This is a problem on machines without a guard digit in
!     add/subtract (Cray XMP, Cray YMP, Cray C 90 and Cray 2).
!     The following code replaces DSIGMA(I) by 2*DSIGMA(I)-DSIGMA(I),
!     which on any of these machines zeros out the bottommost
!     bit of DSIGMA(I) if it is 1; this makes the subsequent
!     subtractions DSIGMA(I)-DSIGMA(J) unproblematic when cancellation
!     occurs. On binary machines with a guard digit (almost all
!     machines) it does not change DSIGMA(I) at all. On hexadecimal
!     and decimal machines with a guard digit, it slightly
!     changes the bottommost bits of DSIGMA(I). It does not account
!     for hexadecimal or decimal machines without guard digits
!     (we know of none). We use a subroutine call to compute
!     2*DSIGMA(I) to prevent optimizing compilers from eliminating
!     this code.
!
      DO i = 1 , K
         Dsigma(i) = SLAMC3(Dsigma(i),Dsigma(i)) - Dsigma(i)
      ENDDO
!
!     Keep a copy of Z.
!
      CALL SCOPY(K,Z,1,Q,1)
!
!     Normalize Z.
!
      rho = SNRM2(K,Z,1)
      CALL SLASCL('G',0,0,rho,ONE,K,1,Z,K,Info)
      rho = rho*rho
!
!     Find the new singular values.
!
      DO j = 1 , K
         CALL SLASD4(K,j,Dsigma,Z,U(1,j),rho,D(j),Vt(1,j),Info)
!
!        If the zero finder fails, report the convergence failure.
!
         IF ( Info/=0 ) RETURN
      ENDDO
!
!     Compute updated Z.
!
      DO i = 1 , K
         Z(i) = U(i,K)*Vt(i,K)
         DO j = 1 , i - 1
            Z(i) = Z(i)                                                 &
     &             *(U(i,j)*Vt(i,j)/(Dsigma(i)-Dsigma(j))/(Dsigma(i)+   &
     &             Dsigma(j)))
         ENDDO
         DO j = i , K - 1
            Z(i) = Z(i)                                                 &
     &             *(U(i,j)*Vt(i,j)/(Dsigma(i)-Dsigma(j+1))/(Dsigma(i)  &
     &             +Dsigma(j+1)))
         ENDDO
         Z(i) = SIGN(SQRT(ABS(Z(i))),Q(i,1))
      ENDDO
!
!     Compute left singular vectors of the modified diagonal matrix,
!     and store related information for the right singular vectors.
!
      DO i = 1 , K
         Vt(1,i) = Z(1)/U(1,i)/Vt(1,i)
         U(1,i) = NEGONE
         DO j = 2 , K
            Vt(j,i) = Z(j)/U(j,i)/Vt(j,i)
            U(j,i) = Dsigma(j)*Vt(j,i)
         ENDDO
         temp = SNRM2(K,U(1,i),1)
         Q(1,i) = U(1,i)/temp
         DO j = 2 , K
            jc = Idxc(j)
            Q(j,i) = U(jc,i)/temp
         ENDDO
      ENDDO
!
!     Update the left singular vector matrix.
!
      IF ( K==2 ) THEN
         CALL SGEMM('N','N',n,K,K,ONE,U2,Ldu2,Q,Ldq,ZERO,U,Ldu)
         GOTO 100
      ENDIF
      IF ( Ctot(1)>0 ) THEN
         CALL SGEMM('N','N',Nl,K,Ctot(1),ONE,U2(1,2),Ldu2,Q(2,1),Ldq,   &
     &              ZERO,U(1,1),Ldu)
         IF ( Ctot(3)>0 ) THEN
            ktemp = 2 + Ctot(1) + Ctot(2)
            CALL SGEMM('N','N',Nl,K,Ctot(3),ONE,U2(1,ktemp),Ldu2,       &
     &                 Q(ktemp,1),Ldq,ONE,U(1,1),Ldu)
         ENDIF
      ELSEIF ( Ctot(3)>0 ) THEN
         ktemp = 2 + Ctot(1) + Ctot(2)
         CALL SGEMM('N','N',Nl,K,Ctot(3),ONE,U2(1,ktemp),Ldu2,Q(ktemp,1)&
     &              ,Ldq,ZERO,U(1,1),Ldu)
      ELSE
         CALL SLACPY('F',Nl,K,U2,Ldu2,U,Ldu)
      ENDIF
      CALL SCOPY(K,Q(1,1),Ldq,U(nlp1,1),Ldu)
      ktemp = 2 + Ctot(1)
      ctemp = Ctot(2) + Ctot(3)
      CALL SGEMM('N','N',Nr,K,ctemp,ONE,U2(nlp2,ktemp),Ldu2,Q(ktemp,1), &
     &           Ldq,ZERO,U(nlp2,1),Ldu)
!
!     Generate the right singular vectors.
!
 100  DO i = 1 , K
         temp = SNRM2(K,Vt(1,i),1)
         Q(i,1) = Vt(1,i)/temp
         DO j = 2 , K
            jc = Idxc(j)
            Q(i,j) = Vt(jc,i)/temp
         ENDDO
      ENDDO
!
!     Update the right singular vector matrix.
!
      IF ( K==2 ) THEN
         CALL SGEMM('N','N',K,m,K,ONE,Q,Ldq,Vt2,Ldvt2,ZERO,Vt,Ldvt)
         RETURN
      ENDIF
      ktemp = 1 + Ctot(1)
      CALL SGEMM('N','N',K,nlp1,ktemp,ONE,Q(1,1),Ldq,Vt2(1,1),Ldvt2,    &
     &           ZERO,Vt(1,1),Ldvt)
      ktemp = 2 + Ctot(1) + Ctot(2)
      IF ( ktemp<=Ldvt2 ) CALL SGEMM('N','N',K,nlp1,Ctot(3),ONE,        &
     &                               Q(1,ktemp),Ldq,Vt2(ktemp,1),Ldvt2, &
     &                               ONE,Vt(1,1),Ldvt)
!
      ktemp = Ctot(1) + 1
      nrp1 = Nr + Sqre
      IF ( ktemp>1 ) THEN
         DO i = 1 , K
            Q(i,ktemp) = Q(i,1)
         ENDDO
         DO i = nlp2 , m
            Vt2(ktemp,i) = Vt2(1,i)
         ENDDO
      ENDIF
      ctemp = 1 + Ctot(2) + Ctot(3)
      CALL SGEMM('N','N',K,nrp1,ctemp,ONE,Q(1,ktemp),Ldq,Vt2(ktemp,nlp2)&
     &           ,Ldvt2,ZERO,Vt(1,nlp2),Ldvt)
!
!
!     End of SLASD3
!
      END SUBROUTINE SLASD3
