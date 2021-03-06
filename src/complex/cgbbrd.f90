!*==cgbbrd.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CGBBRD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGBBRD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgbbrd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgbbrd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgbbrd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGBBRD( VECT, M, N, NCC, KL, KU, AB, LDAB, D, E, Q,
!                          LDQ, PT, LDPT, C, LDC, WORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          VECT
!       INTEGER            INFO, KL, KU, LDAB, LDC, LDPT, LDQ, M, N, NCC
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), E( * ), RWORK( * )
!       COMPLEX            AB( LDAB, * ), C( LDC, * ), PT( LDPT, * ),
!      $                   Q( LDQ, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGBBRD reduces a complex general m-by-n band matrix A to real upper
!> bidiagonal form B by a unitary transformation: Q**H * A * P = B.
!>
!> The routine computes B, and optionally forms Q or P**H, or computes
!> Q**H*C for a given matrix C.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] VECT
!> \verbatim
!>          VECT is CHARACTER*1
!>          Specifies whether or not the matrices Q and P**H are to be
!>          formed.
!>          = 'N': do not form Q or P**H;
!>          = 'Q': form Q only;
!>          = 'P': form P**H only;
!>          = 'B': form both.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NCC
!> \verbatim
!>          NCC is INTEGER
!>          The number of columns of the matrix C.  NCC >= 0.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of subdiagonals of the matrix A. KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of superdiagonals of the matrix A. KU >= 0.
!> \endverbatim
!>
!> \param[in,out] AB
!> \verbatim
!>          AB is COMPLEX array, dimension (LDAB,N)
!>          On entry, the m-by-n band matrix A, stored in rows 1 to
!>          KL+KU+1. The j-th column of A is stored in the j-th column of
!>          the array AB as follows:
!>          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl).
!>          On exit, A is overwritten by values generated during the
!>          reduction.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array A. LDAB >= KL+KU+1.
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is REAL array, dimension (min(M,N))
!>          The diagonal elements of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is REAL array, dimension (min(M,N)-1)
!>          The superdiagonal elements of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is COMPLEX array, dimension (LDQ,M)
!>          If VECT = 'Q' or 'B', the m-by-m unitary matrix Q.
!>          If VECT = 'N' or 'P', the array Q is not referenced.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.
!>          LDQ >= max(1,M) if VECT = 'Q' or 'B'; LDQ >= 1 otherwise.
!> \endverbatim
!>
!> \param[out] PT
!> \verbatim
!>          PT is COMPLEX array, dimension (LDPT,N)
!>          If VECT = 'P' or 'B', the n-by-n unitary matrix P'.
!>          If VECT = 'N' or 'Q', the array PT is not referenced.
!> \endverbatim
!>
!> \param[in] LDPT
!> \verbatim
!>          LDPT is INTEGER
!>          The leading dimension of the array PT.
!>          LDPT >= max(1,N) if VECT = 'P' or 'B'; LDPT >= 1 otherwise.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC,NCC)
!>          On entry, an m-by-ncc matrix C.
!>          On exit, C is overwritten by Q**H*C.
!>          C is not referenced if NCC = 0.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C.
!>          LDC >= max(1,M) if NCC > 0; LDC >= 1 if NCC = 0.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (max(M,N))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (max(M,N))
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
!> \ingroup complexGBcomputational
!
!  =====================================================================
      SUBROUTINE CGBBRD(Vect,M,N,Ncc,Kl,Ku,Ab,Ldab,D,E,Q,Ldq,Pt,Ldpt,C, &
     &                  Ldc,Work,Rwork,Info)
      IMPLICIT NONE
!*--CGBBRD197
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Vect
      INTEGER Info , Kl , Ku , Ldab , Ldc , Ldpt , Ldq , M , N , Ncc
!     ..
!     .. Array Arguments ..
      REAL D(*) , E(*) , Rwork(*)
      COMPLEX Ab(Ldab,*) , C(Ldc,*) , Pt(Ldpt,*) , Q(Ldq,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E+0)
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      LOGICAL wantb , wantc , wantpt , wantq
      INTEGER i , inca , j , j1 , j2 , kb , kb1 , kk , klm , klu1 ,     &
     &        kun , l , minmn , ml , ml0 , mu , mu0 , nr , nrt
      REAL abst , rc
      COMPLEX ra , rb , rs , t
!     ..
!     .. External Subroutines ..
      EXTERNAL CLARGV , CLARTG , CLARTV , CLASET , CROT , CSCAL , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , CONJG , MAX , MIN
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      wantb = LSAME(Vect,'B')
      wantq = LSAME(Vect,'Q') .OR. wantb
      wantpt = LSAME(Vect,'P') .OR. wantb
      wantc = Ncc>0
      klu1 = Kl + Ku + 1
      Info = 0
      IF ( .NOT.wantq .AND. .NOT.wantpt .AND. .NOT.LSAME(Vect,'N') )    &
     &     THEN
         Info = -1
      ELSEIF ( M<0 ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Ncc<0 ) THEN
         Info = -4
      ELSEIF ( Kl<0 ) THEN
         Info = -5
      ELSEIF ( Ku<0 ) THEN
         Info = -6
      ELSEIF ( Ldab<klu1 ) THEN
         Info = -8
      ELSEIF ( Ldq<1 .OR. wantq .AND. Ldq<MAX(1,M) ) THEN
         Info = -12
      ELSEIF ( Ldpt<1 .OR. wantpt .AND. Ldpt<MAX(1,N) ) THEN
         Info = -14
      ELSEIF ( Ldc<1 .OR. wantc .AND. Ldc<MAX(1,M) ) THEN
         Info = -16
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGBBRD',-Info)
         RETURN
      ENDIF
!
!     Initialize Q and P**H to the unit matrix, if needed
!
      IF ( wantq ) CALL CLASET('Full',M,M,CZERO,CONE,Q,Ldq)
      IF ( wantpt ) CALL CLASET('Full',N,N,CZERO,CONE,Pt,Ldpt)
!
!     Quick return if possible.
!
      IF ( M==0 .OR. N==0 ) RETURN
!
      minmn = MIN(M,N)
!
      IF ( Kl+Ku>1 ) THEN
!
!        Reduce to upper bidiagonal form if KU > 0; if KU = 0, reduce
!        first to lower bidiagonal form and then transform to upper
!        bidiagonal
!
         IF ( Ku>0 ) THEN
            ml0 = 1
            mu0 = 2
         ELSE
            ml0 = 2
            mu0 = 1
         ENDIF
!
!        Wherever possible, plane rotations are generated and applied in
!        vector operations of length NR over the index set J1:J2:KLU1.
!
!        The complex sines of the plane rotations are stored in WORK,
!        and the real cosines in RWORK.
!
         klm = MIN(M-1,Kl)
         kun = MIN(N-1,Ku)
         kb = klm + kun
         kb1 = kb + 1
         inca = kb1*Ldab
         nr = 0
         j1 = klm + 2
         j2 = 1 - kun
!
         DO i = 1 , minmn
!
!           Reduce i-th column and i-th row of matrix to bidiagonal form
!
            ml = klm + 1
            mu = kun + 1
            DO kk = 1 , kb
               j1 = j1 + kb
               j2 = j2 + kb
!
!              generate plane rotations to annihilate nonzero elements
!              which have been created below the band
!
               IF ( nr>0 ) CALL CLARGV(nr,Ab(klu1,j1-klm-1),inca,       &
     &                                 Work(j1),kb1,Rwork(j1),kb1)
!
!              apply plane rotations from the left
!
               DO l = 1 , kb
                  IF ( j2-klm+l-1>N ) THEN
                     nrt = nr - 1
                  ELSE
                     nrt = nr
                  ENDIF
                  IF ( nrt>0 ) CALL CLARTV(nrt,Ab(klu1-l,j1-klm+l-1),   &
     &                 inca,Ab(klu1-l+1,j1-klm+l-1),inca,Rwork(j1),     &
     &                 Work(j1),kb1)
               ENDDO
!
               IF ( ml>ml0 ) THEN
                  IF ( ml<=M-i+1 ) THEN
!
!                    generate plane rotation to annihilate a(i+ml-1,i)
!                    within the band, and apply rotation from the left
!
                     CALL CLARTG(Ab(Ku+ml-1,i),Ab(Ku+ml,i),Rwork(i+ml-1)&
     &                           ,Work(i+ml-1),ra)
                     Ab(Ku+ml-1,i) = ra
                     IF ( i<N )                                         &
     &                    CALL CROT(MIN(Ku+ml-2,N-i),Ab(Ku+ml-2,i+1),   &
     &                    Ldab-1,Ab(Ku+ml-1,i+1),Ldab-1,Rwork(i+ml-1),  &
     &                    Work(i+ml-1))
                  ENDIF
                  nr = nr + 1
                  j1 = j1 - kb1
               ENDIF
!
               IF ( wantq ) THEN
!
!                 accumulate product of plane rotations in Q
!
                  DO j = j1 , j2 , kb1
                     CALL CROT(M,Q(1,j-1),1,Q(1,j),1,Rwork(j),          &
     &                         CONJG(Work(j)))
                  ENDDO
               ENDIF
!
               IF ( wantc ) THEN
!
!                 apply plane rotations to C
!
                  DO j = j1 , j2 , kb1
                     CALL CROT(Ncc,C(j-1,1),Ldc,C(j,1),Ldc,Rwork(j),    &
     &                         Work(j))
                  ENDDO
               ENDIF
!
               IF ( j2+kun>N ) THEN
!
!                 adjust J2 to keep within the bounds of the matrix
!
                  nr = nr - 1
                  j2 = j2 - kb1
               ENDIF
!
               DO j = j1 , j2 , kb1
!
!                 create nonzero element a(j-1,j+ku) above the band
!                 and store it in WORK(n+1:2*n)
!
                  Work(j+kun) = Work(j)*Ab(1,j+kun)
                  Ab(1,j+kun) = Rwork(j)*Ab(1,j+kun)
               ENDDO
!
!              generate plane rotations to annihilate nonzero elements
!              which have been generated above the band
!
               IF ( nr>0 ) CALL CLARGV(nr,Ab(1,j1+kun-1),inca,          &
     &                                 Work(j1+kun),kb1,Rwork(j1+kun),  &
     &                                 kb1)
!
!              apply plane rotations from the right
!
               DO l = 1 , kb
                  IF ( j2+l-1>M ) THEN
                     nrt = nr - 1
                  ELSE
                     nrt = nr
                  ENDIF
                  IF ( nrt>0 ) CALL CLARTV(nrt,Ab(l+1,j1+kun-1),inca,   &
     &                 Ab(l,j1+kun),inca,Rwork(j1+kun),Work(j1+kun),kb1)
               ENDDO
!
               IF ( ml==ml0 .AND. mu>mu0 ) THEN
                  IF ( mu<=N-i+1 ) THEN
!
!                    generate plane rotation to annihilate a(i,i+mu-1)
!                    within the band, and apply rotation from the right
!
                     CALL CLARTG(Ab(Ku-mu+3,i+mu-2),Ab(Ku-mu+2,i+mu-1), &
     &                           Rwork(i+mu-1),Work(i+mu-1),ra)
                     Ab(Ku-mu+3,i+mu-2) = ra
                     CALL CROT(MIN(Kl+mu-2,M-i),Ab(Ku-mu+4,i+mu-2),1,   &
     &                         Ab(Ku-mu+3,i+mu-1),1,Rwork(i+mu-1),      &
     &                         Work(i+mu-1))
                  ENDIF
                  nr = nr + 1
                  j1 = j1 - kb1
               ENDIF
!
               IF ( wantpt ) THEN
!
!                 accumulate product of plane rotations in P**H
!
                  DO j = j1 , j2 , kb1
                     CALL CROT(N,Pt(j+kun-1,1),Ldpt,Pt(j+kun,1),Ldpt,   &
     &                         Rwork(j+kun),CONJG(Work(j+kun)))
                  ENDDO
               ENDIF
!
               IF ( j2+kb>M ) THEN
!
!                 adjust J2 to keep within the bounds of the matrix
!
                  nr = nr - 1
                  j2 = j2 - kb1
               ENDIF
!
               DO j = j1 , j2 , kb1
!
!                 create nonzero element a(j+kl+ku,j+ku-1) below the
!                 band and store it in WORK(1:n)
!
                  Work(j+kb) = Work(j+kun)*Ab(klu1,j+kun)
                  Ab(klu1,j+kun) = Rwork(j+kun)*Ab(klu1,j+kun)
               ENDDO
!
               IF ( ml>ml0 ) THEN
                  ml = ml - 1
               ELSE
                  mu = mu - 1
               ENDIF
            ENDDO
         ENDDO
      ENDIF
!
      IF ( Ku==0 .AND. Kl>0 ) THEN
!
!        A has been reduced to complex lower bidiagonal form
!
!        Transform lower bidiagonal form to upper bidiagonal by applying
!        plane rotations from the left, overwriting superdiagonal
!        elements on subdiagonal elements
!
         DO i = 1 , MIN(M-1,N)
            CALL CLARTG(Ab(1,i),Ab(2,i),rc,rs,ra)
            Ab(1,i) = ra
            IF ( i<N ) THEN
               Ab(2,i) = rs*Ab(1,i+1)
               Ab(1,i+1) = rc*Ab(1,i+1)
            ENDIF
            IF ( wantq ) CALL CROT(M,Q(1,i),1,Q(1,i+1),1,rc,CONJG(rs))
            IF ( wantc ) CALL CROT(Ncc,C(i,1),Ldc,C(i+1,1),Ldc,rc,rs)
         ENDDO
!
!        A has been reduced to complex upper bidiagonal form or is
!        diagonal
!
      ELSEIF ( Ku>0 .AND. M<N ) THEN
!
!           Annihilate a(m,m+1) by applying plane rotations from the
!           right
!
         rb = Ab(Ku,M+1)
         DO i = M , 1 , -1
            CALL CLARTG(Ab(Ku+1,i),rb,rc,rs,ra)
            Ab(Ku+1,i) = ra
            IF ( i>1 ) THEN
               rb = -CONJG(rs)*Ab(Ku,i)
               Ab(Ku,i) = rc*Ab(Ku,i)
            ENDIF
            IF ( wantpt ) CALL CROT(N,Pt(i,1),Ldpt,Pt(M+1,1),Ldpt,rc,   &
     &                              CONJG(rs))
         ENDDO
      ENDIF
!
!     Make diagonal and superdiagonal elements real, storing them in D
!     and E
!
      t = Ab(Ku+1,1)
      DO i = 1 , minmn
         abst = ABS(t)
         D(i) = abst
         IF ( abst/=ZERO ) THEN
            t = t/abst
         ELSE
            t = CONE
         ENDIF
         IF ( wantq ) CALL CSCAL(M,t,Q(1,i),1)
         IF ( wantc ) CALL CSCAL(Ncc,CONJG(t),C(i,1),Ldc)
         IF ( i<minmn ) THEN
            IF ( Ku==0 .AND. Kl==0 ) THEN
               E(i) = ZERO
               t = Ab(1,i+1)
            ELSE
               IF ( Ku==0 ) THEN
                  t = Ab(2,i)*CONJG(t)
               ELSE
                  t = Ab(Ku,i+1)*CONJG(t)
               ENDIF
               abst = ABS(t)
               E(i) = abst
               IF ( abst/=ZERO ) THEN
                  t = t/abst
               ELSE
                  t = CONE
               ENDIF
               IF ( wantpt ) CALL CSCAL(N,t,Pt(i+1,1),Ldpt)
               t = Ab(Ku+1,i+1)*CONJG(t)
            ENDIF
         ENDIF
      ENDDO
!
!     End of CGBBRD
!
      END SUBROUTINE CGBBRD
