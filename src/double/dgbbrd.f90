!*==dgbbrd.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DGBBRD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGBBRD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbbrd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbbrd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbbrd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGBBRD( VECT, M, N, NCC, KL, KU, AB, LDAB, D, E, Q,
!                          LDQ, PT, LDPT, C, LDC, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          VECT
!       INTEGER            INFO, KL, KU, LDAB, LDC, LDPT, LDQ, M, N, NCC
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AB( LDAB, * ), C( LDC, * ), D( * ), E( * ),
!      $                   PT( LDPT, * ), Q( LDQ, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGBBRD reduces a real general m-by-n band matrix A to upper
!> bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
!>
!> The routine computes B, and optionally forms Q or P**T, or computes
!> Q**T*C for a given matrix C.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] VECT
!> \verbatim
!>          VECT is CHARACTER*1
!>          Specifies whether or not the matrices Q and P**T are to be
!>          formed.
!>          = 'N': do not form Q or P**T;
!>          = 'Q': form Q only;
!>          = 'P': form P**T only;
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
!>          AB is DOUBLE PRECISION array, dimension (LDAB,N)
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
!>          D is DOUBLE PRECISION array, dimension (min(M,N))
!>          The diagonal elements of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (min(M,N)-1)
!>          The superdiagonal elements of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ,M)
!>          If VECT = 'Q' or 'B', the m-by-m orthogonal matrix Q.
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
!>          PT is DOUBLE PRECISION array, dimension (LDPT,N)
!>          If VECT = 'P' or 'B', the n-by-n orthogonal matrix P'.
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
!>          C is DOUBLE PRECISION array, dimension (LDC,NCC)
!>          On entry, an m-by-ncc matrix C.
!>          On exit, C is overwritten by Q**T*C.
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
!>          WORK is DOUBLE PRECISION array, dimension (2*max(M,N))
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
!> \ingroup doubleGBcomputational
!
!  =====================================================================
      SUBROUTINE DGBBRD(Vect,M,N,Ncc,Kl,Ku,Ab,Ldab,D,E,Q,Ldq,Pt,Ldpt,C, &
     &                  Ldc,Work,Info)
      IMPLICIT NONE
!*--DGBBRD191
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
      DOUBLE PRECISION Ab(Ldab,*) , C(Ldc,*) , D(*) , E(*) , Pt(Ldpt,*) &
     &                 , Q(Ldq,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL wantb , wantc , wantpt , wantq
      INTEGER i , inca , j , j1 , j2 , kb , kb1 , kk , klm , klu1 ,     &
     &        kun , l , minmn , ml , ml0 , mn , mu , mu0 , nr , nrt
      DOUBLE PRECISION ra , rb , rc , rs
!     ..
!     .. External Subroutines ..
      EXTERNAL DLARGV , DLARTG , DLARTV , DLASET , DROT , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
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
         CALL XERBLA('DGBBRD',-Info)
         RETURN
      ENDIF
!
!     Initialize Q and P**T to the unit matrix, if needed
!
      IF ( wantq ) CALL DLASET('Full',M,M,ZERO,ONE,Q,Ldq)
      IF ( wantpt ) CALL DLASET('Full',N,N,ZERO,ONE,Pt,Ldpt)
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
!        The sines of the plane rotations are stored in WORK(1:max(m,n))
!        and the cosines in WORK(max(m,n)+1:2*max(m,n)).
!
         mn = MAX(M,N)
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
               IF ( nr>0 ) CALL DLARGV(nr,Ab(klu1,j1-klm-1),inca,       &
     &                                 Work(j1),kb1,Work(mn+j1),kb1)
!
!              apply plane rotations from the left
!
               DO l = 1 , kb
                  IF ( j2-klm+l-1>N ) THEN
                     nrt = nr - 1
                  ELSE
                     nrt = nr
                  ENDIF
                  IF ( nrt>0 ) CALL DLARTV(nrt,Ab(klu1-l,j1-klm+l-1),   &
     &                 inca,Ab(klu1-l+1,j1-klm+l-1),inca,Work(mn+j1),   &
     &                 Work(j1),kb1)
               ENDDO
!
               IF ( ml>ml0 ) THEN
                  IF ( ml<=M-i+1 ) THEN
!
!                    generate plane rotation to annihilate a(i+ml-1,i)
!                    within the band, and apply rotation from the left
!
                     CALL DLARTG(Ab(Ku+ml-1,i),Ab(Ku+ml,i),             &
     &                           Work(mn+i+ml-1),Work(i+ml-1),ra)
                     Ab(Ku+ml-1,i) = ra
                     IF ( i<N )                                         &
     &                    CALL DROT(MIN(Ku+ml-2,N-i),Ab(Ku+ml-2,i+1),   &
     &                    Ldab-1,Ab(Ku+ml-1,i+1),Ldab-1,Work(mn+i+ml-1),&
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
                     CALL DROT(M,Q(1,j-1),1,Q(1,j),1,Work(mn+j),Work(j))
                  ENDDO
               ENDIF
!
               IF ( wantc ) THEN
!
!                 apply plane rotations to C
!
                  DO j = j1 , j2 , kb1
                     CALL DROT(Ncc,C(j-1,1),Ldc,C(j,1),Ldc,Work(mn+j),  &
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
                  Ab(1,j+kun) = Work(mn+j)*Ab(1,j+kun)
               ENDDO
!
!              generate plane rotations to annihilate nonzero elements
!              which have been generated above the band
!
               IF ( nr>0 ) CALL DLARGV(nr,Ab(1,j1+kun-1),inca,          &
     &                                 Work(j1+kun),kb1,Work(mn+j1+kun),&
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
                  IF ( nrt>0 ) CALL DLARTV(nrt,Ab(l+1,j1+kun-1),inca,   &
     &                 Ab(l,j1+kun),inca,Work(mn+j1+kun),Work(j1+kun),  &
     &                 kb1)
               ENDDO
!
               IF ( ml==ml0 .AND. mu>mu0 ) THEN
                  IF ( mu<=N-i+1 ) THEN
!
!                    generate plane rotation to annihilate a(i,i+mu-1)
!                    within the band, and apply rotation from the right
!
                     CALL DLARTG(Ab(Ku-mu+3,i+mu-2),Ab(Ku-mu+2,i+mu-1), &
     &                           Work(mn+i+mu-1),Work(i+mu-1),ra)
                     Ab(Ku-mu+3,i+mu-2) = ra
                     CALL DROT(MIN(Kl+mu-2,M-i),Ab(Ku-mu+4,i+mu-2),1,   &
     &                         Ab(Ku-mu+3,i+mu-1),1,Work(mn+i+mu-1),    &
     &                         Work(i+mu-1))
                  ENDIF
                  nr = nr + 1
                  j1 = j1 - kb1
               ENDIF
!
               IF ( wantpt ) THEN
!
!                 accumulate product of plane rotations in P**T
!
                  DO j = j1 , j2 , kb1
                     CALL DROT(N,Pt(j+kun-1,1),Ldpt,Pt(j+kun,1),Ldpt,   &
     &                         Work(mn+j+kun),Work(j+kun))
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
                  Ab(klu1,j+kun) = Work(mn+j+kun)*Ab(klu1,j+kun)
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
!        A has been reduced to lower bidiagonal form
!
!        Transform lower bidiagonal form to upper bidiagonal by applying
!        plane rotations from the left, storing diagonal elements in D
!        and off-diagonal elements in E
!
         DO i = 1 , MIN(M-1,N)
            CALL DLARTG(Ab(1,i),Ab(2,i),rc,rs,ra)
            D(i) = ra
            IF ( i<N ) THEN
               E(i) = rs*Ab(1,i+1)
               Ab(1,i+1) = rc*Ab(1,i+1)
            ENDIF
            IF ( wantq ) CALL DROT(M,Q(1,i),1,Q(1,i+1),1,rc,rs)
            IF ( wantc ) CALL DROT(Ncc,C(i,1),Ldc,C(i+1,1),Ldc,rc,rs)
         ENDDO
         IF ( M<=N ) D(M) = Ab(1,M)
      ELSEIF ( Ku>0 ) THEN
!
!        A has been reduced to upper bidiagonal form
!
         IF ( M<N ) THEN
!
!           Annihilate a(m,m+1) by applying plane rotations from the
!           right, storing diagonal elements in D and off-diagonal
!           elements in E
!
            rb = Ab(Ku,M+1)
            DO i = M , 1 , -1
               CALL DLARTG(Ab(Ku+1,i),rb,rc,rs,ra)
               D(i) = ra
               IF ( i>1 ) THEN
                  rb = -rs*Ab(Ku,i)
                  E(i-1) = rc*Ab(Ku,i)
               ENDIF
               IF ( wantpt ) CALL DROT(N,Pt(i,1),Ldpt,Pt(M+1,1),Ldpt,rc,&
     &                                 rs)
            ENDDO
         ELSE
!
!           Copy off-diagonal elements to E and diagonal elements to D
!
            DO i = 1 , minmn - 1
               E(i) = Ab(Ku,i+1)
            ENDDO
            DO i = 1 , minmn
               D(i) = Ab(Ku+1,i)
            ENDDO
         ENDIF
      ELSE
!
!        A is diagonal. Set elements of E to zero and copy diagonal
!        elements to D.
!
         DO i = 1 , minmn - 1
            E(i) = ZERO
         ENDDO
         DO i = 1 , minmn
            D(i) = Ab(1,i)
         ENDDO
      ENDIF
!
!     End of DGBBRD
!
      END SUBROUTINE DGBBRD
