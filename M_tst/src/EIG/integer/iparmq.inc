!*==iparmq.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
      INTEGER FUNCTION IPARMQ(Ispec,Name,Opts,N,Ilo,Ihi,Lwork)
      IMPLICIT NONE
!*--IPARMQ299
!
      INTEGER INMIN , INWIN , INIBL , ISHFTS , IACC22
      PARAMETER (INMIN=12,INWIN=13,INIBL=14,ISHFTS=15,IACC22=16)
      INTEGER NMIN , K22MIN , KACMIN , NIBBLE , KNWSWP
      PARAMETER (NMIN=11,K22MIN=14,KACMIN=14,NIBBLE=14,KNWSWP=500)
      REAL TWO
      PARAMETER (TWO=2.0)
!     ..
!     .. Scalar Arguments ..
      INTEGER Ihi , Ilo , Ispec , Lwork , N
      CHARACTER Name*(*) , Opts*(*)
!     ..
!     .. Local Scalars ..
      INTEGER nh , ns
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC LOG , MAX , MOD , NINT , REAL
!     ..
!     .. Executable Statements ..
      IF ( (Ispec==ISHFTS) .OR. (Ispec==INWIN) .OR. (Ispec==IACC22) )   &
     &     THEN
!
!        ==== Set the number simultaneous shifts ====
!
         nh = Ihi - Ilo + 1
         ns = 2
         IF ( nh>=30 ) ns = 4
         IF ( nh>=60 ) ns = 10
         IF ( nh>=150 ) ns = MAX(10,nh/NINT(LOG(REAL(nh))/LOG(TWO)))
         IF ( nh>=590 ) ns = 64
         IF ( nh>=3000 ) ns = 128
         IF ( nh>=6000 ) ns = 256
         ns = MAX(2,ns-MOD(ns,2))
      ENDIF
!
      IF ( Ispec==INMIN ) THEN
!
!
!        ===== Matrices of order smaller than NMIN get sent
!        .     to LAHQR, the classic double shift algorithm.
!        .     This must be at least 11. ====
!
         IPARMQ = NMIN
!
      ELSEIF ( Ispec==INIBL ) THEN
!
!        ==== INIBL: skip a multi-shift qr iteration and
!        .    whenever aggressive early deflation finds
!        .    at least (NIBBLE*(window size)/100) deflations. ====
!
         IPARMQ = NIBBLE
!
      ELSEIF ( Ispec==ISHFTS ) THEN
!
!        ==== NSHFTS: The number of simultaneous shifts =====
!
         IPARMQ = ns
!
      ELSEIF ( Ispec==INWIN ) THEN
!
!        ==== NW: deflation window size.  ====
!
         IF ( nh<=KNWSWP ) THEN
            IPARMQ = ns
         ELSE
            IPARMQ = 3*ns/2
         ENDIF
!
      ELSEIF ( Ispec==IACC22 ) THEN
!
!        ==== IACC22: Whether to accumulate reflections
!        .     before updating the far-from-diagonal elements
!        .     and whether to use 2-by-2 block structure while
!        .     doing it.  A small amount of work could be saved
!        .     by making this choice dependent also upon the
!        .     NH=IHI-ILO+1.
!
         IPARMQ = 0
         IF ( ns>=KACMIN ) IPARMQ = 1
         IF ( ns>=K22MIN ) IPARMQ = 2
!
      ELSE
!        ===== invalid value of ispec =====
         IPARMQ = -1
!
      ENDIF
!
!     ==== End of IPARMQ ====
!
      END FUNCTION IPARMQ
