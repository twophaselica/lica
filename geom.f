!*******************************************************************
        FUNCTION FUNCBODY(X_IBM,Y_IBM,Z_IBM)
!*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      USE IBM_VAR

!       Function for an immersed body
!       Required condition to use secant method :
!          1. FUNCBODY(X,Y,Z)<0  inside the body
!             FUNCBODY(X,Y,Z)=0  at the body
!             FUNCBODY(X,Y,Z)>0  outside the body
!          2. FUNCBODY Must be continuous function
!          *It is sufficient to satisfy only 1. to finding forcing
!           point but 2. is critical to use secant method
!       Ex.
!          FUNCBODY=X**2+Y**2+Z**2-0.5**2       ! for sphere
!          FUNCBODY=X**2+Y**2-0.5**2            ! for infinite cylinder

!        FUNCBODY=-(X**2.+Y**2.-(0.5*YL)**2.)    ! for PIPE


        IF ( IBMON .EQ. 1 ) THEN

       !FOR PIPE
       FUNCBODY=-(X_IBM**2.+Y_IBM**2.-(0.5*YL)**2.)

!!       !FOR FUEL ROD
!        IF (X_IBM .LT. 0.) THEN
!        	 XX=X_IBM+XL !ALL POINT IS LOCATED IN FIRST QUADRANT USING PERIODIC CONDITION.
!        ELSE
!           XX=X_IBM
!        ENDIF
!        IF (Y_IBM .LT. 0.) THEN
!        	YY=Y_IBM+YL
!        ELSE
!        	YY=Y_IBM
!        ENDIF
!        FUNCBODY=(XX-0.5*XL)**2+(YY-0.5*YL)**2-(YL*5/14) !RADIUS IS YL*5/14 NOW.
        ELSE
          FUNCBODY=100 !FOR NO IBM CASE.

        ENDIF

        RETURN
        END

!******************************************************************
!       SECTION FOR GEOMETRY
!      - READ GRID FILE & COMPUTE VARIABLES ABOUT GEOMETRY
!      - GEOMETRY VARIABLES ARE LISTED IN 'GEOM.H'
!      - DATA READ FROM FILE ARE BLOCK 'DIM','COORD','SCALES',
!        'GEOMINPUT'
!      - COMPUTED VARIABLES ARE BLOCK 'GEOMETRY','INDEX','FIX',
!        'POSITION','VAR','VVAR'
!
!                               Revised by DONGJOO KIM
!                               Turbulence & Flow Control Lab.
!                               School of Mech. & Aerospace Eng.
!                               Seoul National Univ.
!                               9/27/1999
!
!
!      -REVISED FOR TWO-PHASE FLOW
!
!                               09/03/2015
!******************************************************************

      SUBROUTINE GEOM(gridfile)

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      USE LVS_COUPLING

      IMPLICIT NONE
      INTEGER*8 I,J,K
      
      CHARACTER*8 gridfile

      !FLOW_GEOM_VAR
      ALLOCATE(IPV(0:M1+1),JPV(0:M2+1),KPV(0:M3+1)
     &                ,IMV(-1:M1),JMV(-1:M2),KMV(-1:M3))
      ALLOCATE(FIXIL(M1M),FIXIU(M1M),FIXJL(M2M),FIXJU(M2M)
     &            ,FIXKL(M3M),FIXKU(M3M))
      ALLOCATE( X(0:M1+1),Y(0:M2+1),Z(0:M3+1))
      ALLOCATE( XP(0:M1+1),YP(0:M2+1),ZP(0:M3+1))
      ALLOCATE( SDX(0:M1),SDY(0:M2),SDZ(0:M3),VDX(M1),VDY(M2),VDZ(M3))
      ALLOCATE( SSDX(0:M1),SSDY(0:M2),SSDZ(0:M3)
     &             ,VVDX(M1),VVDY(M2),VVDZ(M3))

      !LVS_GEOM
      ALLOCATE( XF(M1L:M1U),YF(M2L:M2U),ZF(M3L:M3U))
      ALLOCATE( XPF(M1L:M1U),YPF(M2L:M2U),ZPF(M3L:M3U))
      ALLOCATE( IPF(0:M1F+1),IMF(-1:M1F),JPF(0:M2F+1),JMF(-1:M2F)
     &,KPF(0:M3F+1),KMF(-1:M3F))

      !LVS_COUPLING
      ALLOCATE( ICOU_VEL(M1F),JCOU_VEL(M2F),KCOU_VEL(M3F))
      ALLOCATE( ICOUMP_VEL(M1F),JCOUMP_VEL(M2F),KCOUMP_VEL(M3F))
      ALLOCATE( ICOU1(-1:M1+2),JCOU1(-1:M2+2),KCOU1(-1:M3+2),
     &          ICOUMP1(-1:M1+2),JCOUMP1(-1:M2+2),KCOUMP1(-1:M3+2))
      ALLOCATE( ICOU2(-1:M1+2),JCOU2(-1:M2+2),KCOU2(-1:M3+2),
     &          ICOUMP2(-1:M1+2),JCOUMP2(-1:M2+2),KCOUMP2(-1:M3+2))

      OPEN(11,FILE=gridfile)
      CALL READGRID
      CLOSE(11)

       CALL GEOMETRIC
       CALL GEOMETRIC_LVS
       CALL INDICES
       CALL INDICES_LVS
!      CALL GEOM_COUPLING   !LAST LINE IN THIS SUBROUTINE.

      IF (IPX .EQ. 1) THEN
      FIXIL(1)=0.
      FIXIU(N1M)=0.
      SDX(0)=SDX(N1M)
      SDX(N1)=SDX(1)
      SSDX(0)=SSDX(N1M)
      SSDX(N1)=SSDX(1)
      VDX(1)=0.5*(SDX(1)+SDX(N1M))
      VDX(N1)=0.5*(SDX(1)+SDX(N1M))
      X(0)=X(1)-SDX(1)
      XP(0)=XP(1)-VDX(1)
      XP(N1)=XP(N1M)+VDX(N1)
      DO 48 I=1,N1
   48 VVDX(I)=1./VDX(I)
      IMV(1)=N1M
      IPV(N1M)=1
      !FOR LEVEL-SET
        IPF(N1FM)=1
        IMF(1)=N1FM
      ENDIF
      IF (IPX.EQ.1) THEN
         IBG=1
      ELSE
         IBG=2
      ENDIF
      
      IF (IPY .EQ. 1) THEN
      FIXJL(1)=0.
      FIXJU(N2M)=0.
      SDY(0)=SDY(N2M)
      SDY(N2)=SDY(1)
      SSDY(0)=SSDY(N2M)
      SSDY(N2)=SSDY(1)
      VDY(1)=0.5*(SDY(1)+SDY(N2M))
      VDY(N2)=0.5*(SDY(1)+SDY(N2M))
      Y(0)=Y(1)-SDY(1)
      YP(0)=YP(1)-VDY(1)
      YP(N2)=YP(N2M)+VDY(N2)
      DO 49 J=1,N2
   49 VVDY(J)=1./VDY(J)
      JMV(1)=N2M
      JPV(N2M)=1
      !FOR LEVEL-SET
        JPF(N2FM)=1
        JMF(1)=N2FM
      ENDIF
      IF (IPY.EQ.1) THEN
         JBG=1
      ELSE
         JBG=2
      ENDIF
      
      IF (IPZ .EQ. 1) THEN
      FIXKL(1)=0.
      FIXKU(N3M)=0.
      SDZ(0)=SDZ(N3M)
      SDZ(N3)=SDZ(1)
      SSDZ(0)=SSDZ(N3M)
      SSDZ(N3)=SSDZ(1)
      VDZ(1)=0.5*(SDZ(1)+SDZ(N3M))
      VDZ(N3)=0.5*(SDZ(1)+SDZ(N3M))
      Z(0)=Z(1)-SDZ(1)
      ZP(0)=ZP(1)-VDZ(1)
      ZP(N3)=ZP(N3M)+VDZ(N3)
      DO 47 K=1,N3
   47 VVDZ(K)=1./VDZ(K)
      KMV(1)=N3M
      KPV(N3M)=1
      !FOR LEVEL-SET
        KPF(N3FM)=1       
        KMF(1)=N3FM
      ENDIF
      IF (IPZ.EQ.1) THEN
         KBG=1
      ELSE
         KBG=2
      ENDIF

      CALL GEOM_COUPLING !need ipf,imf,kpf,kmf

      RETURN
      END

!******************************************************************
       SUBROUTINE READGRID
!******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_GEOM_VAR
      
      IMPLICIT NONE
      REAL*8 XLG,YLG,ZLG
      INTEGER*8   I,J,K
      
      !FOR NAVIER-STOKES
      READ(11,*) N1,N2,N3
      READ(11,*) XL,YL,ZL

      READ(11,*) (X(I),I=1,N1)
      READ(11,*) (Y(J),J=1,N2)
      READ(11,*) (Z(K),K=1,N3)

      IF ((N1.GT.M1) .OR. (N2.GT.M2) .OR. (N3.GT.M3)) THEN
      write(*,*) n1,m1,n2,m2,n3,m3
         PRINT*, 'NSE-ARRAY SIZE CAN NOT HANDLE THIS GRID.'
         STOP
      END IF
      
      !FOR LEVEL-SET(G-GRID)
      READ(11,*) N1F,N2F,N3F
      READ(11,*) XLG,YLG,ZLG

      READ(11,*) (XF(I),I=1,N1F)
      READ(11,*) (YF(J),J=1,N2F)
      READ(11,*) (ZF(K),K=1,N3F)

      IF ((N1F.GT.M1F) .OR. (N2F.GT.M2F) .OR. (N3F.GT.M3F)) THEN
         PRINT*, 'LVS-ARRAY SIZE CAN NOT HANDLE THIS GRID.'
         STOP
      END IF

      N1M=N1-1
      N2M=N2-1
      N3M=N3-1
     
      N1FM=N1F-1
      N2FM=N2F-1
      N3FM=N3F-1

      N1F_BD=DFLOAT(N1FM)/DFLOAT(N1M)
      N2F_BD=DFLOAT(N2FM)/DFLOAT(N2M)
      N3F_BD=DFLOAT(N3FM)/DFLOAT(N3M)

      X(0)=X(1) !X(0),Y(0),Z(0) IS ADDED FOR GRID_COUPLING
      Y(0)=Y(1) !WHEN USING PERIODIC CONDITION
      Z(0)=Z(1)      

      !GRID_INFORMATION_WRITE
      WRITE(*,101) N1,N2,N3,N1M*N2*N3M       
      WRITE(*,102) N1FM,N2FM,N3FM,N1FM*N2FM*N3FM
      WRITE(*,103) XL,YL,ZL      
 101  FORMAT('N1 = ',I4,' N2 = ',I4,' N3 = ',I4,' NSE_TOT = ',I10)
 102  FORMAT('N1F =',I4,' N2F =',I4,' N3F =',I4,' REFINE_TOT = ',I10) 
 103  FORMAT('XL = ',F10.5,' YL = ',F10.5,' ZL = ',F10.5) 

      RETURN
      END
      
!******************************************************************
       SUBROUTINE GEOMETRIC
!******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE
      INTEGER*8   I,J,K
      
      DO 10 I=1,N1M
   10 SDX(I)=X(I+1)-X(I)
      DO 20 I=2,N1M
   20 VDX(I)=0.5*(SDX(I)+SDX(I-1))
      VDX(1)=0.5*SDX(1)             ! LEFT BOUNDED COND.
      VDX(N1)=0.5*SDX(N1M)          ! RIGHT BOUNDED COND.

      DO 30 I=1,N1M
   30 SSDX(I)=1./SDX(I)
      DO 40 I=1,N1
   40 VVDX(I)=1./VDX(I)

      DO 15 J=1,N2M
   15 SDY(J)=Y(J+1)-Y(J)
      DO 25 J=2,N2M
   25 VDY(J)=0.5*(SDY(J)+SDY(J-1))
      VDY(1)=0.5*SDY(1)             ! LOWER WALL BOUNDED COND.
      VDY(N2)=0.5*SDY(N2M)          ! UPPER WALL BOUNDED COND.

      DO 35 J=1,N2M
   35 SSDY(J)=1./SDY(J)
      DO 45 J=1,N2
   45 VVDY(J)=1./VDY(J)

      IF (N3M.EQ.1) THEN
      ZL=1.
      SDZ(1)=ZL/N3M
      SSDZ(1)=1./SDZ(1)
      VDZ(1)=SDZ(1)
      VVDZ(1)=SSDZ(1)

      ELSE
      DO 17 K=1,N3M
   17 SDZ(K)=Z(K+1)-Z(K)
      DO 27 K=2,N3M
   27 VDZ(K)=0.5*(SDZ(K)+SDZ(K-1))
      VDZ(1)=0.5*SDZ(1)             ! LOWER WALL BOUNDED COND.
      VDZ(N3)=0.5*SDZ(N3M)          ! UPPER WALL BOUNDED COND.

      DO 37 K=1,N3M
   37 SSDZ(K)=1./SDZ(K)
      DO 47 K=1,N3
   47 VVDZ(K)=1./VDZ(K)
      ENDIF

C-----set by zero
      SDX(0)=0.
      SDX(M1)=0.
      SDY(0)=0.
      SDY(M2)=0.
      SDZ(0)=0.
      SDZ(M3)=0.

      DO 5 I=1,N1M
    5 XP(I)=X(I)+0.5*SDX(I)
      XP(N1)=X(N1)
      XP(0)=X(1)

      DO 6 J=1,N2M
    6 YP(J)=Y(J)+0.5*SDY(J)
      YP(N2)=Y(N2)
      YP(0)=Y(1)

      DO 7 K=1,N3M
    7 ZP(K)=Z(K)+0.5*SDZ(K)
      ZP(N3)=Z(N3)
      ZP(0)=Z(1)

      RETURN
      END
      
!******************************************************************
       SUBROUTINE GEOMETRIC_LVS
!******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_GEOM_VAR
      
      IMPLICIT NONE
      INTEGER*8   I,J,K


        SDXF=XF(2)-XF(1)
        SDYF=YF(2)-YF(1)   
        SDZF=ZF(2)-ZF(1)

        SSDXF=1./SDXF
        SSDYF=1./SDYF
        SSDZF=1./SDZF
        
        !X-DIR
        DO I=1,N1FM
        	XPF(I)=XF(I)+0.5*SDXF
        ENDDO
        DO I=0,2*N1F_BD
          XF(-I)=XF(1)-SDXF*(I+1)
          XF(N1F+1+I)=XF(N1F)+SDXF*(I+1)
        ENDDO
        DO I=0,2*N1F_BD
          XPF(-I)=XF(-I)+0.5*SDXF
          XPF(N1F+I)=XF(N1F+I)+0.5*SDXF
        ENDDO

        !Y-DIR
        DO J=1,N2FM
        	YPF(J)=YF(J)+0.5*SDYF
        ENDDO
        DO J=0,2*N2F_BD
          YF(-J)=YF(1)-SDYF*(J+1)
          YF(N2F+1+J)=YF(N2F)+SDYF*(J+1)
        ENDDO
        DO J=0,2*N2F_BD
          YPF(-J)=YF(-J)+0.5*SDYF
          YPF(N2F+J)=YF(N2F+J)+0.5*SDYF
        ENDDO

        !Z-DIR
        DO K=1,N3FM
        	ZPF(K)=ZF(K)+0.5*SDZF
        ENDDO
        DO K=0,2*N3F_BD
          ZF(-K)=ZF(1)-SDZF*(K+1)
          ZF(N3F+1+K)=ZF(N3F)+SDZF*(K+1)
        ENDDO
        DO K=0,2*N3F_BD
          ZPF(-K)=ZF(-K)+0.5*SDZF
          ZPF(N3F+K)=ZF(N3F+K)+0.5*SDZF
        ENDDO

      RETURN
      END
      
!******************************************************************
       SUBROUTINE INDICES
!******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
       
      IMPLICIT NONE
      INTEGER*8   I,J,K
      
      !NAVIER-STOKE
!-----STREAMWISE DIRECTION
      DO 10 I=0,N1
      IPV(I)=I+1
      IMV(I)=I-1
   10 CONTINUE
      IPV(N1+1)=N1+2
      IMV(-1)=-2
      
!-----NORMAL DIRECTION
      DO 20 J=0,N2
      JPV(J)=J+1
      JMV(J)=J-1
   20 CONTINUE
      JPV(N2+1)=N2+2
      JMV(-1)=-2
      
!-----SPANWISE DIRECTION
      DO 30 K=0,N3
      KPV(K)=K+1
      KMV(K)=K-1
   30 CONTINUE
      KPV(N3+1)=N3+2
      KMV(-1)=-2
      
      !FOR BOUNDARY CONDITION.
      DO 40 I=1,N1M
      FIXIL(I)=0.
      FIXIU(I)=0.
   40 CONTINUE
      FIXIL(1)=1.
      FIXIU(N1M)=1.

      DO 50 J=1,N2M
      FIXJL(J)=0.
      FIXJU(J)=0.
   50 CONTINUE
      FIXJL(1)=1.
      FIXJU(N2M)=1.

      DO 60 K=1,N3M
      FIXKL(K)=0.
      FIXKU(K)=0.
   60 CONTINUE
      FIXKL(1)=1.
      FIXKU(N3M)=1.


      RETURN
      END

!******************************************************************
       SUBROUTINE INDICES_LVS
!******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_GEOM_VAR
            
      IMPLICIT NONE
      INTEGER*8 I,J,K
      
!------ Axial direction (X)
        DO 40 I=0,N1F
        IPF(I)=I+1
        IMF(I)=I-1
  40    CONTINUE
        IPF(N1F+1)=N1F+2 !FOR DIRICHLET AND NEUMANN CONDTIION IN TRANSPORT, REINITIALIZATION.
        IMF(-1)=-2

!------ Radial direction (R)
        DO 50 J=0,N2F
        JPF(J)=J+1
        JMF(J)=J-1
  50    CONTINUE
        JPF(N2F+1)=N2F+2
        JMF(-1)=-2
        
!------ Azimuthal direction (T)
        DO 60 K=0,N3F
        KPF(K)=K+1
        KMF(K)=K-1
  60    CONTINUE
        KPF(N3F+1)=N3F+2
        KMF(-1)=-2

      RETURN
      END

! =====================================================================
       SUBROUTINE GEOM_COUPLING
! =====================================================================
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_GEOM_VAR
      USE LVS_COUPLING
      
            
      IMPLICIT NONE
      INTEGER*8 I,J,K,II,JJ,KK,IF,JF,KF
      REAL*8    XX,YY,ZZ,X_CRI,Y_CRI,Z_CRI
      
      !FOR FIND FLOW-SOLVER GRID USING REFINED GRID
      !FOR VELOCITY COUPLING

      !X-DIR
      IF ( N1M .EQ. 1 ) THEN
      DO I=1,N1FM
        ICOU_VEL(I)=1
        ICOUMP_VEL(I)=1
       ENDDO
      ELSE

      DO I=1,N1F
         XX=XF(I)
       DO II=1,N1M
        IF ( ABS(XX-XP(II)) .LE. 1.E-10 ) THEN
         ICOU_VEL(I)=II
         GOTO 11
        ELSE IF ( XP(1) .GT. XX ) THEN
         ICOU_VEL(I)=1
         GOTO 11
        ELSE IF ( XP(II) .GT. XX )  THEN
          ICOU_VEL(I)=II-1
         GOTO 11
        ELSE IF ( XP(N1M) .LT. XX ) THEN
         ICOU_VEL(I)=N1M
         GOTO 11
         ENDIF
       ENDDO

          WRITE(*,*) 'X-DIR FLOW_SOLVER->REFINED COUPLING PROBLEM!'
 11   ENDDO
 
      DO I=1,N1F
         XX=XPF(I)
       DO II=1,N1M
        IF ( ABS(XX-XP(II)) .LE. 1.E-10 ) THEN
         ICOUMP_VEL(I)=II
         GOTO 10
        ELSE IF ( XP(1) .GT. XX ) THEN
         ICOUMP_VEL(I)=1
         GOTO 10
        ELSE IF ( XP(II) .GT. XX )  THEN
          ICOUMP_VEL(I)=II-1
         GOTO 10
        ELSE IF ( XP(N1M) .LT. XX ) THEN
         ICOUMP_VEL(I)=N1M
         GOTO 10
         ENDIF
       ENDDO

          WRITE(*,*) 'X-DIR FLOW_SOLVER->REFINED COUPLING PROBLEM!'
 10   ENDDO
      ENDIF
      
      !Y-DIR
      DO J=1,N2F
         YY=YF(J)
       DO JJ=1,N2M
        IF ( ABS(YY-YP(JJ)) .LE. 1.E-10 ) THEN
         JCOU_VEL(J)=JJ
         GOTO 21
        ELSE IF ( YP(1) .GT. YY ) THEN
         JCOU_VEL(J)=1
         GOTO 21
        ELSE IF ( YP(JJ) .GT. YY ) THEN
          JCOU_VEL(J)=JJ-1
         GOTO 21
        ELSE IF ( YP(N2M) .LT. YY ) THEN
          JCOU_VEL(J)=N2M
           GOTO 21
         ENDIF
       ENDDO

          WRITE(*,*) 'Y-DIR FLOW_SOLVER->REFINED COUPLING PROBLEM!'
 21   ENDDO

      DO J=1,N2F
         YY=YPF(J)
       DO JJ=1,N2M
        IF ( ABS(YY-YP(JJ)) .LE. 1.E-10 ) THEN
         JCOUMP_VEL(J)=JJ
         GOTO 20
        ELSE IF ( YP(1) .GT. YY ) THEN
         JCOUMP_VEL(J)=1
         GOTO 20
        ELSE IF ( YP(JJ) .GT. YY ) THEN
          JCOUMP_VEL(J)=JJ-1
         GOTO 20
        ELSE IF ( YP(N2M) .LT. YY ) THEN
           JCOUMP_VEL(J)=N2M
           GOTO 20
         ENDIF
       ENDDO

          WRITE(*,*) 'Y-DIR FLOW_SOLVER->REFINED COUPLING PROBLEM!'
 20   ENDDO
 
      !Z-DIR
      IF ( N3M .EQ. 1 ) THEN
      DO K=1,N3F
        KCOU_VEL(K)=1
        KCOUMP_VEL(K)=1
       ENDDO
      ELSE
      DO K=1,N3F
         ZZ=ZF(K)
       DO KK=1,N3M
        IF ( ABS(ZZ-ZP(KK)) .LE. 1.E-10 ) THEN
         KCOU_VEL(K)=KK
         GOTO 31
        ELSE IF ( ZP(1) .GT. ZZ ) THEN
         KCOU_VEL(K)=1
         GOTO 31
        ELSE IF ( ZP(KK) .GT. ZZ ) THEN
          KCOU_VEL(K)=KK-1
         GOTO 31
        ELSE IF ( ZP(N3M) .LT. ZZ ) THEN
         KCOU_VEL(K)=N3M
         GOTO 31
         ENDIF
       ENDDO

          WRITE(*,*) 'Z-DIR FLOW_SOLVER->REFINED COUPLING PROBLEM!'
 31   ENDDO
 
      DO K=1,N3F
         ZZ=ZPF(K)
       DO KK=1,N3M
        IF ( ABS(ZZ-ZP(KK)) .LE. 1.E-10 ) THEN
         KCOUMP_VEL(K)=KK
         GOTO 30
        ELSE IF ( ZP(1) .GT. ZZ ) THEN
         KCOUMP_VEL(K)=1
         GOTO 30
        ELSE IF ( ZP(KK) .GT. ZZ ) THEN
          KCOUMP_VEL(K)=KK-1
         GOTO 30
        ELSE IF ( ZP(N3M) .LT. ZZ ) THEN
         KCOUMP_VEL(K)=N3M
         GOTO 30
         ENDIF
       ENDDO

          WRITE(*,*) 'Z-DIR FLOW_SOLVER->REFINED COUPLING PROBLEM!'
 30   ENDDO
      ENDIF

!       do i=1,n1f
!      write(*,*) i,icou_vel(i),icoump_vel(i)
!       enddo
!       do i=1,n2f
!      write(*,*) i,jcou_vel(i),jcoump_vel(i)
!       enddo
!       do i=1,n3f
!      write(*,*) i,kcou_vel(i),kcoump_vel(i)
!       enddo
       
      !FOR FIND REFINED GRID USING FLOW-SOLVER GRID
      !just before the flow solver grid
      
      !X-DIR
      DO I=1,N1
       X_CRI=X(I)
        DO IF=1,N1F
         IF ( ABS(XF(IF)-X_CRI) .LE. 1.E-10 ) THEN
          ICOU1(I)=IF-1
          ICOU2(I)=IF
           GOTO 100
         ELSE IF( XF(IF) .GT. X_CRI ) THEN
          ICOU1(I)=IF-1
          ICOU2(I)=IF-1
          GOTO 100
         ENDIF
        ENDDO
           WRITE(*,*) 'CAL ICOU HAVE PROBLEM!!'
      !IF YOU RUN FLOWFIELD.F, THIS IS CAUSED BY THE 
      !DIFFERENCE OF "N1F_BD, N2F_BD, N3F_BD" OF TWO FIELDS
      !%IF YOU USE SPARE GRID FOR SIMULATE PIPE USING IBM,
      !SIMULATION COMES HERE DEPENDING ON THE NUMBER OF REFINEMENTS.
 100  ENDDO
         ICOU1(0)=ICOU1(1)-(ICOU1(N1M)-ICOU1(N1M-1))
         ICOU1(-1)=ICOU1(0)-(ICOU1(N1M-1)-ICOU1(N1M-2))
         ICOU1(N1+1)=ICOU1(N1)+(ICOU1(2)-ICOU1(1))
         ICOU1(N1+2)=ICOU1(N1+1)+(ICOU1(3)-ICOU1(2))
         ICOU2(0)=ICOU2(1)-(ICOU2(N1M)-ICOU2(N1M-1))
         ICOU2(-1)=ICOU2(0)-(ICOU2(N1M-1)-ICOU2(N1M-2))
         ICOU2(N1+1)=ICOU2(N1)+(ICOU2(2)-ICOU2(1))
         ICOU2(N1+2)=ICOU2(N1+1)+(ICOU2(3)-ICOU2(2))

      DO I=1,N1M
       X_CRI=XP(I)
        DO IF=1,N1F
        IF ( ABS(XF(IF)-X_CRI) .LE. 1.E-10 ) THEN
          ICOUMP1(I)=IF-1
          ICOUMP2(I)=IF
           GOTO 101
         ELSE IF( XF(IF) .GT. X_CRI ) THEN
          ICOUMP1(I)=IF-1
          ICOUMP2(I)=IF-1
          GOTO 101
         ENDIF
        ENDDO
           WRITE(*,*) 'CAL ICOUMP HAVE PROBLEM!!'
      !IF YOU RUN FLOWFIELD.F, THIS IS CAUSED BY THE 
      !DIFFERENCE OF "N1F_BD, N2F_BD, N3F_BD" OF TWO FIELDS
      !%IF YOU USE SPARE GRID FOR SIMULATE PIPE USING IBM,
      !SIMULATION COMES HERE DEPENDING ON THE NUMBER OF REFINEMENTS.
 101  ENDDO
         ICOUMP1(0)=ICOUMP1(1)-(ICOUMP1(N1M)-ICOUMP1(N1M-1))
         ICOUMP1(-1)=ICOUMP1(0)-(ICOUMP1(N1M-1)-ICOUMP1(N1M-2))
         ICOUMP1(N1)=ICOUMP1(N1M)+(ICOUMP1(N1M)-ICOUMP1(N1M-1)) !APPROXIMATION
         ICOUMP1(N1+1)=ICOUMP1(N1)+(ICOUMP1(2)-ICOUMP1(1))
         ICOUMP2(0)=ICOUMP2(1)-(ICOUMP2(N1M)-ICOUMP2(N1M-1))
         ICOUMP2(-1)=ICOUMP2(0)-(ICOUMP2(N1M-1)-ICOUMP2(N1M-2))
         ICOUMP2(N1)=ICOUMP2(N1M)+(ICOUMP2(N1M)-ICOUMP2(N1M-1)) !APPROXIMATION
         ICOUMP2(N1+1)=ICOUMP2(N1)+(ICOUMP2(2)-ICOUMP2(1))

      !Y-DIR
      DO J=1,N2
       Y_CRI=Y(J)
        DO JF=1,N2F
        IF ( ABS(YF(JF)-Y_CRI) .LE. 1.E-10 ) THEN
          JCOU1(J)=JF-1
          JCOU2(J)=JF
           GOTO 200
         ELSE IF( YF(JF) .GT. Y_CRI ) THEN
          JCOU1(J)=JF-1
          JCOU2(J)=JF-1
          GOTO 200
         ENDIF
        ENDDO
           WRITE(*,*) 'CAL JCOU HAVE PROBLEM!!'
      !IF YOU RUN FLOWFIELD.F, THIS IS CAUSED BY THE 
      !DIFFERENCE OF "N1F_BD, N2F_BD, N3F_BD" OF TWO FIELDS
      !%IF YOU USE SPARE GRID FOR SIMULATE PIPE USING IBM,
      !SIMULATION COMES HERE DEPENDING ON THE NUMBER OF REFINEMENTS.
 200  ENDDO
         JCOU1(0)=JCOU1(1)-(JCOU1(N2M)-JCOU1(N2M-1))
         JCOU1(-1)=JCOU1(0)-(JCOU1(N2M-1)-JCOU1(N2M-2))
         JCOU1(N2+1)=JCOU1(N2)+(JCOU1(2)-JCOU1(1))
         JCOU1(N2+2)=JCOU1(N2+1)+(JCOU1(3)-JCOU1(2))
         JCOU2(0)=JCOU2(1)-(JCOU2(N2M)-JCOU2(N2M-1))
         JCOU2(-1)=JCOU2(0)-(JCOU2(N2M-1)-JCOU2(N2M-2))
         JCOU2(N2+1)=JCOU2(N2)+(JCOU2(2)-JCOU2(1))
         JCOU2(N2+2)=JCOU2(N2+1)+(JCOU2(3)-JCOU2(2))

      DO J=1,N2M
       Y_CRI=YP(J)
        DO JF=1,N2F
        IF ( ABS(YF(JF)-Y_CRI) .LE. 1.E-10 ) THEN
          JCOUMP1(J)=JF-1
          JCOUMP2(J)=JF
           GOTO 201
         ELSE IF( YF(JF) .GT. Y_CRI ) THEN
          JCOUMP1(J)=JF-1
          JCOUMP2(J)=JF-1
          GOTO 201

         ENDIF
        ENDDO
          WRITE(*,*) 'CAL JCOUMP HAVE PROBLEM!!'
      !IF YOU RUN FLOWFIELD.F, THIS IS CAUSED BY THE 
      !DIFFERENCE OF "N1F_BD, N2F_BD, N3F_BD" OF TWO FIELDS
      !%IF YOU USE SPARE GRID FOR SIMULATE PIPE USING IBM,
      !SIMULATION COMES HERE DEPENDING ON THE NUMBER OF REFINEMENTS.
 201  ENDDO
         JCOUMP1(0)=JCOUMP1(1)-(JCOUMP1(N2M)-JCOUMP1(N2M-1))
         JCOUMP1(-1)=JCOUMP1(0)-(JCOUMP1(N2M-1)-JCOUMP1(N2M-2))
         JCOUMP1(N2)=JCOUMP1(N2M)+(JCOUMP1(N2M)-JCOUMP1(N2M-1)) !APPROXIMATION
         JCOUMP1(N2+1)=JCOUMP1(N2)+(JCOUMP1(2)-JCOUMP1(1))
         JCOUMP2(0)=JCOUMP2(1)-(JCOUMP2(N2M)-JCOUMP2(N2M-1))
         JCOUMP2(-1)=JCOUMP2(0)-(JCOUMP2(N2M-1)-JCOUMP2(N2M-2))
         JCOUMP2(N2)=JCOUMP2(N2M)+(JCOUMP2(N2M)-JCOUMP2(N2M-1)) !APPROXIMATION
         JCOUMP2(N2+1)=JCOUMP2(N2)+(JCOUMP2(2)-JCOUMP2(1))
         
       !Z-DIR
      DO 300 K=1,N3
       Z_CRI=Z(K)
        DO KF=1,N3F
        IF ( ABS(ZF(KF)-Z_CRI) .LE. 1.E-10 ) THEN
          KCOU1(K)=KF-1
          KCOU2(K)=KF
           GOTO 300
         ELSE IF( ZF(KF) .GT. Z_CRI) THEN
          KCOU1(K)=KF-1
          KCOU2(K)=KF-1
          GOTO 300
         ENDIF
        ENDDO
           WRITE(*,*) 'CAL KCOU HAVE PROBLEM!!', K
      !IF YOU RUN FLOWFIELD.F, THIS IS CAUSED BY THE 
      !DIFFERENCE OF "N1F_BD, N2F_BD, N3F_BD" OF TWO FIELDS
      !%IF YOU USE SPARE GRID FOR SIMULATE PIPE USING IBM,
      !SIMULATION COMES HERE DEPENDING ON THE NUMBER OF REFINEMENTS.
 300  ENDDO
         KCOU1(0)=KCOU1(1)-(KCOU1(N3M)-KCOU1(N3M-1))
         KCOU1(-1)=KCOU1(0)-(KCOU1(N3M-1)-KCOU1(N3M-2))
         KCOU1(N3+1)=KCOU1(N3)+(KCOU1(2)-KCOU1(1))
         KCOU1(N3+2)=KCOU1(N3+1)+(KCOU1(3)-KCOU1(2))
         KCOU2(0)=KCOU2(1)-(KCOU2(N3M)-KCOU2(N3M-1))
         KCOU2(-1)=KCOU2(0)-(KCOU2(N3M-1)-KCOU2(N3M-2))
         KCOU2(N3+1)=KCOU2(N3)+(KCOU2(2)-KCOU2(1))
         KCOU2(N3+2)=KCOU2(N3+1)+(KCOU2(3)-KCOU2(2))
         
      DO 301 K=1,N3M
       Z_CRI=ZP(K)
        DO KF=1,N3F
        IF ( ABS(ZF(KF)-Z_CRI) .LE. 1.E-10 ) THEN
          KCOUMP1(K)=KF-1
          KCOUMP2(K)=KF
           GOTO 301
         ELSE IF( ZF(KF) .GT. Z_CRI) THEN
          KCOUMP1(K)=KF-1
          KCOUMP2(K)=KF-1
          GOTO 301
         ENDIF
        ENDDO
           WRITE(*,*) 'CAL KCOUMP HAVE PROBLEM!!'
      !IF YOU RUN FLOWFIELD.F, THIS IS CAUSED BY THE 
      !DIFFERENCE OF "N1F_BD, N2F_BD, N3F_BD" OF TWO FIELDS
      !%IF YOU USE SPARE GRID FOR SIMULATE PIPE USING IBM,
      !SIMULATION COMES HERE DEPENDING ON THE NUMBER OF REFINEMENTS.
 301  ENDDO
         KCOUMP1(0)=KCOUMP1(1)-(KCOUMP1(N3M)-KCOUMP1(N3M-1))
         KCOUMP1(-1)=KCOUMP1(0)-(KCOUMP1(N3M-1)-KCOUMP1(N3M-2))
         KCOUMP1(N3)=KCOUMP1(N3M)+(KCOUMP1(N3M)-KCOUMP1(N3M-1)) !APPROXIMATION
         KCOUMP1(N3+1)=KCOUMP1(N3)+(KCOUMP1(2)-KCOUMP1(1))
         KCOUMP2(0)=KCOUMP2(1)-(KCOUMP2(N3M)-KCOUMP2(N3M-1))
         KCOUMP2(-1)=KCOUMP2(0)-(KCOUMP2(N3M-1)-KCOUMP2(N3M-2))
         KCOUMP2(N3)=KCOUMP2(N3M)+(KCOUMP2(N3M)-KCOUMP2(N3M-1)) !APPROXIMATION
         KCOUMP2(N3+1)=KCOUMP2(N3)+(KCOUMP2(2)-KCOUMP2(1))
      
!      DO I=1,N2
!       WRITE(*,678) I,ICOU1(I),JCOU1(I),KCOU1(I)
!      ENDDO
!      DO I=1,N2
!       WRITE(*,678) I,ICOUMP1(I),JCOUMP1(I),KCOUMP1(I)
!      ENDDO
!      DO I=1,N2
!       WRITE(*,678) I,ICOUMP2(I),JCOUMP2(I),KCOUMP2(I)
!      ENDDO
! 678   FORMAT(4I5)  
 
        RETURN
        END   