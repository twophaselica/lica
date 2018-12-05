c-----compile(MACHINE)
cf90 lica.f geom.f lica_sgs.f pcg.f two_phase.f -g -traceback -fpe0 -openmp -o 1


cf90 modules.f lica.f geom.f lica_sgs.f pcg.f two_phase.f -g -traceback -fpe0 -openmp -o 1


c-----file_transfer
C supercomputer -> machine
c  scp -r x1185kky@gaiad.ksc.re.kr:/gpfs2/x1185kky/141118_ret=1024/ /home/kykim/0bubbly_flow
c  scp x1185kky@gaiad.ksc.re.kr:/gpfs2/x1185kky/141103_caf_0.3/fld072000 /home/kykim/0bubbly_flow/0.3

c macine -> supercomputer
c  scp -r /home/kykim/0bubbly_flow x1185kky@gaiad.ksc.re.kr:/gpfs2/x1185kky/
c  scp /home/kykim/0bubbly_flow/fld129000 x1185kky@gaiad.ksc.re.kr:/gpfs2/x1185kky

c machine -> machine
c  scp -r -P 41118 /home/kykim/150518_les_eddy/ kykim@karman18.snu.ac.kr:/home/kykim/0bubbly_flow
c  scp -P 41118 /home/kykim/141103_caf_0.3/fld283000 kykim@karman18.snu.ac.kr:/home/kykim/0bubbly_flow

C supercomputer -> supercomputer
c scp -r x1185kky@gaiad.ksc.re.kr:/gpfs2/x1185kky/141103_caf_0.3 /gpfs2/x1185kky/141012_caf_0.3


!***********************************************************************
!               Two-phase Large eddy simulation code
!   incorporated with Immersed boundary method in a CArtesian coordinate
!***********************************************************************
!
!     MAIN PROPERTIES
!     1. Cartesian coordinate
!     2. IBM
!     3. DNS/LES
!     4. Two-phase flow
!
!     Users and/or Contributors
!
!     Haecheon Choi
!     Sungwon Kang : Backward facing step
!     Dongju Kim   : Flow over sphere with IBM as a first try.
!                    Later he moved on cylindrical coordinate
!                    to simulate flow over sphere
!     Woongjae Jang: Finite cylinder
!     Jungil Lee   : GM model, Circular cylinder
!     Kiyoung Kim  : Two-phase flow
!
!
!     2015.10.16 -'COMMON' IS REMOVED AND MODULE IS USED  BY KIYOUNG KIM     
!***********************************************************************

      PROGRAM MAIN
      USE FLOW_VAR
      USE MG_OPTION
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      USE LVS_COUPLING

      USE FLD_AVG

      USE IBM_VAR

      USE HEAT_VAR
      
      USE TIME_VAR
      
      USE LES_VAR

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      REAL PSI_XN(0:M1,0:M2,0:M3)
      REAL PSI_YN(0:M1,0:M2,0:M3)
      REAL PSI_ZN(0:M1,0:M2,0:M3)
      REAL PSI_CN(0:M1,0:M2,0:M3)

      REAL VOL_TOT_ORI(MLVS)

      REAL QVOL_ORI
      REAL VOL1 !TOTAL FLUID VOLUME

      !FOR INITIAL MEANPGR
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: DENF_Z

      !RK3_OLD
      REAL, ALLOCATABLE :: RK3XO(:,:,:),RK3YO(:,:,:),RK3ZO(:,:,:) 
      
      !RK3_OLD_ENERGY
      REAL, ALLOCATABLE :: ANPSO(:,:,:) 

      !TRACE
      INTEGER NTEST,NTPRINT,NTEST2
      INTEGER ITR(10000),JTR(10000),KTR(10000)
      INTEGER ITR2(10000),JTR2(10000),KTR2(10000)

      CHARACTER*30 gridfile,fileprevel
      CHARACTER*10 DUMMY
      
      INTEGER IREAD

      call timestamp

      CALL real_time(TOTAL_TIME_B)

C=====read input file
      OPEN(10,FILE='lica.in')
      READ(10,*) DUMMY
      READ(10,*) IRESET,IREAD,IAVG,IPZERO,IDTOLD
      READ(10,*) DUMMY
      READ(10,*) NTST,NPRINT,NPRIAVG,NPIN
      READ(10,*) DUMMY
      READ(10,*) IMPL,IDTOPT,DT,CFLMAX
      READ(10,*) DUMMY
      READ(10,*) RESID_POI,NLEV,IWC,NBLI,MGITR,IOLDV,IMGSOR,WWSOR
      READ(10,*) DUMMY
      READ(10,*) IPS
      READ(10,*) DUMMY
      READ(10,*) ILES,ISGSMOD,DYNMON,CSGS,CSGSPS,IFILTER
      READ(10,*) DUMMY
      READ(10,*) IBMON,MASSON
      READ(10,*) DUMMY
      READ(10,*) IPOISS
      READ(10,*) DUMMY
      READ(10,*) ICH,IPX,IPY,IPZ
      READ(10,*) DUMMY
      READ(10,*) DUMMY
      READ(10,*) ITRACKING,MCLS,MGLOBAL
      READ(10,*) DUMMY
      READ(10,*) RE_AIR,VISR,DENR
      READ(10,*) DUMMY
      READ(10,*) SURF_J,FR,FK
      READ(10,*) DUMMY
      READ(10,*) PRM,SCR,TCR
      READ(10,*) DUMMY
      READ(10,*) DUMMY
      READ(10,'(a)') gridfile
      READ(10,'(a)') fileprevel
      READ(10,*) DUMMY
      READ(10,*) NTEST,NTPRINT

      WRITE(*,101) IRESET,IREAD,IAVG,IDTOLD
      WRITE(*,103) NTST,NPRINT,NPRIAVG,NPIN
      WRITE(*,105) IMPL,IDTOPT,DT,CFLMAX
      WRITE(*,109) RESID_POI,NLEV,IWC,NBLI,MGITR,IOLDV,IMGSOR,WWSOR
      WRITE(*,110) IPS
      WRITE(*,111) ILES,ISGSMOD,DYNMON,CSGS,CSGSPS,IFILTER
      WRITE(*,113) IBMON,MASSON
      WRITE(*,115) IPOISS
      WRITE(*,117) ICH,IPX,IPY,IPZ
      WRITE(*,119) ITRACKING,MCLS,MGLOBAL
      WRITE(*,121) RE_AIR,VISR,DENR
      WRITE(*,123) SURF_J,FR,FK
      WRITE(*,125) gridfile
      IF (IREAD .EQ. 1) WRITE(*,127) fileprevel
     
 101  FORMAT('IRESET=',I8,'  IREAD=',I9,'  IAVG=',I10,
     &       '  IPZERO=',I8,'  IDTOLD=',F7.2)
 103  FORMAT('NTST=',I10,'  NPRINT=',I8,'  NPRIAVG=',I7,'  NPIN=',I10)
 105  FORMAT('IMPL=',I10,'  IDTOPT=',I8,'  DT=',ES12.4,
     &       '  CFLMAX=',F8.2)
 107  FORMAT('IUD=',I11,'  IEND=',I10,'  IENDM=',I9,'  ALPZ=',F10.4)
 109  FORMAT('RESID=',ES9.2,'  NLEV=',I10,'  IWC=',I11,'  NBLI=',I10
     &      ,'  MGITR=',I9,'  IOLDV=',I9,'  IMGSOR=',I8,'  WWSOR=',F9.3)
 110  FORMAT('IPS=',I11)
 111  FORMAT('ILES=',I10,'  ISGSMOD=',I7,'  DYNMON=',L8,
     &       '  CSGS=',F10.4,'  CSGSPS=',F8.4,'  IFILTER=',I7)
 113  FORMAT('IBMON=',I9,'  MASSON=',I8)
 115  FORMAT('IPOISS=',I9)
 117  FORMAT('ICH=',I9,'  IPX=',I8,'  IPY=',I8,'  IPZ=',I8)
 119  FORMAT('ITRACKING=',I9,'  MCLS=',I8,'  MGLOBAL=',I8)
 121  FORMAT('RE_AIR=',F12.5,'  VISR=',F12.5,'  DENR=',F12.5)
 123  FORMAT('SUR_J=',F12.5,'  FR=',F12.5,'  FK=',F12.5)
 125  FORMAT('GRID=   ',A30)
 127  FORMAT('FIELD_READ=   ',A30) 

 100  FORMAT('INITIAL FIELD TREATMNT:                READING DONE')

 300  FORMAT(I15,'   LOCAL RESIDUE FOR MASS FLOW CONSV.:',ES18.5,
     &     ' AT ',ES11.4)
 310  FORMAT(I15,'   CFL# EXCEED GIVEN CFL LIMIT :',ES18.5,
     &     ' AT TIME',F12.5)
 320  FORMAT('--------------------------',I6,'  TIME=',F10.5,
     &       '  DT=',F12.8)

      
C------------------MEMORY ALLOCATAION & INITIALIZATION-----------------C
!$OMP PARALLEL DO private(I,J)
       DO K=0,N3
       DO J=0,N2
       DO I=0,N1
        U(I,J,K)=0.
        V(I,J,K)=0.
        W(I,J,K)=0.
        PSI_XN(I,J,K)=0.
        PSI_YN(I,J,K)=0.
        PSI_ZN(I,J,K)=0.
        PSI_CN(I,J,K)=0.
       ENDDO
       ENDDO
       ENDDO

!$OMP PARALLEL DO
       DO LLVS=1,MLVS
       VOL_TOT_ORI(LLVS)=0.
       ENDDO

!$OMP PARALLEL DO
       DO L=1,10000
        ITR(L)=0 ;JTR(L)=0 ;KTR(L)=0 ;ITR2(L)=0 ;JTR2(L)=0 ;KTR2(L)=0
       ENDDO

C-----MODULE
      !FLOW_VAR
      ALLOCATE(UO(0:M1,0:M2,0:M3),VO(0:M1,0:M2,0:M3),WO(0:M1,0:M2,0:M3))
      ALLOCATE(AIU(M1),CIU(M1),AIVW(M1),CIVW(M1)
     &        ,AJV(M2),CJV(M2),AJUW(M2),CJUW(M2)
     &        ,AKW(M3),CKW(M3),AKUV(M3),CKUV(M3))

!$OMP PARALLEL DO private(I,J)
       DO K=0,N3
       DO J=0,N2
       DO I=0,N1
        UO(I,J,K)=0. ;VO(I,J,K)=0.; WO(I,J,K)=0.
       ENDDO ;ENDDO ;ENDDO

!$OMP PARALLEL DO
       DO I=1,N1
        AIU(I)=0. ;CIU(I)=0. ;AIVW(I)=0. ;CIVW(I)=0.
       ENDDO
!$OMP PARALLEL DO
       DO J=1,N2
        AJV(J)=0. ;CJV(J)=0. ;AJUW(J)=0. ;CJUW(J)=0.
       ENDDO
!$OMP PARALLEL DO
       DO K=1,N3
        AKW(K)=0. ;CKW(K)=0. ;AKUV(K)=0. ;CKUV(K)=0.
       ENDDO
       
      !LVS_VAR
      ALLOCATE( ALPHI(MF_BAND,M_MAX,MLVS))
      ALLOCATE( NUMA(M_MAX,MLVS))
      ALLOCATE( I_B(MF_BAND,M_MAX,MLVS),J_B(MF_BAND,M_MAX,MLVS)
     &                           ,K_B(MF_BAND,M_MAX,MLVS))
      ALLOCATE( PSIF(M1L:M1U,M2L:M2U,M3L:M3U))

      ALLOCATE(MASK_BUB(M1M,M2M,M3M))
      ALLOCATE(MASK_GLOBAL(0:M1F,0:M2F,0:M3F))

       DO LLVS=1,MLVS
       DO N=1,M_MAX
        NUMA(N,LLVS)=0
!$OMP PARALLEL DO
       DO NN=1,MF_BAND
        ALPHI(NN,N,LLVS)=0.
        I_B(NN,N,LLVS)=0 ;J_B(NN,N,LLVS)=0 ;K_B(NN,N,LLVS)=0
       ENDDO
       ENDDO
       ENDDO

!$OMP PARALLEL DO private(I,J)
       DO K=M3L,M3U
       DO J=M2L,M2U
       DO I=M1L,M1U
        PSIF(I,J,K)=0.
       ENDDO
       ENDDO
       ENDDO

!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
        MASK_BUB(I,J,K)=0
       ENDDO
       ENDDO
       ENDDO

!$OMP PARALLEL DO private(I,J)
       DO K=0,N3F
       DO J=0,N2F
       DO I=0,N1F
        MASK_GLOBAL(I,J,K)=0
       ENDDO
       ENDDO
       ENDDO
       
      !IBM_VAR
      ALLOCATE(QMASS(0:M1,0:M2,0:M3))

!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
        QMASS(I,J,K)=0.
       ENDDO
       ENDDO
       ENDDO

      !HEAT_VAR
      IF (IPS .EQ. 1) THEN
      ALLOCATE(T(0:M1,0:M2,0:M3))
      ALLOCATE(TALPH(0:M1,0:M2,0:M3))
      ALLOCATE(RHSPS(M1M,M2M,M3M))

!$OMP PARALLEL DO private(I,J)
       DO K=0,N3 ;DO J=0,N2 ;DO I=0,N1
        T(I,J,K)=0.
        TALPH(I,J,K)=0.
       ENDDO ;ENDDO ;ENDDO

!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M ;DO J=1,N2M ;DO I=1,N1M
        RHSPS(I,J,K)=0.
       ENDDO ;ENDDO ;ENDDO
      ENDIF
      
      !TIME_VAR
      ALLOCATE(SGSTIME_B(3),SGSTIME_E(3)
     &   ,ANSETIME_B(100),ANSETIME_E(100),POISSTIME_B(3),POISSTIME_E(3)
     &   ,CONTINUITY_B(3),CONTINUITY_E(3),ALVS_B(3),ALVS_E(3))
     
       DO L=1,3
        SGSTIME_B(L)=0. ;SGSTIME_E(L)=0.
        POISSTIME_B(L)=0. ;POISSTIME_E(L)=0.
        CONTINUITY_B(L)=0. ;CONTINUITY_E(L)=0.
        ALVS_B(L)=0. ;ALVS_E(L)=0.
       ENDDO
      
       DO L=1,100
        ANSETIME_B(L)=0. ;ANSETIME_E(L)=0
       ENDDO
     
      !LES_VAR
      IF (ILES .EQ. 1) THEN
      ALLOCATE(CFX1 (M1M,-1:1),CFX2 (M1M,-2:2),
     &                 CFZP1(M3M,-1:1),CFZP2(M3M,-2:2))
      ALLOCATE(TNU(0:M1,0:M2,0:M3))
!      ALLOCATE(TNUX(0:M1,0:M2,0:M3)
!     &,TNUY(0:M1,0:M2,0:M3),TNUZ(0:M1,0:M2,0:M3))

      DO L=-1,1
!$OMP PARALLEL DO
      DO I=1,N1M
      CFX1 (I,L)=0.
      CFX2 (I,L)=0.
      ENDDO

!$OMP PARALLEL DO
      DO K=1,N3M
      CFZP1(K,L)=0.
      CFZP2(K,L)=0.
      ENDDO
      ENDDO

!$OMP PARALLEL DO private(I,J)
       DO K=0,N3 ;DO J=0,N2 ;DO I=0,N1
      TNU(I,J,K)=0.
!     TNUX(I,J,K)=0. ;TNUY(I,J,K)=0. ;TNUZ(I,J,K)=0.
      ENDDO ;ENDDO ;ENDDO
      ENDIF

C------------------MEMORY ALLOCATAION & INITIALIZATION-----------------C

       IF (ICH .EQ. 0 ) THEN
         WRITE(*,*) 'NO-MEAN_PRESSURE GRADIENT'
       ELSE IF (ICH .EQ. 1 ) THEN
         WRITE(*,*) 'MASS-FLOW RATE CONST'
       ELSE IF (ICH .EQ. 2 ) THEN
         WRITE(*,*) 'PRESSURE GRADIENT CONST'
       ENDIF

      CALL GEOM(gridfile)

!!!!!!INITIALIZATION
      SGSTIME_B=0.
      SGSTIME_E=0.
      ANSETIME_B=0.
      ANSETIME_E=0.
      POISSTIME_B=0.
      POISSTIME_E=0.
      CONTINUITY_B=0.
      CONTINUITY_E=0.
      ALVS_B=0.
      ALVS_E=0.
!!!!!!INITIALIZATION

      vol1=0.
!$OMP PARALLEL DO private(I,J)
!$OMP&reduction(+:VOL1)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
       IF ( FUNCBODY(XP(I),YP(J),Z(K)) .GE. 0. )
     &                          vol1=vol1+SDX(I)*SDY(J)*VDZ(K)
      ENDDO
      ENDDO
      ENDDO

      PI=ACOS(-1.)
C============================FOR TWO-PHASE=============================C
      RE=RE_AIR
      REI=1./RE

      DENM=1.        !NON-DIMENSIONALIZED VIS
      DENP=DENR
      DEN_DIFF=DENP-DENM

      VISM=1.*REI    !NON-DIMENSIONALIZED VIS
      VISP=VISR*REI
      VIS_DIFF=VISP-VISM

      REM=RE_AIR
      REP=REM*DENR/VISR
      
      DENSCM=DENM*1.   !SPECIFIC HEAT CAPACITY
      DENSCP=DENSCM*(DENR*SCR)
      DENSC_DIFF=DENSCP-DENSCM
      
      TCM=1.   !THERMAL CONDUCTIVITY
      TCP=TCR
      TC_DIFF=TCP-TCM
      
      PRM=PRM !PRANDTL NUMBER
      PRP=PRM*VISR*SCR/TCR
      
      PRA=PRM
      PRAI=1./PRM
      IF ( SURF_J .EQ. 0. ) THEN
       SURF=0.
      ELSE
       SURF=1./SURF_J  !SURF_J IS WEBER NUMBER
      ENDIF

      IF (ICH .EQ. 2) THEN
       PMI_CONST=-FK
      ELSE
        PMI_CONST=0.
      ENDIF

      IF(FR .EQ. 0.) THEN
        GRAVITY=0.
      ELSE
        GRAVITY=1./FR**2
      ENDIF

      TIME_ST=SURF*PI/(DENM+DENP)
      TIME_GV=DEN_DIFF/(DENM+DENP)*GRAVITY/PI

      PRINT*,'MEAN PRESSURE GRADIENT :',FK
      PRINT*,'SURF_J :',SURF_J
      PRINT*,'FR :',FR
      PRINT*,'DENSITY_RATIO :',DENR
      PRINT*,'VISCOSITY_RATIO :',VISR
      PRINT*,'RE_AIR :',REM,'RE_WATER :',REP
      PRINT*,'PR_AIR :',PRM,'PR_WATER :',PRP
C============================FOR TWO-PHASE=============================C
      IF (NTEST     .NE. 0) THEN
      	CALL TRACEINIT(NTEST,NTPRINT,NTEST2,ITR,JTR,KTR,ITR2,JTR2,KTR2)
      ELSE
        NTEST2=0
      ENDIF
      IF (IBMON.NE.0) THEN
      CALL IBMINIT
      IF (IPS .EQ. 1) CALL IBMINIT_PS
      IF (ILES .EQ. 1) CALL NUTZEROREAD
      ENDIF

      IF (ILES .EQ. 1) CALL SGSFILTERINIT

C=====read or make initial field
      IF (IREAD.NE.0) THEN
        OPEN(12,FILE=fileprevel)
        CALL PREFLD(IHIST,fileprevel,IDTOLD,U,V,W,P,IPZERO,VOL_TOT_ORI
     &,QVOL_ORI)

        CLOSE(12)
      ELSE

       CALL MAKEFLD(IHIST,NTII,U,V,W,P)
      ENDIF

      IF (ITRACKING .EQ. 0) THEN

      DO LLVS=1,NLVS
       DO N=1,N_MAX
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,MF_BAND
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
          ALPHI(NN,N,LLVS)=0.
        ENDDO
       ENDDO
      ENDDO

!$OMP PARALLEL DO private(I,J)
       DO K=0,N3
       DO J=0,N2
       DO I=0,N1
         PSI_XN(I,J,K)=0.
         PSI_YN(I,J,K)=0.
         PSI_ZN(I,J,K)=0.
         PSI_CN(I,J,K)=0.
       ENDDO
       ENDDO
       ENDDO
        WRITE(*,*) 'LVS IS NOT SOLVED, ALL LVS FUNCGION IS ZERO'
      ALLOCATE(SUR_X(M1M,M2M,M3M),SUR_Y(M1M,M2M,M3M),SUR_Z(M1M,M2M,M3M))
!$OMP PARALLEL DO private(I,J)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
       SUR_X(I,J,K)=0.
       SUR_Y(I,J,K)=0.
       SUR_Z(I,J,K)=0.
      ENDDO
      ENDDO
      ENDDO
      
       ELSE

      NLVS=MLVS
      N_MAX=M_MAX

      !TIME_VARING_PROPERTIES

      !INPUT
      ITIME_VARING_PROPERTIES=0

      TIME_INTERVAL=0.5
      DENP_OLD=831.667
      VISP_OLD=55.56/1466.
      VISM_OLD=1./1466.
      SURF_OLD=1./0.41485
      !INPUT

      IF (ITIME_VARING_PROPERTIES .EQ. 1) THEN

      TIME_OLD=TIME
      T_TAR=TIME+TIME_INTERVAL

      DENP_NEW=DENP
      VISP_NEW=VISP
      VISM_NEW=VISM
      SURF_NEW=SURF
      
      DENP=DENP_OLD
      DEN_DIFF=DENP_OLD-DENM
      VISP=VISP_OLD
      VISM=VISM_OLD
      VIS_DIFF=VISP_OLD-VISM_OLD
      SURF=SURF_OLD
      ENDIF !IF (ITIME_VARING_PROPERTIES .EQ. 1) THEN

        CALL LVSINIT(ITRACKING,U,V,W,PSI_CN,VOL_TOT_ORI)
        CALL GRID_COUPLING(1,PSI_XN,PSI_YN,PSI_ZN,PSI_CN)
       ENDIF

        IF ( ICH .EQ. 0 ) THEN
          PMI=0.
        ELSE
          IF (IRESET .EQ. 1 ) THEN
            OPEN(96,FILE='0MEANP.dat')
          ELSE
            OPEN(96,FILE='0MEANP.dat',POSITION='APPEND')
          ENDIF
        ENDIF

       ALLOCATE (DENF_Z(0:M1,0:M2,0:M3))
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
        DENF_Z(I,J,K)=DENM+DEN_DIFF*PSI_ZN(I,J,K)
       ENDDO
       ENDDO
       ENDDO

       IF(ICH .NE. 0 ) CALL MEANPGR(W,PSI_ZN,PSI_CN,DENF_Z,VOL1) !after LVSINIT FOR AVERAGED DENSITY

       IF(ICH.EQ.1) THEN
        IF (IRESET .EQ. 1) THEN !ELSE READ AT THE FIELD
         CALL QVOLCALC(QVOL_ORI,1,W,PSI_ZN,DENF_Z,VOL1)
         WRITE(*,*) 'TOTAL MASS FLOW RATE IS CONSTANT'
         WRITE(*,*) 'INITIAL MASS FLOW RATE = ',QVOL_ORI
        ENDIF
       ENDIF
       DEALLOCATE (DENF_Z)

C=====initialize others
      IF (IRESET .EQ. 1) THEN
         IHIST=0
         TIME=0.
      END IF
      NTIME=0

      CALL RK3COEF
      IF(ILES .EQ. 1) CALL SGSINIT
!      IF(IPOISS .EQ. 0) CALL POISINIT
      CALL LHSINIT
      IF (IAVG .EQ. 1) CALL FIELD_AVG_INIT_RT(IHIST)
      IF (IAVG .EQ. 2) CALL FIELD_AVG_INIT_XYZ(IHIST)

      CFLM=CFLMAX

      NV=101                  ! INITIAL FILE NUMBER(BINARY)
      NAV=2001
      L1=3001
      L2=5001

      IF (IREAD.EQ.0) THEN

      OPEN(13,FILE='0ftrhist.dat')
      IF (ILES  .EQ. 1) OPEN(89,FILE='0ftrles.dat')
      IF (ILES  .EQ. 1) OPEN(95,FILE='0ftrles2.dat')
      IF (NTEST .NE. 0) OPEN(91,FILE='0ftru.dat')
      IF (NTEST .NE. 0) OPEN(92,FILE='0ftrv.dat')
      IF (NTEST .NE. 0) OPEN(93,FILE='0ftrw.dat')
      IF (NTEST .NE. 0) OPEN(94,FILE='0ftrp.dat')
      IF (IPS   .EQ. 1) OPEN(97,FILE='0ftrps.dat')
      OPEN(87,FILE='0ftrtime.dat')

      ELSE IF (IREAD.EQ.1) THEN
      OPEN(13,FILE='0ftrhist.dat',POSITION='APPEND')
      IF (ILES  .EQ. 1) OPEN(89,FILE='0ftrles.dat',POSITION='APPEND')
      IF (ILES  .EQ. 1) OPEN(95,FILE='0ftrles2.dat',POSITION='APPEND')
      IF (NTEST .NE. 0) OPEN(91,FILE='0ftru.dat',POSITION='APPEND')
      IF (NTEST .NE. 0) OPEN(92,FILE='0ftrv.dat',POSITION='APPEND')
      IF (NTEST .NE. 0) OPEN(93,FILE='0ftrw.dat',POSITION='APPEND')
      IF (NTEST .NE. 0) OPEN(94,FILE='0ftrp.dat',POSITION='APPEND')
      IF (IPS   .EQ. 1) OPEN(97,FILE='0ftrps.dat',POSITION='APPEND')
      OPEN(87,FILE='0ftrtime.dat',POSITION='APPEND')
      IF (NTEST2.NE. 0) OPEN(81,FILE='0ftru.dat',POSITION='APPEND')
      IF (NTEST2.NE. 0) OPEN(82,FILE='0ftrv.dat',POSITION='APPEND')
      IF (NTEST2.NE. 0) OPEN(83,FILE='0ftrw.dat',POSITION='APPEND')
      IF (NTEST2.NE. 0) OPEN(84,FILE='0ftrp.dat',POSITION='APPEND')

      ENDIF

C=====start time dependent calculation
      DO 3000 M=1,NTST              ! TOTAL ITERATION

      CALL real_time(TIME_BEGIN)

      IF (ITIME_VARING_PROPERTIES .EQ. 1) THEN
      IF (TIME .LE. T_TAR) THEN
      DENP=DENP_OLD+(DENP_NEW-DENP_OLD)/(T_TAR-TIME_OLD)*(TIME-TIME_OLD)
      VISP=VISP_OLD+(VISP_NEW-VISP_OLD)/(T_TAR-TIME_OLD)*(TIME-TIME_OLD)
      VISM=VISM_OLD+(VISM_NEW-VISM_OLD)/(T_TAR-TIME_OLD)*(TIME-TIME_OLD)
      SURF=SURF_OLD+(SURF_NEW-SURF_OLD)/(T_TAR-TIME_OLD)*(TIME-TIME_OLD)
      ELSE
       DENP=DENP_NEW
       VISP=VISP_NEW
       VISM=VISM_NEW
       SURF=SURF_NEW
       
      DEN_DIFF=DENP_NEW-DENM
      VIS_DIFF=VISP_NEW-VISM_NEW
      ENDIF
      ENDIF

      NTIME=NTIME+1

      IF (IBMON.NE.0) THEN
         FCVAVG=0.
      ENDIF

C-----determine time step (DT)
      DT_OLD=DT

      !this is coded refering 'capillary wave' at wikipedia
      IF( TIME_ST .EQ. 0. ) THEN
      DT_SOURCE=1000000.  !VERY SMALL NUMBER FOR AVODING DEVIDED BY ZERO
      ELSE
      DT_SOURCE=1000.
       IF (N3M .EQ. 1 ) THEN
      K=1
!$OMP PARALLEL DO private(I,DX1,DX2,DT_TMP1,DT_TMP2)
!$OMP&firstprivate(K)
!$OMP&reduction(MIN:DT_SOURCE)
      DO J=1,N2M
      DO I=1,N1M
        IF ( ABS(PSI_CN(I,J,K)-0.5) .LT. 0.5 ) THEN
           DX1=MIN(SDX(I),SDY(J))
           DT_TMP1=0.5*DX1/SQRT(TIME_ST/DX1+TIME_GV*DX1)
           DX2=MAX(SDX(I),SDY(J))
           DT_TMP2=0.5*DX2/SQRT(TIME_ST/DX2+TIME_GV*DX2)
           DT_SOURCE=MIN(DT_SOURCE,DT_TMP1,DT_TMP2)
        ENDIF
      ENDDO
      ENDDO
       ELSE IF (N1M .EQ. 1 ) THEN
      I=1
!$OMP PARALLEL DO private(J,DX1,DX2,DT_TMP1,DT_TMP2)
!$OMP&firstprivate(I)
!$OMP&reduction(MIN:DT_SOURCE)
      DO K=1,N3M
      DO J=1,N2M
        IF ( ABS(PSI_CN(I,J,K)-0.5) .LT. 0.5 ) THEN
           DX1=MIN(SDY(J),SDZ(K))
           DT_TMP1=0.5*DX1/SQRT(TIME_ST/DX1+TIME_GV*DX1)
           DX2=MAX(SDY(J),SDZ(K))
           DT_TMP2=0.5*DX2/SQRT(TIME_ST/DX2+TIME_GV*DX2)
           DT_SOURCE=MIN(DT_SOURCE,DT_TMP1,DT_TMP2)
        ENDIF
      ENDDO
      ENDDO
       ELSE
!$OMP PARALLEL DO private(I,J,DX1,DX2,DT_TMP1,DT_TMP2)
!$OMP&reduction(MIN:DT_SOURCE)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
        IF ( ABS(PSI_CN(I,J,K)-0.5) .LT. 0.5 ) THEN
             DX1=MIN(SDX(I),SDY(J),SDZ(K))
           DT_TMP1=0.5*DX1/SQRT(TIME_ST/DX1+TIME_GV*DX1)
             DX2=MAX(SDX(I),SDY(J),SDZ(K))
           DT_TMP2=0.5*DX2/SQRT(TIME_ST/DX2+TIME_GV*DX2)
           DT_SOURCE=MIN(DT_SOURCE,DT_TMP1,DT_TMP2)
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDIF

      ENDIF

      CALL CFL(IMPL,CFLM,U,V,W)

      CFLCONVMAX=CFLMAX
      CFLSURFMAX=CFLMAX

      CFL_SURF=DT/DT_SOURCE
      IF (IDTOPT .EQ. 1) THEN
C-----DT_DETERMINE 1ND TYPE
        CFL_SUM=0.5*(CFLM+SQRT(CFLM**2+4*CFL_SURF))
        DT=DT*CFLCONVMAX/CFL_SUM
!C-----DT_DETERMINE 2ND TYPE
!        CFL_SUM=CFLM+CFL_SURF
!        DT=DT*AMIN1(CFLCONVMAX/CFLM,CFLSURFMAX/CFL_SURF)
      ENDIF

      IF ( DT/DT_OLD .GT. 2. ) THEN
        DT=DT_OLD*1.2    !CHECK LATTER, SMOOTH STARTING WHEN INITIAL VEL. IS ZERO.
        WRITE(*,*) 'DT CHANGE IS LARGE ALGORITHM IS ACTIVATED'
      ENDIF

      CFL_MAX=MAX(CFLM,CFL_SURF)
      IF (CFL_MAX .GT.(CFLMAX*1.1)) THEN
      PRINT*,' '
!      WRITE(*,310) NTIME,CFLM,TIME
      WRITE(*,320) NTIME,TIME,DT
      WRITE(*,*) 'CFL NUMBER IS EXCEEDED!!'
      ELSE
      PRINT*,' '
      WRITE(*,320) NTIME,TIME,DT
      ENDIF

      WRITE(*,153) CFLM,CFL_SURF,CFL_SUM
  153 FORMAT('CFL_CON= ',F7.4,'  CFL_SURF= ',F7.4,'  CFLSUM= ',F7.4)
C-----MAIN SOLVER-------------------------------------------------------
      IF ( IMPL .EQ. 0 ) THEN
        ITER_NS=3

      ALLOCATE(RK3XO(M1M,M2M,M3M),RK3YO(M1M,M2M,M3M),RK3ZO(M1M,M2M,M3M))
      ALLOCATE(ANPSO(M1M,M2M,M3M))
!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
        RK3XO(I,J,K)=0.
        RK3YO(I,J,K)=0.
        RK3ZO(I,J,K)=0.
        ANPSO(I,J,K)=0.
       ENDDO
       ENDDO
       ENDDO

      DO 2000 MSUB=1,3                    ! SUB-ITERATION(K=1,2,3)
      ALPHA_RK3=0.5*(GAMMA(MSUB)+RO(MSUB))
      ACOEF=ALPHA_RK3*DT
      ACOEFI=1./ACOEF
      DTCONST=2.*ALPHA_RK3*DT
      DTCONSTI=1./DTCONST

      CALL real_time(SGSTIME_B(MSUB))
      IF(ILES .EQ. 1) CALL SGSCALC(U,V,W,PSI_CN,IFILTER,IREAD)
!      IF (ILES .EQ. 1) CALL RHSSGS !we just modify VIS
      CALL real_time(SGSTIME_E(MSUB))

      CALL RHSNLHS(ITRACKING,IPOISS,U,V,W,P,PSI_XN,PSI_YN,PSI_ZN,PSI_CN
     &,RK3XO,RK3YO,RK3ZO,ANPSO,QVOL_ORI,VOL1)
!        write(*,*) 'no navier-stokes'

!        T=8.
!!$OMP PARALLEL DO private(I,J,xx,yy)
!        DO K=1,N3M
!        DO J=1,N2M
!        DO I=1,N1M
!        XX=X(I)
!        YY=YP(J)
!         U(I,J,K)=-2.*SIN(PI*XX)**2*SIN(PI*YY)*COS(PI*YY)*COS(PI*TIME/T)
!        XX=XP(I)
!        YY=Y(J)
!         V(I,J,K)=2.*SIN(PI*XX)*COS(PI*XX)*SIN(PI*YY)**2*COS(PI*TIME/T)
!         W(I,J,K)=0.
!        ENDDO
!        ENDDO
!        ENDDO
!
!!$OMP PARALLEL DO private(I,J)
!        DO K=1,N3M
!        DO J=1,N2M
!        DO I=1,N1M
!          U(I,J,K)=2.*(SIN(PI*X(I))**2)*SIN(2*PI*YP(J))
!     &                                *SIN(2*PI*ZP(K))*COS(PI*TIME/3.)
!          V(I,J,K)=-SIN(2*PI*XP(I))*(SIN(PI*Y(J))**2)
!     &                                *SIN(2*PI*ZP(K))*COS(PI*TIME/3.)
!          W(I,J,K)=-SIN(2*PI*XP(I))*SIN(2.*PI*YP(J))
!     &                             *(SIN(PI*Z(K))**2)*COS(PI*TIME/3.)
!       ENDDO
!       ENDDO
!       ENDDO

      CALL real_time(ALVS_B(MSUB))
       IF (ITRACKING.NE.0) THEN
      WRITE(*,*)'>>>>>>>>TRACKING CACULATION START>>>>>>>>'  
          CALL TRANSPORT(MCLS,MGLOBAL,U,V,W,PSI_CN,VOL_TOT_ORI)
          CALL GRID_COUPLING(0,PSI_XN,PSI_YN,PSI_ZN,PSI_CN)
 !        write(*,*) '##################transport is removed!!'
      WRITE(*,*)'<<<<<<<<<TRACKING CACULATION END<<<<<<<<<'    
       ENDIF    
      CALL real_time(ALVS_E(MSUB))

 2000  CONTINUE
      DEALLOCATE(RK3XO,RK3YO,RK3ZO)
      DEALLOCATE(ANPSO)
       ELSE IF (IMPL .EQ. 1 ) THEN
!       MSUB=1
!       ALPHA_RK3=0.5    !BECAUSE 2*ALPHI=1
!      DTCONST=DT
!      DTCONSTI=1./DTCONST
!
!      WRITE(*,154) CFLM,CFL_SURF,CFLM_GRAV,CFL_SUM
!  154 FORMAT('CFL_CON= ',F7.4,'  CFL_SURF= ',F7.4,'  CFL_GRAV= ',F7.4
!     &                                              ,'  CFLSUM= ',F7.4)
!
!      IF (IPS.EQ.1) CALL RHSPSCD1(U,V,W,T)
!      IF (IPS.EQ.1 .AND. IBMON.EQ.2) CALL RHS_IBM_PS
!
!       CALL RHSNLHS_CN2(ITER_NS)
!!            write(*,*) 'no rhsnlhs'
!
!      CALL real_time(ALVS_B(1))
!       IF (ITRACKING.NE.0) THEN
!      WRITE(*,*)'>>>>>>>>TRACKING CACULATION START>>>>>>>>'  
!          CALL TRANSPORT(MCLS,MGLOBAL,U,V,W,PSI_CN,VOL_TOT_ORI)
!!        write(*,*) '##################transport is removed!!'
!      WRITE(*,*)'<<<<<<<<<TRACKING CACULATION END<<<<<<<<<'    
!       ENDIF 
!      CALL real_time(ALVS_E(1))

!      IF (IPS.EQ.1) CALL LHSPS(U,V,W)

      ENDIF

!       PHI_MAX=-10000.
!       PHI_MIN=10000.
!      DO K=1,N3M
!      DO J=1,N2M
!      DO I=1,N1M
!       PHI_MAX=MAX(PHI_MAX,P(I,J,K))
!       PHI_MIN=MIN(PHI_MIN,P(I,J,K))
!      ENDDO
!      ENDDO
!      ENDDO
!        WRITE(*,965) TIME,PHI_MAX,PHI_MIN,ABS(PHI_MAX-PHI_MIN)-36.5
!      CLOSE(923)
! 965   FORMAT(4es15.5)

c-----MAXIMUM_VEL
      IF (DENR .EQ. 1.) THEN
      DEN_DIFFI=0.
      ELSE
      DEN_DIFFI=1./DEN_DIFF
      ENDIF
       U_MAX=0.
       E_SUM=0.
       VOL=0.
       U_BUB_AVG=0.
       VOL_BUB=0.
!$OMP PARALLEL DO private(I,J,DVOL,DEN_CN,UU,VV,WW,VEL_KIM)
!$OMP&reduction(MAX:U_MAX)
!$OMP&reduction(+:E_SUM,VOL,U_BUB_AVG,VOL_BUB)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
        IF (FUNCBODY(XP(I),YP(J),ZP(K)) .GT. 0.) THEN
           DVOL=SDX(I)*SDY(J)*SDZ(K)
           DEN_CN=DENM+DEN_DIFF*PSI_CN(I,J,K)
           UU=0.5*(U(I,J,K)+U(IPV(I),J,K))
           VV=0.5*(V(I,J,K)+V(I,JPV(J),K))
           WW=0.5*(W(I,J,k)+W(I,J,KPV(K)))
           VEL_KIM=UU**2+VV**2+WW**2

           U_MAX=MAX(U_MAX,sqrt(VEL_KIM))
           E_SUM=E_SUM+0.5*DEN_CN*VEL_KIM*DVOL
           VOL=VOL+DVOL

           !BUBBLE_BULK_VELOCITY
           U_BUB_AVG=U_BUB_AVG+W(I,J,K)*DVOL*(1.-PSI_CN(I,J,K)) 
           VOL_BUB=VOL_BUB+DVOL*(1.-PSI_CN(I,J,K))
        ENDIF
      ENDDO
      ENDDO
      ENDDO
         E_AVG=E_SUM/VOL
         IF(VOL_BUB.NE.0.) U_BUB_AVG=U_BUB_AVG/VOL_BUB
         write(*,966) TIME,U_MAX,E_AVG
 966  format ('time',f10.5,' u_max',es10.3,' e_avg',es10.3)

      OPEN(923,FILE='0VEL,ENER.DAT',POSITION='APPEND')
       WRITE(923,924) TIME,U_MAX,U_BUB_AVG,E_AVG !TIME IS N TIME, NOT N+1 TIME. 
 924   FORMAT(4ES15.5)                           !SO I USE PSI_CN(I,J,K)
      CLOSE(923)
c-----MAXIMUM_VEL
C-----------------------------------------------------------------------

      TIME=TIME+DT

      CALL CONVRGE1(DVMAX,U,V,W)          ! CONVERGENCY CHECK for FLOW
      IF (DVMAX.GT.RESID_POI) THEN
         WRITE(*,300) NTIME,DVMAX,TIME
      ENDIF

!      CALL MASSCHECK(QMMAX)

      IF (MOD(NTIME,NPIN).EQ.0) THEN
         IHIST=IHIST+1
         CALL WRITEHISTORY(CFLM,DVMAX,QMMAX)
      ENDIF

      IF ((NTEST.NE.0) .AND. (MOD(M,NTPRINT).EQ.0)) THEN
        CALL TRACER(U,V,W,P,NTEST,NTPRINT,NTEST2,ITR,JTR,KTR
     &,ITR2,JTR2,KTR2)
      ENDIF
       
      IF (MOD(NTIME,NPRINT).EQ.0) THEN
         CALL WRITEFIELD(NV,IHIST,U,V,W,P,VOL_TOT_ORI,QVOL_ORI)  ! binary file
      ENDIF


      IF (IAVG .EQ. 1) CALL FIELD_AVG_RT(IHIST,NAV,U,V,W,P,PSI_CN)
      IF (IAVG .EQ. 2) CALL FIELD_AVG_XYZ(IHIST,NAV,U,V,W,P,PSI_CN)
      
      CALL real_time(TIME_END)

      FTRTIME1=(CONTINUITY_E(1)-CONTINUITY_B(1))
     &        +(CONTINUITY_E(2)-CONTINUITY_B(2))
     &        +(CONTINUITY_E(3)-CONTINUITY_B(3))
      FTRTIME2=0.
      FTRTIME3=0.
      DO I=1,ITER_NS
      FTRTIME2=FTRTIME2+(ANSETIME_E(I)-ANSETIME_B(I))
      FTRTIME3=FTRTIME3+(POISSTIME_E(I)-POISSTIME_B(I))
      ENDDO
      FTRTIME4=(ALVS_E(1)-ALVS_B(1))
     &        +(ALVS_E(2)-ALVS_B(2))
     &        +(ALVS_E(3)-ALVS_B(3))

      FTRTIME5=TIME_END-TIME_BEGIN

!      WRITE(*,201) FTRTIME1
!      WRITE(*,202) FTRTIME2
!      WRITE(*,203) FTRTIME3
!      WRITE(*,204) FTRTIME4
!      WRITE(*,205) FTRTIME5
! 201  FORMAT('TIME FOR CONTI   : ',F12.3,' SECONDS')
! 202  FORMAT('TIME FOR NSE     : ',F12.3,' SECONDS')
! 203  FORMAT('TIME FOR POISSON : ',F12.3,' SECONDS')
! 204  FORMAT('TIME OF LVS      : ',F12.3,' SECONDS')
! 205  FORMAT('TIME OF OPERATION: ',F12.3,' SECONDS')

      WRITE(*,206) FTRTIME1,FTRTIME2,FTRTIME3,FTRTIME4,FTRTIME5
 206   FORMAT('CONTI',F6.2,'  NSE',F6.2,'  POI',F6.2,'  LVS',F6.2,
     &'  TOT',F6.2)

      WRITE(87,211)TIME,FTRTIME1,FTRTIME2,FTRTIME3,FTRTIME4,FTRTIME5

 3000 CONTINUE

      CLOSE(100)
      CLOSE(13)
      IF (ILES  .EQ. 1) CLOSE(89)
      IF (NTEST .NE. 0) CLOSE(91)
      IF (NTEST .NE. 0) CLOSE(92)
      IF (NTEST .NE. 0) CLOSE(93)
      IF (NTEST .NE. 0) CLOSE(94)
      IF (ILES  .EQ. 1) CLOSE(95)
      IF (ICH   .NE. 0) CLOSE(96)
      IF (IPS   .EQ. 1) CLOSE(97)
      CLOSE(87)

      IF (MOD(NTIME,NPRINT).NE.0) THEN
         CALL WRITEFIELD(NV,IHIST,U,V,W,P,VOL_TOT_ORI,QVOL_ORI)
      ENDIF

      CALL real_time(TOTAL_TIME_E)
      WRITE(*,*) ' '
      WRITE(*,210) TOTAL_TIME_E-TOTAL_TIME_B
 210  FORMAT('TOTAL TIME      : ',F12.2,' SECONDS')
 211  FORMAT(F13.5,5F12.4)

      call timestamp

      STOP
      END

C===============================================================
C     SUBROUTINES FOR INPUT AND OUTPUT
C===============================================================

C*******************************************************************
C     SUBROUTINE MAKEFLD
C*******************************************************************
C     THIS SUBROUTINE IS ON lica_**.f

C*******************************************************************
      SUBROUTINE PREFLD(IHIST,fileprevel,IDTOLD,U,V,W,P,IPZERO
     &,VOL_TOT_ORI,QVOL_ORI)
C*******************************************************************
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR

      USE HEAT_VAR

      CHARACTER*30 fileprevel

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      REAL VOL_TOT_ORI(MLVS)
      
      REAL QVOL_ORI

!     dum for future use
      READ(12,1001) N1,N2,N3
      READ(12,1002) IHIST,M,TIME,DT_O
      READ(12,1003) IPSS,IIPX,IIPY,IIPZ,PRAA
      READ(12,1004) QVOL_ORI
      READ(12,1005) ((( U(I,J,K) ,I=1,N1),J=0,N2),K=0,N3)
      READ(12,1005) ((( V(I,J,K) ,I=0,N1),J=1,N2),K=0,N3)
      READ(12,1005) ((( W(I,J,K) ,I=0,N1),J=0,N2),K=1,N3)
      READ(12,1005) ((( P(I,J,K) ,I=1,N1M),J=1,N2M),K=1,N3M)

      !LVS
      READ(12,1011) N1F,N2F,N3F
      READ(12,1012) DENR,VISR,FR,SURF_J,FK,RE_AIR

      READ(12,1013) NLVS,N_MAX
      READ(12,1014) (VOL_TOT_ORI(LLVS), LLVS=1,NLVS)
      READ(12,1015) ((NUMA(N,LLVS) ,N=1,N_MAX),LLVS=1,NLVS)
      DO LLVS=1,NLVS
       DO N=1,N_MAX
        DO NN=1,NUMA(N,LLVS)
        READ(12,1016) I,J,K,ALPHI_TMP
         I_B(NN,N,LLVS)=I
         J_B(NN,N,LLVS)=J
         K_B(NN,N,LLVS)=K
         ALPHI(NN,N,LLVS)=ALPHI_TMP
        ENDDO
       ENDDO
      ENDDO
      
      IF (IPS.EQ.1) THEN
      IF (IPSS.EQ.0) THEN
       CALL MAKEFLD_TEMP
      ELSE IF (IPSS.EQ.1) THEN
       READ(12,1006) PRM,SCR,TCR
       READ(12,1005) ((( T(I,J,K) ,I=0,N1),J=0,N2),K=0,N3)
      ENDIF
      ENDIF

 1001   FORMAT(3I8)
 1002   FORMAT(2I8,2ES20.12)
 1003   FORMAT(4I8,1ES20.12)
 1004   FORMAT(1ES20.12)
 1005   FORMAT(5ES20.12)
 1006   FORMAT(3ES20.12)
 
 1011   FORMAT(3I8)
 1012   FORMAT(6ES20.12)
 1013   FORMAT(2I8)
 1014   FORMAT(8ES20.12)
 1015   FORMAT(8I8)
 1016   FORMAT(3I8,1ES20.12)


      IF (IPZERO .EQ. 1) P=0.

      WRITE(*,*)' '
      WRITE(*,100)
      WRITE(*,101)
      WRITE(*,102) fileprevel
      WRITE(*,104) N1,N2,N3
      WRITE(*,105) IHIST,M,TIME,DT_O
      WRITE(*,106) IPSS,PRAA

  100 FORMAT('----------- INITIAL FIELD INFORMATION -----------')
  101 FORMAT('INITIAL FIELD      : READING DONE')
  102 FORMAT('INITIAL FIELD NAME : ',A30)
  104 FORMAT('N1=',I12,'  N2=',I12,'  N3=',I12)
  105 FORMAT('IHIST=',I9,'  M=',I13,'  TIME=',F10.5,'  DT=',F12.8)
  106 FORMAT('IPS=',I11,'  PRA=',F11.3)

       IF (IDTOLD .EQ. 1) THEN
        DT=DT_O
        WRITE(*,*) 'DT_OLD IS USED'
       ELSE
        WRITE(*,*) 'DT_NEW IS USED'
       ENDIF

      DO LLVS=1,NLVS
       NUMA_MAX=0
       DO N=1,N_MAX
         NUMA_MAX=max(NUMA_MAX,NUMA(N,LLVS))
       ENDDO
        WRITE(*,*) 'NUMA_MAX/MF_BAND=',FLOAT(NUMA_MAX)/FLOAT(MF_BAND)
      ENDDO

      RETURN
      END

C*******************************************************************
      SUBROUTINE WRITEFIELD(NV,IHIST,U,V,W,P,VOL_TOT_ORI,QVOL_ORI)
C*******************************************************************
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR

      USE HEAT_VAR
      
      USE IBM_VAR
      
      IMPLICIT NONE

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      REAL VOL_TOT_ORI(MLVS)
      REAL QVOL_ORI
      
      INTEGER NV,IHIST
      
      INTEGER IDUM
      REAL DUM
      INTEGER I,J,K,N,NN,LLVS,L
      CHARACTER*9 tname
      CHARACTER*3 tfn1
      INTEGER IDG1,IDG2,IDG3,IDG4,IDG5,IDG6

      idum=0
      dum=0.

      tfn1='fld'
      idg1=ihist/100000
      idg2=(ihist-idg1*100000)/10000
      idg3=(ihist-idg1*100000-idg2*10000)/1000
      idg4=(ihist-idg1*100000-idg2*10000-idg3*1000)/100
      idg5=(ihist-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
      idg6=ihist-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10
      tname=tfn1//char(idg1+48)//char(idg2+48)//
     &      char(idg3+48)//char(idg4+48)//char(idg5+48)//char(idg6+48)

      OPEN(NV,FILE=tname)
ccc      IVER=-2
ccc      WRITE(NV) IVER
      WRITE(NV,1001) N1,N2,N3
      WRITE(NV,1002) IHIST,M,TIME,DT
      WRITE(NV,1003) IPS,IPX,IPY,IPZ,PRA
      WRITE(NV,1004) QVOL_ORI
      WRITE(NV,1005) ((( U(I,J,K) ,I=1,N1),J=0,N2),K=0,N3)
      WRITE(NV,1005) ((( V(I,J,K) ,I=0,N1),J=1,N2),K=0,N3)
      WRITE(NV,1005) ((( W(I,J,K) ,I=0,N1),J=0,N2),K=1,N3)
      WRITE(NV,1005) ((( P(I,J,K) ,I=1,N1M),J=1,N2M),K=1,N3M)

      !LVS
      WRITE(NV,1011) N1F,N2F,N3F
      WRITE(NV,1012) DENR,VISR,FR,SURF_J,FK,RE_AIR

      WRITE(NV,1013) NLVS,N_MAX
      WRITE(NV,1014) (VOL_TOT_ORI(LLVS), LLVS=1,NLVS)
      WRITE(NV,1015) ((NUMA(N,LLVS) ,N=1,N_MAX),LLVS=1,NLVS)
      DO LLVS=1,NLVS
       DO N=1,N_MAX
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
        WRITE(NV,1016) I,J,K,ALPHI(NN,N,LLVS)
        ENDDO
       ENDDO
      ENDDO

      IF (IPS.EQ.1) THEN
      WRITE(NV,1006) PRM,SCR,TCR
      WRITE(NV,1005) ((( T(I,J,K) ,I=0,N1),J=0,N2),K=0,N3)
      ENDIF

      CLOSE(NV)

 1001   FORMAT(3I8)
 1002   FORMAT(2I8,2ES20.12)
 1003   FORMAT(4I8,1ES20.12)
 1004   FORMAT(1ES20.12)
 1005   FORMAT(5ES20.12)
 1006   FORMAT(3ES20.12)
 
 1011   FORMAT(3I8)
 1012   FORMAT(6ES20.12)
 1013   FORMAT(2I8)
 1014   FORMAT(8ES20.12)
 1015   FORMAT(8I8)
 1016   FORMAT(3I8,1ES20.12)

      NV=NV+1

      RETURN
      END

C*******************************************************************
      SUBROUTINE WRITEHISTORY(CFLM,DVMAX,QMMAX)
C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE FLOW_VAR

      IMPLICIT NONE

      REAL CFLM,DVMAX,QMMAX

      WRITE(13,130) TIME,DT,CFLM,DVMAX,QMMAX
 130  FORMAT(F13.5,4ES15.7)

      IF (ICH.NE.0) THEN
      WRITE(96,140) TIME,PMI_DUDY,PMI
 140  FORMAT(F13.5,3ES20.12)
      ENDIF

      RETURN
      END
      
C*******************************************************************
      SUBROUTINE WRITEHISTORY_HEAT(W,PSI_C)
C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE FLOW_VAR

      USE TWO_PHASE_PROPERTY
      USE HEAT_VAR

      IMPLICIT NONE

      REAL W(0:M1,0:M2,0:M3)
      REAL PSI_C(0:M1,0:M2,0:M3)

      INTEGER I,J,K
      REAL AA,UU,FUNCBODY,DVOL,UU_TMP,DENSCF,AMEAN_TEMP,ANUSSET

!      PSAVG=0.
!      DO J=1,N2M
!      TMP=0.
!      DO K=1,N3M
!      DO I=1,N1M
!      TMP=TMP+T(I,J,K)
!      ENDDO
!      ENDDO
!      PSAVG=PSAVG+TMP*SDY(J)
!      ENDDO
!      PSAVG=PSAVG/2./FLOAT(N1M*N3M)
!
!      DPSDYW=0.
!      DO K=1,N3M
!      DO I=1,N1M
!      DPSDYW=DPSDYW+T(I,1,K)-T(I,0,K)
!      ENDDO
!      ENDDO
!      DPSDYW=DPSDYW/FLOAT(N1M*N3M)*VVDY(1)
!
!      WRITE(97,140)TIME,PSAVG,DPSDYW
      
      AA=0.
      UU=0.
!$OMP PARALLEL DO private(I,J,DVOL,UU_TMP,DENSCF)
!$OMP&reduction(+:AA,UU)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
       	IF (FUNCBODY(XP(I),YP(J),ZP(K)) .GT. 0.) THEN
       	 DVOL=SDX(I)*SDY(J)*SDZ(K)

         UU_TMP=0.5*(W(I,J,KPV(K))+W(I,J,K))
         DENSCF=DENSCM+DENSC_DIFF*PSI_C(I,J,K)

         AA=AA+DENSCF*UU_TMP*T(I,J,K)*DVOL
         UU=UU+DENSCF*UU_TMP*DVOL
        endif
      ENDDO
      ENDDO
      ENDDO

       IF (AA .NE. 0. .AND. UU.NE. 0.) THEN
        AMEAN_TEMP=AA/UU
        ANUSSET=2./AMEAN_TEMP
        WRITE(*,*) 'THETA_MEAN=',AMEAN_TEMP,'NUSSET=',ANUSSET
       ELSE
        AMEAN_TEMP=0.
        ANUSSET=0.
       ENDIF
       
      WRITE(97,140)TIME,ANUSSET,AMEAN_TEMP
 140  FORMAT(F13.5,3ES20.12)
      

      RETURN
      END

C*******************************************************************
      SUBROUTINE TRACEINIT(NTEST,NTPRINT,NTEST2,ITR,JTR,KTR
     &,ITR2,JTR2,KTR2)
C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE

      INTEGER NTEST,NTPRINT,NTEST2
      INTEGER ITR(10000),JTR(10000),KTR(10000)
      INTEGER ITR2(10000),JTR2(10000),KTR2(10000)
      
      INTEGER N

      DO 10 N=1,NTEST
      READ(10,*) ITR(N),JTR(N),KTR(N)
   10 CONTINUE
      READ(10,*) NTEST2
      IF (NTEST2.NE.0) THEN
      DO 11 N=1,NTEST2
      READ(10,*) ITR2(N),JTR2(N),KTR2(N)
   11 CONTINUE
      ENDIF

      PRINT*, '===== TRACE POSITION ====='
      DO 15 N=1,NTEST
      WRITE(*,16) N,ITR(N),JTR(N),KTR(N),
     &            XP(ITR(N)),YP(JTR(N)),ZP(KTR(N))
   15 CONTINUE
      IF (NTEST2.NE.0) THEN
      PRINT*, '===== TRACE POSITION 2 ====='
      DO 17 N=1,NTEST2
      WRITE(*,16) N,ITR2(N),JTR2(N),KTR2(N),
     &            XP(ITR2(N)),YP(JTR2(N)),ZP(KTR2(N))
   17 CONTINUE
      ENDIF
   16 FORMAT(I5,3I6,3F12.5)

      RETURN
      END

C*******************************************************************
      SUBROUTINE TRACER(U,V,W,P,NTEST,NTPRINT,NTEST2,ITR,JTR,KTR
     &,ITR2,JTR2,KTR2)
C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE FLOW_VAR
      
      IMPLICIT NONE
      
      INTEGER N

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      INTEGER NTEST,NTPRINT,NTEST2
      INTEGER ITR(10000),JTR(10000),KTR(10000)
      INTEGER ITR2(10000),JTR2(10000),KTR2(10000)

      WRITE(91,101)TIME,(0.5*(U(ITR(N),JTR(N),KTR(N))
     &                       +U(ITR(N)+1,JTR(N),KTR(N))),N=1,NTEST)
      WRITE(92,101)TIME,(0.5*(V(ITR(N),JTR(N),KTR(N))
     &                       +V(ITR(N),JTR(N)+1,KTR(N))),N=1,NTEST)
      WRITE(93,101)TIME,(0.5*(W(ITR(N),JTR(N),KTR(N))
     &                       +W(ITR(N),JTR(N),KPV(KTR(N)))),N=1,NTEST)
      WRITE(94,101)TIME,(P(ITR(N),JTR(N),KTR(N)),N=1,NTEST)
!      WRITE(91,101)TIME,(U(ITR(N),JTR(N),KTR(N)),N=1,NTEST)
!      WRITE(92,101)TIME,(V(ITR(N),JTR(N),KTR(N)),N=1,NTEST)
!      WRITE(93,101)TIME,(W(ITR(N),JTR(N),KTR(N)),N=1,NTEST)
!      WRITE(94,101)TIME,(P(ITR(N),JTR(N),KTR(N)),N=1,NTEST)
      IF (NTEST2.NE.0) THEN
      WRITE(81,101)TIME,(0.5*(U(ITR2(N),JTR2(N),KTR2(N))
     &                       +U(ITR2(N)+1,JTR2(N),KTR2(N))),N=1,NTEST2)
      WRITE(82,101)TIME,(0.5*(V(ITR2(N),JTR2(N),KTR2(N))
     &                       +V(ITR2(N),JTR2(N)+1,KTR2(N))),N=1,NTEST2)
      WRITE(83,101)TIME,(0.5*(W(ITR2(N),JTR2(N),KTR2(N))
     &                    +W(ITR2(N),JTR2(N),KPV(KTR2(N)))),N=1,NTEST2)
      WRITE(84,101)TIME,(P(ITR2(N),JTR2(N),KTR2(N)),N=1,NTEST2)
      ENDIF
  101 FORMAT(F15.7,10000ES15.7)

      RETURN
      END

C*******************************************************************
      SUBROUTINE RK3COEF
C*******************************************************************
      USE FLOW_VAR

      GAMMA(1)=8./15.
      GAMMA(2)=5./12.
      GAMMA(3)=3./4.
      RO(1)=0.
      RO(2)=-17./60.
      RO(3)=-5./12.

      RETURN
      END

C******************************************************************
      SUBROUTINE CFL(IMPL,CFLM,U,V,W)
C******************************************************************
C     CALCULATE THE MAXIMUM CFL NUMBER OF FLOW FIELD
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      INTEGER IMPL
      REAL CFLM
      
      INTEGER I,J,K
      REAL CFLI1,CFLI2,CFLI3,CFLI

      CFLM=0.

       IF ( N3M .EQ. 1 ) THEN
!$OMP PARALLEL DO private(I,CFLI1,CFLI2,CFLI)
!$OMP&reduction(MAX:CFLM)
      DO 12 J=1,N2M
      DO 12 I=1,N1M
        CFLI1=ABS(U(I,J,1)+U(IPV(I),J,1))*SSDX(I)
        CFLI2=ABS(V(I,J,1)+V(I,JPV(J),1))*SSDY(J)
        CFLI=0.5*(CFLI1+CFLI2)*DT
        CFLM=MAX(CFLI,CFLM)
 12    CONTINUE
          WRITE(*,*) '2D CFL IS CAL.'

        ELSE
!$OMP PARALLEL DO private(I,J,CFLI1,CFLI2,CFLI3,CFLI)
!$OMP&reduction(MAX:CFLM)
      DO 13 K=1,N3M
      DO 13 J=1,N2M
      DO 13 I=1,N1M
         CFLI1=ABS(U(I,J,K)+U(IPV(I),J,K))*SSDX(I)
         CFLI2=ABS(V(I,J,K)+V(I,JPV(J),K))*SSDY(J)
         CFLI3=ABS(W(I,J,K)+W(I,J,KPV(K)))*SSDZ(K)
         CFLI=0.5*(CFLI1+CFLI2+CFLI3)*DT
         CFLM=MAX(CFLI,CFLM)
!         IF (CFLI .GE. CFLM) THEN
!              write(*,*) i,j,k,cflm
!         ENDIF
   13 CONTINUE
       ENDIF

      RETURN
      END

C*******************************************************************
      SUBROUTINE MAKEFLD(IHIST,NTII,U,V,W,P)
C*******************************************************************
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      USE HEAT_VAR
      
      IMPLICIT NONE

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)
      
      INTEGER IHIST,NTII      
      
      INTEGER I,J,K,ICHANNEL
      REAL PI,RCOEF,RE_WATER,EPS_PTR,MFD_TURB,UMAX,RR,Y_PLUS,W1,W2
      REAL RAN1,RAN_NUM1,IDUM,EVMM,EVM,EVM_DIVIDE
      REAL FUNCBODY,PBODY
      REAL UU,VV,WW,U_BULK
      REAL VOL,DVOL
      REAL DTDX,COEF,AA,BB,CC

      PI=ACOS(-1.)

        U=0.
        V=0.
        W=0.

        IHIST=0
        NTII=0
        TIME=0.

C-----INITIAL FIELD OPTION--------C
      RCOEF=(2./YL)**2   !for the case that pipe_radius is not 0.5
      RE_WATER=RE_AIR*DENR/VISR
      EPS_PTR=1.
      MFD_TURB=0 !0:LAMINAR PROFILE, 1:TURBULENT PROFILE

!      ICHANNEL=0

      IF (MFD_TURB .EQ. 0) THEN
       UMAX=2.     !for the case that MAXIMUN VEL OF LAMINAR PROFILE is not 1.
      ELSE IF (MFD_TURB .EQ. 1) THEN
       UMAX=1.
      ENDIF
C-----INITIAL FIELD OPTION--------C


!     INITIALIZE BY LAMINAR PROFILE (ONLY TO DEFINE AT THE CORNER)
!-----ADD PERTURBATION
!     GENERATED FIELD ADDED BY PERTURBATION DOES NOT SATISFY PERIODICITY
!     HOWEVER, PERIODICITY WILL BE SATISFIED IN 1 TIME STEP MARCHING

!     U: LAMINAR PROFILE + PERTURBATION
c!$OMP PARALLEL DO private(I,J,PBODY,RAN_NUM1,RR,UU,Y_PLUS,W1,W2) !no_omp for random_number_generation
      DO K=1,N3M
      DO J=1,N2M
      DO I=IBG,N1M
        PBODY=FUNCBODY(X(I),YP(J),ZP(K))
       IF ( PBODY .GT. 0. ) THEN
      RAN_NUM1=2.*RAN1(IDUM)-1.
       IF (MFD_TURB .EQ. 0) THEN
       	IF (ICHANNEL .EQ. 1) THEN
        RR=YP(J)**2
       	ELSE
        RR=(XP(I)**2+YP(J)**2)
        ENDIF
         UU=1.-RCOEF*RR
       ELSE IF (MFD_TURB .EQ. 1) THEN
         Y_PLUS=PBODY*RE_WATER
        IF (Y_PLUS .LE. 1.) THEN
           UU=Y_PLUS            !LAW OF WALL
        ELSE IF (Y_PLUS .GE. 80.) THEN
           UU=2.5*LOG(Y_PLUS)+5.5      !LOG LAW(KIM ET AL. 1987)
        ELSE !BUFFER LAYER
          W1=Y_PLUS                 
          W2=2.5*LOG(Y_PLUS)+5.5
         IF (W1 .LT. W2) THEN
          UU=W1 
         ELSE
          UU=W2 
         ENDIF
        ENDIF
       ENDIF
      U(I,J,K)=EPS_PTR*RAN_NUM1*UU
       ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP PARALLEL DO private(I,J,EVM,EVM_DIVIDE,EVMM)
      DO K=1,N3M
      EVM=0.
      EVM_DIVIDE=0.
      DO J=1,N2M
      DO I=IBG,N1M
       IF ( FUNCBODY(X(I),YP(J),ZP(K)) .GT. 0. ) THEN
      EVM=EVM+U(I,J,K)*SDX(I)*SDY(J)
      EVM_DIVIDE=EVM_DIVIDE+SDX(I)*SDY(J)
       ENDIF
      ENDDO
      ENDDO

      IF (EVM .EQ. 0.) THEN
       EVMM=0.
      ELSE
       EVMM=EVM/EVM_DIVIDE !/YL/ZL
      ENDIF

      DO J=1,N2M
      DO I=IBG,N1M
       IF ( FUNCBODY(X(I),YP(J),ZP(K)) .GT. 0. ) THEN
      U(I,J,K)=UMAX*( U(I,J,K)-EVMM )
       ENDIF
      ENDDO
      ENDDO
      ENDDO

!     V: LAMINAR PROFILE + PERTURBATION
c!$OMP PARALLEL DO private(I,J,PBODY,RAN_NUM1,RR,UU,Y_PLUS,W1,W2) !no_omp for random_number_generation
      DO K=1,N3M
      DO J=JBG,N2M
      DO I=1,N1M
        PBODY=FUNCBODY(XP(I),Y(J),ZP(K))
       IF ( PBODY .GT. 0. ) THEN
      RAN_NUM1=2.*RAN1(IDUM)-1.
       IF (MFD_TURB .EQ. 0) THEN
       	IF (ICHANNEL .EQ. 1) THEN
        RR=YP(J)**2
       	ELSE
        RR=(XP(I)**2+YP(J)**2)
        ENDIF
         UU=1.-RCOEF*RR
       ELSE IF (MFD_TURB .EQ. 1) THEN
         Y_PLUS=PBODY*RE_WATER
        IF (Y_PLUS .LE. 1.) THEN
           UU=Y_PLUS            !LAW OF WALL
        ELSE IF (Y_PLUS .GE. 80.) THEN
           UU=2.5*LOG(Y_PLUS)+5.5      !LOG LAW(KIM ET AL. 1987)
        ELSE !BUFFER LAYER
          W1=Y_PLUS                 
          W2=2.5*LOG(Y_PLUS)+5.5
         IF (W1 .LT. W2) THEN
          UU=W1 
         ELSE
          UU=W2 
         ENDIF
        ENDIF
       ENDIF
      V(I,J,K)=EPS_PTR*RAN_NUM1*UU
       ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP PARALLEL DO private(I,J,EVM,EVM_DIVIDE,EVMM)
      DO K=1,N3M
      EVM=0.
      EVM_DIVIDE=0.
      DO J=JBG,N2M
      DO I=1,N1M
       IF ( FUNCBODY(XP(I),Y(J),ZP(K)) .GT. 0. ) THEN
      EVM=EVM+V(I,J,K)*SDX(I)*SDY(J)
      EVM_DIVIDE=EVM_DIVIDE+SDX(I)*SDY(J)
       ENDIF
      ENDDO
      ENDDO
      
      IF (EVM .EQ. 0.) THEN
       EVMM=0.
      ELSE
       EVMM=EVM/EVM_DIVIDE !/YL/ZL
      ENDIF

      DO J=JBG,N2M
      DO I=1,N1M
       IF ( FUNCBODY(XP(I),Y(J),ZP(K)) .GT. 0. ) THEN
      V(I,J,K)=UMAX*( V(I,J,K)-EVMM )
       ENDIF
      ENDDO
      ENDDO
      ENDDO

!     W: LAMINAR PROFILE + PERTURBATION
c!$OMP PARALLEL DO private(I,J,PBODY,RAN_NUM1,RR,UU,Y_PLUS,W1,W2) !no_omp for random_number_generation
      DO K=KBG,N3M
      DO J=1,N2M
      DO I=1,N1M
        PBODY=FUNCBODY(XP(I),YP(J),Z(K))
       IF ( PBODY .GT. 0. ) THEN
      RAN_NUM1=2.*RAN1(IDUM)-1.
       IF (MFD_TURB .EQ. 0) THEN
       	IF (ICHANNEL .EQ. 1) THEN
        RR=YP(J)**2
       	ELSE
        RR=(XP(I)**2+YP(J)**2)
        ENDIF
         UU=1.-RCOEF*RR
       ELSE IF (MFD_TURB .EQ. 1) THEN
         Y_PLUS=PBODY*RE_WATER
        IF (Y_PLUS .LE. 1.) THEN
           UU=Y_PLUS            !LAW OF WALL
        ELSE IF (Y_PLUS .GE. 80.) THEN
           UU=2.5*LOG(Y_PLUS)+5.5      !LOG LAW(KIM ET AL. 1987)
        ELSE !BUFFER LAYER
          W1=Y_PLUS                 
          W2=2.5*LOG(Y_PLUS)+5.5
         IF (W1 .LT. W2) THEN
          UU=W1 
         ELSE
          UU=W2 
         ENDIF
        ENDIF
       ENDIF
       W(I,J,K)=EPS_PTR*RAN_NUM1*UU
       ENDIF
      ENDDO
      ENDDO
      ENDDO
!$OMP PARALLEL DO private(I,J,EVM,EVM_DIVIDE,EVMM)
!$OMP&private(PBODY,RR,UU,Y_PLUS,W1,W2)
      DO K=KBG,N3M
      EVM=0.
      EVM_DIVIDE=0.
      DO J=1,N2M
      DO I=1,N1M
       IF ( FUNCBODY(XP(I),YP(J),Z(K)) .GT. 0. ) THEN
      EVM=EVM+W(I,J,K)*SDX(I)*SDY(J)
      EVM_DIVIDE=EVM_DIVIDE+SDX(I)*SDY(J)
       ENDIF
      ENDDO
      ENDDO

      IF (EVM .EQ. 0.) THEN
       EVMM=0.
      ELSE
       EVMM=EVM/EVM_DIVIDE !/YL/ZL
      ENDIF

      DO J=1,N2M
      DO I=1,N1M
        PBODY=FUNCBODY(XP(I),YP(J),Z(K))
       IF ( PBODY .GT. 0. ) THEN
       IF (MFD_TURB .EQ. 0) THEN
       	IF (ICHANNEL .EQ. 1) THEN
        RR=YP(J)**2
       	ELSE
        RR=(XP(I)**2+YP(J)**2)
        ENDIF
         UU=1.-RCOEF*RR
       ELSE IF (MFD_TURB .EQ. 1) THEN
         Y_PLUS=PBODY*RE_WATER
        IF (Y_PLUS .LE. 1.) THEN
           UU=Y_PLUS            !LAW OF WALL
        ELSE IF (Y_PLUS .GE. 80.) THEN
           UU=2.5*LOG(Y_PLUS)+5.5      !LOG LAW(KIM ET AL. 1987)
        ELSE !BUFFER LAYER
          W1=Y_PLUS                 
          W2=2.5*LOG(Y_PLUS)+5.5
         IF (W1 .LT. W2) THEN
          UU=W1 
         ELSE
          UU=W2 
         ENDIF
        ENDIF
       ENDIF
        W(I,J,K)=UMAX*( UU+W(I,J,K)-EVMM )
       ENDIF
      ENDDO
      ENDDO
      ENDDO

c-----velocity scaling to make bulk_velocity=1
       WW=0.
       VOL=0.
!$OMP PARALLEL DO private(I,J,DVOL)
!$OMP&reduction(+:WW,VOL)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
       IF ( FUNCBODY(XP(I),YP(J),Z(K)) .GT. 0. ) THEN
         DVOL=SDX(I)*SDY(J)*VDZ(K)
         WW=WW+W(I,J,K)*DVOL
         VOL=VOL+DVOL
       ENDIF
      ENDDO
      ENDDO
      ENDDO
       WW=WW/VOL
       IF (WW .EQ. 0.) THEN
       	WRITE(*,*) 'ZERO INITIAL STREAMWISE VELOCITY'
       	GOTO 100
       ENDIF
       WRITE(*,*) 'BULK_VELOCITY(before_scaling)=',WW

!$OMP PARALLEL DO private(I,J)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
       W(I,J,K)=W(I,J,K)*(1./WW)
      ENDDO
      ENDDO
      ENDDO

       WW=0.
       VOL=0.
!$OMP PARALLEL DO private(I,J,DVOL)
!$OMP&reduction(+:WW,VOL)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
       IF ( FUNCBODY(XP(I),YP(J),Z(K)) .GT. 0. ) THEN
         DVOL=SDX(I)*SDY(J)*VDZ(K)
         WW=WW+W(I,J,K)*DVOL
         VOL=VOL+DVOL
       ENDIF
      ENDDO
      ENDDO
      ENDDO
       WW=WW/VOL
       U_BULK=WW
       WRITE(*,*) 'BULK_VELOCITY(after_scaling)=',WW

 100   CONTINUE
c-----velocity scaling to make bulk_velocity=1

!BOUNDARY_CONDITION
        IF (IPX .EQ. 0) THEN
!$OMP PARALLEL DO private(J)
          DO K=1,N3M
          DO J=1,N2M
          U(1,J,K)=U(2,J,K)
          V(0,J,K)=V(1,J,K)
          W(0,J,K)=W(1,J,K)
          U(N1,J,K)=U(N1M,J,K)
          V(N1,J,K)=V(N1M,J,K)
          W(N1,J,K)=W(N1M,J,K)
         ENDDO
         ENDDO
        ENDIF
        IF (IPY .EQ. 0) THEN
!$OMP PARALLEL DO private(I)
          DO K=1,N3M
          DO I=1,N1M
          U(I,0,K)=U(I,1,K)
          V(I,1,K)=V(I,2,K)
          W(I,0,K)=W(I,1,K)
          U(I,N2,K)=U(I,N2M,K)
          V(I,N2,K)=V(I,N2M,K)
          W(I,N2,K)=W(I,N2M,K) 
         ENDDO
         ENDDO
        ENDIF
        IF (IPZ .EQ. 0) THEN
!$OMP PARALLEL DO private(I)
          DO J=1,N2M
          DO I=1,N1M
          U(I,J,0)=U(I,J,1)
          V(I,J,0)=V(I,J,1)
          W(I,J,1)=W(I,J,2)
          U(I,J,N3)=U(I,J,N3M)
          V(I,J,N3)=V(I,J,N3M)
          W(I,J,N3)=W(I,J,N3M) 
         ENDDO
         ENDDO
        ENDIF

!$OMP PARALLEL DO private(I,J)
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M
         P(I,J,K)=0.  
        ENDDO
        ENDDO
        ENDDO
        
         IF (IPS.EQ.1) CALL MAKEFLD_TEMP

        WRITE(*,*) 'INITIAL FIELD (U,V,W, AND P) IS MADE'
        WRITE(*,*)' '

C=====SAVE
       K=1
      !POSITION IS A LITTLE DIFFERENT FOR CONVINIENCE
       OPEN(146,FILE='0INITIAL_VEL_2D.DAT')
      IF (IPS .EQ. 0) THEN
       WRITE(146,*) 'VARIABLES="X","Y","U","V","W"'
      WRITE(146,*) 'ZONE I=',N1M,',J=',N2M,',F=POINT'
      DO J=1,N2M
      DO I=1,N1M
         UU=0.5*(U(I,J,K)+U(IPV(I),J,K))
         VV=0.5*(V(I,J,K)+V(I,JPV(J),K))
         WW=0.5*(W(I,J,K)+W(I,J,KPV(K)))
          WRITE(146,147) XP(I),YP(J),UU,VV,WW
      ENDDO
      ENDDO
      ELSE IF (IPS .EQ. 1) THEN
       WRITE(146,*) 'VARIABLES="X","Y","U","V","W","T"'
      WRITE(146,*) 'ZONE I=',N1M,',J=',N2M,',F=POINT'
      DO J=1,N2M
      DO I=1,N1M
         UU=0.5*(U(I,J,K)+U(IPV(I),J,K))
         VV=0.5*(V(I,J,K)+V(I,JPV(J),K))
         WW=0.5*(W(I,J,K)+W(I,J,KPV(K)))
          WRITE(146,147) XP(I),YP(J),UU,VV,WW,T(I,J,K)
      ENDDO
      ENDDO
      ENDIF !IF (IPS .EQ. 0) THEN
       CLOSE(146)
 147  FORMAT(6F15.8)
 
 
!       OPEN(146,FILE='0INITIAL_VEL.DAT')
!       WRITE(146,*) 'VARIABLES="X","Y","Z","U","V","W"'
!      WRITE(146,*) 'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'
!      DO K=1,N3M
!      DO J=1,N2M
!      DO I=1,N1M
!         UU=0.5*(U(I,J,K)+U(IPV(I),J,K))
!         VV=0.5*(V(I,J,K)+V(I,JPV(J),K))
!         WW=0.5*(W(I,J,K)+W(I,J,KPV(K)))
!        WRITE(146,147) XP(I),YP(J),ZP(K),UU,VV,WW
!      ENDDO
!      ENDDO
!      ENDDO
!       CLOSE(146)
! 147  FORMAT(6F15.8)

        RETURN
        END
        
C*******************************************************************
      SUBROUTINE MAKEFLD_TEMP
C*******************************************************************
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE HEAT_VAR

      IMPLICIT NONE

      INTEGER I,J,K
      REAL FUNCBODY,RR
      REAL DTDX,U_BULK,COEF,AA,BB,CC

       T=0.

       DTDX=2.
       U_BULK=1.
       COEF=2.*U_BULK*(0.5*YL)**2/(TCM/DENSCM)*DTDX !NOMALIZED COEFFICIENT
       AA=3./16.
       BB=1./16.
       CC=1./4.

       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
       	IF (FUNCBODY(XP(I),YP(J),ZP(K)) .GT. 0.) THEN
       	 RR=SQRT(XP(I)**2+YP(J)**2)
         T(I,J,K)=COEF*(AA+BB*RR**4-CC*RR**2)	
        ENDIF
       ENDDO
       ENDDO
       ENDDO

!      DO K=1,N3M
!      DO I=1,N1M
!      T(I,N2,K)=1.
!      T(I,0 ,K)=-1.
!      ENDDO
!      ENDDO
!
!      EVM=0.
!      DO K=1,N3M
!      DO J=1,N2M
!      DO I=1,N1M
!      RAN_NUM1=2.*RAN1(IDUM)-1.
!      T(I,J,K)=EPS_PTR*RAN_NUM1*(1.-YP(J)**2.)
!      EVM=EVM+T(I,J,K)
!      ENDDO
!      ENDDO
!      ENDDO
!      EVM=EVM/FLOAT(N1M*N2M*N3M)
!      DO K=1,N3M
!      DO J=1,N2M
!      DO I=1,N1M
!      T(I,J,K)=YP(J)+T(I,J,K)-EVM
!      ENDDO
!      ENDDO
!      ENDDO
!
!!     Z PERIODICITY
!      IF (IPZ .EQ. 1) THEN
!      DO J=1,N2M
!      DO I=1,N1M
!         T(I,J,0) =T(I,J,N3M)
!         T(I,J,N3)=T(I,J,1)
!      ENDDO
!      ENDDO
!      ENDIF
!
!!     X PERIODICITY
!      IF (IPX .EQ. 1) THEN
!      DO K=1,N3M
!      DO J=1,N2M
!         T(0 ,J,K)=T(N1M,J,K)
!         T(N1,J,K)=T(1  ,J,K)
!      ENDDO
!      ENDDO
!      ENDIF

        RETURN
        END
        
        
C------------------------------------------------
C RANDOM NUMBER GENERATOR

      FUNCTION RAN1(IDUM)
      INTEGER IDUM,IA,IM,IQ,IR,NTAB,NDIV
      REAL RAN1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2E-7,RNMX=1.-EPS)
      INTEGER J,K,IV(NTAB),IY
      SAVE IV,IY
      DATA IV /NTAB*0/, IY /0/
      IF (IDUM.LE.0.OR.IY.EQ.0) THEN
        IDUM=MAX(-IDUM,1)
        DO 11 J=NTAB+8,1,-1
          K=IDUM/IQ
          IDUM=IA*(IDUM-K*IQ)-IR*K
          IF (IDUM.LT.0) IDUM=IDUM+IM
          IF (J.LE.NTAB) IV(J)=IDUM
11      CONTINUE
        IY=IV(1)
      ENDIF
      K=IDUM/IQ
      IDUM=IA*(IDUM-K*IQ)-IR*K
      IF (IDUM.LT.0) IDUM=IDUM+IM

      J=1+IY/NDIV
      IY=IV(J)
      IV(J)=IDUM
      RAN1=MIN(AM*IY,RNMX)

      RETURN
      END


C******************************************************************
C     SECTION FOR SOLVING UNCOUPLED MOMENTUM EQS. OF NAVIER-STOKES EQS.
C     - RHSNLHS-EI.F
C     - USES 2ND. ORDER CRANK-NICHOLSON METHOD FOR VISCOUS TERMS
C     - USES 3RD. ORDER RUNGE-KUTTA METHOD FOR CONVECTIVE TERMS
C     - ADI SCHEME IS USED TO SOLVE ELLIPTIC MOMENTUM EQS. IN Z->Y->X
C       ORDER
C     - CAN HANDLE OPEN( 2 TYPES ) OR CLOSED TOP ACCORDING TO ITOPEN.
C     - FOR THESE SUBROUTINES, 'COMMON/WRK2~4/~ ' VARIABLES WORK
C       AS TEMPORARY STORAGE. THEY STORE NUt AT EXCEPT CELL CENTER AND
C       SHOULD NOT BE MODIFIED FOR WHOLE SGS & RHSNLHS PROCESS
C     - TO REDUCE REQUIRED MEMORY FOR DNS, ERASE COMMON STATEMENTS
C       HAVING TNU,TNUX,TNUY,TNUZ WHICH STORE SGS VISCOSITY FOR LES
C       IN WHOLE RHSNLHS ROUTINES
C     - MUST BE ACCOMPANIED BY LHS-EI.F OR LHS-EI-L.F(IN CASE OF LES)
C
C                                  LAST REVISED BY JUNGWOO KIM
C                                            1999.08.04.   11:00
c     -code for Interpolation scehme(bilinear interploation)
C******************************************************************

!***********************************************************************
!       FULLY IMPLICIT & FULLY EXPLICIT DIRECT NUMERICAL SIMULATION CODE
!        -LARGE EDDY SIMULATION
!        -IMMERSED BOUNDARY METHOD
!        -LEVEL SET METHOD    
!
!       ORININAL INCOMPRESSIBLE NAVIER-STOKES EQUATION 
!                                IS CODED BY JUNIL LEE(SEPTEMBER 2007)
!
!                                                  2013.10. KIYOUNH KIM
!***********************************************************************
C*******************************************************************
      SUBROUTINE RHSNLHS(ITRACKING,IPOISS,U,V,W,P,PSI_XN,PSI_YN,PSI_ZN
     &,PSI_CN,RK3XO,RK3YO,RK3ZO,ANPSO,QVOL_ORI,VOL1)
C*******************************************************************
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE IBM_VAR
      
      USE HEAT_VAR
      
      USE TIME_VAR
      
      IMPLICIT NONE

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      REAL PSI_XN(0:M1,0:M2,0:M3)
      REAL PSI_YN(0:M1,0:M2,0:M3)
      REAL PSI_ZN(0:M1,0:M2,0:M3)
      REAL PSI_CN(0:M1,0:M2,0:M3)

      REAL RK3XO(M1M,M2M,M3M),RK3YO(M1M,M2M,M3M),RK3ZO(M1M,M2M,M3M)
      REAL ANPSO(M1M,M2M,M3M)

      REAL DENF_X(0:M1,0:M2,0:M3)
      REAL DENF_Y(0:M1,0:M2,0:M3)
      REAL DENF_Z(0:M1,0:M2,0:M3)
      REAL DENF_C(0:M1,0:M2,0:M3)
      
      REAL QVOL_ORI
      REAL VOL1
      
      REAL DF_X(0:M1,0:M2,0:M3,3),VF_X(0:M1,0:M2,0:M3,3)
      REAL DF_Y(0:M1,0:M2,0:M3,3),VF_Y(0:M1,0:M2,0:M3,3)
      REAL DF_Z(0:M1,0:M2,0:M3,3),VF_Z(0:M1,0:M2,0:M3,3)
      REAL DF_C(0:M1,0:M2,0:M3,3)

      REAL RHS1(M1M,M2M,M3M,3)
      REAL DENF_XI(0:M1,0:M2,0:M3),DENF_YI(0:M1,0:M2,0:M3)
     &                                         ,DENF_ZI(0:M1,0:M2,0:M3) 
     
      REAL PSI_X(0:M1,0:M2,0:M3),PSI_Y(0:M1,0:M2,0:M3)
     &                     ,PSI_Z(0:M1,0:M2,0:M3),PSI_C(0:M1,0:M2,0:M3)

      INTEGER IPOISS
      INTEGER ITRACKING

      INTEGER I,J,K
      REAL QM1

      REAL WRITE_ON
      REAL RHS11,RHS22,RHS33,UU,VV,WW
      REAL DENSC,DENSCF
      

      write_on=0

!      IF(ICH.EQ.1)CALL QVOLCALC(QVOL1,1,W,PSI_ZN,DENF_Z,VOL1)
      CALL REAL_TIME(CONTINUITY_B(MSUB))
      CALL CAL_CONTINUITY(ITRACKING,DF_X,DF_Y,DF_Z,DF_C,VF_X,VF_Y,VF_Z
     &,U,V,W,PSI_XN,PSI_YN,PSI_ZN,PSI_CN,PSI_X,PSI_Y,PSI_Z,PSI_C)

!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M !not start at IBG. WE need I=1 for POISSON_SOLVER even though we don`t use periodic condition.
        DENF_X(I,J,K)=DENM+DEN_DIFF*PSI_X(I,J,K)
        DENF_XI(I,J,K)=1./DENF_X(I,J,K)
       ENDDO
       ENDDO
       ENDDO
       IF (IPX .NE. 1) THEN
!$OMP PARALLEL DO private(J)
       DO K=1,N3M
       DO J=1,N2M
        DENF_X(N1,J,K)=DENF_X(N1M,J,K)
        DENF_XI(N1,J,K)=DENF_XI(N1M,J,K) !FOR VARIABLE POISSON EQN.
       ENDDO
       ENDDO
       ENDIF

!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
        DENF_Y(I,J,K)=DENM+DEN_DIFF*PSI_Y(I,J,K)
        DENF_YI(I,J,K)=1./DENF_Y(I,J,K)
       ENDDO
       ENDDO
       ENDDO
       IF (IPY .NE. 1) THEN
!$OMP PARALLEL DO private(I)
       DO K=1,N3M
       DO I=1,N1M
        DENF_Y(I,N2,K)=DENF_Y(I,N2M,K)
        DENF_YI(I,N2,K)=DENF_YI(I,N2M,K) !FOR VARIABLE POISSON EQN.
       ENDDO
       ENDDO
       ENDIF

       IF (N3M .NE. 1) THEN
!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
         DENF_Z(I,J,K)=DENM+DEN_DIFF*PSI_Z(I,J,K)
         DENF_ZI(I,J,K)=1./DENF_Z(I,J,K)
        ENDDO
        ENDDO
        ENDDO
       IF (IPZ .NE. 1) THEN
!$OMP PARALLEL DO private(I)
       DO J=1,N2M
       DO I=1,N1M
        DENF_Z(I,J,N3)=DENF_Z(I,J,N3M)
        DENF_ZI(I,J,N3)=DENF_ZI(I,J,N3M) !FOR VARIABLE POISSON EQN.
       ENDDO
       ENDDO
       ENDIF
       ENDIF

!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
        DENF_C(I,J,K)=DENM+DEN_DIFF*PSI_C(I,J,K)
       ENDDO
       ENDDO
       ENDDO
      CALL REAL_TIME(CONTINUITY_E(MSUB))

       !ENERGY_EQUATION
        IF (IPS .EQ. 1) THEN

      	 CALL CAL_RHSPS(U,V,W,ANPSO,PSI_CN,DF_C)
         IF (IBMON.EQ.1) CALL RHS_IBM_PS(PSI_CN,PSI_C)
!$OMP PARALLEL DO private(I,J,DENSC,DENSCF)
         DO K=1,N3M
         DO J=1,N2M
         DO I=1,N1M
         	DENSC=DENSCM+DENSC_DIFF*PSI_CN(I,J,K)
         	DENSCF=DENSCM+DENSC_DIFF*PSI_C(I,J,K)

          RHSPS(I,J,K)=( (DENSC-DENSCF)*T(I,J,K)+RHSPS(I,J,K)	)/DENSCF
         ENDDO
         ENDDO
         ENDDO

         CALL LHSPS(U,V,W,PSI_CN,PSI_C,DF_C)

       ENDIF
        
C-----MOMENTUM EQUATION SOLVER
      !CAL RHS
      CALL REAL_TIME(ANSETIME_B(MSUB))
       CALL RHSX(ITRACKING,RHS1,U,V,W,P,DF_X,VF_X,PSI_XN,RK3XO)
       CALL RHSY(ITRACKING,RHS1,U,V,W,P,DF_Y,VF_Y,PSI_YN,RK3YO)
      IF (N3M .NE. 1) CALL RHSZ(ITRACKING,RHS1,U,V,W,P,DF_Z,VF_Z,PSI_ZN,
     &DENF_Z,RK3ZO)

       IF (IBMON.NE.0) THEN
          CALL RHS_IBMX(RHS1,DENF_X,U,V,W)
          CALL RHS_IBMY(RHS1,DENF_Y,U,V,W)
          IF (N3M .NE. 1)CALL RHS_IBMZ(RHS1,DENF_Z,U,V,W)
       ENDIF
       IF (ITRACKING .NE. 0) DEALLOCATE(SUR_X,SUR_Y,SUR_Z)

      !FOR DELTA-FORMULATION. 
!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M
       DO J=1,N2M
       DO I=IBG,N1M
          RHS1(I,J,K,1)=RHS1(I,J,K,1)*DENF_XI(I,J,K)-U(I,J,K)
         ENDDO
         ENDDO
         ENDDO

!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M
       DO J=JBG,N2M
       DO I=1,N1M
          RHS1(I,J,K,2)=RHS1(I,J,K,2)*DENF_YI(I,J,K)-V(I,J,K)
         ENDDO
         ENDDO
         ENDDO

      IF (N3M .NE. 1) THEN
!$OMP PARALLEL DO private(I,J)
       DO K=KBG,N3M
       DO J=1,N2M
       DO I=1,N1M
          RHS1(I,J,K,3)=RHS1(I,J,K,3)*DENF_ZI(I,J,K)-W(I,J,K)
         ENDDO
         ENDDO
         ENDDO
      ENDIF
 
      !CAL LHS
       CALL LHSU(RHS1,U(0,0,0),V(0,0,0),W(0,0,0),DF_X,VF_X,DENF_XI)
       CALL LHSV(RHS1,U(0,0,0),V(0,0,0),W(0,0,0),DF_Y,VF_Y,DENF_YI)
       IF (N3M .NE. 1) CALL LHSW(RHS1,U(0,0,0),V(0,0,0),W(0,0,0)
     &,DF_Z,VF_Z,DENF_ZI)
 
!        DO K=1,N3M
!        DO J=1,N2M
!        DO I=IBG,N1M
!         U(I,J,K)=RHS1(I,J,K,1)*DENF_XI(I,J,K)
!        ENDDO
!        ENDDO
!        ENDDO
!        DO K=1,N3M
!        DO J=JBG,N2M
!        DO I=1,N1M
!         V(I,J,K)=RHS1(I,J,K,2)*DENF_YI(I,J,K)
!        ENDDO
!        ENDDO
!        ENDDO
!        DO K=KBG,N3M
!        DO J=1,N2M
!        DO I=1,N1M
!         W(I,J,K)=RHS1(I,J,K,3)*DENF_ZI(I,J,K)
!        ENDDO
!        ENDDO
!        ENDDO
!        write(*,*) 'explicit Euler is used for NS equation'
 
      !!BOUNDARY UPDATE BEFORE FSM_CHOI FOR 2ND-ORDER ACCURACY.
        CALL BC(RHS1,U,V,W)
        CALL FSM_CHOI(DENF_XI,DENF_YI,DENF_ZI,U,V,W,P,PSI_ZN,DENF_Z
     &,QVOL_ORI,VOL1)

      CALL REAL_TIME(ANSETIME_E(MSUB))

      CALL real_time(POISSTIME_B(MSUB))

       CALL PRESSURE_CORRECTION(IPOISS,U(0,0,0),V(0,0,0),W(0,0,0)
     &,P(0,0,0),DENF_Y,DENF_C,DENF_XI,DENF_YI,DENF_ZI)
       CALL QVOLCALC(QM1,3,W,PSI_ZN,DENF_Z,VOL1) !FOR WRITE THE MASS FLOW RATE

       IF(ICH .NE. 0 ) CALL MEANPGR(W,PSI_ZN,PSI_CN,DENF_Z,VOL1)
      CALL real_time(POISSTIME_E(MSUB))
      
      IF (IPS .EQ. 1) THEN
      	IF (MSUB .EQ. 3) CALL WRITEHISTORY_HEAT(W,PSI_C)
      ENDIF

      if ( write_on .eq. 1 ) then
!!C=====SAVE
      IF ( MSUB .EQ. 1 ) THEN
       OPEN(144,FILE='0RHS1.DAT')
       WRITE(144,*) 'VARIABLES="X","Y","Z","u","v","w"'
      WRITE(144,*) 'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'

      OPEN(145,FILE='0u_hat1.DAT')
       WRITE(145,*) 'VARIABLES="X","Y","Z","u","v","w"'
      WRITE(145,*) 'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'  

       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
         RHS11=0.5*(RHS1(I,J,K,1)+RHS1(IPV(I),J,K,1))
         RHS22=0.5*(RHS1(I,J,K,2)+RHS1(I,JPV(J),K,2))
         RHS33=0.5*(RHS1(I,J,K,3)+RHS1(I,J,KPV(K),3))
!         RHS11=RHS1(I,J,K,1)
!         RHS22=RHS1(I,J,K,2)
!         RHS33=RHS1(I,J,K,3)  
       WRITE(144,148) XP(I),YP(J),ZP(K),RHS11,RHS22,RHS33
         UU=0.5*(U(I,J,K)+U(IPV(I),J,K))
         VV=0.5*(V(I,J,K)+V(I,JPV(J),K))
         WW=0.5*(W(I,J,K)+W(I,J,KPV(K)))
!         UU=U(I,J,K)
!         VV=V(I,J,K)
!         WW=W(I,J,K)
       WRITE(145,147) XP(I),YP(J),ZP(K),UU,VV,WW
        ENDDO
        ENDDO
        ENDDO

       close(144)
       close(145)

      ENDIF


!!C=====SAVE
      IF ( MSUB .EQ. 2 ) THEN
       OPEN(144,FILE='0RHS2.DAT')
       WRITE(144,*) 'VARIABLES="X","Y","Z","u","v","w"'
      WRITE(144,*) 'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'

      OPEN(145,FILE='0u_hat2.DAT')
       WRITE(145,*) 'VARIABLES="X","Y","Z","u","v","w"'
      WRITE(145,*) 'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'  

       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
         RHS11=0.5*(RHS1(I,J,K,1)+RHS1(IPV(I),J,K,1))
         RHS22=0.5*(RHS1(I,J,K,2)+RHS1(I,JPV(J),K,2))
         RHS33=0.5*(RHS1(I,J,K,3)+RHS1(I,J,KPV(K),3))
!         RHS11=RHS1(I,J,K,1)
!         RHS22=RHS1(I,J,K,2)
!         RHS33=RHS1(I,J,K,3)  
       WRITE(144,148) XP(I),YP(J),ZP(K),RHS11,RHS22,RHS33
         UU=0.5*(U(I,J,K)+U(IPV(I),J,K))
         VV=0.5*(V(I,J,K)+V(I,JPV(J),K))
         WW=0.5*(W(I,J,K)+W(I,J,KPV(K)))
!         UU=U(I,J,K)
!         VV=V(I,J,K)
!         WW=W(I,J,K)
       WRITE(145,147) XP(I),YP(J),ZP(K),UU,VV,WW
        ENDDO
        ENDDO
        ENDDO

       close(144)
       close(145)

      ENDIF

!!C=====SAVE
      IF ( MSUB .EQ. 3 ) THEN
       OPEN(144,FILE='0RHS3.DAT')
       WRITE(144,*) 'VARIABLES="X","Y","Z","u","v","w"'
      WRITE(144,*) 'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'
      
      OPEN(145,FILE='0u_hat3.DAT')
       WRITE(145,*) 'VARIABLES="X","Y","Z","u","v","w"'
      WRITE(145,*) 'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'  

       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
         RHS11=0.5*(RHS1(I,J,K,1)+RHS1(IPV(I),J,K,1))
         RHS22=0.5*(RHS1(I,J,K,2)+RHS1(I,JPV(J),K,2))
         RHS33=0.5*(RHS1(I,J,K,3)+RHS1(I,J,KPV(K),3))
!         RHS11=RHS1(I,J,K,1)
!         RHS22=RHS1(I,J,K,2)
!         RHS33=RHS1(I,J,K,3)  
       WRITE(144,148) XP(I),YP(J),ZP(K),RHS11,RHS22,RHS33
         UU=0.5*(U(I,J,K)+U(IPV(I),J,K))
         VV=0.5*(V(I,J,K)+V(I,JPV(J),K))
         WW=0.5*(W(I,J,K)+W(I,J,KPV(K)))
!         UU=U(I,J,K)
!         VV=V(I,J,K)
!         WW=W(I,J,K)
       WRITE(145,147) XP(I),YP(J),ZP(K),UU,VV,WW

        ENDDO
        ENDDO
        ENDDO

       close(144)
       close(145)

      ENDIF

 148  FORMAT(3F18.14,3ES20.12)
 147  FORMAT(3F18.14,3ES20.12) 
       endif

      RETURN
      END

C******************************************************************
      SUBROUTINE RHSX(ITRACKING,RHS1,U,V,W,P,DF_X,VF_X,PSI_XN,RK3XO)
C******************************************************************
C     CALCULATION OF RHS(NAVIER-STOKES)
      USE PARAM_VAR
      USE FLOW_VAR

      USE FLOW_GEOM_VAR

      USE TWO_PHASE_PROPERTY

      IMPLICIT NONE

      REAL RHS1(M1M,M2M,M3M,3)
      REAL RK3XO(M1M,M2M,M3M)

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      REAL PSI_XN(0:M1,0:M2,0:M3)
      REAL DF_X(0:M1,0:M2,0:M3,3),VF_X(0:M1,0:M2,0:M3,3)
      
      INTEGER ITRACKING
      
      INTEGER I,J,K,IP1,IM1,JP1,JM1,KP1,KM1
      REAL UE,UW,UN,US,UT,UB,VN,VS,WT,WB
      REAL ANXX,ANXY,ANXZ,ANX
      REAL ALX1,ALX2,ALX3,ALX4,ALX5,ALX6,ALX7,ALX8,ALX9,ALX10
      REAL ALXX,ALXY,ALXZ,TALXX,TALXY,TALXZ,ALX_SUM1,ALX_SUM2
      REAL PGR_X,RK3X,CN2X,DEN_XN,BUOYANCY

C-----RHS1 calculation for u momentum -----------------
!$OMP PARALLEL DO private(I,J,IP1,IM1,JP1,JM1,KP1,KM1)
!$OMP&private(UE,UW,UN,US,UT,UB,VN,VS,WT,WB)
!$OMP&private(ANXX,ANXY,ANXZ,ANX,ALX1,ALX2,ALX3,ALX4,ALX5,ALX6,ALX7)
!$OMP&private(ALX8,ALX9,ALX10,ALXX,ALXY,ALXZ,TALXX,TALXY,TALXZ)
!$OMP&private(ALX_SUM1,ALX_SUM2,PGR_X,RK3X,CN2X,DEN_XN,BUOYANCY)
       DO 1000 K=1,N3M
         KP1=KPV(K)
         KM1=KMV(K)
       DO 1000 J=1,N2M
         JP1=JPV(J)
         JM1=JMV(J)
       DO 1000 I=IBG,N1M
         IP1=IPV(I)
         IM1=IMV(I)

c-----Cdvection
      UE=0.5*(U(IP1,J,K)+U(I,J,K))
      UW=0.5*(U(I,J,K)+U(IM1,J,K))

      UN=(0.5*(SDY(JP1)*U(I,J,K)+SDY(J)*U(I,JP1,K))*VVDY(JP1))
     &*(1.-FIXJU(J))+U(I,N2,K)*FIXJU(J)
      US=(0.5*(SDY(J)*U(I,JM1,K)+SDY(JM1)*U(I,J,K))*VVDY(J))
     &*(1.-FIXJL(J))+U(I,0,K)*FIXJL(J)

      UT=(0.5*(SDZ(KP1)*U(I,J,K)+SDZ(K)*U(I,J,KP1))*VVDZ(KP1))
     &*(1.-FIXKU(K))+U(I,J,N3)*FIXKU(K)
      UB=(0.5*(SDZ(K)*U(I,J,KM1)+SDZ(KM1)*U(I,J,K))*VVDZ(K))
     &*(1.-FIXKL(K))+U(I,J,0)*FIXKL(K)

      VN=0.5*(SDX(IM1)*V(I,JP1,K)+SDX(I)*V(IM1,JP1,K))*VVDX(I)
      VS=0.5*(SDX(IM1)*V(I,J,K)+SDX(I)*V(IM1,J,K))*VVDX(I)

      WT=0.5*(SDX(IM1)*W(I,J,KP1)+SDX(I)*W(IM1,J,KP1))*VVDX(I)
      WB=0.5*(SDX(IM1)*W(I,J,K)+SDX(I)*W(IM1,J,K))*VVDX(I)

!      IF (ABS(PSI_XN(I,J,K)-0.5) .NE. 0.5 ) THEN  !UPWIND FOR UNIFORM MESH
!!      IF (UE .GE. 0.) THEN
!!        UE2=U(I,J,K)
!!      ELSE
!!        UE2=U(IP1,J,K)
!!      ENDIF
!!      IF (UW .GE. 0.) THEN
!!        UW2=U(IM1,J,K)
!!      ELSE
!!        UW2=U(I,J,K)
!!      ENDIF
!!      IF (VN .GE. 0.) THEN
!!        UN2=U(I,J,K)
!!      ELSE
!!        UN2=U(I,JP1,K)
!!      ENDIF
!!      IF (VS .GE. 0.) THEN
!!        US2=U(I,JM1,K)
!!      ELSE
!!        US2=U(I,J,K)
!!      ENDIF
!!      IF (WT .GE. 0.) THEN
!!        UT2=U(I,J,K)
!!      ELSE
!!        UT2=U(I,J,KP1)
!!      ENDIF
!!      IF (WB .GE. 0.) THEN
!!        UB2=U(I,J,KM1)
!!      ELSE
!!        UB2=U(I,J,K)
!!      ENDIF
!      IF (UE .GE. 0.) THEN
!        UE2=1.5*U(I,J,K)-0.5*U(IM1,J,K)
!      ELSE
!        UE2=1.5*U(IP1,J,K)-0.5*U(IPV(IP1),J,K)
!      ENDIF
!      IF (UW .GE. 0.) THEN
!        UW2=1.5*U(IM1,J,K)-0.5*U(IMV(IM1),J,K)
!      ELSE
!        UW2=1.5*U(I,J,K)-0.5*U(IP1,J,K)
!      ENDIF
!      IF (VN .GE. 0.) THEN
!        UN2=1.5*U(I,J,K)-0.5*U(I,JM1,K)
!      ELSE
!        UN2=1.5*U(I,JP1,K)-0.5*U(I,JPV(JP1),K)
!      ENDIF
!      IF (VS .GE. 0.) THEN
!        US2=1.5*U(I,JM1,K)-0.5*U(I,JMV(JM1),K)
!      ELSE
!        US2=1.5*U(I,J,K)-0.5*U(I,JP1,K)
!      ENDIF
!      IF (WT .GE. 0.) THEN
!        UT2=1.5*U(I,J,K)-0.5*U(I,J,KM1)
!      ELSE
!        UT2=1.5*U(I,J,KP1)-0.5*U(I,J,KPV(KP1))
!      ENDIF
!      IF (WB .GE. 0.) THEN
!        UB2=1.5*U(I,J,KM1)-0.5*U(I,J,KMV(KM1))
!      ELSE
!        UB2=1.5*U(I,J,K)-0.5*U(I,J,KP1)
!      ENDIF
!      ANXX=(DF_X(I,J,K,1)*UE2*UE-DF_X(IM1,J,K,1)*UW2*UW)*VVDX(I)
!      ANXY=(DF_X(I,J,K,2)*UN2*VN-DF_X(I,JM1,K,2)*US2*VS)*SSDY(J)
!      ANXZ=(DF_X(I,J,K,3)*UT2*WT-DF_X(I,J,KM1,3)*UB2*WB)*SSDZ(K)
!
!      ELSE !2ND-ORDER CENTRAL DIFFERENCE
      ANXX=(DF_X(I,J,K,1)*UE**2-DF_X(IM1,J,K,1)*UW**2)*VVDX(I)
      ANXY=(DF_X(I,J,K,2)*UN*VN-DF_X(I,JM1,K,2)*US*VS)*SSDY(J)
      ANXZ=(DF_X(I,J,K,3)*UT*WT-DF_X(I,J,KM1,3)*UB*WB)*SSDZ(K)
!      ENDIF
      ANX=-ANXX-ANXY-ANXZ

      !DIFFUSION
      ALX1=(U(IP1,J,K)-U(I,J,K))*SSDX(I)
      ALX2=(U(I,J,K)-U(IM1,J,K))*SSDX(IM1)
      ALX3=(U(I,JP1,K)-U(I,J,K))*VVDY(JP1)
      ALX4=(U(I,J,K)-U(I,JM1,K))*VVDY(J)
      ALX5=(U(I,J,KP1)-U(I,J,K))*VVDZ(KP1)
      ALX6=(U(I,J,K)-U(I,J,KM1))*VVDZ(K)
      ALXX=(VF_X(I,J,K,1)*ALX1-VF_X(IM1,J,K,1)*ALX2)*VVDX(I)
      ALXY=(VF_X(I,J,K,2)*ALX3-VF_X(I,JM1,K,2)*ALX4)*SSDY(J)
      ALXZ=(VF_X(I,J,K,3)*ALX5-VF_X(I,J,KM1,3)*ALX6)*SSDZ(K)

      ALX7=(V(I,JP1,K)-V(IM1,JP1,K))
      ALX8=(V(I,J,K)-V(IM1,J,K))
      ALX9=(W(I,J,KP1)-W(IM1,J,KP1))
      ALX10=(W(I,J,K)-W(IM1,J,K))
      TALXX=ALXX
      TALXY=(VF_X(I,J,K,2)*ALX7-VF_X(I,JM1,K,2)*ALX8)*SSDY(J)*VVDX(I)
      TALXZ=(VF_X(I,J,K,3)*ALX9-VF_X(I,J,KM1,3)*ALX10)*SSDZ(K)*VVDX(I)

      ALX_SUM1=ALXX+ALXY+ALXZ+TALXX
      ALX_SUM2=TALXY+TALXZ

      !OTHER TERM
      PGR_X=(P(I,J,K)-P(IM1,J,K))*VVDX(I)

      !BUOYANCY
      !DEN_XN=(DENM+DEN_DIFF*PSI_XN(I,J,K))
      BUOYANCY=0.!0.5*(DEN_XN+DENF_X(I,J,K) )*GRAVITY

      RK3X=ANX+ALX_SUM2
      CN2X=ALX_SUM1
!      RK3X=ANX+ALX_SUM1+ALX_SUM2
!      CN2X=0.

      RHS1(I,J,K,1)=(DENM+DEN_DIFF*PSI_XN(I,J,K))*U(I,J,K)
     &        +( GAMMA(MSUB)*RK3X+RO(MSUB)*RK3XO(I,J,K)
     &        +2.*ALPHA_RK3*CN2X
     &        +2.*ALPHA_RK3*(SUR_X(I,J,K)-PGR_X)
     &        -2.*ALPHA_RK3*BUOYANCY
     &        )*DT
      RK3XO(I,J,K)=RK3X

 1000  CONTINUE

      RETURN
      END

C******************************************************************
      SUBROUTINE RHSY(ITRACKING,RHS1,U,V,W,P,DF_Y,VF_Y,PSI_YN,RK3YO)
C******************************************************************
C     CALCULATION OF RHS(NAVIER-STOKES)
      USE PARAM_VAR
      USE FLOW_VAR

      USE FLOW_GEOM_VAR

      USE TWO_PHASE_PROPERTY
      
      IMPLICIT NONE

      REAL RHS1(M1M,M2M,M3M,3)
      REAL RK3YO(M1M,M2M,M3M)

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      REAL PSI_YN(0:M1,0:M2,0:M3)
      REAL DF_Y(0:M1,0:M2,0:M3,3),VF_Y(0:M1,0:M2,0:M3,3)

      INTEGER ITRACKING
      
      INTEGER I,J,K,IP1,IM1,JP1,JM1,KP1,KM1
      REAL VE,VW,VN,VS,VT,VB,UE,UW,WT,WB
      REAL ANYX,ANYY,ANYZ,ANY
      REAL ALY1,ALY2,ALY3,ALY4,ALY5,ALY6,ALY7,ALY8,ALY9,ALY10
      REAL ALYX,ALYY,ALYZ,TALYX,TALYY,TALYZ,ALY_SUM1,ALY_SUM2
      REAL PGR_Y,RK3Y,CN2Y,DEN_YN,BUOYANCY
      
C-----RHS1 calculation for v momentum -----------------
!$OMP PARALLEL DO private(I,J,IP1,IM1,JP1,JM1,KP1,KM1)
!$OMP&private(VE,VW,VN,VS,VT,VB,UE,UW,WT,WB)
!$OMP&private(ANYX,ANYY,ANYZ,ANY,ALY1,ALY2,ALY3,ALY4,ALY5,ALY6,ALY7)
!$OMP&private(ALY8,ALY9,ALY10,ALYX,ALYY,ALYZ,TALYX,TALYY,TALYZ)
!$OMP&private(ALY_SUM1,ALY_SUM2,PGR_Y,RK3Y,CN2Y,DEN_YN,BUOYANCY)
       DO 1000 K=1,N3M
         KP1=KPV(K)
         KM1=KMV(K)
       DO 1000 J=JBG,N2M
         JP1=JPV(J)
         JM1=JMV(J)
       DO 1000 I=1,N1M
         IP1=IPV(I)
         IM1=IMV(I)

C-----CD2 FOR CONVECTION
      UE=0.5*(SDY(JM1)*U(IP1,J,K)+SDY(J)*U(IP1,JM1,K))*VVDY(J)
      UW=0.5*(SDY(JM1)*U(I,J,K)+SDY(J)*U(I,JM1,K))*VVDY(J)

      VE=(0.5*(SDX(IP1)*V(I,J,K)+SDX(I)*V(IP1,J,K))*VVDX(IP1))
     &*(1.-FIXIU(I))+V(N1,J,K)*FIXIU(I)
      VW=(0.5*(SDX(I)*V(IM1,J,K)+SDX(IM1)*V(I,J,K))*VVDX(I))
     &*(1.-FIXIL(I))+V(0,J,K)*FIXIL(I)

      VN=0.5*(V(I,JP1,K)+V(I,J,K))
      VS=0.5*(V(I,J,K)+V(I,JM1,K))

      VT=(0.5*(SDZ(KP1)*V(I,J,K)+SDZ(K)*V(I,J,KP1))*VVDZ(KP1))
     &*(1.-FIXKU(K))+V(I,J,N3)*FIXKU(K)
      VB=(0.5*(SDZ(K)*V(I,J,KM1)+SDZ(KM1)*V(I,J,K))*VVDZ(K))
     &*(1.-FIXKL(K))+V(I,J,0)*FIXKL(K)

      WT=0.5*(SDY(JM1)*W(I,J,KP1)+SDY(J)*W(I,JM1,KP1))*VVDY(J)
      WB=0.5*(SDY(JM1)*W(I,J,K)+SDY(J)*W(I,JM1,K))*VVDY(J)

!      IF (ABS(PSI_YN(I,J,K)-0.5) .NE. 0.5 ) THEN   !UPWIND FOR UNIFORM MESH
!!      IF (UE .GE. 0.) THEN
!!        VE2=V(I,J,K)
!!      ELSE
!!        VE2=V(IP1,J,K)
!!      ENDIF
!!      IF (UW .GE. 0.) THEN
!!        VW2=V(IM1,J,K)
!!      ELSE
!!        VW2=V(I,J,K)
!!      ENDIF
!!      IF (VN .GE. 0.) THEN
!!        VN2=V(I,J,K)
!!      ELSE
!!        VN2=V(I,JP1,K)
!!      ENDIF
!!      IF (VS .GE. 0.) THEN
!!        VS2=V(I,JM1,K)
!!      ELSE
!!        VS2=V(I,J,K)
!!      ENDIF
!!      IF (WT .GE. 0.) THEN
!!        VT2=V(I,J,K)
!!      ELSE
!!        VT2=V(I,J,KP1)
!!      ENDIF
!!      IF (WB .GE. 0.) THEN
!!        VB2=V(I,J,KM1)
!!      ELSE
!!        VB2=V(I,J,K)
!!      ENDIF
!      IF (UE .GE. 0.) THEN
!        VE2=1.5*V(I,J,K)-0.5*V(IM1,J,K)
!      ELSE
!        VE2=1.5*V(IP1,J,K)-0.5*V(IPV(IP1),J,K)
!      ENDIF
!      IF (UW .GE. 0.) THEN
!        VW2=1.5*V(IM1,J,K)-0.5*V(IMV(IM1),J,K)
!      ELSE
!        VW2=1.5*V(I,J,K)-0.5*V(IP1,J,K)
!      ENDIF
!      IF (VN .GE. 0.) THEN
!        VN2=1.5*V(I,J,K)-0.5*V(I,JM1,K)
!      ELSE
!        VN2=1.5*V(I,JP1,K)-0.5*V(I,JPV(JP1),K)
!      ENDIF
!      IF (VS .GE. 0.) THEN
!        VS2=1.5*V(I,JM1,K)-0.5*V(I,JMV(J),K)
!      ELSE
!        VS2=1.5*V(I,J,K)-0.5*V(I,JP1,K)
!      ENDIF
!      IF (WT .GE. 0.) THEN
!        VT2=1.5*V(I,J,K)-0.5*V(I,J,KM1)
!      ELSE
!        VT2=1.5*V(I,J,KP1)-0.5*V(I,J,KPV(KP1))
!      ENDIF
!      IF (WB .GE. 0.) THEN
!        VB2=1.5*V(I,J,KM1)-0.5*V(I,J,KMV(KM1))
!      ELSE
!        VB2=1.5*V(I,J,K)-0.5*V(I,J,KP1)
!      ENDIF
!      ANYX=(DF_Y(I,J,K,1)*UE*VE2-DF_Y(IM1,J,K,1)*UW*VW2)*SSDX(I)
!      ANYY=(DF_Y(I,J,K,2)*VN*VN2-DF_Y(I,JM1,K,2)*VS*VS2)*VVDY(J)
!      ANYZ=(DF_Y(I,J,K,3)*WT*VT2-DF_Y(I,J,KM1,3)*WB*VB2)*SSDZ(K)
!
!      ELSE !2ND-ORDER CENTRAL DIFFERENCE
      ANYX=(DF_Y(I,J,K,1)*UE*VE-DF_Y(IM1,J,K,1)*UW*VW)*SSDX(I)
      ANYY=(DF_Y(I,J,K,2)*VN**2-DF_Y(I,JM1,K,2)*VS**2)*VVDY(J)
      ANYZ=(DF_Y(I,J,K,3)*VT*WT-DF_Y(I,J,KM1,3)*VB*WB)*SSDZ(K)
!      ENDIF
      ANY=-ANYX-ANYY-ANYZ

      !DIFFUSION
      ALY1=(V(IP1,J,K)-V(I,J,K))*VVDX(IP1)
      ALY2=(V(I,J,K)-V(IM1,J,K))*VVDX(I)
      ALY3=(V(I,JP1,K)-V(I,J,K))*SSDY(J)
      ALY4=(V(I,J,K)-V(I,JM1,K))*SSDY(JM1)
      ALY5=(V(I,J,KP1)-V(I,J,K))*VVDZ(KP1)
      ALY6=(V(I,J,K)-V(I,J,KM1))*VVDZ(K)
      ALYX=(VF_Y(I,J,K,1)*ALY1-VF_Y(IM1,J,K,1)*ALY2)*SSDX(I)
      ALYY=(VF_Y(I,J,K,2)*ALY3-VF_Y(I,JM1,K,2)*ALY4)*VVDY(J)
      ALYZ=(VF_Y(I,J,K,3)*ALY5-VF_Y(I,J,KM1,3)*ALY6)*SSDZ(K)

      ALY7=(U(IP1,J,K)-U(IP1,JM1,K))
      ALY8=(U(I,J,K)-U(I,JM1,K))
      ALY9=(W(I,J,KP1)-W(I,JM1,KP1))
      ALY10=(W(I,J,K)-W(I,JM1,K))
      TALYX=(VF_Y(I,J,K,1)*ALY7-VF_Y(IM1,J,K,1)*ALY8)*SSDX(I)*VVDY(J)
      TALYY=ALYY
      TALYZ=(VF_Y(I,J,K,3)*ALY9-VF_Y(I,J,KM1,3)*ALY10)*SSDZ(K)*VVDY(J)

      ALY_SUM1=ALYX+ALYY+ALYZ+TALYY
      ALY_SUM2=TALYX+TALYZ

      !OTHER TERM
       PGR_Y=(P(I,J,K)-P(I,JM1,K))*VVDY(J)

      !BUOYANCY
       !DEN_YN=DENM+DEN_DIFF*PSI_YN(I,J,K)
      BUOYANCY=0.!0.5*(DEN_YN+DENF_Y(I,J,K) )*GRAVITY

      RK3Y=ANY+ALY_SUM2
      CN2Y=ALY_SUM1
!      RK3Y=ANY+ALY_SUM2+ALY_SUM1
!      CN2Y=0.

      RHS1(I,J,K,2)=(DENM+DEN_DIFF*PSI_YN(I,J,K))*V(I,J,K)
     &       +( GAMMA(MSUB)*RK3Y+RO(MSUB)*RK3YO(I,J,K)
     &         +2.*ALPHA_RK3*CN2Y
     &         +2.*ALPHA_RK3*(SUR_Y(I,J,K)-PGR_Y)
     &         -2.*ALPHA_RK3*BUOYANCY
     &        )*DT
      RK3YO(I,J,K)=RK3Y

 1000  CONTINUE

      RETURN
      END

C******************************************************************
      SUBROUTINE RHSZ(ITRACKING,RHS1,U,V,W,P,DF_Z,VF_Z,PSI_ZN,DENF_Z
     &,RK3ZO)
C******************************************************************
C     CALCULATION OF RHS(NAVIER-STOKES)
      USE PARAM_VAR
      USE FLOW_VAR

      USE FLOW_GEOM_VAR

      USE TWO_PHASE_PROPERTY
      
      IMPLICIT NONE

      REAL RHS1(M1M,M2M,M3M,3)
      REAL RK3ZO(M1M,M2M,M3M)

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      REAL DENF_Z(0:M1,0:M2,0:M3)
      REAL PSI_ZN(0:M1,0:M2,0:M3)
      REAL DF_Z(0:M1,0:M2,0:M3,3),VF_Z(0:M1,0:M2,0:M3,3)

      INTEGER ITRACKING

      INTEGER I,J,K,IP1,IM1,JP1,JM1,KP1,KM1
      REAL UE,UW,VN,VS,WE,WW,WN,WS,WT,WB
      REAL ANZX,ANZY,ANZZ,ANZ
      REAL ALZ1,ALZ2,ALZ3,ALZ4,ALZ5,ALZ6,ALZ7,ALZ8,ALZ9,ALZ10
      REAL ALZX,ALZY,ALZZ,TALZX,TALZY,TALZZ,ALZ_SUM1,ALZ_SUM2
      REAL PGR_Z,RK3Z,CN2Z,DEN_ZN,BUOYANCY

      REAL FUNCBODY

C-----RHS1 calculation for w momentum -----------------
!$OMP PARALLEL DO private(I,J,IP1,IM1,JP1,JM1,KP1,KM1)
!$OMP&private(UE,UW,VN,VS,WE,WW,WN,WS,WT,WB)
!$OMP&private(ANZX,ANZY,ANZZ,ANZ,ALZ1,ALZ2,ALZ3,ALZ4,ALZ5,ALZ6,ALZ7)
!$OMP&private(ALZ8,ALZ9,ALZ10,ALZX,ALZY,ALZZ,TALZX,TALZY,TALZZ)
!$OMP&private(ALZ_SUM1,ALZ_SUM2,PGR_Z,RK3Z,CN2Z,DEN_ZN,BUOYANCY)
       DO 1000 K=KBG,N3M
         KP1=KPV(K)
         KM1=KMV(K)
       DO 1000 J=1,N2M
         JP1=JPV(J)
         JM1=JMV(J)
       DO 1000 I=1,N1M
         IP1=IPV(I)
         IM1=IMV(I)

      !CONVECTION
      UE=0.5*(SDZ(KM1)*U(IP1,J,K)+SDZ(K)*U(IP1,J,KM1))*VVDZ(K)
      UW=0.5*(SDZ(KM1)*U(I,J,K)+SDZ(K)*U(I,J,KM1))*VVDZ(K)

      VN=0.5*(SDZ(KM1)*V(I,JP1,K)+SDZ(K)*V(I,JP1,KM1))*VVDZ(K)
      VS=0.5*(SDZ(KM1)*V(I,J,K)+SDZ(K)*V(I,J,KM1))*VVDZ(K)

      WE=(0.5*(SDX(IP1)*W(I,J,K)+SDX(I)*W(IP1,J,K))*VVDX(IP1))
     &*(1.-FIXIU(I))+FIXIU(I)*W(N1,J,K)
      WW=(0.5*(SDX(I)*W(IM1,J,K)+SDX(IM1)*W(I,J,K))*VVDX(I))
     &*(1.-FIXIL(I))+FIXIL(I)*W(0,J,K)

      WN=(0.5*(SDY(JP1)*W(I,J,K)+SDY(J)*W(I,JP1,K))*VVDY(JP1))
     & *(1.-FIXJU(J))+FIXJU(J)*W(I,N2,K)
      WS=(0.5*(SDY(J)*W(I,JM1,K)+SDY(JM1)*W(I,J,K))*VVDY(J))
     & *(1.-FIXJL(J))+FIXJL(J)*W(I,0,K)

      WT=0.5*(W(I,J,KP1)+W(I,J,K))
      WB=0.5*(W(I,J,K)+W(I,J,KM1))

!      IF (ABS(PSI_ZN(I,J,K)-0.5) .NE. 0.5 ) THEN !UPWIND FOR UNIFORM MESH
!!      IF (UE .GE. 0.) THEN
!!        WE2=W(I,J,K)
!!      ELSE
!!        WE2=W(IP1,J,K)
!!      ENDIF
!!      IF (UW .GE. 0.) THEN
!!        WW2=W(IM1,J,K)
!!      ELSE
!!        WW2=W(I,J,K)
!!      ENDIF
!!      IF (VN .GE. 0.) THEN
!!        WN2=W(I,J,K)
!!      ELSE
!!        WN2=W(I,JP1,K)
!!      ENDIF
!!      IF (VS .GE. 0.) THEN
!!        WS2=W(I,JM1,K)
!!      ELSE
!!        WS2=W(I,J,K)
!!      ENDIF
!!      IF (WT .GE. 0.) THEN
!!        WT2=W(I,J,K)
!!      ELSE
!!        WT2=W(I,J,KP1)
!!      ENDIF
!!      IF (WB .GE. 0.) THEN
!!        WB2=W(I,J,KM1)
!!      ELSE
!!        WB2=W(I,J,K)
!!      ENDIF
!      IF (UE .GE. 0.) THEN
!        WE2=1.5*W(I,J,K)-0.5*W(I,IM1,K)
!      ELSE
!        WE2=1.5*W(IP1,J,K)-0.5*W(IPV(IP1),J,K)
!      ENDIF
!      IF (UW .GE. 0.) THEN
!        WW2=1.5*W(IM1,J,K)-0.5*W(IMV(IM1),J,K)
!      ELSE
!        WW2=1.5*W(I,J,K)-0.5*W(IP1,J,K)
!      ENDIF
!      IF (VN .GE. 0.) THEN
!        WN2=1.5*W(I,J,K)-0.5*W(I,JM1,K)
!      ELSE
!        WN2=1.5*W(I,JP1,K)-0.5*W(I,JPV(JP1),K)
!      ENDIF
!      IF (VS .GE. 0.) THEN
!        WS2=1.5*W(I,JM1,K)-0.5*W(I,JMV(JM1),K)
!      ELSE
!        WS2=1.5*W(I,J,K)-0.5*W(I,JP1,K)
!      ENDIF
!      IF (WT .GE. 0.) THEN
!        WT2=1.5*W(I,J,K)-0.5*W(I,J,KM1)
!      ELSE
!        WT2=1.5*W(I,J,KP1)-0.5*W(I,J,KPV(KP1))
!      ENDIF
!      IF (WB .GE. 0.) THEN
!        WB2=1.5*W(I,J,KM1)-0.5*W(I,J,KMV(KM1))
!      ELSE
!        WB2=1.5*W(I,J,K)-0.5*W(I,J,KP1)
!      ENDIF
!
!      ANZX=(DF_Z(I,J,K,1)*UE*WE2-DF_Z(IM1,J,K,1)*UW*WW2)*SSDX(I)
!      ANZY=(DF_Z(I,J,K,2)*VN*WN2-DF_Z(I,JM1,K,2)*VS*WS2)*SSDY(J)
!      ANZZ=(DF_Z(I,J,K,3)*WT*WT2-DF_Z(I,J,KM1,3)*WB*WB2)*VVDZ(K)
!
!      ELSE !2ND-ORDER CENTRAL DIFFERENCE
      ANZX=(DF_Z(I,J,K,1)*UE*WE-DF_Z(IM1,J,K,1)*UW*WW)*SSDX(I)
      ANZY=(DF_Z(I,J,K,2)*VN*WN-DF_Z(I,JM1,K,2)*VS*WS)*SSDY(J)
      ANZZ=(DF_Z(I,J,K,3)*WT**2-DF_Z(I,J,KM1,3)*WB**2)*VVDZ(K)
!      ENDIF
      ANZ=-ANZX-ANZY-ANZZ                     

      !DIFFUSION
      ALZ1=(W(IP1,J,K)-W(I,J,K))*VVDX(IP1)
      ALZ2=(W(I,J,K)-W(IM1,J,K))*VVDX(I)
      ALZ3=(W(I,JP1,K)-W(I,J,K))*VVDY(JP1)
      ALZ4=(W(I,J,K)-W(I,JM1,K))*VVDY(J)
      ALZ5=(W(I,J,KP1)-W(I,J,K))*SSDZ(K)
      ALZ6=(W(I,J,K)-W(I,J,KM1))*SSDZ(KM1)
      ALZX=(VF_Z(I,J,K,1)*ALZ1-VF_Z(IM1,J,K,1)*ALZ2)*SSDX(I)
      ALZY=(VF_Z(I,J,K,2)*ALZ3-VF_Z(I,JM1,K,2)*ALZ4)*SSDY(J)
      ALZZ=(VF_Z(I,J,K,3)*ALZ5-VF_Z(I,J,KM1,3)*ALZ6)*VVDZ(K)

      ALZ7=(U(IP1,J,K)-U(IP1,J,KM1))
      ALZ8=(U(I,J,K)-U(I,J,KM1))
      ALZ9=(V(I,JP1,K)-V(I,JP1,KM1))
      ALZ10=(V(I,J,K)-V(I,J,KM1))
      TALZX=(VF_Z(I,J,K,1)*ALZ7-VF_Z(IM1,J,K,1)*ALZ8)*SSDX(I)*VVDZ(K)
      TALZY=(VF_Z(I,J,K,2)*ALZ9-VF_Z(I,JM1,K,2)*ALZ10)*SSDY(J)*VVDZ(K)
      TALZZ=ALZZ

      ALZ_SUM1=ALZX+ALZY+ALZZ+TALZZ
      ALZ_SUM2=TALZX+TALZY

      !OTHER TERM
       PGR_Z=(P(I,J,K)-P(I,J,KM1))*VVDZ(K)

      !BUOYANCY
      DEN_ZN=(DENM+DEN_DIFF*PSI_ZN(I,J,K))
      BUOYANCY=0.5*(DEN_ZN+DENF_Z(I,J,K) )*GRAVITY

      RK3Z=ANZ+ALZ_SUM2
      CN2Z=ALZ_SUM1
!      RK3Z=ANZ+ALZ_SUM1+ALZ_SUM2
!      CN2Z=0.

       IF ( FUNCBODY(XP(I),YP(J),Z(K)) .GT. 0. ) THEN

!      if (i.eq.48 .and. j.eq.52 .and. k.eq.43) then
!        write(*,*) rhs1(i,j,k,3)
!     &      ,DEN_ZN*W(I,J,K)+
!     &       ( GAMMA(MSUB)*RK3Z+RO(MSUB)*RK3ZO(I,J,K)
!     &         +2.*ALPHA_RK3*CN2Z
!     &         +2.*ALPHA_RK3*(SUR_Z(I,J,K)-PGR_Z)
!     &         -2.*ALPHA_RK3*BUOYANCY
!     &         -2.*ALPHA_RK3*(PMI+PMI_GRAV) )*DT
!     
!       WRITE(*,*)   RO(MSUB),RK3ZO(I,J,K)
!
!      endif

      RHS1(I,J,K,3)=DEN_ZN*W(I,J,K)+
     &       ( GAMMA(MSUB)*RK3Z+RO(MSUB)*RK3ZO(I,J,K)
     &         +2.*ALPHA_RK3*CN2Z
     &         +2.*ALPHA_RK3*(SUR_Z(I,J,K)-PGR_Z)
     &         -2.*ALPHA_RK3*BUOYANCY
     &         -2.*ALPHA_RK3*(PMI+PMI_GRAV) )*DT

     
       ELSE
      RHS1(I,J,K,3)=DEN_ZN*W(I,J,K)+
     &       ( GAMMA(MSUB)*RK3Z+RO(MSUB)*RK3ZO(I,J,K)
     &         +2.*ALPHA_RK3*CN2Z
     &         +2.*ALPHA_RK3*(SUR_Z(I,J,K)-PGR_Z)
     &       )*DT
       ENDIF

      RK3ZO(I,J,K)=RK3Z

 1000  CONTINUE
 
      RETURN
      END

C******************************************************************
      SUBROUTINE LHSU(RHS1,U,V,W,DF_X,VF_X,DENF_XI)
C******************************************************************
      USE PARAM_VAR
      
      USE FLOW_VAR
      USE FLOW_GEOM_VAR

      IMPLICIT NONE

      REAL RHS1(M1M,M2M,M3M,3)
      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL DF_X(0:M1,0:M2,0:M3,3),VF_X(0:M1,0:M2,0:M3,3)
      REAL DENF_XI(0:M1,0:M2,0:M3)

      real, dimension (:,:), allocatable :: AI,BI,CI,GI
      real, dimension (:,:), allocatable :: AJ,BJ,CJ,GJ
      real, dimension (:,:), allocatable :: AK,BK,CK,GK

      INTEGER I,J,K
      REAL DI1,DI2,DI3,DJ1,DJ2,DJ3,DK1,DK2,DK3
      
C=====ADI STARTS

      IF (N3M.EQ.1) GOTO 100

C-----Z-DIRECTION
!$OMP PARALLEL private(I,K,DK1,DK2,DK3)
!$OMP&private(AK,CK,BK,GK)
      allocate(AK(M1,M3),BK(M1,M3),CK(M1,M3),GK(M1,M3))
!$OMP DO
      DO 131 J=1,N2M

      DO 141 K=1,N3M
      DO 141 I=IBG,N1M

!      DK1=ACOEF*DENF_XI(I,J,K)*VF_X(I,J,KMV(K),3)*AKUV(K)*(1.-FIXKL(K)) !NEWMANN
!      DK3=ACOEF*DENF_XI(I,J,K)*VF_X(I,J,K,3)*CKUV(K)*(1.-FIXKU(K))
      DK1=ACOEF*DENF_XI(I,J,K)*VF_X(I,J,KMV(K),3)*AKUV(K) !PERIODIC
      DK3=ACOEF*DENF_XI(I,J,K)*VF_X(I,J,K,3)*CKUV(K)
      DK2=-DK1-DK3

      AK(I,K)=DK1
      CK(I,K)=DK3
      BK(I,K)=1.+DK2
      GK(I,K)=RHS1(I,J,K,1)

      !BOUNDARY CONDITION IS NOT NEEDED BECAUSE OF periodicity.

  141 CONTINUE
      IF (IPZ .EQ. 1) THEN
       CALL TRDIAG3P(AK,BK,CK,GK,1,N3M,IBG,N1M)
      ELSE
       CALL TRDIAG3(AK,BK,CK,GK,GK,1,N3M,IBG,N1M)
      ENDIF

      DO 151 K=1,N3M
      DO 151 I=IBG,N1M
      RHS1(I,J,K,1)=GK(I,K)
  151 CONTINUE

  131 CONTINUE
!$OMP END DO
      deallocate(AK,BK,CK,GK)
!$OMP END PARALLEL
      
  100 CONTINUE

!$OMP PARALLEL private(I,J,DJ1,DJ2,DJ3,DI1,DI2,DI3)
!$OMP&private(AJ,BJ,CJ,GJ,AI,BI,CI,GI)
      allocate(AI(M2,M1),BI(M2,M1),CI(M2,M1),GI(M2,M1))
      allocate(AJ(M1,M2),BJ(M1,M2),CJ(M1,M2),GJ(M1,M2))
!$OMP DO
      DO 40 K=1,N3M

C-----Y-DIRECTION
      DO 91 J=1,N2M
      DO 91 I=IBG,N1M

!      DJ1=ACOEF*DENF_XI(I,J,K)*VF_X(I,JMV(J),K,2)*AJUW(J)*(1.-FIXJL(J)) !NEUMANN
!      DJ3=ACOEF*DENF_XI(I,J,K)*VF_X(I,J,K,2)*CJUW(J)*(1.-FIXJU(J))
      DJ1=ACOEF*DENF_XI(I,J,K)*VF_X(I,JMV(J),K,2)*AJUW(J)  !PERIODIC
      DJ3=ACOEF*DENF_XI(I,J,K)*VF_X(I,J,K,2)*CJUW(J)
      DJ2=-DJ1-DJ3

      AJ(I,J)=DJ1
      CJ(I,J)=DJ3
      BJ(I,J)=1.+DJ2
      GJ(I,J)=RHS1(I,J,K,1)
   91 CONTINUE

      !DEL_U(I,0,K) AND DEL_U(I,N2,K) IS ZERO WITH CONVERGENCE
!c     boundary condition (USE PREVIOUS STEP VALUE)
!      DO I=1,N1M
!      GJ(I,1)=GJ(I,1)-AJ(I,1)*(U(I,0,K)-UR(I,0,K))
!      GJ(I,N2M)=GJ(I,N2M)-CJ(I,N2M)*(U(I,N2,K)-UR(I,N2,K))
!      ENDDO

      IF (IPY .EQ. 1) THEN
      CALL TRDIAG2P(AJ,BJ,CJ,GJ,1,N2M,IBG,N1M)
      ELSE
      CALL TRDIAG2(AJ,BJ,CJ,GJ,GJ,1,N2M,IBG,N1M)
      ENDIF

C-----X-DIRECTION
      DO 51 J=1,N2M
      DO 51 I=IBG,N1M

      DI1=2.*ACOEF*DENF_XI(I,J,K)*VF_X(IMV(I),J,K,1)*AIU(I)
      DI3=2.*ACOEF*DENF_XI(I,J,K)*VF_X(I,J,K,1)*CIU(I)
      DI2=-DI1-DI3

      AI(J,I)=DI1
      CI(J,I)=DI3
      BI(J,I)=1.+DI2
      GI(J,I)=GJ(I,J)
   51 CONTINUE

!c     boundary condition (USE PREVIOUS STEP VALUE)
!      DO J=1,N2M
!      GI(J,1)=GI(J,1)-AI(J,1)*(U(0,J,K)-UR(0,J,K))
!      GI(J,N1M)=GI(J,N1M)-CI(J,N1M)*(U(N1,J,K)-UR(N1,J,K))
!      ENDDO

      IF (IPX .EQ. 1) THEN
       CALL TRDIAG1P(AI,BI,CI,GI,IBG,N1M,1,N2M)
      ELSE
       CALL TRDIAG1(AI,BI,CI,GI,GI,IBG,N1M,1,N2M)
      ENDIF

      DO 61 J=1,N2M
      DO 61 I=IBG,N1M
        U(I,J,K)=GI(J,I)+U(I,J,K)
   61 CONTINUE

   40 CONTINUE
!$OMP END DO
      deallocate(AI,BI,CI,GI)
      deallocate(AJ,BJ,CJ,GJ)
!$OMP END PARALLEL

      RETURN
      END

C******************************************************************
      SUBROUTINE LHSV(RHS1,U,V,W,DF_Y,VF_Y,DENF_YI)
C******************************************************************
      USE PARAM_VAR
      
      USE FLOW_VAR
      USE FLOW_GEOM_VAR

      IMPLICIT NONE

      REAL RHS1(M1M,M2M,M3M,3)
      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL DF_Y(0:M1,0:M2,0:M3,3),VF_Y(0:M1,0:M2,0:M3,3)
      REAL DENF_YI(0:M1,0:M2,0:M3)

      real, dimension (:,:), allocatable :: AI,BI,CI,GI
      real, dimension (:,:), allocatable :: AJ,BJ,CJ,GJ
      real, dimension (:,:), allocatable :: AK,BK,CK,GK

      INTEGER I,J,K
      REAL DI1,DI2,DI3,DJ1,DJ2,DJ3,DK1,DK2,DK3
      
C=====ADI STARTS

      IF (N3M.EQ.1) GOTO 100
C-----Z-DIRECTION
!$OMP PARALLEL private(I,K,DK1,DK2,DK3)
!$OMP&private(AK,CK,BK,GK)
      allocate(AK(M1,M3),BK(M1,M3),CK(M1,M3),GK(M1,M3))
!$OMP DO
      DO 131 J=JBG,N2M

      DO 141 K=1,N3M
      DO 141 I=1,N1M

!      DK1=ACOEF*DENF_YI(I,J,K)*VF_Y(I,J,KMV(K),3)*AKUV(K)*(1.-FIXKL(K))  !NEUMANN
!      DK3=ACOEF*DENF_YI(I,J,K)*VF_Y(I,J,K,3)*CKUV(K)*(1.-FIXKU(K))
      DK1=ACOEF*DENF_YI(I,J,K)*VF_Y(I,J,KMV(K),3)*AKUV(K)  !PERIODIC
      DK3=ACOEF*DENF_YI(I,J,K)*VF_Y(I,J,K,3)*CKUV(K)
      DK2=-DK1-DK3

      AK(I,K)=DK1
      CK(I,K)=DK3
      BK(I,K)=1.+DK2
      GK(I,K)=RHS1(I,J,K,2)

      !BOUNDARY CONDITION IS NOT NEEDED BECAUSE OF periodicity.

  141 CONTINUE
      IF (IPZ .EQ. 1) THEN
       CALL TRDIAG3P(AK,BK,CK,GK,1,N3M,1,N1M)
      ELSE
       CALL TRDIAG3(AK,BK,CK,GK,GK,1,N3M,1,N1M)
      ENDIF

      DO 151 K=1,N3M
      DO 151 I=1,N1M
      RHS1(I,J,K,2)=GK(I,K)
  151 CONTINUE

  131 CONTINUE
!$OMP END DO
      deallocate(AK,BK,CK,GK)
!$OMP END PARALLEL

  100 CONTINUE

!$OMP PARALLEL private(I,J,DJ1,DJ2,DJ3,DI1,DI2,DI3)
!$OMP&private(AJ,BJ,CJ,GJ,AI,BI,CI,GI)
      allocate(AI(M2,M1),BI(M2,M1),CI(M2,M1),GI(M2,M1))
      allocate(AJ(M1,M2),BJ(M1,M2),CJ(M1,M2),GJ(M1,M2))
!$OMP DO
      DO 40 K=1,N3M

C-----Y-DIRECTION
      DO 91 J=JBG,N2M
      DO 91 I=1,N1M

      DJ1=2.*ACOEF*DENF_YI(I,J,K)*VF_Y(I,JMV(J),K,2)*AJV(J)
      DJ3=2.*ACOEF*DENF_YI(I,J,K)*VF_Y(I,J,K,2)*CJV(J)
      DJ2=-DJ1-DJ3

      AJ(I,J)=DJ1
      CJ(I,J)=DJ3
      BJ(I,J)=1.+DJ2
      GJ(I,J)=RHS1(I,J,K,2)
   91 CONTINUE

!c     boundary condition (USE PREVIOUS STEP VALUE)
!      DO I=1,N1M
!      GJ(I,1)=GJ(I,1)-AJ(I,1)*(V(I,0,K)-VR(I,0,K))
!      GJ(I,N2M)=GJ(I,N2M)-CJ(I,N2M)*(V(I,N2,K)-VR(I,N2,K))
!      ENDDO
      IF (IPY .EQ. 1) THEN
       CALL TRDIAG2P(AJ,BJ,CJ,GJ,JBG,N2M,1,N1M)
      ELSE
       CALL TRDIAG2(AJ,BJ,CJ,GJ,GJ,JBG,N2M,1,N1M)
      ENDIF

C-----X-DIRECTION
      DO 51 J=JBG,N2M
      DO 51 I=1,N1M

!      DI1=ACOEF*DENF_YI(I,J,K)*VF_Y(IMV(I),J,K,1)*AIVW(I)*(1.-FIXIL(I)) !NEUMANN
!      DI3=ACOEF*DENF_YI(I,J,K)*VF_Y(I,J,K,1)*CIVW(I)*(1.-FIXIU(I))
      DI1=ACOEF*DENF_YI(I,J,K)*VF_Y(IMV(I),J,K,1)*AIVW(I)  !PERIODIC
      DI3=ACOEF*DENF_YI(I,J,K)*VF_Y(I,J,K,1)*CIVW(I)
      DI2=-DI1-DI3

      AI(J,I)=DI1
      CI(J,I)=DI3
      BI(J,I)=1.+DI2
      GI(J,I)=GJ(I,J)
   51 CONTINUE

!c     boundary condition (USE PREVIOUS STEP VALUE)
!      DO J=1,N2M
!      GI(J,1)=GI(J,1)-AI(J,1)*(V(0,J,K)-VR(0,J,K))
!      GI(J,N1M)=GI(J,N1M)-CI(J,N1M)*(V(N1,J,K)-VR(N1,J,K))
!      ENDDO
      IF (IPX .EQ. 1) THEN
       CALL TRDIAG1P(AI,BI,CI,GI,1,N1M,JBG,N2M)
      ELSE
       CALL TRDIAG1(AI,BI,CI,GI,GI,1,N1M,JBG,N2M)
      ENDIF

      DO 61 J=JBG,N2M
      DO 61 I=1,N1M
      V(I,J,K)=GI(J,I)+V(I,J,K)
   61 CONTINUE

   40 CONTINUE
!$OMP END DO
      deallocate(AI,BI,CI,GI)
      deallocate(AJ,BJ,CJ,GJ)
!$OMP END PARALLEL

      RETURN
      END

C******************************************************************
      SUBROUTINE LHSW(RHS1,U,V,W,DF_Z,VF_Z,DENF_ZI)
C******************************************************************
      USE PARAM_VAR
      
      USE FLOW_VAR
      USE FLOW_GEOM_VAR

      IMPLICIT NONE

      REAL RHS1(M1M,M2M,M3M,3)
      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL DF_Z(0:M1,0:M2,0:M3,3),VF_Z(0:M1,0:M2,0:M3,3)
      REAL DENF_ZI(0:M1,0:M2,0:M3)

      real, dimension (:,:), allocatable :: AI,BI,CI,GI
      real, dimension (:,:), allocatable :: AJ,BJ,CJ,GJ
      real, dimension (:,:), allocatable :: AK,BK,CK,GK
      
      INTEGER I,J,K
      REAL DI1,DI2,DI3,DJ1,DJ2,DJ3,DK1,DK2,DK3
 
C-----Z-DIRECTION
!$OMP PARALLEL private(I,K,DK1,DK2,DK3)
!$OMP&private(AK,CK,BK,GK)
      allocate(AK(M1,M3),BK(M1,M3),CK(M1,M3),GK(M1,M3))
!$OMP DO
      DO 131 J=1,N2M

      DO 141 K=KBG,N3M
      DO 141 I=1,N1M

      DK1=2.*ACOEF*DENF_ZI(I,J,K)*VF_Z(I,J,KMV(K),3)*AKW(K)
      DK3=2.*ACOEF*DENF_ZI(I,J,K)*VF_Z(I,J,K,3)*CKW(K)
      DK2=-DK1-DK3

      AK(I,K)=DK1
      CK(I,K)=DK3
      BK(I,K)=1.+DK2
      GK(I,K)=RHS1(I,J,K,3)
  141 CONTINUE
      IF (IPZ .EQ. 1) THEN
       CALL TRDIAG3P(AK,BK,CK,GK,KBG,N3M,1,N1M)
      ELSE
       CALL TRDIAG3(AK,BK,CK,GK,GK,KBG,N3M,1,N1M)
      ENDIF

      DO 151 K=KBG,N3M
      DO 151 I=1,N1M
      RHS1(I,J,K,3)=GK(I,K)
  151 CONTINUE

  131 CONTINUE
!$OMP END DO
      deallocate(AK,BK,CK,GK)
!$OMP END PARALLEL

!$OMP PARALLEL private(I,J,DJ1,DJ2,DJ3,DI1,DI2,DI3)
!$OMP&private(AJ,BJ,CJ,GJ,AI,BI,CI,GI)
      allocate(AI(M2,M1),BI(M2,M1),CI(M2,M1),GI(M2,M1))
      allocate(AJ(M1,M2),BJ(M1,M2),CJ(M1,M2),GJ(M1,M2))
!$OMP DO
      DO 40 K=KBG,N3M

C-----Y-DIRECTION
      DO 91 J=1,N2M
      DO 91 I=1,N1M

!      DJ1=ACOEF*DENF_ZI(I,J,K)*VF_Z(I,JMV(J),K,2)*AJUW(J)*(1.-FIXJL(J))  !NEUMANN
!      DJ3=ACOEF*DENF_ZI(I,J,K)*VF_Z(I,J,K,2)*CJUW(J)*(1.-FIXJU(J))
      DJ1=ACOEF*DENF_ZI(I,J,K)*VF_Z(I,JMV(J),K,2)*AJUW(J) !PERIODIC
      DJ3=ACOEF*DENF_ZI(I,J,K)*VF_Z(I,J,K,2)*CJUW(J)
      DJ2=-DJ1-DJ3

      AJ(I,J)=DJ1
      CJ(I,J)=DJ3
      BJ(I,J)=1.+DJ2
      GJ(I,J)=RHS1(I,J,K,3)
   91 CONTINUE

!c     boundary condition (USE PREVIOUS STEP VALUE)
!      DO I=1,N1M
!      GJ(I,1)=GJ(I,1)-AJ(I,1)*(W(I,0,K)-WR(I,0,K))
!      GJ(I,N2M)=GJ(I,N2M)-CJ(I,N2M)*(W(I,N2,K)-WR(I,N2,K))
!      ENDDO
      IF (IPY .EQ. 1) THEN
       CALL TRDIAG2P(AJ,BJ,CJ,GJ,1,N2M,1,N1M)
      ELSE
       CALL TRDIAG2(AJ,BJ,CJ,GJ,GJ,1,N2M,1,N1M)
      ENDIF

C-----X-DIRECTION
      DO 51 J=1,N2M
      DO 51 I=1,N1M

!      DI1=ACOEF*VF_Z(IMV(I),J,K,1)*DENF_ZI(I,J,K)*AIVW(I)*(1.-FIXIL(I))  !NEUMANN
!      DI3=ACOEF*VF_Z(I,J,K,1)*DENF_ZI(I,J,K)*CIVW(I)*(1.-FIXIU(I))
      DI1=ACOEF*VF_Z(IMV(I),J,K,1)*DENF_ZI(I,J,K)*AIVW(I) !PERIODIC
      DI3=ACOEF*VF_Z(I,J,K,1)*DENF_ZI(I,J,K)*CIVW(I)
      DI2=-DI1-DI3

      AI(J,I)=DI1
      CI(J,I)=DI3
      BI(J,I)=1.+DI2
      GI(J,I)=GJ(I,J)
   51 CONTINUE

!c     boundary condition (USE PREVIOUS STEP VALUE)
!      DO J=1,N2M
!      GI(J,1)=GI(J,1)-AI(J,1)*(W(0,J,K)-WR(0,J,K))
!      GI(J,N1M)=GI(J,N1M)-CI(J,N1M)*(W(N1,J,K)-WR(N1,J,K))
!      ENDDO
      IF (IPX .EQ. 1) THEN
       CALL TRDIAG1P(AI,BI,CI,GI,1,N1M,1,N2M)
      ELSE
       CALL TRDIAG1(AI,BI,CI,GI,GI,1,N1M,1,N2M)
      ENDIF

      DO 61 J=1,N2M
      DO 61 I=1,N1M
      W(I,J,K)=GI(J,I)+W(I,J,K)
   61 CONTINUE

   40 CONTINUE
!$OMP END DO
      deallocate(AI,BI,CI,GI)
      deallocate(AJ,BJ,CJ,GJ)
!$OMP END PARALLEL
 
      RETURN
      END

C******************************************************************
      SUBROUTINE BC(RHS1,U,V,W)
C******************************************************************
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL RHS1(M1M,M2M,M3M,3)
      
      INTEGER I,J,K
      REAL UBC_AVG,Q_NP,AA
      
!       CALL BC_FLUX(Q_N,AA,U,V,W)        !CALCULATE MASS FLUX AT THE BOUNDARY AT N-STEP.

      !APPLY NEUMANN BOUNDARY CONDITION    
!BOUNDARY_CONDITION
        IF (IPX .EQ. 0) THEN  !NEUMANN CONDITION WITH INTERMEDIATE VELOCITIES.
!$OMP PARALLEL DO private(J)
          DO K=1,N3M
          DO J=1,N2M
          U(1,J,K)=U(2,J,K)
          U(N1,J,K)=U(N1M,J,K)
         ENDDO
         ENDDO
        ENDIF
        IF (IPY .EQ. 0) THEN
!$OMP PARALLEL DO private(I)
          DO K=1,N3M
          DO I=1,N1M
          V(I,1,K)=V(I,2,K)
          V(I,N2,K)=V(I,N2M,K)
         ENDDO
         ENDDO
        ENDIF
        IF (IPZ .EQ. 0) THEN
!$OMP PARALLEL DO private(I)
          DO J=1,N2M
          DO I=1,N1M
          W(I,J,1)=W(I,J,2)
          W(I,J,N3)=W(I,J,N3M)
         ENDDO
         ENDDO
        ENDIF

       CALL BC_FLUX(Q_NP,AA,U,V,W) !CALCULATE MASS FLUX AT THE BOUNDARY AT (N+1)-STEP.
       IF (AA .NE. 0.) UBC_AVG=-Q_NP/AA

      !VELOCITY CORRECTION FOR MASS CONSERVATION.
       IF (IPX .EQ. 0) THEN
!$OMP PARALLEL DO private(J)
        DO K=1,N3M
        DO J=1,N2M
         U(1,J,K)=U(1,J,K)-UBC_AVG
         U(N1,J,K)=U(N1,J,K)+UBC_AVG
        ENDDO
        ENDDO
       ENDIF
       IF (IPY .EQ. 0) THEN
!$OMP PARALLEL DO private(I)
        DO K=1,N3M
        DO I=1,N1M
         V(I,1,K)=V(I,1,K)-UBC_AVG
         V(I,N2,K)=V(I,N2,K)+UBC_AVG
        ENDDO
        ENDDO
       ENDIF
       IF (IPZ .EQ. 0) THEN
!$OMP PARALLEL DO private(I)
        DO J=1,N2M
        DO I=1,N1M
         W(I,J,1)=W(I,J,1)-UBC_AVG
         W(I,J,N3)=W(I,J,N3)+UBC_AVG
        ENDDO
        ENDDO
       ENDIF

      !BOUNDARY CONDITION IS TREATED EXPLICITLY.
      IF (IPX .EQ. 0) THEN
!$OMP PARALLEL DO private(J)
      DO K=1,N3M
      DO J=1,N2M
      RHS1(2,J,K,1)=RHS1(2,J,K,1)
     &           -ACOEF*CIU(2)*(U(2,J,K)-U(1,J,K))  !CIU HAS '-'SIGN.
      RHS1(N1M,J,K,1)=RHS1(N1M,J,K,1)
     &           -ACOEF*CIU(N1M)*(U(N1M,J,K)-U(N1,J,K))  !CIU HAS '-'SIGN.
      ENDDO
      ENDDO 
      ENDIF
      IF (IPY .EQ. 0) THEN
!$OMP PARALLEL DO private(I)
      DO K=1,N3M
      DO I=1,N1M
      RHS1(I,2,K,2)=RHS1(I,2,K,2)
     &           -ACOEF*CJV(2)*(V(I,2,K)-V(I,1,K))  !CJV HAS '-'SIGN.
      RHS1(I,N2M,K,2)=RHS1(I,N2M,K,2)
     &           -ACOEF*CJV(N2M)*(V(I,N2M,K)-V(I,N2,K))  !CJV HAS '-'SIGN.
      ENDDO
      ENDDO
      ENDIF
      IF (IPZ .EQ. 0) THEN
!$OMP PARALLEL DO private(I)
      DO J=1,N2M
      DO I=1,N1M
      RHS1(I,J,2,3)=RHS1(I,J,2,3)
     &           -ACOEF*CKW(2)*(W(I,J,2)-W(I,J,1))  !CKW HAS '-'SIGN.
      RHS1(I,J,N3M,3)=RHS1(I,J,N3M,3)
     &           -ACOEF*CKW(N3M)*(W(I,J,N3M)-W(I,J,N3))  !CKW HAS '-'SIGN.
      ENDDO
      ENDDO
      ENDIF

      RETURN
      END

C******************************************************************
      SUBROUTINE BC_FLUX(QQ,AA,U,V,W)
C******************************************************************
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE IBM_VAR
      
      IMPLICIT NONE

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL QQ,AA
      
      INTEGER I,J,K
      REAL QQ1,QQ2,QQ3,Q_IBM

      !OUTWARD IS POSITIVE.
        QQ1=0.
        QQ2=0.
        QQ3=0.
        AA=0.
       IF (IPX .EQ. 0) THEN
         AA=AA+2.*YL*ZL
!$OMP PARALLEL DO private(J)
!$OMP&reduction(+:QQ1)
        DO K=1,N3M
        DO J=1,N2M
         QQ1=QQ1+(U(N1,J,K)-U(1,J,K))*SDY(J)*SDZ(K)
        ENDDO
        ENDDO
       ENDIF
       IF (IPY .EQ. 0) THEN
         AA=AA+2.*XL*ZL
!$OMP PARALLEL DO private(I)
!$OMP&reduction(+:QQ2)
        DO K=1,N3M
        DO I=1,N1M
         QQ2=QQ2+(V(I,N2,K)-V(I,1,K))*SDX(I)*SDZ(K)
        ENDDO
        ENDDO
       ENDIF
       IF (IPZ .EQ. 0) THEN
         AA=AA+2.*XL*YL
!$OMP PARALLEL DO private(I)
!$OMP&reduction(+:QQ3)
        DO J=1,N2M
        DO I=1,N1M
         QQ3=QQ3+(W(I,J,N3)-W(I,J,1))*SDX(I)*SDY(J)
        ENDDO
        ENDDO
       ENDIF

      Q_IBM=0.
      IF (IBMON.EQ.1) THEN
C-----MASS FLUX OF ACTUATION
!$OMP PARALLEL DO private(I,J)
!$OMP&reduction(+:Q_IBM)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
      Q_IBM=Q_IBM+(-QMASS(I,J,K)*SDX(I)*SDY(J)*SDZ(K))  !OUTWARD IS POSITIVE.
      ENDDO
      ENDDO
      ENDDO
      ENDIF

      QQ=QQ1+QQ2+QQ3+Q_IBM

      RETURN
      END


C*******************************************************************
      SUBROUTINE LHSINIT
C*******************************************************************
C     CALCULATE COEFFICIENTS FOR LHS
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      USE FLOW_VAR
      
      IMPLICIT NONE
      
      INTEGER IC,IM,IP,JC,JM,JP,KC,KM,KP

!$OMP PARALLEL DO private(IM)
      DO 10 IC=IBG,N1M
      IM=IMV(IC)
      AIU(IC)=-VVDX(IC)*SSDX(IM)                 !   :i-1
      CIU(IC)=-VVDX(IC)*SSDX(IC)                 !   :i+1
   10 CONTINUE

!$OMP PARALLEL DO private(IP)
      DO 20 IC=1,N1M
      IP=IPV(IC)
      AIVW(IC)=-VVDX(IC)*SSDX(IC)
      CIVW(IC)=-VVDX(IP)*SSDX(IC)
   20 CONTINUE

!$OMP PARALLEL DO private(JM)
      DO 30 JC=JBG,N2M
      JM=JMV(JC)
      AJV(JC)=-VVDY(JC)*SSDY(JM)
      CJV(JC)=-VVDY(JC)*SSDY(JC)
   30 CONTINUE

!$OMP PARALLEL DO private(JP)
      DO 40 JC=1,N2M
      JP=JPV(JC)
      AJUW(JC)=-VVDY(JC)*SSDY(JC)
      CJUW(JC)=-VVDY(JP)*SSDY(JC)
   40 CONTINUE

!$OMP PARALLEL DO private(KM)
      DO 50 KC=KBG,N3M
      KM=KMV(KC)
      AKW(KC)=-VVDZ(KC)*SSDZ(KM)
      CKW(KC)=-VVDZ(KC)*SSDZ(KC)
   50 CONTINUE

!$OMP PARALLEL DO private(KP)
      DO 60 KC=1,N3M
      KP=KPV(KC)
      AKUV(KC)=-VVDZ(KC)*SSDZ(KC)
      CKUV(KC)=-VVDZ(KP)*SSDZ(KC)
   60 CONTINUE

      RETURN
      END
      
C=======================CONTINUITY EQUATION============================C
C      CAL. TWO-PHASE CONTINUITY EQUATION                              C
C                                                                      C
C      FLUX IS CALCULATED AT REFINED LEVEL-SET GRID                    C
C      MASS-REDISTRIBUTION ALGORITHM IS IMPLEMETAION                   C
C      (D.KIM STANFORD THESIS 2011)                                    C
C                                                                      C
C                                            KIYOUNG KIM 2015.3.18     C
C=======================CONTINUITY EQUATION============================C
C =====================================================================
      SUBROUTINE CAL_CONTINUITY(ITRACKING,DF_X,DF_Y,DF_Z,DF_C
     &,VF_X,VF_Y,VF_Z,U,V,W,PSI_XN,PSI_YN,PSI_ZN,PSI_CN
     &,PSI_X,PSI_Y,PSI_Z,PSI_C)
C =====================================================================
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      USE LES_VAR
      
      IMPLICIT NONE
      
      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL PSI_XN(0:M1,0:M2,0:M3)
      REAL PSI_YN(0:M1,0:M2,0:M3)
      REAL PSI_ZN(0:M1,0:M2,0:M3)
      REAL PSI_CN(0:M1,0:M2,0:M3)
      
      REAL DF_X(0:M1,0:M2,0:M3,3),VF_X(0:M1,0:M2,0:M3,3)
      REAL DF_Y(0:M1,0:M2,0:M3,3),VF_Y(0:M1,0:M2,0:M3,3)
      REAL DF_Z(0:M1,0:M2,0:M3,3),VF_Z(0:M1,0:M2,0:M3,3)
      REAL DF_C(0:M1,0:M2,0:M3,3)
      INTEGER REGION(-1:M1+1,-1:M2+1,-1:M3+1)

      REAL PSI_X(0:M1,0:M2,0:M3),PSI_Y(0:M1,0:M2,0:M3)
     &                     ,PSI_Z(0:M1,0:M2,0:M3),PSI_C(0:M1,0:M2,0:M3)
     
      INTEGER ITRACKING
      
      INTEGER I,J,K,II,JJ,KK,I2,J2,K2
      REAL DT_CON
      
      IF ((DENM .NE. DENP) .OR. (VISM .NE. VISP) ) THEN
C=====!DEVIDE CAL. REGION FOR EFFICIENT CALCULATION
      !BUT IF INTERFACE IS EXACTLY AT CELL-INTERFACE, THIS CAN BE PROBLEM
      !BECAUSE PSI_CN(I,J,K) CAN NOT HAVE THE VALUE BETWEEN 0 AND 1.
      !SO, INITIAL FIELD TAHT EXACTLY LAID CELL-FACE SHOULD BE AVOIDED.
      !OTHERWISE, THIS ALGORITHM HAVE TO BE REVISED!

!$OMP PARALLEL DO private(I,J)
      DO K=-1,N3+1
      DO J=-1,N2+1
      DO I=-1,N1+1
      REGION(I,J,K)=0
      ENDDO
      ENDDO
      ENDDO

C-----DEVIDE REGION WITH GRADIENT
      IF ( N3M .EQ. 1 ) THEN  !2D
      K=1
!$OMP PARALLEL DO private(I,II,JJ,I2,J2)
       DO J=1,N2M
       DO I=1,N1M
        IF (   ( PSI_CN(IPV(I),J,K)-PSI_CN(I,J,K) .NE. 0.)      !THIS SHOULD BE EVALUATED FIRST THAN
     &    .OR. ( PSI_CN(I,J,K)-PSI_CN(IMV(I),J,K) .NE. 0.)  !BELOW 'ELSE IF' SENTANCE!!
     &    .OR. ( PSI_CN(I,JPV(J),K)-PSI_CN(I,J,K) .NE. 0.)
     &    .OR. ( PSI_CN(I,J,K)-PSI_CN(I,JMV(J),K) .NE. 0.) ) THEN
         DO J2=J-2,J+2
            JJ=J2
             IF (IPY .EQ. 1) THEN
              IF (JJ .LE. 0) JJ=JJ+N2M
              IF (JJ .GE. N2) JJ=JJ-N2M
             ENDIF
         DO I2=I-2,I+2
            II=I2
             IF (IPX .EQ. 1) THEN
              IF (II .LE. 0) II=II+N1M
              IF (II .GE. N1) II=II-N1M
             ENDIF
          REGION(II,JJ,K)=2
         ENDDO
         ENDDO
        ELSE IF (PSI_CN(I,J,K) .EQ. 1.) THEN 
          IF (REGION(I,J,K) .NE. 2) REGION(I,J,K)=1
       ENDIF
         ENDDO
         ENDDO
      ELSE !3D

!$OMP PARALLEL DO private(I,J,II,JJ,KK,I2,J2,K2)
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M

        IF (   ( PSI_CN(IPV(I),J,K)-PSI_CN(I,J,K) .NE. 0.)      !THIS SHOULD BE EVALUATED FIRST THAN
     &    .OR. ( PSI_CN(I,J,K)-PSI_CN(IMV(I),J,K) .NE. 0.)      !BELOW 'ELSE IF' SENTANCE!!
     &    .OR. ( PSI_CN(I,JPV(J),K)-PSI_CN(I,J,K) .NE. 0.)
     &    .OR. ( PSI_CN(I,J,K)-PSI_CN(I,JMV(J),K) .NE. 0.)
     &    .OR. ( PSI_CN(I,J,KPV(K))-PSI_CN(I,J,K) .NE. 0.)
     &    .OR. ( PSI_CN(I,J,K)-PSI_CN(I,J,KMV(K)) .NE. 0.) ) THEN
         DO K2=K-2,K+2
            KK=K2
             IF (IPZ .EQ. 1) THEN
              IF (KK .LE. 0) KK=KK+N3M
              IF (KK .GE. N3) KK=KK-N3M
             ENDIF
         DO J2=J-2,J+2
            JJ=J2
             IF (IPY .EQ. 1) THEN
              IF (JJ .LE. 0) JJ=JJ+N2M
              IF (JJ .GE. N2) JJ=JJ-N2M
             ENDIF
         DO I2=I-2,I+2
            II=I2
             IF (IPX .EQ. 1) THEN
              IF (II .LE. 0) II=II+N1M
              IF (II .GE. N1) II=II-N1M
             ENDIF
          REGION(II,JJ,KK)=2
         ENDDO
         ENDDO
         ENDDO
        ELSE IF (PSI_CN(I,J,K) .EQ. 1.) THEN 
          IF (REGION(I,J,K) .NE. 2)   REGION(I,J,K)=1
       ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDIF !      IF ( N3M .EQ. 1 ) THEN  !2D
      
      ENDIF !      IF ((DENM .NE. DENP) .OR. (VISM .NE. VISP) ) THEN

!C-----CONTINUITY SOLVER
      IF ( DENR .EQ. 1. ) THEN  !NEED NOT TO SOLVE CONTINUITY EQN.

!       IF ( M .EQ. 1 ) THEN !NEED ONLY INITIAL STEP
!$OMP PARALLEL DO private(I,J)
       DO K=0,N3
       DO J=0,N2
       DO I=0,N1
         PSI_X(I,J,K)=DENM
         PSI_Y(I,J,K)=DENM
         PSI_Z(I,J,K)=DENM
         PSI_C(I,J,K)=DENM
       ENDDO
       ENDDO
       ENDDO

!$OMP PARALLEL DO private(I,J)
       DO K=0,N3M
       DO J=0,N2M
       DO I=0,N1M
          DF_X(I,J,K,1)=DENM
          DF_X(I,J,K,2)=DENM
          DF_X(I,J,K,3)=DENM
          
          DF_Y(I,J,K,1)=DENM
          DF_Y(I,J,K,2)=DENM
          DF_Y(I,J,K,3)=DENM
          
          DF_Z(I,J,K,1)=DENM
          DF_Z(I,J,K,2)=DENM
          DF_Z(I,J,K,3)=DENM
          
          DF_C(I,J,K,1)=DENM
          DF_C(I,J,K,2)=DENM
          DF_C(I,J,K,3)=DENM
       ENDDO
       ENDDO
       ENDDO

        IF ( VISM .EQ. VISP ) THEN

!$OMP PARALLEL DO private(I,J)
        DO K=0,N3M
        DO J=0,N2M
        DO I=0,N1M
          VF_X(I,J,K,1)=VISM
          VF_X(I,J,K,2)=VISM
          VF_X(I,J,K,3)=VISM
          
          VF_Y(I,J,K,1)=VISM
          VF_Y(I,J,K,2)=VISM
          VF_Y(I,J,K,3)=VISM
          
          VF_Z(I,J,K,1)=VISM
          VF_Z(I,J,K,2)=VISM
          VF_Z(I,J,K,3)=VISM
        ENDDO
        ENDDO
        ENDDO

        ENDIF
         
!       ENDIF

      IF ( VISM .NE. VISP ) THEN
         WRITE(*,*) 'DENSITY IS SAME, BUT VISCOSITY IS DIFF. STOP!'
         STOP
!        CALL VIS_FOR SINGLE_DENSITY(REGION,VF_X,VF_Y,VF_Z) 
      ENDIF

       ELSE !denr .ne. 1

         DT_CON=0.05*DTCONST
      CALL CONTI_X(ITRACKING,REGION,DF_X,VF_X,PSI_X,DT_CON,U,V,W,PSI_XN)
      CALL CONTI_Y(ITRACKING,REGION,DF_Y,VF_Y,PSI_Y,DT_CON,U,V,W,PSI_YN)
      IF(N3M.NE.1)CALL CONTI_Z(ITRACKING,REGION,DF_Z,VF_Z,PSI_Z,DT_CON
     &,U,V,W,PSI_ZN)
       CALL CONTI_C(ITRACKING,REGION,DF_C,PSI_C,DT_CON,U,V,W
     &,PSI_XN,PSI_YN,PSI_ZN,PSI_CN) !FOR MULTIGRID AT POISSON
      ENDIF

       IF (ILES .EQ. 1) THEN
!$OMP PARALLEL DO private(I,J)
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M
          VF_X(I,J,K,1)=VF_X(I,J,K,1)+TNU(I,J,K)
          VF_X(I,J,K,2)=VF_X(I,J,K,2)+TNU(I,J,K)
          VF_X(I,J,K,3)=VF_X(I,J,K,3)+TNU(I,J,K)

          VF_Y(I,J,K,1)=VF_Y(I,J,K,1)+TNU(I,J,K)
          VF_Y(I,J,K,2)=VF_Y(I,J,K,2)+TNU(I,J,K)
          VF_Y(I,J,K,3)=VF_Y(I,J,K,3)+TNU(I,J,K)

          VF_Z(I,J,K,1)=VF_Z(I,J,K,1)+TNU(I,J,K)
          VF_Z(I,J,K,2)=VF_Z(I,J,K,2)+TNU(I,J,K)
          VF_Z(I,J,K,3)=VF_Z(I,J,K,3)+TNU(I,J,K)
        ENDDO
        ENDDO
        ENDDO
       ENDIF
         
      RETURN
      END

C =====================================================================
      SUBROUTINE VIS_FOR SINGLE_DENSITY(REGION,VF_X,VF_Y,VF_Z) 
C =====================================================================
!      USE TWO_PHASE_PROPERTY
!
!      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
     
!      INTEGER REGION(-1:M1+1,-1:M2+1,-1:M3+1)
!      REAL DF_X(0:M1,0:M2,0:M3,3),VF_X(0:M1,0:M2,0:M3,3)
!      REAL DF_Y(0:M1,0:M2,0:M3,3),VF_Y(0:M1,0:M2,0:M3,3)
!      REAL DF_Z(0:M1,0:M2,0:M3,3),VF_Z(0:M1,0:M2,0:M3,3)
!
!C-----X-DIR
!!$OMP PARALLEL DO private(I,J,UU,VV,WW,PSI_SUM)
!       DO K=KBG,N3M
!       DO J=JBG,N2M
!       DO I=IBG,N1M
!
!       UU=0.5*(U(I+1,J,K)+U(I,J,K))
!       VV=0.5*(SDX(I-1)*V(I,J+1,K)+SDX(I)*V(I-1,J+1,K))*VVDX(I)
!       WW=0.5*(SDX(I-1)*W(I,J,K+1)+SDX(I)*W(I-1,J,K+1))*VVDX(I)
!
!       CALL PSIFACE_X(I,J,K,UU,PSI_SUM,1)
!       VF_X(I,J,K,1)=VISM+(VISP-VISM)*PSI_SUM
!       CALL PSIFACE_X(I,J,K,VV,PSI_SUM,2)
!       VF_X(I,J,K,2)=VISM+(VISP-VISM)*PSI_SUM
!       CALL PSIFACE_X(I,J,K,WW,PSI_SUM,3)
!       VF_X(I,J,K,3)=VISM+(VISP-VISM)*PSI_SUM
!
!      ENDDO
!      ENDDO
!      ENDDO
!
!      I=IBG-1
!!$OMP PARALLEL DO private(J,UU,PSI_SUM)
!       DO K=KBG,N3M
!       DO J=JBG,N2M
!       UU=0.5*(U(I+1,J,K)+U(I,J,K))
!       CALL PSIFACE_X(I,J,K,UU,PSI_SUM,1)
!       VF_X(I,J,K,1)=VISM+(VISP-VISM)*PSI_SUM
!      ENDDO
!      ENDDO
!      
!      J=0
!!$OMP PARALLEL DO private(I,VV,PSI_SUM)
!       DO K=KBG,N3M
!       DO I=IBG,N1M
!       VV=0.5*(SDX(I-1)*V(I,J+1,K)+SDX(I)*V(I-1,J+1,K))*VVDX(I)
!       CALL PSIFACE_X(I,J,K,VV,PSI_SUM,2)
!       VF_X(I,J,K,2)=VISM+(VISP-VISM)*PSI_SUM
!       ENDDO
!       ENDDO
!      
!      K=KMV(1)
!!$OMP PARALLEL DO private(I,WW,PSI_SUM)
!       DO J=JBG,N2M
!       DO I=IBG,N1M
!       WW=0.5*(SDX(I-1)*W(I,J,K+1)+SDX(I)*W(I-1,J,K+1))*VVDX(I)
!       CALL PSIFACE_X(I,J,K,WW,PSI_SUM,3)
!       VF_X(I,J,K,3)=VISM+(VISP-VISM)*PSI_SUM
!       ENDDO
!       ENDDO  
!      
!C-----Y-DIR
!!$OMP PARALLEL DO private(I,J,UU,VV,WW,PSI_SUM)
!       DO K=KBG,N3M
!       DO J=JBG,N2M
!       DO I=IBG,N1M
!      UU=0.5*(SDY(J-1)*U(I+1,J,K)+SDY(J)*U(I+1,J-1,K))*VVDY(J)
!      VV=0.5*(V(I,J+1,K)+V(I,J,K))
!      WW=0.5*(SDY(J-1)*W(I,J,K+1)+SDY(J)*W(I,J-1,K+1))*VVDY(J)
!       CALL PSIFACE_Y(I,J,K,UU,PSI_SUM,1)
!       VF_Y(I,J,K,1)=VISM+(VISP-VISM)*PSI_SUM
!       CALL PSIFACE_Y(I,J,K,VV,PSI_SUM,2)
!       VF_Y(I,J,K,2)=VISM+(VISP-VISM)*PSI_SUM
!       CALL PSIFACE_Y(I,J,K,WW,PSI_SUM,3)
!       VF_Y(I,J,K,3)=VISM+(VISP-VISM)*PSI_SUM
!      ENDDO
!      ENDDO
!      ENDDO
!      
!      I=0.
!!$OMP PARALLEL DO private(J,UU,PSI_SUM)
!       DO K=KBG,N3M
!       DO J=JBG,N2M
!       UU=0.5*(SDY(J-1)*U(I+1,J,K)+SDY(J)*U(I+1,J-1,K))*VVDY(J)
!       CALL PSIFACE_Y(I,J,K,UU,PSI_SUM,1)
!       VF_Y(I,J,K,1)=VISM+(VISP-VISM)*PSI_SUM
!      ENDDO
!      ENDDO
!      
!      J=1
!!$OMP PARALLEL DO private(I,VV,PSI_SUM)
!       DO K=KBG,N3M
!       DO I=IBG,N1M
!       VV=0.5*(V(I,J+1,K)+V(I,J,K))
!       CALL PSIFACE_Y(I,J,K,VV,PSI_SUM,2)
!       VF_Y(I,J,K,2)=VISM+(VISP-VISM)*PSI_SUM
!      ENDDO
!      ENDDO
!      
!      K=KMV(1)
!!$OMP PARALLEL DO private(I,WW,PSI_SUM)
!       DO J=JBG,N2M
!       DO I=IBG,N1M
!       WW=0.5*(SDY(J-1)*W(I,J,K+1)+SDY(J)*W(I,J-1,K+1))*VVDY(J)
!       CALL PSIFACE_Y(I,J,K,WW,PSI_SUM,3)
!       VF_Y(I,J,K,3)=VISM+(VISP-VISM)*PSI_SUM
!      ENDDO
!      ENDDO  
!      
!C-----Z-DIR
!!$OMP PARALLEL DO private(I,J,UU,VV,WW,PSI_SUM)
!       DO K=KBG,N3M
!       DO J=JBG,N2M
!       DO I=IBG,N1M
!      UU=0.5*(SDZ(K-1)*U(I+1,J,K)+SDZ(K)*U(I+1,J,K-1))*VVDZ(K)
!      VV=0.5*(SDZ(K-1)*V(I,J+1,K)+SDZ(K)*V(I,J+1,K-1))*VVDZ(K)
!      WW=0.5*(W(I,J,K+1)+W(I,J,K))
!
!       CALL PSIFACE_Z(I,J,K,UU,PSI_SUM,1)
!       VF_Z(I,J,K,1)=VISM+(VISP-VISM)*PSI_SUM
!       CALL PSIFACE_Z(I,J,K,VV,PSI_SUM,2)
!       VF_Z(I,J,K,2)=VISM+(VISP-VISM)*PSI_SUM
!       CALL PSIFACE_Z(I,J,K,WW,PSI_SUM,3)
!       VF_Z(I,J,K,3)=VISM+(VISP-VISM)*PSI_SUM
!      ENDDO
!      ENDDO
!      ENDDO
!
!      I=0
!!$OMP PARALLEL DO private(J,UU,PSI_SUM)
!       DO K=KBG,N3M
!       DO J=JBG,N2M
!      UU=0.5*(SDZ(K-1)*U(I+1,J,K)+SDZ(K)*U(I+1,J,K-1))*VVDZ(K)
!       CALL PSIFACE_Z(I,J,K,UU,PSI_SUM,1)
!       VF_Z(I,J,K,1)=VISM+(VISP-VISM)*PSI_SUM
!      ENDDO
!      ENDDO
!
!      J=0
!!$OMP PARALLEL DO private(I,VV,PSI_SUM)
!       DO K=KBG,N3M
!       DO I=IBG,N1M
!      VV=0.5*(SDZ(K-1)*V(I,J+1,K)+SDZ(K)*V(I,J+1,K-1))*VVDZ(K)
!       CALL PSIFACE_Z(I,J,K,VV,PSI_SUM,2)
!       VF_Z(I,J,K,2)=VISM+(VISP-VISM)*PSI_SUM
!      ENDDO
!      ENDDO
!
!      K=KMV(KBG)
!!$OMP PARALLEL DO private(I,WW,PSI_SUM)
!       DO J=JBG,N2M
!       DO I=IBG,N1M
!      WW=0.5*(W(I,J,K+1)+W(I,J,K))
!       CALL PSIFACE_Z(I,J,K,WW,PSI_SUM,3)
!       VF_Z(I,J,K,3)=VISM+(VISP-VISM)*PSI_SUM
!      ENDDO
!      ENDDO  

      RETURN
      END

C =====================================================================
      SUBROUTINE CONTI_X(ITRACKING,REGION,DF_X,VF_X,PSI_X,DT_CON,U,V,W
     &,PSI_XN)
C =====================================================================
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY
      
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE
      
      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL PSI_XN(0:M1,0:M2,0:M3)

      REAL DFX(0:M1,0:M2,0:M3),DFY(0:M1,0:M2,0:M3),DFZ(0:M1,0:M2,0:M3)
      REAL EPS(0:M1,0:M2,0:M3),EPS_OLD(0:M1,0:M2,0:M3)
      INTEGER PSI_ITERATION,STEADY_EPS
      
      INTEGER REGION(-1:M1+1,-1:M2+1,-1:M3+1)
      REAL DF_X(0:M1,0:M2,0:M3,3),VF_X(0:M1,0:M2,0:M3,3)
      REAL PSI_X(0:M1,0:M2,0:M3)
      REAL TAUC(3,2)
      
      INTEGER ITRACKING
      REAL DT_CON
      
      INTEGER I,J,K,IP1,IM1,JP1,JM1,KP1,KM1
      REAL UE,UW,VN,VS,WT,WB
      REAL PSI_SUM,FLUX_X,FLUX_Y,FLUX_Z,DEN_XN,DEN_X
      
!$OMP PARALLEL DO private(I,J)
       DO K=0,N3
       DO J=0,N2
       DO I=0,N1
        DF_X(I,J,K,1)=0.    !INITIALIZE for  IF ( DF_X(I,J,K,1) .EQ. 0. ) THEN
        DF_X(I,J,K,2)=0.
        DF_X(I,J,K,3)=0.
       ENDDO
       ENDDO
       ENDDO

      !DEFINE PSI_X (defined at cell-face)
      !MODIFIED VERSION OF
      !'A LEVEL SET BASED METHOD FOR CALCULATING FLUX DENSITIES IN
      !TWO-PHASE FLOWS - Raessi(CTR,2008)'
      !SOLVE WHOLE DOMAIN FOR MASS CONSERCATION.
      !IF SOLVE IN THE REGION 2, MASS CAN BE NOT CONSERVED BECAUSE OF NUMERICAL ERROR.
!$OMP PARALLEL DO private(I,J,IP1,IM1,JP1,JM1,KP1,KM1,UE,UW,VN,VS,WT,WB)
!$OMP&private(PSI_SUM,FLUX_X,FLUX_Y,FLUX_Z,DEN_XN,DEN_X)
       DO K=1,N3M
        KP1=KPV(K)
        KM1=KMV(K)
       DO J=1,N2M
         JP1=JPV(J)
         JM1=JMV(J)
       DO I=IBG,N1M
         IP1=IPV(I)
         IM1=IMV(I)

      IF ( REGION(I,J,K) .EQ. 0 ) THEN
       VF_X(I,J,K,1)=VISM
       VF_X(IM1,J,K,1)=VISM
       VF_X(I,J,K,2)=VISM
       VF_X(I,JM1,K,2)=VISM
       VF_X(I,J,K,3)=VISM
       VF_X(I,J,KM1,3)=VISM

       DF_X(I,J,K,1)=DENM
       DF_X(IM1,J,K,1)=DENM
       DF_X(I,J,K,2)=DENM
       DF_X(I,JM1,K,2)=DENM
       DF_X(I,J,K,3)=DENM
       DF_X(I,J,KM1,3)=DENM

       PSI_X(I,J,K)=PSI_XN(I,J,K)

      ELSE IF (REGION(I,J,K) .EQ. 1 ) THEN
       VF_X(I,J,K,1)=VISP
       VF_X(IM1,J,K,1)=VISP
       VF_X(I,J,K,2)=VISP
       VF_X(I,JM1,K,2)=VISP
       VF_X(I,J,K,3)=VISP
       VF_X(I,J,KM1,3)=VISP
       
       DF_X(I,J,K,1)=DENP
       DF_X(IM1,J,K,1)=DENP
       DF_X(I,J,K,2)=DENP
       DF_X(I,JM1,K,2)=DENP
       DF_X(I,J,K,3)=DENP
       DF_X(I,J,KM1,3)=DENP

       PSI_X(I,J,K)=PSI_XN(I,J,K)
        
      ELSE

       UE=0.5*(U(IP1,J,K)+U(I,J,K))
       UW=0.5*(U(I,J,K)+U(IM1,J,K))
       VN=0.5*(SDX(IM1)*V(I,JP1,K)+SDX(I)*V(IM1,JP1,K))*VVDX(I)
       VS=0.5*(SDX(IM1)*V(I,J,K)+SDX(I)*V(IM1,J,K))*VVDX(I)
       WT=0.5*(SDX(IM1)*W(I,J,KP1)+SDX(I)*W(IM1,J,KP1))*VVDX(I)
       WB=0.5*(SDX(IM1)*W(I,J,K)+SDX(I)*W(IM1,J,K))*VVDX(I)
       
!      PSI1=0.5*(PSI_XN(I,J,K)+PSI_XN(IP1,J,K))
!      PSI2=0.5*(PSI_XN(I,J,K)+PSI_XN(IM1,J,K))
!      PSI3=0.5*(PSI_XN(I,J,K)*SDY(JP1)+PSI_XN(I,JP1,K)*SDY(J))*VVDY(JP1)
!      PSI4=0.5*(PSI_XN(I,JM1,K)*SDY(J)+PSI_XN(I,J,K)*SDY(JM1))*VVDY(J)
!      PSI5=0.5*(PSI_XN(I,J,K)*SDZ(KP1)+PSI_XN(I,J,KP1)*SDZ(K))*VVDZ(KP1)
!      PSI6=0.5*(PSI_XN(I,J,KM1)*SDZ(K)+PSI_XN(I,J,K)*SDZ(KM1))*VVDZ(K)
!
!       DF_X(I,J,K,1)=DENM+DEN_DIFF*PSI1
!       DF_X(IM1,J,K,1)=DENM+DEN_DIFF*PSI2
!       DF_X(I,J,K,2)=DENM+DEN_DIFF*PSI3
!       DF_X(I,JM1,K,2)=DENM+DEN_DIFF*PSI4
!       DF_X(I,J,K,3)=DENM+DEN_DIFF*PSI5
!       DF_X(I,J,KM1,3)=DENM+DEN_DIFF*PSI6
!
!       VF_X(I,J,K,1)=VISM+VIS_DIFF*PSI1
!       VF_X(IM1,J,K,1)=VISM+VIS_DIFF*PSI2
!       VF_X(I,J,K,2)=VISM+VIS_DIFF*PSI3
!       VF_X(I,JM1,K,2)=VISM+VIS_DIFF*PSI4
!       VF_X(I,J,K,3)=VISM+VIS_DIFF*PSI5
!       VF_X(I,J,KM1,3)=VISM+VIS_DIFF*PSI6

      IF ( DF_X(I,J,K,1) .EQ. 0. ) THEN
       IF (ABS(UE) .LE. 1.E-8) THEN
         PSI_SUM=0.5*(PSI_XN(I,J,K)+PSI_XN(IP1,J,K))
       ELSE
        CALL PSIFACE_X(I,J,K,UE,PSI_SUM,1)
       ENDIF
        VF_X(I,J,K,1)=VISM+(VISP-VISM)*PSI_SUM
        DF_X(I,J,K,1)=DENM+DEN_DIFF*PSI_SUM
      ENDIF
      IF ( DF_X(IM1,J,K,1) .EQ. 0. ) THEN
       IF (ABS(UW) .LE. 1.E-8) THEN
         PSI_SUM=0.5*(PSI_XN(I,J,K)+PSI_XN(IM1,J,K))
       ELSE
        CALL PSIFACE_X(IM1,J,K,UW,PSI_SUM,1)
       ENDIF
        VF_X(IM1,J,K,1)=VISM+(VISP-VISM)*PSI_SUM
        DF_X(IM1,J,K,1)=DENM+DEN_DIFF*PSI_SUM
      ENDIF
      IF ( DF_X(I,J,K,2) .EQ. 0. ) THEN
       IF (ABS(VN) .LE. 1.E-8) THEN
         PSI_SUM=0.5*(PSI_XN(I,J,K)*SDY(JP1)
     &                             +PSI_XN(I,JP1,K)*SDY(J))*VVDY(JP1)
       ELSE
        CALL PSIFACE_X(I,J,K,VN,PSI_SUM,2)
       ENDIF
        VF_X(I,J,K,2)=VISM+(VISP-VISM)*PSI_SUM
        DF_X(I,J,K,2)=DENM+DEN_DIFF*PSI_SUM
       ENDIF
      IF ( DF_X(I,JM1,K,2) .EQ. 0. ) THEN
       IF (ABS(VS) .LE. 1.E-8) THEN
         PSI_SUM=0.5*(PSI_XN(I,JM1,K)*SDY(J)
     &                             +PSI_XN(I,J,K)*SDY(JM1))*VVDY(J)
       ELSE
        CALL PSIFACE_X(I,JM1,K,VS,PSI_SUM,2)
       ENDIF
        VF_X(I,JM1,K,2)=VISM+(VISP-VISM)*PSI_SUM
        DF_X(I,JM1,K,2)=DENM+DEN_DIFF*PSI_SUM
      ENDIF
      IF ( DF_X(I,J,K,3) .EQ. 0. ) THEN
       IF (ABS(WT) .LE. 1.E-8) THEN
         PSI_SUM=0.5*(PSI_XN(I,J,K)*SDZ(KP1)
     &                             +PSI_XN(I,J,KP1)*SDZ(K))*VVDZ(KP1)
       ELSE
        CALL PSIFACE_X(I,J,K,WT,PSI_SUM,3)
       ENDIF
        VF_X(I,J,K,3)=VISM+(VISP-VISM)*PSI_SUM
        DF_X(I,J,K,3)=DENM+DEN_DIFF*PSI_SUM
      ENDIF
      IF ( DF_X(I,J,KM1,3) .EQ. 0. ) THEN
       IF (ABS(WB) .LE. 1.E-8) THEN
         PSI_SUM=0.5*(PSI_XN(I,J,KM1)*SDZ(K)
     &                             +PSI_XN(I,J,K)*SDZ(KM1))*VVDZ(K)
       ELSE
       CALL PSIFACE_X(I,J,KM1,WB,PSI_SUM,3)        !K-1 NOT KMV(K)
       ENDIF
        VF_X(I,J,KM1,3)=VISM+(VISP-VISM)*PSI_SUM
        DF_X(I,J,KM1,3)=DENM+DEN_DIFF*PSI_SUM 
      ENDIF
            
        FLUX_X=( DF_X(I,J,K,1)*UE-DF_X(IM1,J,K,1)*UW )*VVDX(I)
        FLUX_Y=( DF_X(I,J,K,2)*VN-DF_X(I,JM1,K,2)*VS )*SSDY(J)
        FLUX_Z=( DF_X(I,J,K,3)*WT-DF_X(I,J,KM1,3)*WB )*SSDZ(K)

       DEN_XN=DENM+DEN_DIFF*PSI_XN(I,J,K)
       DEN_X=DEN_XN-DTCONST*( FLUX_X+FLUX_Y+FLUX_Z )
       PSI_X(I,J,K)=(DEN_X-DENM)/DEN_DIFF
       
      ENDIF
       ENDDO
       ENDDO
       ENDDO
       
       !INSTEAD OF MASS REDISTRIBUTION, SIMPLY REMOVE THE ERROR.
       !SHOULD REVISITE HERE.
!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M
       DO J=1,N2M
       DO I=IBG,N1M
         IF (PSI_X(I,J,K) .GT. 1.) THEN
           PSI_X(I,J,K)=1.
         ELSE IF (PSI_X(I,J,K) .LT. 0.) THEN
           PSI_X(I,J,K)=0.
         ENDIF
       ENDDO
       ENDDO
       ENDDO
       
C-----CONTINUITY BOUNDARY CONDITION : NEUMAN+PERIODIC CONDITION.
      IF ( IPX .NE. 1 ) THEN
!$OMP PARALLEL DO private(J)
       DO K=1,N3M
       DO J=1,N2M
       PSI_X(1,J,K)=PSI_X(2,J,K)
       PSI_X(N1,J,K)=PSI_X(N1M,J,K)
       ENDDO
       ENDDO
      ENDIF
      IF ( IPY .NE. 1 ) THEN
!$OMP PARALLEL DO private(I)
       DO K=1,N3M
       DO I=IBG,N1M
       PSI_X(I,0,K)=PSI_X(I,1,K)
       PSI_X(I,N2,K)=PSI_X(I,N2M,K)
       ENDDO
       ENDDO
      ENDIF

      IF ( IPZ .NE. 1 ) THEN
!$OMP PARALLEL DO private(I)
       DO J=1,N2M
       DO I=IBG,N1M
       PSI_X(I,J,0)=PSI_X(I,J,1)
       PSI_X(I,J,N3)=PSI_X(I,J,N3M)
      ENDDO
      ENDDO
      ENDIF
C-----CONTINUITY BOUNDARY CONDITION : NEUMAN+PERIODIC CONDITION.

!        CRI1=1.E-4
!        CRI2=1.-CRI1
!c-----REDISTRIBUTION ALGORITHM
!       DO 5000 PSI_ITERATION=1,20
!        EPSMAX=0.
!!$OMP PARALLEL DO private(I,J)
!       DO K=0,N3
!       DO J=0,N2
!       DO I=0,N1
!        EPS(I,J,K)=0.
!       ENDDO
!       ENDDO
!       ENDDO
!
!      !A MASS-CONSERVING LEVEL-SET METHOD FOR MODELLING OF MULTI-PHASE 
!      !FLOWSPIJL ET AL. (INTERNATIONAL JOURNAL FOR NUMERICAL METHODS 
!      !IN FLUIDS,2004)
!      !(NUMERICAL ERRORS ARE IN GENERAL OF THE ORDER 1.E-4)
!
!      !DEFINE EPS (defined at cell-face)
!      
!!$OMP PARALLEL DO private(I,J)
!       DO K=1,N3M
!       DO J=1,N2M
!       DO I=IBG,N1M
!      IF((PSI_X(I,J,K) .GT. 0.).AND.(PSI_X(I,J,K) .LE. CRI1)) THEN
!        EPS(I,J,K)=-PSI_X(I,J,K)
!      ELSE IF((PSI_X(I,J,K) .GE. CRI2).AND.(PSI_X(I,J,K) .LT. 1.))THEN
!        EPS(I,J,K)=1.-PSI_X(I,J,K)
!      ELSE IF(PSI_X(I,J,K) .GT. 1.) THEN
!         EPS(I,J,K)=1.-PSI_X(I,J,K)
!      ELSE IF(PSI_X(I,J,K) .LT. 0.) THEN
!         EPS(I,J,K)=-PSI_X(I,J,K)
!      ENDIF
!
!       ENDDO
!       ENDDO
!       ENDDO
!       
!!$OMP PARALLEL DO private(I,J)
!       DO K=1,N3M
!       DO J=1,N2M
!       DO I=IBG,N1M
!        IF ( ABS(EPS(I,J,K)) .GT. 1.E-8 ) THEN
!          EPS_OLD(I,J,K)=EPS(I,J,K)
!        ENDIF
!       ENDDO
!       ENDDO
!       ENDDO
!
!!$OMP PARALLEL DO private(I,J)
!       DO K=0,N3
!       DO J=0,N2
!       DO I=0,N1
!       DFX(I,J,K)=0.
!       DFY(I,J,K)=0.
!       DFZ(I,J,K)=0.
!       ENDDO
!       ENDDO
!       ENDDO
!       
!       !define the unit velocity toward phase-interface
!!$OMP PARALLEL DO private(I,J,JP1,JM1,KP1,KM1,AA,BB,CC,AABBCC,AABBCCI)
!!$OMP&private(PSI_XE,PSI_XW,PSI_XM,PSI_XS,PSI_XT,PSI_XB)
!       DO K=1,N3M
!       DO J=1,N2M
!       DO I=IBG,N1M  
!           
!        IF ( ABS(EPS(I,J,K)) .GT. 1.E-8 ) THEN
!      PSI_XE=0.5*(PSI_X(IPV(I),J,K)+PSI_X(I,J,K))
!      PSI_XW=0.5*(PSI_X(I,J,K)+PSI_X(IMV(I),J,K))
!        AA=( PSI_XE-PSI_XW )*VVDX(I)
!
!           JP1=JPV(J)
!           JM1=JMV(J)
!      PSI_Xm=0.5*(PSI_X(I,JP1,K)*SDY(J)+PSI_X(I,J,K)*SDY(JP1))*VVDY(JP1)
!      PSI_XS=0.5*(PSI_X(I,J,K)*SDY(JM1)+PSI_X(I,JM1,K)*SDY(J))*VVDY(J)
!        BB=( PSI_Xm-PSI_XS )*SSDY(J)  !PSI_Xm BECAUSE OF PSI_XN
!
!           KP1=KPV(K)
!           KM1=KMV(K)
!      PSI_XT=0.5*(PSI_X(I,J,KP1)*SDZ(K)+PSI_X(I,J,K)*SDZ(KP1))*VVDZ(KP1)
!      PSI_XB=0.5*(PSI_X(I,J,K)*SDZ(KM1)+PSI_X(I,J,KM1)*SDZ(K))*VVDZ(K)
!        CC=( PSI_XT-PSI_XB )*SSDZ(K)
!
!        AABBCC=SQRT(AA**2+BB**2+CC**2)
!        IF (AABBCC .GE. 1.E-4 ) THEN
!          AABBCCI=1./AABBCC
!         IF( PSI_X(I,J,K) .GE. 0.5) THEN
!          DFX(I,J,K)=-AA*AABBCCI
!          DFY(I,J,K)=-BB*AABBCCI
!          DFZ(I,J,K)=-CC*AABBCCI
!         ELSE
!          DFX(I,J,K)=AA*AABBCCI
!          DFY(I,J,K)=BB*AABBCCI
!          DFZ(I,J,K)=CC*AABBCCI
!         ENDIF
!        ENDIF
!       ENDIF
!        
!       ENDDO
!       ENDDO
!       ENDDO
!
!      IF ( IPX .NE. 1 ) THEN
!!$OMP PARALLEL DO private(J)
!       DO K=1,N3M
!       DO J=1,N2M
!       DFX(1,J,K)=0.!DFX(2,J,K)    !NO FLUX CROSS THE WALL. IS THIS RIGHT?
!       DFX(N1,J,K)=0.!DFX(N1M,J,K)  !THEN HERE SHOULD BE ZERO.
!      ENDDO
!      ENDDO
!      ENDIF
!      IF ( IPY .NE. 1 ) THEN
!!$OMP PARALLEL DO private(I)
!       DO K=1,N3M
!       DO I=IBG,N1M
!        DFY(I,0,K)=0.!DFY(I,1,K)
!        DFY(I,N2,K)=0.!DFY(I,N2M,K)
!      ENDDO
!      ENDDO
!      ENDIF
!      IF ( IPZ .NE. 1 ) THEN
!!$OMP PARALLEL DO private(I)
!       DO J=1,N2M
!       DO I=IBG,N1M
!       DFZ(I,J,0)=0.!DFZ(I,J,1)
!       DFZ(I,J,N3)=0.!DFZ(I,J,N3M)
!      ENDDO
!      ENDDO
!      ENDIF
!
!      DO 2500 STEADY_EPS=1,5  !FIND STEADY SOL.
!
!        EPS_DIF=0.
!       !cal transport eqn.
!!$OMP PARALLEL DO private(I,J,JP1,JM1,KP1,KM1)
!!$OMP&private(DFXE,DFXW,DFYN,DFYS,DFZT,DFZB,EPSX,EPSY,EPSZ)
!!$OMP&reduction(MAX:EPS_DIF)
!       DO K=1,N3M
!       DO J=1,N2M
!       DO I=IBG,N1M
!           
!        IF ( ABS(EPS(I,J,K)) .GT. 1.E-8 ) THEN
!       DFXE=0.5*(EPS(IPV(I),J,K)+EPS(I,J,K))
!       DFXW=0.5*(EPS(I,J,K)+EPS(IMV(I),J,K))
!        EPSX=( DFX(I,J,K)*DFXE-DFX(IMV(I),J,K)*DFXW )*VVDX(I)
!
!           JP1=JPV(J)
!           JM1=JMV(J)
!       DFYN=0.5*(EPS(I,JP1,K)*SDY(J)+EPS(I,J,K)*SDY(JP1))*VVDY(JP1)
!       DFYS=0.5*(EPS(I,J,K)*SDY(JM1)+EPS(I,JM1,K)*SDY(J))*VVDY(J)
!        EPSY=( DFY(I,J,K)*DFYN-DFY(I,JM1,K)*DFYS )*SSDY(J)
!
!           KP1=KPV(K)
!           KM1=KMV(K)
!       DFZT=0.5*(EPS(I,J,KP1)*SDZ(K)+EPS(I,J,K)*SDZ(KP1))*VVDZ(KP1)
!       DFZB=0.5*(EPS(I,J,K)*SDZ(KM1)+EPS(I,J,KM1)*SDZ(K))*VVDZ(K)
!        EPSZ=( DFZ(I,J,K)*DFZT-DFZ(I,J,KM1)*DFZB )*SSDZ(K)
!
!        EPS(I,J,K)=EPS_OLD(I,J,K)-DT_CON*( EPSX+EPSY+EPSZ )
!        EPS_DIF=AMAX1(EPS_DIF,ABS(EPS(I,J,K)-EPS_OLD(I,J,K)))
!       ENDIF
!
!       ENDDO
!       ENDDO
!       ENDDO
!       
!!$OMP PARALLEL DO private(I,J)
!       DO K=1,N3M
!       DO J=1,N2M
!       DO I=IBG,N1M
!        IF ( ABS(EPS(I,J,K)) .GT. 1.E-8 ) THEN
!          EPS_OLD(I,J,K)=EPS(I,J,K)
!        ENDIF
!       ENDDO
!       ENDDO
!       ENDDO
!       
!!         write(*,124) PSI_ITERATION,STEADY_EPS,EPS_DIF
!! 124  format(2I4,1ES12.5)
!          IF ( EPS_DIF .LE. 1.E-8 ) GOTO 100
!
! 2500 ENDDO
!c------------------UPDATE PSI_X
!
! 100   CONTINUE
!!$OMP PARALLEL DO private(I,J)
!       DO K=1,N3M
!       DO J=1,N2M
!       DO I=IBG,N1M
!        IF ( ABS(EPS(I,J,K)) .GT. 1.E-8 ) THEN
!         PSI_X(I,J,K)=PSI_X(I,J,K)+EPS(I,J,K)
!         EPSMAX=MAX(EPSMAX,ABS(EPS(I,J,K)))
!        ENDIF
!       ENDDO
!       ENDDO
!       ENDDO
!
!C-----CONTINUITY BOUNDARY CONDITION : NEUMAN+PERIODIC CONDITION.
!      IF ( IPX .NE. 1 ) THEN
!!$OMP PARALLEL DO private(J)
!       DO K=1,N3M
!       DO J=1,N2M
!       PSI_X(1,J,K)=PSI_X(2,J,K)
!       PSI_X(N1,J,K)=PSI_X(N1M,J,K)
!      ENDDO
!      ENDDO
!      ENDIF
!      IF ( IPY .NE. 1 ) THEN
!!$OMP PARALLEL DO private(I)
!       DO K=1,N3M
!       DO I=IBG,N1M
!       PSI_X(I,0,K)=PSI_X(I,1,K)
!       PSI_X(I,N2,K)=PSI_X(I,N2M,K)
!      ENDDO
!      ENDDO
!      ENDIF
!      IF ( IPZ .NE. 1 ) THEN
!!$OMP PARALLEL DO private(I)
!       DO J=1,N2M
!       DO I=IBG,N1M
!       PSI_X(I,J,0)=PSI_X(I,J,1)
!       PSI_X(I,J,N3)=PSI_X(I,J,N3M)
!      ENDDO
!      ENDDO
!      ENDIF
!C-----CONTINUITY BOUNDARY CONDITION : NEUMAN+PERIODIC CONDITION.
!
!!       write(*,123) PSI_ITERATION,STEADY_EPS,EPSMAX
!! 123  format('eps_MAX ',2I4,1ES12.5)
!         IF ( EPSMAX .LE. 1.E-8) GOTO 101
! 5000  ENDDO
!
!       WRITE(*,*) 'X-CONTINUITY ITERATION LIMIT EXCEEDED!!!'
! 101   CONTINUE
!       
!!       write(*,124) PSI_ITERATION,EPSMAX
!! 124  format('eps_MAX ',I4,1ES12.5)

       RETURN
       END

C =====================================================================
      SUBROUTINE CONTI_Y(ITRACKING,REGION,DF_Y,VF_Y,PSI_Y,DT_CON,U,V,W
     &,PSI_YN)
C =====================================================================
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL PSI_YN(0:M1,0:M2,0:M3)

      REAL DFX(0:M1,0:M2,0:M3),DFY(0:M1,0:M2,0:M3),DFZ(0:M1,0:M2,0:M3)
      REAL EPS(0:M1,0:M2,0:M3),EPS_OLD(0:M1,0:M2,0:M3)
      INTEGER PSI_ITERATION,STEADY_EPS

      INTEGER REGION(-1:M1+1,-1:M2+1,-1:M3+1)
      REAL DF_Y(0:M1,0:M2,0:M3,3),VF_Y(0:M1,0:M2,0:M3,3)
      REAL PSI_Y(0:M1,0:M2,0:M3)
      REAL TAUC(3,2)

      INTEGER ITRACKING
      REAL DT_CON
      
      INTEGER I,J,K,IP1,IM1,JP1,JM1,KP1,KM1
      REAL UE,UW,VN,VS,WT,WB
      REAL PSI_SUM,FLUX_X,FLUX_Y,FLUX_Z,DEN_YN,DEN_Y
      
!$OMP PARALLEL DO private(I,J)
       DO K=0,N3
       DO J=0,N2
       DO I=0,N1
        DF_Y(I,J,K,1)=0.
        DF_Y(I,J,K,2)=0.
        DF_Y(I,J,K,3)=0.
       ENDDO
       ENDDO
       ENDDO

      !DEFINE PSI_Y
!$OMP PARALLEL DO private(I,J,IP1,IM1,JP1,JM1,KP1,KM1,UE,UW,VN,VS,WT,WB)
!$OMP&private(PSI_SUM,FLUX_X,FLUX_Y,FLUX_Z,DEN_YN,DEN_Y)
       DO K=1,N3M
        KP1=KPV(K)
        KM1=KMV(K)
       DO J=JBG,N2M
         JP1=JPV(J)
         JM1=JMV(J)
       DO I=1,N1M
         IP1=IPV(I)
         IM1=IMV(I)

      IF ( REGION(I,J,K) .EQ. 0 ) THEN
       VF_Y(I,J,K,1)=VISM
       VF_Y(IM1,J,K,1)=VISM
       VF_Y(I,J,K,2)=VISM
       VF_Y(I,JM1,K,2)=VISM
       VF_Y(I,J,K,3)=VISM
       VF_Y(I,J,KM1,3)=VISM

       DF_Y(I,J,K,1)=DENM
       DF_Y(IM1,J,K,1)=DENM
       DF_Y(I,J,K,2)=DENM
       DF_Y(I,JM1,K,2)=DENM
       DF_Y(I,J,K,3)=DENM
       DF_Y(I,J,KM1,3)=DENM

       PSI_Y(I,J,K)=PSI_YN(I,J,K)

      ELSE IF (REGION(I,J,K) .EQ. 1 ) THEN
       VF_Y(I,J,K,1)=VISP
       VF_Y(IM1,J,K,1)=VISP
       VF_Y(I,J,K,2)=VISP
       VF_Y(I,JM1,K,2)=VISP
       VF_Y(I,J,K,3)=VISP
       VF_Y(I,J,KM1,3)=VISP

       DF_Y(I,J,K,1)=DENP
       DF_Y(IM1,J,K,1)=DENP
       DF_Y(I,J,K,2)=DENP
       DF_Y(I,JM1,K,2)=DENP
       DF_Y(I,J,K,3)=DENP
       DF_Y(I,J,KM1,3)=DENP

       PSI_Y(I,J,K)=PSI_YN(I,J,K)

      ELSE 

      UE=0.5*(SDY(JM1)*U(IP1,J,K)+SDY(J)*U(IP1,JM1,K))*VVDY(J)
      UW=0.5*(SDY(JM1)*U(I,J,K)+SDY(J)*U(I,JM1,K))*VVDY(J)
      VN=0.5*(V(I,JP1,K)+V(I,J,K))
      VS=0.5*(V(I,J,K)+V(I,JM1,K))
      WT=0.5*(SDY(JM1)*W(I,J,KP1)+SDY(J)*W(I,JM1,KP1))*VVDY(J)
      WB=0.5*(SDY(JM1)*W(I,J,K)+SDY(J)*W(I,JM1,K))*VVDY(J)

!      PSI1=0.5*(PSI_YN(I,J,K)*SDX(IP1)+PSI_YN(IP1,J,K)*SDX(I))*VVDX(IP1)
!      PSI2=0.5*(PSI_YN(IM1,J,K)*SDX(I)+PSI_YN(I,J,K)*SDX(IM1))*VVDX(I)
!      PSI3=0.5*(PSI_YN(I,J,K)+PSI_YN(I,JP1,K))
!      PSI4=0.5*(PSI_YN(I,J,K)+PSI_YN(I,JM1,K))
!      PSI5=0.5*(PSI_YN(I,J,K)*SDZ(KP1)+PSI_YN(I,J,KP1)*SDZ(K))*VVDZ(KP1)
!      PSI6=0.5*(PSI_YN(I,J,KM1)*SDZ(K)+PSI_YN(I,J,K)*SDZ(KM1))*VVDZ(K)
!        
!      DF_Y(I,J,K,1)=DENM+DEN_DIFF*PSI1
!      DF_Y(IM1,J,K,1)=DENM+DEN_DIFF*PSI2
!      DF_Y(I,J,K,2)=DENM+DEN_DIFF*PSI3
!      DF_Y(I,JM1,K,2)=DENM+DEN_DIFF*PSI4
!      DF_Y(I,J,K,3)=DENM+DEN_DIFF*PSI5
!      DF_Y(I,J,KM1,3)=DENM+DEN_DIFF*PSI6
!
!      VF_Y(I,J,K,1)=VISM+VIS_DIFF*PSI1
!      VF_Y(IM1,J,K,1)=VISM+VIS_DIFF*PSI2
!      VF_Y(I,J,K,2)=VISM+VIS_DIFF*PSI3
!      VF_Y(I,JM1,K,2)=VISM+VIS_DIFF*PSI4
!      VF_Y(I,J,K,3)=VISM+VIS_DIFF*PSI5
!      VF_Y(I,J,KM1,3)=VISM+VIS_DIFF*PSI6

      IF ( DF_Y(I,J,K,1) .EQ. 0. ) THEN
       IF (ABS(UE) .LE. 1.E-8) THEN
         PSI_SUM=0.5*(PSI_YN(I,J,K)*SDX(IP1)
     &                              +PSI_YN(IP1,J,K)*SDX(I))*VVDX(IP1)
       ELSE
        CALL PSIFACE_Y(I,J,K,UE,PSI_SUM,1)
       ENDIF
        VF_Y(I,J,K,1)=VISM+(VISP-VISM)*PSI_SUM
        DF_Y(I,J,K,1)=DENM+DEN_DIFF*PSI_SUM
      ENDIF
      IF ( DF_Y(IM1,J,K,1) .EQ. 0. ) THEN
       IF (ABS(UW) .LE. 1.E-8) THEN
         PSI_SUM=0.5*(PSI_YN(IM1,J,K)*SDX(I)
     &                              +PSI_YN(I,J,K)*SDX(IM1))*VVDX(I)
       ELSE
        CALL PSIFACE_Y(IM1,J,K,UW,PSI_SUM,1)
       ENDIF
        VF_Y(IM1,J,K,1)=VISM+(VISP-VISM)*PSI_SUM
        DF_Y(IM1,J,K,1)=DENM+DEN_DIFF*PSI_SUM
      ENDIF
      IF ( DF_Y(I,J,K,2) .EQ. 0. ) THEN
       IF (ABS(VN) .LE. 1.E-8) THEN
         PSI_SUM=0.5*(V(I,JP1,K)+V(I,J,K))
       ELSE
        CALL PSIFACE_Y(I,J,K,VN,PSI_SUM,2)
       ENDIF
        VF_Y(I,J,K,2)=VISM+(VISP-VISM)*PSI_SUM
        DF_Y(I,J,K,2)=DENM+DEN_DIFF*PSI_SUM
       ENDIF
      IF ( DF_Y(I,JM1,K,2) .EQ. 0. ) THEN
       IF (ABS(VS) .LE. 1.E-8) THEN
         PSI_SUM=0.5*(PSI_YN(I,J,K)+PSI_YN(I,JM1,K))
       ELSE
        CALL PSIFACE_Y(I,JM1,K,VS,PSI_SUM,2)
       ENDIF
        VF_Y(I,JM1,K,2)=VISM+(VISP-VISM)*PSI_SUM
        DF_Y(I,JM1,K,2)=DENM+DEN_DIFF*PSI_SUM
      ENDIF
      IF ( DF_Y(I,J,K,3) .EQ. 0. ) THEN
       IF (ABS(WT) .LE. 1.E-8) THEN
         PSI_SUM=0.5*(PSI_YN(I,J,K)*SDZ(KP1)
     &                              +PSI_YN(I,J,KP1)*SDZ(K))*VVDZ(KP1)
       ELSE
        CALL PSIFACE_Y(I,J,K,WT,PSI_SUM,3)
       ENDIF
        VF_Y(I,J,K,3)=VISM+(VISP-VISM)*PSI_SUM
        DF_Y(I,J,K,3)=DENM+DEN_DIFF*PSI_SUM
      ENDIF
      IF ( DF_Y(I,J,KM1,3) .EQ. 0. ) THEN
       IF (ABS(WB) .LE. 1.E-8) THEN
         PSI_SUM=0.5*(PSI_YN(I,J,KM1)*SDZ(K)
     &                              +PSI_YN(I,J,K)*SDZ(KM1))*VVDZ(K)
       ELSE
        CALL PSIFACE_Y(I,J,KM1,WB,PSI_SUM,3)        !K-1 NOT KMV(K)
       ENDIF
        VF_Y(I,J,KM1,3)=VISM+(VISP-VISM)*PSI_SUM
        DF_Y(I,J,KM1,3)=DENM+DEN_DIFF*PSI_SUM
      ENDIF

       FLUX_X=( DF_Y(I,J,K,1)*UE-DF_Y(IM1,J,K,1)*UW )*SSDX(I)
       FLUX_Y=( DF_Y(I,J,K,2)*VN-DF_Y(I,JM1,K,2)*VS )*VVDY(J)
       FLUX_Z=( DF_Y(I,J,K,3)*WT-DF_Y(I,J,KM1,3)*WB )*SSDZ(K)

       DEN_YN=DENM+DEN_DIFF*PSI_YN(I,J,K)
       DEN_Y=DEN_YN-DTCONST*(FLUX_X+FLUX_Y+FLUX_Z)
       PSI_Y(I,J,K)=(DEN_Y-DENM)/DEN_DIFF

       ENDIF
       ENDDO
       ENDDO
       ENDDO

!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M
       DO J=JBG,N2M
       DO I=1,N1M
         IF (PSI_Y(I,J,K) .GT. 1.) THEN
           PSI_Y(I,J,K)=1.
         ELSE IF (PSI_Y(I,J,K) .LT. 0.) THEN
           PSI_Y(I,J,K)=0.
         ENDIF
       ENDDO
       ENDDO
       ENDDO

C-----CONTINUITY BOUNDARY CONDITION : NEUMAN+PERIODIC CONDITION.
      IF ( IPX .NE. 1 ) THEN
!$OMP PARALLEL DO private(J)
       DO K=1,N3M
       DO J=JBG,N2M
       PSI_Y(0,J,K)=PSI_Y(1,J,K)
       PSI_Y(N1,J,K)=PSI_Y(N1M,J,K)
      ENDDO
      ENDDO
      ENDIF
      IF ( IPY .NE. 1 ) THEN
!$OMP PARALLEL DO private(I)
       DO K=1,N3M
       DO I=1,N1M
       PSI_Y(I,1,K)=PSI_Y(I,2,K)
       PSI_Y(I,N2,K)=PSI_Y(I,N2M,K)
      ENDDO
      ENDDO
      ENDIF
      IF ( IPZ .NE. 1 ) THEN
!$OMP PARALLEL DO private(I)
       DO J=JBG,N2M
       DO I=1,N1M
       PSI_Y(I,J,0)=PSI_Y(I,J,1)
       PSI_Y(I,J,N3)=PSI_Y(I,J,N3M)
      ENDDO
      ENDDO
      ENDIF
C-----CONTINUITY BOUNDARY CONDITION : NEUMAN+PERIODIC CONDITION.
!
!        CRI1=1.E-4
!        CRI2=1.-CRI1
!C-----REDISTRIBUTION ALGORITHM
!       DO 5000 PSI_ITERATION=1,20
!        EPSMAX=0.
!!$OMP PARALLEL DO private(I,J)
!       DO K=0,N3
!       DO J=0,N2
!       DO I=0,N1
!        EPS(I,J,K)=0.
!       ENDDO
!       ENDDO
!       ENDDO
!       
!!$OMP PARALLEL DO private(I,J)
!       DO K=1,N3M
!       DO J=JBG,N2M
!       DO I=1,N1M
!      IF((PSI_Y(I,J,K) .GT. 0.).AND.(PSI_Y(I,J,K) .LE. CRI1)) THEN
!        EPS(I,J,K)=-PSI_Y(I,J,K)
!      ELSE IF((PSI_Y(I,J,K) .GE. CRI2).AND.(PSI_Y(I,J,K) .LT. 1.))THEN  
!        EPS(I,J,K)=1.-PSI_Y(I,J,K)
!      ELSE IF(PSI_Y(I,J,K) .GT. 1.) THEN
!         EPS(I,J,K)=1.-PSI_Y(I,J,K)
!      ELSE IF(PSI_Y(I,J,K) .LT. 0.) THEN
!         EPS(I,J,K)=-PSI_Y(I,J,K) 
!      ENDIF
!       ENDDO
!       ENDDO
!       ENDDO
!
!!$OMP PARALLEL DO private(I,J)
!       DO K=1,N3M
!       DO J=JBG,N2M
!       DO I=1,N1M
!        IF ( ABS(EPS(I,J,K)) .GT. 1.E-8 ) THEN
!          EPS_OLD(I,J,K)=EPS(I,J,K)
!        ENDIF
!       ENDDO
!       ENDDO
!       ENDDO
!       
!!$OMP PARALLEL DO private(I,J)
!       DO K=0,N3
!       DO J=0,N2
!       DO I=0,N1
!       DFX(I,J,K)=0.
!       DFY(I,J,K)=0.
!       DFZ(I,J,K)=0.
!       ENDDO
!       ENDDO
!       ENDDO
!       
!       !define the unit velocity toward phase-interface
!!$OMP PARALLEL DO private(I,J,IP1,IM1,KP1,KM1,AA,BB,CC,AABBCC,AABBCCI)
!!$OMP&private(PSI_YE,PSI_YW,PSI_YM,PSI_YS,PSI_YT,PSI_YB)
!       DO K=1,N3M
!       DO J=JBG,N2M
!       DO I=1,N1M  
!
!        IF ( ABS(EPS(I,J,K)) .GT. 1.E-8 ) THEN
!           IP1=IPV(I)
!           IM1=IMV(I)
!      PSI_YE=0.5*(PSI_Y(IP1,J,K)*SDX(I)+PSI_Y(I,J,K)*SDX(IP1))*VVDX(IP1)
!      PSI_YW=0.5*(PSI_Y(I,J,K)*SDX(IM1)+PSI_Y(IM1,J,K)*SDX(I))*VVDX(I)
!        AA=( PSI_YE-PSI_YW )*SSDX(I)
!        
!      PSI_Ym=0.5*(PSI_Y(I,JPV(J),K)+PSI_Y(I,J,K))
!      PSI_YS=0.5*(PSI_Y(I,J,K)+PSI_Y(I,JMV(J),K))
!        BB=( PSI_Ym-PSI_YS )*VVDY(J)  !PSI_Ym BECAUSE OF PSI_YN
!
!           KP1=KPV(K)
!           KM1=KMV(K)
!      PSI_YT=0.5*(PSI_Y(I,J,KP1)*SDZ(K)+PSI_Y(I,J,K)*SDZ(KP1))*VVDZ(KP1)
!      PSI_YB=0.5*(PSI_Y(I,J,K)*SDZ(KM1)+PSI_Y(I,J,KM1)*SDZ(K))*VVDZ(K)
!        CC=( PSI_YT-PSI_YB )*SSDZ(K)
!        
!        AABBCC=SQRT(AA**2+BB**2+CC**2)
!        IF (AABBCC .GE. 1.E-4 ) THEN
!          AABBCCI=1./AABBCC
!         IF( PSI_Y(I,J,K) .GE. 0.5) THEN
!          DFX(I,J,K)=-AA*AABBCCI
!          DFY(I,J,K)=-BB*AABBCCI
!          DFZ(I,J,K)=-CC*AABBCCI
!         ELSE
!          DFX(I,J,K)=AA*AABBCCI
!          DFY(I,J,K)=BB*AABBCCI
!          DFZ(I,J,K)=CC*AABBCCI
!         ENDIF
!        ENDIF
!
!       ENDIF
!
!       ENDDO
!       ENDDO
!       ENDDO
!
!      IF ( IPX .NE. 1 ) THEN
!!$OMP PARALLEL DO private(J)
!      DO K=1,N3M
!      DO J=JBG,N2M
!       DFX(0,J,K)=0.!DFX(1,J,K)
!       DFX(N1,J,K)=0.!DFX(N1M,J,K)
!      ENDDO
!      ENDDO
!      ENDIF
!      IF ( IPY .NE. 1 ) THEN
!!$OMP PARALLEL DO private(I)
!      DO K=1,N3M
!      DO I=1,N1M
!       DFY(I,1,K)=0.!DFY(I,2,K)
!       DFY(I,N2,K)=0.!DFY(I,N2M,K)
!      ENDDO
!      ENDDO
!      ENDIF
!      IF ( IPZ .NE. 1 ) THEN
!!$OMP PARALLEL DO private(I)
!      DO J=JBG,N2M
!      DO I=1,N1M
!       DFZ(I,J,0)=0.!DFZ(I,J,1)
!       DFZ(I,J,N3)=0.!DFZ(I,J,N3M)
!      ENDDO
!      ENDDO
!      ENDIF
!
!      DO 2500 STEADY_EPS=1,5  !FIND STEADY SOL.
!
!        EPS_DIF=0.
!       !cal transport eqn.
!!$OMP PARALLEL DO private(I,J,IP1,IM1,KP1,KM1)
!!$OMP&private(DFXE,DFXW,DFYN,DFYS,DFZT,DFZB,EPSX,EPSY,EPSZ)
!!$OMP&reduction(MAX:EPS_DIF)
!       DO K=1,N3M
!       DO J=JBG,N2M
!       DO I=1,N1M
!
!        IF ( ABS(EPS(I,J,K)) .GT. 1.E-8 ) THEN
!           IP1=IPV(I)
!           IM1=IMV(I)
!       DFXE=0.5*(EPS(IP1,J,K)*SDX(I)+EPS(I,J,K)*SDX(IP1))*VVDX(IP1)
!       DFXW=0.5*(EPS(I,J,K)*SDX(IM1)+EPS(IM1,J,K)*SDX(I))*VVDX(I)
!        EPSX=( DFX(I,J,K)*DFXE-DFX(IM1,J,K)*DFXW )*SSDX(I)
!        
!       DFYN=0.5*(EPS(I,JPV(J),K)+EPS(I,J,K))
!       DFYS=0.5*(EPS(I,J,K)+EPS(I,JMV(J),K))
!        EPSY=( DFY(I,J,K)*DFYN-DFY(I,JMV(J),K)*DFYS )*SSDY(J)
!
!           KP1=KPV(K)
!           KM1=KMV(K)
!       DFZT=0.5*(EPS(I,J,KP1)*SDZ(K)+EPS(I,J,K)*SDZ(KP1))*VVDZ(KP1)
!       DFZB=0.5*(EPS(I,J,K)*SDZ(KM1)+EPS(I,J,KM1)*SDZ(K))*VVDZ(K)
!        EPSZ=( DFZ(I,J,K)*DFZT-DFZ(I,J,KM1)*DFZB )*SSDZ(K)
!
!        EPS(I,J,K)=EPS_OLD(I,J,K)-DT_CON*( EPSX+EPSY+EPSZ )
!        EPS_DIF=AMAX1(EPS_DIF,ABS(EPS(I,J,K)-EPS_OLD(I,J,K)))
!         ENDIF
!
!       ENDDO
!       ENDDO
!       ENDDO
!       
!!$OMP PARALLEL DO private(I,J)
!       DO K=1,N3M
!       DO J=JBG,N2M
!       DO I=1,N1M
!        IF ( ABS(EPS(I,J,K)) .GT. 1.E-8 ) THEN
!          EPS_OLD(I,J,K)=EPS(I,J,K)
!        ENDIF
!       ENDDO
!       ENDDO
!       ENDDO
!
!!         write(*,124) PSI_ITERATION,STEADY_EPS,EPS_DIF
!! 124  format(2I4,1ES12.5)
!          IF ( EPS_DIF .LE. 1.E-8 ) GOTO 200
!        
! 2500  CONTINUE
!
!c------------------UPDATE PSI_R
! 200   CONTINUE
! 
!!$OMP PARALLEL DO private(I,J)
!!$OMP&reduction(MAX:EPSMAX)
!       DO K=1,N3M 
!       DO J=JBG,N2M
!       DO I=1,N1M    
!        IF ( ABS(EPS(I,J,K)) .GT. 1.E-8 ) THEN
!         PSI_Y(I,J,K)=PSI_Y(I,J,K)+EPS(I,J,K)
!         EPSMAX=MAX(EPSMAX,ABS(EPS(I,J,K)))  
!        ENDIF
!       ENDDO
!       ENDDO
!       ENDDO
!
!C-----CONTINUITY BOUNDARY CONDITION : NEUMAN+PERIODIC CONDITION.
!      IF ( IPX .NE. 1 ) THEN
!!$OMP PARALLEL DO private(J)
!       DO K=1,N3M
!       DO J=JBG,N2M
!       PSI_Y(0,J,K)=PSI_Y(1,J,K)
!       PSI_Y(N1,J,K)=PSI_Y(N1M,J,K)
!      ENDDO
!      ENDDO
!      ENDIF
!      IF ( IPY .NE. 1 ) THEN
!!$OMP PARALLEL DO private(I)
!       DO K=1,N3M
!       DO I=1,N1M  
!       PSI_Y(I,1,K)=PSI_Y(I,2,K)
!       PSI_Y(I,N2,K)=PSI_Y(I,N2M,K)
!      ENDDO
!      ENDDO
!      ENDIF
!      IF ( IPZ .NE. 1 ) THEN
!!$OMP PARALLEL DO private(I)
!       DO J=JBG,N2M
!       DO I=1,N1M  
!       PSI_Y(I,J,0)=PSI_Y(I,J,1)
!       PSI_Y(I,J,N3)=PSI_Y(I,J,N3M)
!      ENDDO
!      ENDDO
!      ENDIF
!C-----CONTINUITY BOUNDARY CONDITION : NEUMAN+PERIODIC CONDITION.
!
!!       write(*,123) PSI_ITERATION,STEADY_EPS,EPSMAX
!! 123  format('eps_MAX ',2I4,1ES12.5)      
!          IF ( EPSMAX .LE. 1.E-8) GOTO 201
! 5000  ENDDO
!
!       WRITE(*,*) 'Y-CONTINUITY ITERATION LIMIT EXCEEDED!!!'
!
! 201   CONTINUE
!       
!!       write(*,124) PSI_ITERATION,EPSMAX
!! 124  format('eps_MAX ',I4,1ES12.5)

       RETURN
       END

C =====================================================================
      SUBROUTINE CONTI_Z(ITRACKING,REGION,DF_Z,VF_Z,PSI_Z,DT_CON,U,V,W
     &,PSI_ZN)
C =====================================================================
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL PSI_ZN(0:M1,0:M2,0:M3)

      REAL DFX(0:M1,0:M2,0:M3),DFY(0:M1,0:M2,0:M3),DFZ(0:M1,0:M2,0:M3)
      REAL EPS(0:M1,0:M2,0:M3),EPS_OLD(0:M1,0:M2,0:M3)
      INTEGER PSI_ITERATION,STEADY_EPS
      
      INTEGER REGION(-1:M1+1,-1:M2+1,-1:M3+1)
      REAL DF_Z(0:M1,0:M2,0:M3,3),VF_Z(0:M1,0:M2,0:M3,3) 
      REAL PSI_Z(0:M1,0:M2,0:M3)
      REAL TAUC(3,2)

      INTEGER ITRACKING
      REAL DT_CON
      
      INTEGER I,J,K,IP1,IM1,JP1,JM1,KP1,KM1
      REAL UE,UW,VN,VS,WT,WB
      REAL PSI_SUM,FLUX_X,FLUX_Y,FLUX_Z,DEN_ZN,DEN_Z
      
!$OMP PARALLEL DO private(I,J)
       DO K=0,N3
       DO J=0,N2
       DO I=0,N1
        DF_Z(I,J,K,1)=0.
        DF_Z(I,J,K,2)=0.
        DF_Z(I,J,K,3)=0.
       ENDDO
       ENDDO
       ENDDO

      !DEFINE PSI_Z
!$OMP PARALLEL DO private(I,J,IP1,IM1,JP1,JM1,KP1,KM1,UE,UW,VN,VS,WT,WB)
!$OMP&private(PSI_SUM,FLUX_X,FLUX_Y,FLUX_Z,DEN_ZN,DEN_Z)
       DO K=KBG,N3M
        KP1=KPV(K)
        KM1=KMV(K)
       DO J=1,N2M
         JP1=JPV(J)
         JM1=JMV(J)
       DO I=1,N1M
         IP1=IPV(I)
         IM1=IMV(I)

      IF ( REGION(I,J,K) .EQ. 0 ) THEN
       VF_Z(I,J,K,1)=VISM
       VF_Z(IM1,J,K,1)=VISM
       VF_Z(I,J,K,2)=VISM
       VF_Z(I,JM1,K,2)=VISM
       VF_Z(I,J,K,3)=VISM
       VF_Z(I,J,KM1,3)=VISM

       DF_Z(I,J,K,1)=DENM
       DF_Z(IM1,J,K,1)=DENM
       DF_Z(I,J,K,2)=DENM
       DF_Z(I,JM1,K,2)=DENM
       DF_Z(I,J,K,3)=DENM
       DF_Z(I,J,KM1,3)=DENM

       PSI_Z(I,J,K)=PSI_ZN(I,J,K)

      ELSE IF (REGION(I,J,K) .EQ. 1 ) THEN
       VF_Z(I,J,K,1)=VISP
       VF_Z(IM1,J,K,1)=VISP
       VF_Z(I,J,K,2)=VISP
       VF_Z(I,JM1,K,2)=VISP
       VF_Z(I,J,K,3)=VISP
       VF_Z(I,J,KM1,3)=VISP

       DF_Z(I,J,K,1)=DENP
       DF_Z(IM1,J,K,1)=DENP
       DF_Z(I,J,K,2)=DENP
       DF_Z(I,JM1,K,2)=DENP
       DF_Z(I,J,K,3)=DENP
       DF_Z(I,J,KM1,3)=DENP

       PSI_Z(I,J,K)=PSI_ZN(I,J,K)

      ELSE

      UE=0.5*(SDZ(KM1)*U(IP1,J,K)+SDZ(K)*U(IP1,J,KM1))*VVDZ(K)
      UW=0.5*(SDZ(KM1)*U(I,J,K)+SDZ(K)*U(I,J,KM1))*VVDZ(K)
      VN=0.5*(SDZ(KM1)*V(I,JP1,K)+SDZ(K)*V(I,JP1,KM1))*VVDZ(K)
      VS=0.5*(SDZ(KM1)*V(I,J,K)+SDZ(K)*V(I,J,KM1))*VVDZ(K)
      WT=0.5*(W(I,J,KP1)+W(I,J,K))
      WB=0.5*(W(I,J,K)+W(I,J,KM1))

!      PSI1=0.5*(PSI_ZN(I,J,K)*SDX(IP1)+PSI_ZN(IP1,J,K)*SDX(I))*VVDX(IP1)
!      PSI2=0.5*(PSI_ZN(IM1,J,K)*SDX(I)+PSI_ZN(I,J,K)*SDX(IM1))*VVDX(I)
!      PSI3=0.5*(PSI_ZN(I,J,K)*SDY(JP1)+PSI_ZN(I,JP1,K)*SDY(J))*VVDY(JP1)
!      PSI4=0.5*(PSI_ZN(I,JM1,K)*SDY(J)+PSI_ZN(I,J,K)*SDY(JM1))*VVDY(J)
!      PSI5=0.5*(PSI_ZN(I,J,K)+PSI_ZN(I,J,KP1))
!      PSI6=0.5*(PSI_ZN(I,J,K)+PSI_ZN(I,J,KM1))
!        
!      DF_Z(I,J,K,1)=DENM+DEN_DIFF*PSI1
!      DF_Z(IM1,J,K,1)=DENM+DEN_DIFF*PSI2
!      DF_Z(I,J,K,2)=DENM+DEN_DIFF*PSI3
!      DF_Z(I,JM1,K,2)=DENM+DEN_DIFF*PSI4
!      DF_Z(I,J,K,3)=DENM+DEN_DIFF*PSI5
!      DF_Z(I,J,KM1,3)=DENM+DEN_DIFF*PSI6
!
!      VF_Z(I,J,K,1)=VISM+VIS_DIFF*PSI1
!      VF_Z(IM1,J,K,1)=VISM+VIS_DIFF*PSI2
!      VF_Z(I,J,K,2)=VISM+VIS_DIFF*PSI3
!      VF_Z(I,JM1,K,2)=VISM+VIS_DIFF*PSI4
!      VF_Z(I,J,K,3)=VISM+VIS_DIFF*PSI5
!      VF_Z(I,J,KM1,3)=VISM+VIS_DIFF*PSI6
      
      IF ( DF_Z(I,J,K,1) .EQ. 0. ) THEN
       IF (ABS(UE) .LE. 1.E-8) THEN
         PSI_SUM=0.5*(PSI_ZN(I,J,K)*SDX(IP1)
     &                             +PSI_ZN(IP1,J,K)*SDX(I))*VVDX(IP1)
       ELSE
        CALL PSIFACE_Z(I,J,K,UE,PSI_SUM,1)
       ENDIF
        VF_Z(I,J,K,1)=VISM+(VISP-VISM)*PSI_SUM
        DF_Z(I,J,K,1)=DENM+DEN_DIFF*PSI_SUM
      ENDIF
      IF ( DF_Z(IM1,J,K,1) .EQ. 0. ) THEN
       IF (ABS(UW) .LE. 1.E-8) THEN
         PSI_SUM=0.5*(PSI_ZN(IM1,J,K)*SDX(I)
     &                             +PSI_ZN(I,J,K)*SDX(IM1))*VVDX(I)
       ELSE
        CALL PSIFACE_Z(IM1,J,K,UW,PSI_SUM,1)
       ENDIF
        VF_Z(IM1,J,K,1)=VISM+(VISP-VISM)*PSI_SUM
        DF_Z(IM1,J,K,1)=DENM+DEN_DIFF*PSI_SUM
      ENDIF
      IF ( DF_Z(I,J,K,2) .EQ. 0. ) THEN
       IF (ABS(VN) .LE. 1.E-8) THEN
         PSI_SUM=0.5*(PSI_ZN(I,J,K)*SDY(JP1)
     &                             +PSI_ZN(I,JP1,K)*SDY(J))*VVDY(JP1)
       ELSE
        CALL PSIFACE_Z(I,J,K,VN,PSI_SUM,2)
       ENDIF
        VF_Z(I,J,K,2)=VISM+(VISP-VISM)*PSI_SUM
        DF_Z(I,J,K,2)=DENM+DEN_DIFF*PSI_SUM
       ENDIF
      IF ( DF_Z(I,JM1,K,2) .EQ. 0. ) THEN
       IF (ABS(VS) .LE. 1.E-8) THEN
         PSI_SUM=0.5*(PSI_ZN(I,JM1,K)*SDY(J)
     &                             +PSI_ZN(I,J,K)*SDY(JM1))*VVDY(J)
       ELSE
        CALL PSIFACE_Z(I,JM1,K,VS,PSI_SUM,2)
       ENDIF
        VF_Z(I,JM1,K,2)=VISM+(VISP-VISM)*PSI_SUM
        DF_Z(I,JM1,K,2)=DENM+DEN_DIFF*PSI_SUM
      ENDIF
      IF ( DF_Z(I,J,K,3) .EQ. 0. ) THEN
       IF (ABS(WT) .LE. 1.E-8) THEN
         PSI_SUM=0.5*(PSI_ZN(I,J,K)+PSI_ZN(I,J,KP1))
       ELSE
        CALL PSIFACE_Z(I,J,K,WT,PSI_SUM,3)
       ENDIF
        VF_Z(I,J,K,3)=VISM+(VISP-VISM)*PSI_SUM
        DF_Z(I,J,K,3)=DENM+DEN_DIFF*PSI_SUM
      ENDIF
      IF ( DF_Z(I,J,KM1,3) .EQ. 0. ) THEN
       IF (ABS(WB) .LE. 1.E-8) THEN
         PSI_SUM=0.5*(PSI_ZN(I,J,K)+PSI_ZN(I,J,KM1))
       ELSE
        CALL PSIFACE_Z(I,J,KM1,WB,PSI_SUM,3)        !K-1 NOT KMV(K)
       ENDIF
        VF_Z(I,J,KM1,3)=VISM+(VISP-VISM)*PSI_SUM
        DF_Z(I,J,KM1,3)=DENM+DEN_DIFF*PSI_SUM
      ENDIF

       FLUX_X=( DF_Z(I,J,K,1)*UE-DF_Z(IM1,J,K,1)*UW )*SSDX(I)
       FLUX_Y=( DF_Z(I,J,K,2)*VN-DF_Z(I,JM1,K,2)*VS )*SSDY(J)
       FLUX_Z=( DF_Z(I,J,K,3)*WT-DF_Z(I,J,KM1,3)*WB )*VVDZ(K)


       DEN_ZN=DENM+DEN_DIFF*PSI_ZN(I,J,K)
       DEN_Z=DEN_ZN-DTCONST*( FLUX_X+FLUX_Y+FLUX_Z )
       PSI_Z(I,J,K)=(DEN_Z-DENM)/DEN_DIFF

        ENDIF
       ENDDO
       ENDDO
       ENDDO

!$OMP PARALLEL DO private(I,J)
       DO K=KBG,N3M
       DO J=1,N2M
       DO I=1,N1M
         IF (PSI_Z(I,J,K) .GT. 1.) THEN
           PSI_Z(I,J,K)=1.
         ELSE IF (PSI_Z(I,J,K) .LT. 0.) THEN
           PSI_Z(I,J,K)=0.
         ENDIF
       ENDDO
       ENDDO
       ENDDO

C-----CONTINUITY BOUNDARY CONDITION : NEUMAN+PERIODIC CONDITION.
      IF ( IPX .NE. 1 ) THEN
!$OMP PARALLEL DO private(J)
      DO K=KBG,N3M
      DO J=1,N2M
       PSI_Z(0,J,K)=PSI_Z(1,J,K)
       PSI_Z(N1,J,K)=PSI_Z(N1M,J,K)
      ENDDO
      ENDDO
      ENDIF
      IF ( IPY .NE. 1 ) THEN
!$OMP PARALLEL DO private(I)
      DO K=KBG,N3M
      DO I=1,N1M
       PSI_Z(I,0,K)=PSI_Z(I,1,K)
       PSI_Z(I,N2,K)=PSI_Z(I,N2M,K)
      ENDDO
      ENDDO
      ENDIF
      IF ( IPZ .NE. 1 ) THEN
!$OMP PARALLEL DO private(I)
       DO J=1,N2M
       DO I=1,N1M
       PSI_Z(I,J,1)=PSI_Z(I,J,2)
       PSI_Z(I,J,N3)=PSI_Z(I,J,N3M)
      ENDDO
      ENDDO
      ENDIF
C-----CONTINUITY BOUNDARY CONDITION : NEUMAN+PERIODIC CONDITION.

!        CRI1=1.E-4
!        CRI2=1.-CRI1
!C-----REDISTRIBUTION ALGORITHM
!       DO 5000 PSI_ITERATION=1,20
!        EPSMAX=0.
!!$OMP PARALLEL DO private(I,J)
!       DO K=0,N3
!       DO J=0,N2
!       DO I=0,N1
!        EPS(I,J,K)=0.
!       ENDDO
!       ENDDO
!       ENDDO
!
!!$OMP PARALLEL DO private(I,J)
!       DO K=KBG,N3M
!       DO J=1,N2M
!       DO I=1,N1M
!      IF((PSI_Z(I,J,K) .GT. 0.).AND.(PSI_Z(I,J,K) .LE. CRI1))THEN  
!        EPS(I,J,K)=-PSI_Z(I,J,K) 
!      ELSE IF((PSI_Z(I,J,K) .GE. CRI2).AND.(PSI_Z(I,J,K) .LT. 1.))THEN  
!         EPS(I,J,K)=1.-PSI_Z(I,J,K)
!      ELSE IF(PSI_Z(I,J,K) .GT. 1.) THEN
!         EPS(I,J,K)=1.-PSI_Z(I,J,K) 
!      ELSE IF(PSI_Z(I,J,K) .LT. 0.) THEN
!         EPS(I,J,K)=-PSI_Z(I,J,K) 
!      ENDIF
!       ENDDO
!       ENDDO
!       ENDDO
!
!!$OMP PARALLEL DO private(I,J)
!       DO K=KBG,N3M
!       DO J=1,N2M
!       DO I=1,N1M
!        IF ( ABS(EPS(I,J,K)) .GT. 1.E-8 ) THEN
!          EPS_OLD(I,J,K)=EPS(I,J,K)
!        ENDIF
!       ENDDO
!       ENDDO
!       ENDDO
!       
!!$OMP PARALLEL DO private(I,J)
!       DO K=0,N3
!       DO J=0,N2
!       DO I=0,N1
!        DFX(I,J,K)=0.
!        DFY(I,J,K)=0.
!        DFZ(I,J,K)=0.
!       ENDDO
!       ENDDO
!       ENDDO
!
!       !define the unit velocity toward phase-interface
!!$OMP PARALLEL DO private(I,J,IP1,IM1,JP1,JM1,AA,BB,CC,AABBCC,AABBCCI)
!!$OMP&private(PSI_ZE,PSI_ZW,PSI_ZM,PSI_ZS,PSI_ZT,PSI_ZB)
!       DO K=KBG,N3M
!       DO J=1,N2M
!       DO I=1,N1M  
!           
!        IF ( ABS(EPS(I,J,K)) .GT. 1.E-8 ) THEN
!           IP1=IPV(I)
!           IM1=IMV(I)
!      PSI_ZE=0.5*(PSI_Z(IP1,J,K)*SDX(I)+PSI_Z(I,J,K)*SDX(IP1))*VVDX(IP1)
!      PSI_ZW=0.5*(PSI_Z(I,J,K)*SDX(IM1)+PSI_Z(IM1,J,K)*SDX(I))*VVDX(I)
!        AA=( PSI_ZE-PSI_ZW )*SSDX(I)
!        
!           JP1=JPV(J)
!           JM1=JMV(J)
!      PSI_Zm=0.5*(PSI_Z(I,JP1,K)*SDY(J)+PSI_Z(I,J,K)*SDY(JP1))*VVDY(JP1)
!      PSI_ZS=0.5*(PSI_Z(I,J,K)*SDY(JM1)+PSI_Z(I,JM1,K)*SDY(J))*VVDY(J)
!        BB=( PSI_Zm-PSI_ZS )*SSDY(J)  !PSI_Zm BECAUSE OF PSI_ZN
!
!      PSI_ZT=0.5*(PSI_Z(I,J,KPV(K))+PSI_Z(I,J,K))
!      PSI_ZB=0.5*(PSI_Z(I,J,K)+PSI_Z(I,J,KMV(K)))
!        CC=( PSI_ZT-PSI_ZB )*SSDZ(K)
!
!        AABBCC=SQRT(AA**2+BB**2+CC**2)
!        IF (AABBCC .GE. 1.E-4 ) THEN
!          AABBCCI=1./AABBCC
!         IF( PSI_Z(I,J,K) .GE. 0.5) THEN
!          DFX(I,J,K)=-AA*AABBCCI
!          DFY(I,J,K)=-BB*AABBCCI
!          DFZ(I,J,K)=-CC*AABBCCI
!         ELSE
!          DFX(I,J,K)=AA*AABBCCI
!          DFY(I,J,K)=BB*AABBCCI
!          DFZ(I,J,K)=CC*AABBCCI
!         ENDIF
!        ENDIF
!        
!       ENDIF
!
!       ENDDO
!       ENDDO
!       ENDDO  
!      IF ( IPX .NE. 1 ) THEN
!!$OMP PARALLEL DO private(J)
!      DO K=KBG,N3M
!      DO J=1,N2M
!       DFX(0,J,K)=0.!DFX(1,J,K)
!       DFX(N1,J,K)=0.!DFX(N1M,J,K)
!      ENDDO
!      ENDDO
!      ENDIF
!      IF ( IPY .NE. 1 ) THEN
!!$OMP PARALLEL DO private(I)
!      DO K=KBG,N3M
!      DO I=1,N1M
!       DFY(I,0,K)=0.!DFY(I,1,K)
!       DFY(I,N2,K)=0.!DFY(I,N2M,K)
!      ENDDO
!      ENDDO
!      ENDIF
!      IF ( IPZ .NE. 1 ) THEN
!!$OMP PARALLEL DO private(I)
!      DO J=1,N2M
!      DO I=1,N1M
!       DFZ(I,J,1)=0.!DFZ(I,J,1)
!       DFZ(I,J,N3)=0.!DFZ(I,J,N3M)
!      ENDDO
!      ENDDO
!      ENDIF
!      
!      DO 2500 STEADY_EPS=1,5  !FIND STEADY SOL.
!
!        EPS_DIF=0.
!       !cal transport eqn.
!!$OMP PARALLEL DO private(I,J,IP1,IM1,JP1,JM1,KP1,KM1)
!!$OMP&private(DFXE,DFXW,DFYN,DFYS,DFZT,DFZB,EPSX,EPSY,EPSZ)
!!$OMP&reduction(MAX:EPS_DIF)
!       DO K=KBG,N3M
!           KP1=KPV(K)
!           KM1=KMV(K)
!       DO J=1,N2M
!           JP1=JPV(J)
!           JM1=JMV(J)
!       DO I=1,N1M
!           IP1=IPV(I)
!           IM1=IMV(I)
!
!        IF ( ABS(EPS(I,J,K)) .GT. 1.E-8 ) THEN
!       DFXE=0.5*(EPS(IP1,J,K)*SDX(I)+EPS(I,J,K)*SDX(IP1))*VVDX(IP1)
!       DFXW=0.5*(EPS(I,J,K)*SDX(IM1)+EPS(IM1,J,K)*SDX(I))*VVDX(I)
!        EPSX=( DFX(I,J,K)*DFXE-DFX(IM1,J,K)*DFXW )*SSDX(I)
!
!       DFYN=0.5*(EPS(I,JP1,K)*SDY(J)+EPS(I,J,K)*SDY(JP1))*VVDY(JP1)
!       DFYS=0.5*(EPS(I,J,K)*SDY(JM1)+EPS(I,JM1,K)*SDY(J))*VVDY(J)
!        EPSY=( DFY(I,J,K)*DFYN-DFY(I,JM1,K)*DFYS )*SSDY(J)
!
!       DFZT=0.5*(EPS(I,J,KPV(K))+EPS(I,J,K))
!       DFZB=0.5*(EPS(I,J,K)+EPS(I,J,KMV(K)))
!        EPSZ=( DFZ(I,J,K)*DFZT-DFZ(I,J,KMV(K))*DFZB )*SSDZ(K)
!
!        EPS(I,J,K)=EPS_OLD(I,J,K)-DT_CON*( EPSX+EPSY+EPSZ )
!        EPS_DIF=AMAX1(EPS_DIF,ABS(EPS(I,J,K)-EPS_OLD(I,J,K)))
!       ENDIF
!
!       ENDDO
!       ENDDO
!       ENDDO
!
!!         write(*,124) PSI_ITERATION,STEADY_EPS,EPS_DIF
!! 124  format(2I4,1ES12.5)
!          IF ( EPS_DIF .LE. 1.E-8 ) GOTO 300
!          
!!$OMP PARALLEL DO private(I,J)
!       DO K=KBG,N3M
!       DO J=1,N2M
!       DO I=1,N1M
!        IF ( ABS(EPS(I,J,K)) .GT. 1.E-8 ) THEN
!          EPS_OLD(I,J,K)=EPS(I,J,K)
!        ENDIF
!       ENDDO
!       ENDDO
!       ENDDO
! 2500 ENDDO
!
!c------------------UPDATE PSI_T
! 300   CONTINUE
! 
!!$OMP PARALLEL DO private(I,J)
!!$OMP&reduction(MAX:EPSMAX)
!       DO K=KBG,N3M
!       DO J=1,N2M
!       DO I=1,N1M
!        IF ( ABS(EPS(I,J,K)) .GT. 1.E-8 ) THEN
!         PSI_Z(I,J,K)=PSI_Z(I,J,K)+EPS(I,J,K)
!         EPSMAX=MAX(EPSMAX,ABS(EPS(I,J,K)))
!        ENDIF   
!       ENDDO
!       ENDDO
!       ENDDO 
!
!C-----CONTINUITY BOUNDARY CONDITION : NEUMAN+PERIODIC CONDITION.
!      IF ( IPX .NE. 1 ) THEN
!!$OMP PARALLEL DO private(J)
!      DO K=KBG,N3M
!      DO J=1,N2M
!       PSI_Z(0,J,K)=PSI_Z(1,J,K)
!       PSI_Z(N1,J,K)=PSI_Z(N1M,J,K)
!      ENDDO
!      ENDDO
!      ENDIF
!      IF ( IPY .NE. 1 ) THEN
!!$OMP PARALLEL DO private(I)
!      DO K=KBG,N3M
!      DO I=1,N1M
!       PSI_Z(I,0,K)=PSI_Z(I,1,K)
!       PSI_Z(I,N2,K)=PSI_Z(I,N2M,K)
!      ENDDO
!      ENDDO
!      ENDIF
!      IF ( IPZ .NE. 1 ) THEN
!!$OMP PARALLEL DO private(I)
!      DO J=1,N2M
!      DO I=1,N1M
!       PSI_Z(I,J,0)=PSI_Z(I,J,1)
!       PSI_Z(I,J,N3)=PSI_Z(I,J,N3M)
!      ENDDO
!      ENDDO
!      ENDIF
!        
!C-----CONTINUITY BOUNDARY CONDITION : NEUMAN+PERIODIC CONDITION.
!
!!       write(*,123) PSI_ITERATION,STEADY_EPS,EPSMAX
!! 123  format('eps_MAX ',2I4,1ES12.5)
!          IF ( EPSMAX .LE. 1.E-8) GOTO 301
! 5000  ENDDO
!       WRITE(*,*) 'Z-CONTINUITY ITERATION LIMIT EXCEEDED!!!'
!       
! 301   CONTINUE
!       
!!       write(*,124) PSI_ITERATION,EPSMAX
!! 124  format('eps_MAX ',I4,1ES12.5)

       RETURN
       END

C =====================================================================
      SUBROUTINE CONTI_C(ITRACKING,REGION,DF_C,PSI_C,DT_CON,U,V,W
     &,PSI_XN,PSI_YN,PSI_ZN,PSI_CN)
C =====================================================================
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL PSI_XN(0:M1,0:M2,0:M3)
      REAL PSI_YN(0:M1,0:M2,0:M3)
      REAL PSI_ZN(0:M1,0:M2,0:M3)
      REAL PSI_CN(0:M1,0:M2,0:M3)

      INTEGER REGION(-1:M1+1,-1:M2+1,-1:M3+1)
      REAL DF_C(0:M1,0:M2,0:M3,3)
      REAL PSI_C(0:M1,0:M2,0:M3)

      INTEGER ITRACKING
      REAL DT_CON
      
      INTEGER I,J,K,IP1,IM1,JP1,JM1,KP1,KM1
      REAL UE,UW,VN,VS,WT,WB
      REAL PSI_SUM,FLUX_X,FLUX_Y,FLUX_Z,DEN_CN,DEN_C
      
!$OMP PARALLEL DO private(I,J)
       DO K=0,N3
       DO J=0,N2
       DO I=0,N1
        DF_C(I,J,K,1)=0.
        DF_C(I,J,K,2)=0.
        DF_C(I,J,K,3)=0.
       ENDDO
       ENDDO
       ENDDO

      !DEFINE PSI_C
!$OMP PARALLEL DO private(I,J,IP1,IM1,JP1,JM1,KP1,KM1,UE,UW,VN,VS,WT,WB)
!$OMP&private(PSI_SUM,FLUX_X,FLUX_Y,FLUX_Z,DEN_CN,DEN_C)
       DO K=1,N3M
        KP1=KPV(K)
        KM1=KMV(K)
       DO J=1,N2M
         JP1=JPV(J)
         JM1=JMV(J)
       DO I=1,N1M
         IP1=IPV(I)
         IM1=IMV(I)

      IF ( REGION(I,J,K) .EQ. 0 ) THEN
       PSI_C(I,J,K)=PSI_CN(I,J,K)
       DF_C(I,J,K,1)=DENM
       DF_C(IM1,J,K,1)=DENM
       DF_C(I,J,K,2)=DENM
       DF_C(I,JM1,K,2)=DENM
       DF_C(I,J,K,3)=DENM
       DF_C(I,J,KM1,3)=DENM

      ELSE IF (REGION(I,J,K) .EQ. 1 ) THEN
       PSI_C(I,J,K)=PSI_CN(I,J,K)
       DF_C(I,J,K,1)=DENP
       DF_C(IM1,J,K,1)=DENP
       DF_C(I,J,K,2)=DENP
       DF_C(I,JM1,K,2)=DENP
       DF_C(I,J,K,3)=DENP
       DF_C(I,J,KM1,3)=DENP

      ELSE

      UE=U(IP1,J,K)
      UW=U(I,J,K)
      VN=V(I,JP1,K)
      VS=V(I,J,K)
      WT=W(I,J,KP1)
      WB=W(I,J,K)

      IF ( DF_C(I,J,K,1) .EQ. 0. ) THEN
       IF (ABS(UE) .LE. 1.E-8) THEN
         PSI_SUM=PSI_XN(IP1,J,K)
       ELSE
        CALL PSIFACE_C(I,J,K,UE,PSI_SUM,1)
       ENDIF
        DF_C(I,J,K,1)=DENM+DEN_DIFF*PSI_SUM
      ENDIF
      IF ( DF_C(IM1,J,K,1) .EQ. 0. ) THEN
       IF (ABS(UW) .LE. 1.E-8) THEN
         PSI_SUM=PSI_XN(I,J,K)
       ELSE
        CALL PSIFACE_C(IM1,J,K,UW,PSI_SUM,1)
       ENDIF
        DF_C(IM1,J,K,1)=DENM+DEN_DIFF*PSI_SUM
      ENDIF
      IF ( DF_C(I,J,K,2) .EQ. 0. ) THEN
       IF (ABS(VN) .LE. 1.E-8) THEN
         PSI_SUM=PSI_YN(I,JP1,K)
       ELSE
        CALL PSIFACE_C(I,J,K,VN,PSI_SUM,2)
       ENDIF
        DF_C(I,J,K,2)=DENM+DEN_DIFF*PSI_SUM
       ENDIF
      IF ( DF_C(I,JM1,K,2) .EQ. 0. ) THEN
       IF (ABS(VS) .LE. 1.E-8) THEN
         PSI_SUM=PSI_YN(I,J,K)
       ELSE
        CALL PSIFACE_C(I,JM1,K,VS,PSI_SUM,2)
       ENDIF
        DF_C(I,JM1,K,2)=DENM+DEN_DIFF*PSI_SUM
      ENDIF
      IF ( DF_C(I,J,K,3) .EQ. 0. ) THEN
       IF (ABS(WT) .LE. 1.E-8) THEN
         PSI_SUM=PSI_CN(I,J,KP1)
       ELSE
        CALL PSIFACE_C(I,J,K,WT,PSI_SUM,3)
       ENDIF
        DF_C(I,J,K,3)=DENM+DEN_DIFF*PSI_SUM
      ENDIF
      IF ( DF_C(I,J,KM1,3) .EQ. 0. ) THEN
       IF (ABS(WB) .LE. 1.E-8) THEN
         PSI_SUM=PSI_CN(I,J,K)
       ELSE
        CALL PSIFACE_C(I,J,KM1,WB,PSI_SUM,3)        !K-1 NOT KMV(K)
       ENDIF
        DF_C(I,J,KM1,3)=DENM+DEN_DIFF*PSI_SUM
      ENDIF

       FLUX_X=( DF_C(I,J,K,1)*UE-DF_C(IM1,J,K,1)*UW )*SSDX(I)
       FLUX_Y=( DF_C(I,J,K,2)*VN-DF_C(I,JM1,K,2)*VS )*SSDY(J)
       FLUX_Z=( DF_C(I,J,K,3)*WT-DF_C(I,J,KM1,3)*WB )*SSDZ(K)

       DEN_CN=DENM+DEN_DIFF*PSI_CN(I,J,K)
       DEN_C=DEN_CN-DTCONST*( FLUX_X+FLUX_Y+FLUX_Z )
       PSI_C(I,J,K)=(DEN_C-DENM)/DEN_DIFF

        ENDIF
       ENDDO
       ENDDO
       ENDDO

!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
         IF (PSI_C(I,J,K) .GT. 1.) THEN
           PSI_C(I,J,K)=1.
         ELSE IF (PSI_C(I,J,K) .LT. 0.) THEN
           PSI_C(I,J,K)=0.
         ENDIF
       ENDDO
       ENDDO
       ENDDO

C-----CONTINUITY BOUNDARY CONDITION : NEUMAN+PERIODIC CONDITION.
      IF ( IPX .NE. 1 ) THEN
!$OMP PARALLEL DO private(J)
      DO K=1,N3M
      DO J=1,N2M
       PSI_C(0,J,K)=PSI_C(1,J,K)
       PSI_C(N1,J,K)=PSI_C(N1M,J,K)
      ENDDO
      ENDDO
      ENDIF
      IF ( IPY .NE. 1 ) THEN
!$OMP PARALLEL DO private(I)
      DO K=1,N3M
      DO I=1,N1M
       PSI_C(I,0,K)=PSI_C(I,1,K)
       PSI_C(I,N2,K)=PSI_C(I,N2M,K)
      ENDDO
      ENDDO
      ENDIF
      IF ( IPZ .NE. 1 ) THEN
!$OMP PARALLEL DO private(I)
       DO J=1,N2M
       DO I=1,N1M
       PSI_C(I,J,1)=PSI_C(I,J,2)
       PSI_C(I,J,N3)=PSI_C(I,J,N3M)
      ENDDO
      ENDDO
      ENDIF
C-----CONTINUITY BOUNDARY CONDITION : NEUMAN+PERIODIC CONDITION.

!       OPEN(145,FILE='0PSI_C.DAT')
!       WRITE(145,*) 'VARIABLES="X","Y","Z","PSI_C","DPSI_C"'
!      WRITE(145,*) 'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'
!       DO K=1,N3M
!       DO J=1,N2M
!       DO I=1,N1M
!         UU=0.5*(U(I,J,K)+U(IPV(I),J,K))
!         VV=0.5*(V(I,J,K)+V(I,JPV(J),K))
!         WW=0.5*(W(I,J,K)+W(I,J,KPV(K)))     
!       WRITE(145,147) XP(I),YP(J),ZP(K),PSI_C(I,J,K)
!     &                                  ,PSI_CN(I,J,K)-PSI_C(I,J,K)
!        ENDDO
!        ENDDO
!        ENDDO
!       close(145)
! 147  FORMAT(5F18.14) 

       RETURN
       END

C=======================================================================
      SUBROUTINE PSIFACE_X(I,J,K,UU,PSI_SUM,IDIR)
C=======================================================================  
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE FLOW_VAR

      USE LVS_COUPLING
      
      IMPLICIT NONE
      
      INTEGER I,J,K,IDIR
      REAL UU,PSI_SUM
      
      INTEGER I1,I2,J1,J2,K1,K2
      REAL X1,X2,Y1,Y2,Z1,Z2

       !X-DIR
       IF ( IDIR .EQ. 1) THEN
       IF ( UU .GE. 0.) THEN
        X1=XP(I)-UU*DTCONST
        X2=XP(I)
!              IF ( X1 .LT. X(1) ) X1=X(1)
        I1=ICOUMP2(I-1)      !IF CFL IS GREATER THAN 1, THIS MAY HAVE PROBLEM.
        I2=ICOUMP1(I)
       ELSE
        X1=XP(I)
        X2=XP(I)-UU*DTCONST
!              IF ( X2 .GT. X(N1) ) X2=X(N1)
        I1=ICOUMP2(I)      !IF CFL IS GREATER THAN 1, THIS MAY HAVE PROBLEM.
        I2=ICOUMP1(I+1)
       ENDIF
       Y1=Y(J)
       Y2=Y(J+1)
       Z1=Z(K)
       Z2=Z(K+1)

        J1=JCOU2(J)
        J2=JCOU1(J+1)
        K1=KCOU2(K)
        K2=KCOU1(K+1)

      !Y-DIR
      ELSE IF ( IDIR .EQ. 2) THEN

       X1=XP(I-1)
       X2=XP(I)
       IF ( UU .GE. 0.) THEN
       Y1=Y(J+1)-UU*DTCONST    
       Y2=Y(J+1)
!            IF ( Y1 .LT. Y(1) ) Y1=Y(1)
        J1=JCOU2(J)
        J2=JCOU1(J+1)
       ELSE
       Y1=Y(J+1)
       Y2=Y(J+1)-UU*DTCONST
!            IF ( Y2 .GT. Y(N2) ) Y2=Y(N2)
        J1=JCOU2(J+1)
        J2=JCOU1(J+2)
       ENDIF
       Z1=Z(K)
       Z2=Z(K+1)

        I1=ICOUMP2(I-1)
        I2=ICOUMP1(I)
        K1=KCOU2(K)
        K2=KCOU1(K+1)

       ELSE IF ( IDIR .EQ. 3) THEN

       X1=XP(I-1)
       X2=XP(I)
       Y1=Y(J)
       Y2=Y(J+1)
         IF ( UU .GE. 0. ) THEN
       Z1=Z(K+1)-UU*DTCONST
       Z2=Z(K+1)
!            IF ( Z1 .LT. Z(1) ) Z1=Z(1)
        K1=KCOU2(K)
        K2=KCOU1(K+1)
         ELSE
       Z1=Z(K+1)
       Z2=Z(K+1)-UU*DTCONST
!            IF ( Z2 .GT. Z(N3) ) Z2=Z(N3)
        K1=KCOU2(K+1)
        K2=KCOU1(K+2)
         ENDIF

        I1=ICOUMP2(I-1)
        I2=ICOUMP1(I)
        J1=JCOU2(J)
        J2=JCOU1(J+1)

      ENDIF

       CALL INTER_PSI_FACE(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,PSI_SUM)

!       if (idir .eq. 1) then
!      write(*,*) i,j,k,i1,i2,psi_sum
!      write(*,*) xf(i1),xf(i2)
!      write(*,*) x1,x2
!       endif

       RETURN
       END

C=======================================================================
      SUBROUTINE PSIFACE_Y(I,J,K,UU,PSI_SUM,IDIR)
C=======================================================================
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE FLOW_VAR

      USE LVS_COUPLING
      
      IMPLICIT NONE
      
      INTEGER I,J,K,IDIR
      REAL UU,PSI_SUM
      
      INTEGER I1,I2,J1,J2,K1,K2
      REAL X1,X2,Y1,Y2,Z1,Z2

      !X-DIR
       IF ( IDIR .EQ. 1) THEN
          IF ( UU .GE. 0. ) THEN
            X1=X(I+1)-UU*DTCONST
            X2=X(I+1)
!            IF ( X1 .LT. X(1) ) X1=X(1)
        I1=ICOU2(I)
        I2=ICOU1(I+1)
          ELSE
            X1=X(I+1)
            X2=X(I+1)-UU*DTCONST
 !           IF ( X2 .GT. X(N1) ) X2=X(N1)
        I1=ICOU2(I+1)
        I2=ICOU1(I+2)
          ENDIF
        Y1=YP(J-1)
        Y2=YP(J)
        Z1=Z(K)
        Z2=Z(K+1)

        J1=JCOUMP2(J-1)
        J2=JCOUMP1(J)
        K1=KCOU2(K)
        K2=KCOU1(K+1)

      !Y-DIR
      ELSE IF ( IDIR .EQ. 2) THEN

        X1=X(I)
        X2=X(I+1)
        IF ( UU .GE. 0. ) THEN
         Y1=YP(J)-UU*DTCONST
         Y2=YP(J)
!            IF ( Y1 .LT. Y(1) ) Y1=Y(1)
        J1=JCOUMP2(J-1)
        J2=JCOUMP1(J)
        ELSE
         Y1=YP(J)
         Y2=YP(J)-UU*DTCONST
!            IF ( Y2 .GT. Y(N2) ) Y2=Y(N2)
        J1=JCOUMP2(J)
        J2=JCOUMP1(J+1)
       ENDIF
        Z1=Z(K)
        Z2=Z(K+1)

        I1=ICOU2(I)
        I2=ICOU1(I+1)
        K1=KCOU2(K)
        K2=KCOU1(K+1)

      !Z-DIR
      ELSE IF ( IDIR .EQ. 3) THEN

        X1=X(I)
        X2=X(I+1)
        Y1=YP(J-1)
        Y2=YP(J)
         IF ( UU .GE. 0. ) THEN
       Z1=Z(K+1)-UU*DTCONST
       Z2=Z(K+1)
 !           IF ( Z1 .LT. Z(1) ) Z1=Z(1)
        K1=KCOU2(K)
        K2=KCOU1(K+1)
         ELSE
       Z1=Z(K+1)
       Z2=Z(K+1)-UU*DTCONST
!            IF ( Z2 .GT. Z(N3) ) Z2=Z(N3)
        K1=KCOU2(K+1)
        K2=KCOU1(K+2)
         ENDIF

        I1=ICOU2(I)
        I2=ICOU1(I+1)
        J1=JCOUMP2(J-1)
        J2=JCOUMP1(J)

      ENDIF

       CALL INTER_PSI_FACE(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,PSI_SUM)

       RETURN
       END

C=======================================================================
      SUBROUTINE PSIFACE_Z(I,J,K,UU,PSI_SUM,IDIR)
C=======================================================================  
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE FLOW_VAR

      USE LVS_COUPLING
      
      IMPLICIT NONE
      
      INTEGER I,J,K,IDIR
      REAL UU,PSI_SUM
      
      INTEGER I1,I2,J1,J2,K1,K2
      REAL X1,X2,Y1,Y2,Z1,Z2

      !X-DIR
       IF ( IDIR .EQ. 1) THEN  
        IF ( UU .GE. 0. ) THEN
       X1=X(I+1)-UU*DTCONST
       X2=X(I+1)
!            IF ( X1 .LT. X(1) ) X1=X(1)
        I1=ICOU2(I)
        I2=ICOU1(I+1)
        ELSE
       X1=X(I+1)
       X2=X(I+1)-UU*DTCONST
!            IF ( X2 .GT. X(N1) ) X2=X(N1)
        I1=ICOU2(I+1)
        I2=ICOU1(I+2)
        ENDIF
        Y1=Y(J)
        Y2=Y(J+1)
        Z1=ZP(K-1)
        Z2=ZP(K)

        J1=JCOU2(J)
        J2=JCOU1(J+1)
        K1=KCOUMP2(K-1)
        K2=KCOUMP1(K)

       !Y-DIR
       ELSE IF ( IDIR .EQ. 2) THEN

        X1=X(I)
        X2=X(I+1)
        IF ( UU .GE. 0. ) THEN
       Y1=Y(J+1)-UU*DTCONST
       Y2=Y(J+1)
!            IF ( Y1 .LT. Y(1) ) Y1=Y(1)
        J1=JCOU2(J)
        J2=JCOU1(J+1)
        ELSE
       Y1=Y(J+1)
       Y2=Y(J+1)-UU*DTCONST
!            IF ( Y2 .GT. Y(N2) ) Y2=Y(N2)
        J1=JCOU2(J+1)
        J2=JCOU1(J+2)
        ENDIF
        Z1=ZP(K-1)
        Z2=ZP(K)

        I1=ICOU2(I)
        I2=ICOU1(I+1)
        K1=KCOUMP2(K-1)
        K2=KCOUMP1(K)

      !Z-DIR
      ELSE IF ( IDIR .EQ. 3) THEN

        X1=X(I)
        X2=X(I+1)
        Y1=Y(J)
        Y2=Y(J+1)
         IF ( UU .GE. 0. ) THEN
       Z1=ZP(K)-UU*DTCONST
       Z2=ZP(K)
!            IF ( Z1 .LT. Z(1) ) Z1=Z(1)
        K1=KCOUMP2(K-1)
        K2=KCOUMP1(K)
         ELSE
       Z1=ZP(K)
       Z2=ZP(K)-UU*DTCONST
!            IF ( Z2 .GT. Z(N3) ) Z2=Z(N3)
        K1=KCOUMP2(K)
        K2=KCOUMP1(K+1)
         ENDIF

        I1=ICOU2(I)
        I2=ICOU1(I+1)
        J1=JCOU2(J)
        J2=JCOU1(J+1)

      ENDIF

       CALL INTER_PSI_FACE(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,PSI_SUM)

        RETURN
        END

C=======================================================================
      SUBROUTINE PSIFACE_C(I,J,K,UU,PSI_SUM,IDIR)
C=======================================================================  
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE FLOW_VAR

      USE LVS_COUPLING
      
      IMPLICIT NONE
      
      INTEGER I,J,K,IDIR
      REAL UU,PSI_SUM
      
      INTEGER I1,I2,J1,J2,K1,K2
      REAL X1,X2,Y1,Y2,Z1,Z2
      
       !X-DIR
       IF ( IDIR .EQ. 1) THEN
       IF ( UU .GE. 0.) THEN
        X1=X(I+1)-UU*DTCONST
        X2=X(I+1)
!              IF ( X1 .LT. X(1) ) X1=X(1)
        I1=ICOU2(I)      !IF CFL IS GREATER THAN 1, THIS MAY HAVE PROBLEM.
        I2=ICOU1(I+1)
       ELSE
        X1=X(I+1)
        X2=X(I+1)-UU*DTCONST
!              IF ( X2 .GT. X(N1) ) X2=X(N1)
        I1=ICOU2(I+1)      !IF CFL IS GREATER THAN 1, THIS MAY HAVE PROBLEM.
        I2=ICOU1(I+2)
       ENDIF
       Y1=Y(J)
       Y2=Y(J+1)
       Z1=Z(K)
       Z2=Z(K+1)

        J1=JCOU2(J)
        J2=JCOU1(J+1)
        K1=KCOU2(K)
        K2=KCOU1(K+1)

      !Y-DIR
      ELSE IF ( IDIR .EQ. 2) THEN

       X1=X(I)
       X2=X(I+1)
       IF ( UU .GE. 0.) THEN
       Y1=Y(J+1)-UU*DTCONST    
       Y2=Y(J+1)
!            IF ( Y1 .LT. Y(1) ) Y1=Y(1)
        J1=JCOU2(J)
        J2=JCOU1(J+1)
       ELSE
       Y1=Y(J+1)
       Y2=Y(J+1)-UU*DTCONST
!            IF ( Y2 .GT. Y(N2) ) Y2=Y(N2)
        J1=JCOU2(J+1)
        J2=JCOU1(J+2)
       ENDIF
       Z1=Z(K)
       Z2=Z(K+1)

        I1=ICOU2(I)
        I2=ICOU1(I+1)
        K1=KCOU2(K)
        K2=KCOU1(K+1)

       ELSE IF ( IDIR .EQ. 3) THEN

       X1=X(I)
       X2=X(I+1)
       Y1=Y(J)
       Y2=Y(J+1)
         IF ( UU .GE. 0. ) THEN
       Z1=Z(K+1)-UU*DTCONST
       Z2=Z(K+1)
!            IF ( Z1 .LT. Z(1) ) Z1=Z(1)
        K1=KCOU2(K)
        K2=KCOU1(K+1)
         ELSE
       Z1=Z(K+1)
       Z2=Z(K+1)-UU*DTCONST
!            IF ( Z2 .GT. Z(N3) ) Z2=Z(N3)
        K1=KCOU2(K+1)
        K2=KCOU1(K+2)
         ENDIF

        I1=ICOU2(I)
        I2=ICOU1(I+1)
        J1=JCOU2(J)
        J2=JCOU1(J+1)

      ENDIF

       CALL INTER_PSI_FACE(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,PSI_SUM)

       RETURN
       END

C=======================================================================
      SUBROUTINE INTER_PSI_FACE(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2
     &,PSI_SUM)
C=======================================================================
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      IMPLICIT NONE
      
      INTEGER I1,I2,J1,J2,K1,K2
      REAL X1,X2,Y1,Y2,Z1,Z2
      REAL PSI_SUM
      
      INTEGER II,JJ,KK
      REAL VOLF
      REAL AI_CUT,AJ_CUT

      !SHOULD BE DONE GEOMETRICALLY.

          PSI_SUM=0.
          VOLF=0.

      DO II=I1,I2
       !IIIIIIII0
        IF ( (XF(II) .LE. X1).AND. (XF(II+1) .GE. X2 ) ) THEN  
       DO JJ=J1,J2
        !JJJJJJJJJ0
        IF ( YF(JJ) .LE. Y1 .AND. YF(JJ+1) .GE. Y2 ) THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)
              VOLF=VOLF+1.
              GOTO 40
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)
        VOLF=VOLF+1.
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 40
           ENDIF
         ENDDO
        !JJJJJJJJJ1
        ELSE IF ( (YF(JJ) .LT. Y1).AND. (YF(JJ+1) .GT. Y1 ) ) THEN
            AJ_CUT=(YF(JJ+1)-Y1)*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AJ_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AJ_CUT
              GOTO 40
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AJ_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 40
           ENDIF
         ENDDO
        !JJJJJJJJJ2
       ELSE IF( (YF(JJ) .GE. Y1).AND.(YF(JJ+1) .LE. Y2) )THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)
              VOLF=VOLF+1.
              GOTO 40
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)
        VOLF=VOLF+1.
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 40
           ENDIF
         ENDDO
        !JJJJJJJJJ3
       ELSE IF( (YF(JJ) .LT. Y2) .AND.((YF(JJ+1) .GT. Y2)) )THEN
            AJ_CUT=(Y2-YF(JJ))*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AJ_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AJ_CUT
              GOTO 40
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AJ_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 40
           ENDIF
         ENDDO
        !JJJJJJJJJ4
        ELSE IF( YF(JJ) .GE. Y2 ) THEN
         GOTO 100
       ENDIF
 40    CONTINUE
       ENDDO

       !IIIIIIII1
        ELSE IF ( (XF(II) .LT. X1).AND. (XF(II+1) .GT. X1 ) ) THEN  
            AI_CUT=(XF(II+1)-X1)*SSDXF
       DO JJ=J1,J2
        !JJJJJJJJJ0
        IF ( YF(JJ) .LE. Y1 .AND. YF(JJ+1) .GE. Y2 ) THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AI_CUT
              GOTO 50
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AI_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 50
           ENDIF
         ENDDO
        !JJJJJJJJJ1
        ELSE IF ( (YF(JJ) .LT. Y1).AND. (YF(JJ+1) .GT. Y1 ) ) THEN
            AJ_CUT=(YF(JJ+1)-Y1)*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AI_CUT*AJ_CUT
              GOTO 50 
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AI_CUT*AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 50
           ENDIF
         ENDDO
        !JJJJJJJJJ2
       ELSE IF( (YF(JJ) .GE. Y1).AND.(YF(JJ+1) .LE. Y2) )THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AI_CUT
              GOTO 50
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AI_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 50
           ENDIF
         ENDDO
        !JJJJJJJJJ3
       ELSE IF( (YF(JJ) .LT. Y2) .AND.((YF(JJ+1) .GT. Y2)) )THEN
            AJ_CUT=(Y2-YF(JJ))*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AI_CUT*AJ_CUT
              GOTO 50
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AI_CUT*AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 50
           ENDIF
         ENDDO
        !JJJJJJJJJ4
        ELSE IF( YF(JJ) .GE. Y2 ) THEN
         GOTO 100
       ENDIF
 50    CONTINUE
       ENDDO

       !IIIIIIII2
       ELSE IF( (XF(II) .GE. X1).AND.(XF(II+1) .LE. X2) )THEN
       DO JJ=J1,J2
        !JJJJJJJJJ0
        IF ( YF(JJ) .LE. Y1 .AND. YF(JJ+1) .GE. Y2 ) THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)
              VOLF=VOLF+1.
              GOTO 60
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)
        VOLF=VOLF+1.
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 60
           ENDIF
         ENDDO
        !JJJJJJJJJ1
       ELSE IF ( (YF(JJ) .LT. Y1).AND. (YF(JJ+1) .GT. Y1 ) ) THEN
            AJ_CUT=(YF(JJ+1)-Y1)*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AJ_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AJ_CUT
              GOTO 60
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AJ_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 60
           ENDIF     
         ENDDO
        !JJJJJJJJJ2
       ELSE IF( (YF(JJ) .GE. Y1).AND.(YF(JJ+1) .LE. Y2) )THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)
              VOLF=VOLF+1.
              GOTO 60
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)
        VOLF=VOLF+1.
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 60
           ENDIF  
       ENDDO

        !JJJJJJJJJ3
       ELSE IF( (YF(JJ) .LT. Y2) .AND.((YF(JJ+1) .GT. Y2)) )THEN
            AJ_CUT=(Y2-YF(JJ))*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AJ_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AJ_CUT
              GOTO 60
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AJ_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 60
           ENDIF
         ENDDO
        !JJJJJJJJJ4
           ELSE IF( YF(JJ) .GE. Y2 ) THEN
         GOTO 100
       ENDIF
 60    CONTINUE
       ENDDO

       !IIIIIIII3
       ELSE IF( (XF(II) .LT. X2) .AND.((XF(II+1) .GT. X2)) )THEN
            AI_CUT=(X2-XF(II))*SSDXF
       DO JJ=J1,J2
        !JJJJJJJJJ0
        IF ( YF(JJ) .LE. Y1 .AND. YF(JJ+1) .GE. Y2 ) THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AI_CUT
              GOTO 70
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AI_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 70
           ENDIF
         ENDDO
        !JJJJJJJJJ1
       ELSE IF ( (YF(JJ) .LT. Y1).AND. (YF(JJ+1) .GT. Y1 ) ) THEN
            AJ_CUT=(YF(JJ+1)-Y1)*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AI_CUT*AJ_CUT
              GOTO 70
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AI_CUT*AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 70
           ENDIF
         ENDDO        
        !JJJJJJJJJ2
       ELSE IF( (YF(JJ) .GE. Y1).AND.(YF(JJ+1) .LE. Y2) )THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AI_CUT
              GOTO 70
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AI_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 70
           ENDIF
         ENDDO
        !JJJJJJJJJ3
       ELSE IF( (YF(JJ) .LT. Y2) .AND.((YF(JJ+1) .GT. Y2)) )THEN
            AJ_CUT=(Y2-YF(JJ))*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AI_CUT*AJ_CUT
              GOTO 70
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AI_CUT*AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 70
           ENDIF
         ENDDO
        !JJJJJJJJJ4
           ELSE IF( YF(JJ) .GE. Y2 ) THEN
         GOTO 100
       ENDIF
 70    CONTINUE
       ENDDO

       !IIIIIIII4
       ELSE IF( XF(II) .GE. X2 ) THEN
         GOTO 101
       ENDIF
 100   ENDDO

 101   CONTINUE

         PSI_SUM=PSI_SUM/VOLF

       RETURN
       END

C*******************************************************************
      SUBROUTINE RHS_IBMX(RHS1,DENF_X,U,V,W)
C*******************************************************************
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE IBM_VAR
      
      IMPLICIT NONE

      REAL RHS1(M1M,M2M,M3M,3)

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL DENF_X(0:M1,0:M2,0:M3)

      INTEGER N,II,JJ,KK,I1,J1,K1,I2,J2,K2,L
      REAL UU1,UU2,B1,B2,B3
      REAL UBODY,UTARG,DENFU
      REAL DTI
      
C-----compute target velocities & forcing values at forcing points
      UBODY=0.
      DTI=1./DT
!$OMP PARALLEL DO private(II,JJ,KK,I1,J1,K1,I2,J2,K2)
!$OMP&private(UU1,UU2,B1,B2,B3,DENFU)
      DO N=1,NFC_INTP(1)
         II =IFC(N,1)
         JJ =JFC(N,1)
         KK =KFC(N,1)

         I1=INT(GFI(N,1,0,0,1))
         J1=INT(GFI(N,1,1,0,1))
         K1=INT(GFI(N,1,2,0,1))
         I2=INT(GFI(N,1,0,1,1))
         J2=INT(GFI(N,1,1,1,1))
         K2=INT(GFI(N,1,2,1,1))

      UU1=DENF_X(I1,J1,K1)*U(I1,J1,K1)
      UU2=DENF_X(I2,J2,K2)*U(I2,J2,K2)

       B1=GFI(N,1,0,0,0)*UU1+GFI(N,1,0,1,0)*UU2
       B2=GFI(N,1,1,0,0)*UU1+GFI(N,1,1,1,0)*UU2
       B3=GFI(N,1,2,0,0)*UU1+GFI(N,1,2,1,0)*UU2

       DENFU=B1+B2*X(II)+B3*YP(JJ)

       FCV(N,1)=(DENFU-RHS1(II,JJ,KK,1))*DTI
       RHS1(II,JJ,KK,1)=DENFU

!       aa=B1+B2*X(I1)+B3*YP(J1)
!       bb=B1+B2*X(I2)+B3*YP(J2)
!        if (abs(aa-10.) .ge. 1.e-5)  write(*,*) 'a',ii,jj,aa
!        if (abs(bb-10.) .ge. 1.e-5)  write(*,*) 'b',ii,jj,bb

      ENDDO
! 136  format(3i5,2es25.13)


!*************************       CAUTION       *************************
!     FOLLOWING 3 DO-LOOPS REQUIRES THE OPTION OF
!     '-Wf "-pvctl vwork=stack"'
!     WHEN THE CODE IS COMPILED ON NEC MACHINE.
!     FOR DETAILED EXPLANATION OF THE REASON,
!     YOU CAN CONSULT "NEC PORTING GUIDE"
!     PROVIDED BY NEC SUPERCOMPUTIONG CENTER.
!***********************************************************************

!$OMP PARALLEL DO private(II,JJ,KK,UTARG)
      DO N=NFC_INTP(1)+1,NFC_INTP(1)+NFC_INNER(1)
         II=IFC(N,1)
         JJ=JFC(N,1)
         KK=KFC(N,1)
          UTARG=UBODY
          FCV(N,1)=(UTARG-RHS1(II,JJ,KK,1))*DTI
          RHS1(II,JJ,KK,1)=UTARG
      ENDDO

C-----compute average forcing values during RK3 steps
      L=1
!$OMP PARALLEL DO
      DO N=1,NFC_INTP(L)+NFC_INNER(L)
         FCVAVG(N,L)=FCVAVG(N,L)+FCV(N,L)
      ENDDO

      RETURN
      END

C*******************************************************************
      SUBROUTINE RHS_IBMY(RHS1,DENF_Y,U,V,W)
C*******************************************************************
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE IBM_VAR
      
      IMPLICIT NONE

      REAL RHS1(M1M,M2M,M3M,3)

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL DENF_Y(0:M1,0:M2,0:M3)

      INTEGER N,II,JJ,KK,I1,J1,K1,I2,J2,K2,L
      REAL VV1,VV2,B1,B2,B3
      REAL VBODY,VTARG,DENFV
      REAL DTI
      
C-----compute target velocities & forcing values at forcing points
      VBODY=0.
      DTI=1./DT
!$OMP PARALLEL DO private(II,JJ,KK,I1,J1,K1,I2,J2,K2)
!$OMP&private(VV1,VV2,B1,B2,B3,DENFV)
      DO N=1,NFC_INTP(2)
         II =IFC(N,2)
         JJ =JFC(N,2)
         KK =KFC(N,2)

         I1=INT(GFI(N,2,0,0,1))
         J1=INT(GFI(N,2,1,0,1))
         K1=INT(GFI(N,2,2,0,1))
         I2=INT(GFI(N,2,0,1,1))
         J2=INT(GFI(N,2,1,1,1))
         K2=INT(GFI(N,2,2,1,1))

      VV1=DENF_Y(I1,J1,K1)*V(I1,J1,K1)
      VV2=DENF_Y(I2,J2,K2)*V(I2,J2,K2)

       B1=GFI(N,2,0,0,0)*VV1+GFI(N,2,0,1,0)*VV2
       B2=GFI(N,2,1,0,0)*VV1+GFI(N,2,1,1,0)*VV2
       B3=GFI(N,2,2,0,0)*VV1+GFI(N,2,2,1,0)*VV2

       DENFV=B1+B2*XP(II)+B3*Y(JJ)

       FCV(N,2)=(DENFV-RHS1(II,JJ,KK,2))*DTI
       RHS1(II,JJ,KK,2)=DENFV

!       aa=B1+B2*XP(I1)+B3*Y(J1)
!       bb=B1+B2*XP(I2)+B3*Y(J2)
!        if (abs(aa-10.) .ge. 1.e-5)  write(*,*) 'a',ii,jj,aa
!        if (abs(bb-10.) .ge. 1.e-5)  write(*,*) 'b',ii,jj,bb

      ENDDO

!$OMP PARALLEL DO private(II,JJ,KK,VTARG)
      DO N=NFC_INTP(2)+1,NFC_INTP(2)+NFC_INNER(2)
         II=IFC(N,2)
         JJ=JFC(N,2)
         KK=KFC(N,2)
          VTARG=VBODY
          FCV(N,2)=(VTARG-RHS1(II,JJ,KK,2))*DTI
          RHS1(II,JJ,KK,2)=VTARG
      ENDDO

C-----compute average forcing values during RK3 steps
      L=2
!$OMP PARALLEL DO
      DO N=1,NFC_INTP(L)+NFC_INNER(L)
         FCVAVG(N,L)=FCVAVG(N,L)+FCV(N,L)
      ENDDO

      RETURN
      END

C*******************************************************************
      SUBROUTINE RHS_IBMZ(RHS1,DENF_Z,U,V,W)
C*******************************************************************
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE IBM_VAR
      
      IMPLICIT NONE

      REAL RHS1(M1M,M2M,M3M,3)

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL DENF_Z(0:M1,0:M2,0:M3)
      
      INTEGER N,II,JJ,KK,I1,J1,K1,I2,J2,K2,L
      REAL WW1,WW2,B1,B2,B3
      REAL DENC,WBODY,WTARG,DENFW
      REAL DTI

      DENC=0.05*(DENP-DENM)

C-----compute target velocities & forcing values at forcing points
      DTI=1./DT
      WBODY=0.
!$OMP PARALLEL DO private(II,JJ,KK,I1,J1,K1,I2,J2,K2)
!$OMP&private(WW1,WW2,B1,B2,B3,DENFW)
      DO N=1,NFC_INTP(3)
         II =IFC(N,3)
         JJ =JFC(N,3)
         KK =KFC(N,3)

         I1=INT(GFI(N,3,0,0,1))
         J1=INT(GFI(N,3,1,0,1))
         K1=INT(GFI(N,3,2,0,1))
         I2=INT(GFI(N,3,0,1,1))
         J2=INT(GFI(N,3,1,1,1))
         K2=INT(GFI(N,3,2,1,1))

      IF ( ABS(DENF_Z(II,JJ,KK)-DENM) .LE. DENC   !DENC IS CRITERIA
     &  .OR. ABS(DENF_Z(II,JJ,KK)-DENP) .LE. DENC ) THEN
      WW1=DENF_Z(I1,J1,K1)*W(I1,J1,K1)
      WW2=DENF_Z(I2,J2,K2)*W(I2,J2,K2)

       B1=GFI(N,3,0,0,0)*WW1+GFI(N,3,0,1,0)*WW2
       B2=GFI(N,3,1,0,0)*WW1+GFI(N,3,1,1,0)*WW2
       B3=GFI(N,3,2,0,0)*WW1+GFI(N,3,2,1,0)*WW2

       DENFW=B1+B2*XP(II)+B3*YP(JJ)

       FCV(N,3)=(DENFW-RHS1(II,JJ,KK,3))*DTI
       RHS1(II,JJ,KK,3)=DENFW

       ELSE

       DENFW=DENF_Z(I1,J1,K1)*W(I1,J1,K1)
       FCV(N,3)=(DENFW-RHS1(II,JJ,KK,3))*DTI
       RHS1(II,JJ,KK,3)=DENFW
      ENDIF

!       aa=B1+B2*XP(I1)+B3*YP(J1)
!       bb=B1+B2*XP(I2)+B3*YP(J2)
!        if (abs(aa-6.) .ge. 1.e-5)  write(*,*) 'a',ii,jj,aa
!        if (abs(bb-6.) .ge. 1.e-5)  write(*,*) 'b',ii,jj,bb

      ENDDO

!$OMP PARALLEL DO private(II,JJ,KK,WTARG)
      DO N=NFC_INTP(3)+1,NFC_INTP(3)+NFC_INNER(3)
         II=IFC(N,3)
         JJ=JFC(N,3)
         KK=KFC(N,3)
          WTARG=WBODY
          FCV(N,3)=(WTARG-RHS1(II,JJ,KK,3))*DTI
          RHS1(II,JJ,KK,3)=WTARG
      ENDDO

C-----compute average forcing values during RK3 steps
      L=3
!$OMP PARALLEL DO
      DO N=1,NFC_INTP(L)+NFC_INNER(L)
         FCVAVG(N,L)=FCVAVG(N,L)+FCV(N,L)
      ENDDO

      RETURN
      END

C===========================CORRECTION_STEP============================C
C      TWO-PHASE PRESSURE-CORRECTION ALGORITHM                          C
C                                                                      C
C      CALCULATE CONTANT POISSON EQUATION ITERATIVELY INSTEAD OF        C
C     SOLVING VARIABLE POISSON EQUATION ( D.KIM STANFORD THESIS 2011)  C
C                                                                      C
C                                            KIYOUNG KIM 2012.8.10     C
C===========================CORRECTION_STEP============================C
C*******************************************************************
      SUBROUTINE DIVGS(U,V,W,DIVGSUM)
C*******************************************************************
C     POISSON EQUATION PRE-PROCESSING
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE IBM_VAR

      IMPLICIT NONE

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL DIVGSUM(M1M,M2M,M3M)

      INTEGER I,J,K,II,JJ,KK,IMM,JMM,KMM,N
      REAL DIVG1,DIVG2,DIVG3

!$OMP PARALLEL DO private(I,J,DIVG1,DIVG2,DIVG3)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
       DIVG1=(U(IPV(I),J,K)-U(I,J,K))*SSDX(I)
       DIVG2=(V(I,JPV(J),K)-V(I,J,K))*SSDY(J)
       DIVG3=(W(I,J,KPV(K))-W(I,J,K))*SSDZ(K)
       DIVGSUM(I,J,K)=(DIVG1+DIVG2+DIVG3)
      ENDDO
      ENDDO
      ENDDO

!$OMP PARALLEL DO private(I,J)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
       QMASS(I,J,K)=0.
      ENDDO
      ENDDO
      ENDDO

      if (masson.eq.1) then
!$OMP PARALLEL DO private(II,JJ,KK)
      DO N=1,NFC_INTP(1)+NFC_INNER(1)
        II=IFC(N,1)
        JJ=JFC(N,1)
        KK=KFC(N,1)
        QMASS(II,JJ,KK)=QMASS(II,JJ,KK)-U(II,JJ,KK)*SSDX(II)
        DIVGSUM(II,JJ,KK)=DIVGSUM(II,JJ,KK)+U(II,JJ,KK)*SSDX(II)
      ENDDO
!$OMP PARALLEL DO private(II,JJ,KK,IMM)
      DO N=1,NFC_INTP(1)+NFC_INNER(1)
        II=IFC(N,1)
        JJ=JFC(N,1)
        KK=KFC(N,1)
        IMM=IMV(II)
        QMASS(IMM,JJ,KK)=QMASS(IMM,JJ,KK)+U(II,JJ,KK)*SSDX(IMM)
        DIVGSUM(IMM,JJ,KK)=DIVGSUM(IMM,JJ,KK)-U(II,JJ,KK)*SSDX(IMM)
      ENDDO
!$OMP PARALLEL DO private(II,JJ,KK)
      DO N=1,NFC_INTP(2)+NFC_INNER(2)
        II=IFC(N,2)
        JJ=JFC(N,2)
        KK=KFC(N,2)
        QMASS(II,JJ,KK)=QMASS(II,JJ,KK)-V(II,JJ,KK)*SSDY(JJ)
        DIVGSUM(II,JJ,KK)=DIVGSUM(II,JJ,KK)+V(II,JJ,KK)*SSDY(JJ)
      ENDDO
!$OMP PARALLEL DO private(II,JJ,KK,JMM)
      DO N=1,NFC_INTP(2)+NFC_INNER(2)
        II=IFC(N,2)
        JJ=JFC(N,2)
        KK=KFC(N,2)
        JMM=JMV(JJ)
        QMASS(II,JMM,KK)=QMASS(II,JMM,KK)+V(II,JJ,KK)*SSDY(JMM)
        DIVGSUM(II,JMM,KK)=DIVGSUM(II,JMM,KK)-V(II,JJ,KK)*SSDY(JMM)
      ENDDO
!$OMP PARALLEL DO private(II,JJ,KK)
      DO N=1,NFC_INTP(3)+NFC_INNER(3)
        II=IFC(N,3)
        JJ=JFC(N,3)
        KK=KFC(N,3)
        QMASS(II,JJ,KK)=QMASS(II,JJ,KK)-W(II,JJ,KK)*SSDZ(KK)
        DIVGSUM(II,JJ,KK)=DIVGSUM(II,JJ,KK)+W(II,JJ,KK)*SSDZ(KK)
      ENDDO
!$OMP PARALLEL DO private(II,JJ,KK,KMM)
      DO N=1,NFC_INTP(3)+NFC_INNER(3)
        II=IFC(N,3)
        JJ=JFC(N,3)
        KK=KFC(N,3)
        KMM=KMV(KK)
        QMASS(II,JJ,KMM)=QMASS(II,JJ,KMM)+W(II,JJ,KK)*SSDZ(KMM)
        DIVGSUM(II,JJ,KMM)=DIVGSUM(II,JJ,KMM)-W(II,JJ,KK)*SSDZ(KMM)
      ENDDO

      endif
 
!$OMP PARALLEL DO private(I,J)
      DO 80 K=1,N3M
      DO 80 J=1,N2M
      DO 80 I=1,N1M
      DIVGSUM(I,J,K)=DIVGSUM(I,J,K)*DTCONSTI
   80 CONTINUE
 
!!      IF ( MSUB .EQ. 1 ) THEN
!C=====SAVE
!       OPEN(143,FILE='0DIVG.DAT')
!      WRITE(143,*) 'VARIABLES="X","Y","Z","DIVGX","DIVGY","DIVGZ","SUM"'
!      WRITE(143,*) 'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'
!       DO K=1,N3M
!       DO J=1,N2M
!       DO I=1,N1M
!       DIV1=(U(IPV(I),J,K)-U(I,J,K))*SSDX(I)
!       DIV2=(V(I,JPV(J),K)-V(I,J,K))*SSDY(J)
!       DIV3=(W(I,J,KPV(K))-W(I,J,K))*SSDZ(K)
!       WRITE(143,149) XP(I),YP(J),ZP(K),DIV1,DIV2,DIV3,DIV1+DIV2+DIV3
!        ENDDO
!        ENDDO
!        ENDDO
!       CLOSE(143)
! 149  FORMAT(7F20.10)
!!      ENDIF

      RETURN
      END

C*******************************************************************
      SUBROUTINE CONVRGE1(DVMAX,U,V,W)
C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      USE IBM_VAR
      
      IMPLICIT NONE
      
      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL DVMAX
      
      INTEGER I,J,K
      REAL DVG11,DVG12,DVG13,DIVGMAX_TMP

      DVMAX=0.
!$OMP PARALLEL DO private(I,J,DVG11,DVG12,DVG13,DIVGMAX_TMP)
!$OMP&reduction(MAX:DVMAX)
      DO 13 K=1,N3M
      DO 13 J=1,N2M
      DO 13 I=1,N1M
       DVG11=(U(IPV(I),J,K)-U(I,J,K))*SSDX(I)
       DVG12=(V(I,JPV(J),K)-V(I,J,K))*SSDY(J)
       DVG13=(W(I,J,KPV(K))-W(I,J,K))*SSDZ(K)
       DIVGMAX_TMP=ABS(DVG11+DVG12+DVG13-QMASS(I,J,K))  !check qmass..
       DVMAX=MAX(DIVGMAX_TMP,DVMAX)
!        if (dvmax .eq. divgmax_tmp) then
!          write(*,*) i,j,k,divgmax_tmp
!        endif
   13 CONTINUE

       RETURN
       END

C******************************************************************
      SUBROUTINE MEANPGR(W,PSI_ZN,PSI_CN,DENF_Z,VOL1)
C******************************************************************
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE IBM_VAR
      
       IMPLICIT NONE

      REAL W(0:M1,0:M2,0:M3)

      REAL PSI_ZN(0:M1,0:M2,0:M3)
      REAL PSI_CN(0:M1,0:M2,0:M3)

      REAL DENF_Z(0:M1,0:M2,0:M3)

      REAL VOL1
      
      INTEGER I,J,K,N
      REAL PMI1,PMI2,PMI3,PMI4,PMI5
      REAL VISU,VISL,WGU1,WGL1,WGU2,WGL2
      REAL DEN_SUM,DEN_ZN
      REAL PI
      REAL FUNCBODY

      PMI=0.
      PMI1=0.
      PMI2=0.
      PI=ACOS(-1.)

      !CAL PMI1,PMI2
!$OMP PARALLEL DO private(J,VISU,VISL,WGU1,WGL1)
!$OMP&reduction(+:PMI1)
      DO K=1,N3M
      DO J=1,N2M
      VISU=VISM+VIS_DIFF*PSI_CN(N1M,J,K)
      VISL=VISM+VIS_DIFF*PSI_CN(1,J,K)
!      WGU1=VISU*(W(N1,J,K)-W(N1M,J,K))*VVDX(N1)*VDZ(K)*SDY(J)
!      WGL1=VISL*(W(1,J,K)-W(0,J,K))*VVDX(1)*VDZ(K)*SDY(J)
      WGU1=VISU*(W(IPV(N1M),J,K)-W(N1M,J,K))*VVDX(N1)*VDZ(K)*SDY(J) !PERIODIC CONDITION->ZERO
      WGL1=VISL*(W(1,J,K)-W(IMV(1),J,K))*VVDX(1)*VDZ(K)*SDY(J)
      PMI1=PMI1+(WGU1-WGL1)
      ENDDO
      ENDDO

!$OMP PARALLEL DO private(I,VISU,VISL,WGU2,WGL2)
!$OMP&reduction(+:PMI2)
      DO K=1,N3M
      DO I=1,N1M
      VISU=VISM+VIS_DIFF*PSI_CN(I,N2M,K)
      VISL=VISM+VIS_DIFF*PSI_CN(I,1,K)
!      WGU2=VISU*(W(I,N2,K)-W(I,N2M,K))*VVDY(N2)*VDZ(K)*SDX(I)
!      WGL2=VISL*(W(I,1,K)-W(I,0,K))*VVDY(1)*VDZ(K)*SDX(I)
      WGU2=VISU*(W(I,JPV(N2M),K)-W(I,N2M,K))*VVDY(N2)*VDZ(K)*SDX(I) !PERIODIC CONDITION->ZERO
      WGL2=VISL*(W(I,1,K)-W(I,JMV(1),K))*VVDY(1)*VDZ(K)*SDX(I)
      PMI2=PMI2+(WGU2-WGL2)
      ENDDO
      ENDDO

      !CAL PMI3
      DEN_SUM=0.
!$OMP PARALLEL DO private(I,J,DEN_ZN)
!$OMP&reduction(+:DEN_SUM)
      DO K=KBG,N3M
      DO J=1,N2M
      DO I=1,N1M
       IF ( FUNCBODY(XP(I),YP(J),Z(K)) .GT. 0. ) THEN
        DEN_ZN=DENM+DEN_DIFF*PSI_ZN(I,J,K)
        DEN_SUM=DEN_SUM+0.5*(DEN_ZN+DENF_Z(I,J,K))*SDX(I)*SDY(J)*VDZ(K)
       ENDIF
      ENDDO
      ENDDO
      ENDDO

      DEN_SUM=DEN_SUM/VOL1
      PMI_GRAV=-DEN_SUM*GRAVITY
!        PMI_GRAV=-DENP*GRAVITY

      PMI4=0.
      IF (IBMON.NE.0) THEN
      !CAL PMI4
!$OMP PARALLEL DO
!$OMP&reduction(+:PMI4)
      DO 40 N=1,NFC_INTP(3)+NFC_INNER(3)
      PMI4=PMI4+FCV(N,3)*SDX(IFC(N,3))*SDY(JFC(N,3))*VDZ(KFC(N,3))
 40   CONTINUE
        IF (ALPHA_RK3 .NE. 0.) PMI4=PMI4/(2.*ALPHA_RK3)       !AT INITIAL STEP, ALPHA_RK3=0.
      ENDIF

      IF (ICH .EQ. 1 ) THEN
       PMI=(PMI1+PMI2+PMI4)/VOL1 !+(DEN_SUM-DENP)*GRAVITY
      ELSE
       PMI=PMI_CONST+(DEN_SUM-DENP)*GRAVITY
      ENDIF


      IF (MSUB.EQ.3 .AND. IBMON.NE.0) THEN !FOR SAVE
       PMI5=0.
      DO 50 N=1,NFC_INTP(3)+NFC_INNER(3)
       PMI5=PMI5+FCVAVG(N,3)*SDX(IFC(N,3))*SDY(JFC(N,3))*VDZ(KFC(N,3))
 50   CONTINUE

      ENDIF

       PMI_DUDY=(PMI1+PMI2+PMI5)/VOL1

      RETURN
      END

C******************************************************************
      SUBROUTINE QVOLCALC(QM1,II,W,PSI_ZN,DENF_Z,VOL1)
C******************************************************************
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE

      REAL W(0:M1,0:M2,0:M3)

      REAL PSI_ZN(0:M1,0:M2,0:M3)
      REAL DENF_Z(0:M1,0:M2,0:M3)

      REAL VOL1
      INTEGER II
      
      INTEGER I,J,K
      REAL QM_TOT,QM_BUB,QM_WAT,QM1,DVOL
      REAL DEN,DEN_DIFFI,PSI
      REAL FUNCBODY

      !ATTENTION
      !HERE MASS FLUX IS NOT DENF_Z(I,J,K)*W(I,J,K)*PSI, BUT DENP*W(I,J,K)*PSI.
      !DENF_Z(I,J,K) AND PSI SHOULD NOT BE USED AT THE SAME TIME. (VOID FRACTION EFFECTS ARE DOUBLED)

      IF ( II .EQ. 1 ) THEN

      QM_TOT=0.
!$OMP PARALLEL DO private(I,J,DEN)
!$OMP&reduction(+:QM_TOT)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
       IF ( FUNCBODY(XP(I),YP(J),Z(K)) .GT. 0. ) THEN

        !whole domain
        DEN=DENM+DEN_DIFF*PSI_ZN(I,J,K)
        QM_TOT=QM_TOT+DEN*W(I,J,K)*SDX(I)*SDY(J)*VDZ(K)

!        !water only
!        IF(PSI_ZN(I,J,K) .NE. 0.) THEN
!         QM1=QM1+DENP*W(I,J,K)*SDX(I)*SDY(J)*VDZ(K)*PSI_ZN(I,J,K) 
!        ENDIF

       ENDIF
      ENDDO
      ENDDO
      ENDDO
       QM1=QM_TOT/VOL1

      ELSE IF( II .EQ. 2 ) THEN

      QM_TOT=0.
      QM_BUB=0.
      QM_WAT=0.
!$OMP PARALLEL DO private(I,J,DVOL)
!$OMP&reduction(+:QM_TOT)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
       IF ( FUNCBODY(XP(I),YP(J),Z(K)) .GT. 0. ) THEN
         DVOL=SDX(I)*SDY(J)*VDZ(K)
        QM_TOT=QM_TOT+DENF_Z(I,J,K)*W(I,J,K)*DVOL

!        PSI=(DENF_Z(I,J,K)-DENM)*DEN_DIFFI
!        QM_BUB=QM_BUB+DENM*W(I,J,K)*DVOL*(1.-PSI)
!         QM_WAT=QM_WAT+DENP*W(I,J,K)*DVOL*PSI
       ENDIF
      ENDDO
      ENDDO
      ENDDO

        QM_TOT=QM_TOT/VOL1
!        QM_BUB=QM_BUB/VOL1
!        QM_WAT=QM_WAT/VOL1

        QM1=QM_TOT !whole domain   
       !QM1=QM_WAT !water only

      ELSE IF( II .EQ. 3 ) THEN

      IF (DENR .EQ. 1.) THEN
      DEN_DIFFI=0.
      ELSE
      DEN_DIFFI=1./DEN_DIFF
      ENDIF
      QM_TOT=0.
      QM_BUB=0.
      QM_WAT=0.
!$OMP PARALLEL DO private(I,J,DVOL,PSI)
!$OMP&reduction(+:QM_TOT,QM_BUB,QM_WAT)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
       IF ( FUNCBODY(XP(I),YP(J),Z(K)) .GT. 0. ) THEN
         DVOL=SDX(I)*SDY(J)*VDZ(K)
        QM_TOT=QM_TOT+DENF_Z(I,J,K)*W(I,J,K)*DVOL

        PSI=(DENF_Z(I,J,K)-DENM)*DEN_DIFFI
        QM_BUB=QM_BUB+DENM*W(I,J,K)*DVOL*(1.-PSI)
         QM_WAT=QM_WAT+DENP*W(I,J,K)*DVOL*PSI
       ENDIF
      ENDDO
      ENDDO
      ENDDO

        QM_TOT=QM_TOT/VOL1
        QM_BUB=QM_BUB/VOL1
        QM_WAT=QM_WAT/VOL1

        QM1=QM_TOT !whole domain   
       !QM1=QM_WAT !water only

      IF (MSUB .EQ. 3) THEN
       OPEN(100,FILE='0FLUX.DAT',POSITION='APPEND')
        WRITE(100,140) TIME,QM_BUB,QM_WAT,QM_TOT
 140    FORMAT(F13.5,4ES20.12)
       CLOSE(100)
      ENDIF

      ENDIF

      RETURN
      END
 
C =====================================================================      
      SUBROUTINE FSM_CHOI(DENF_XI,DENF_YI,DENF_ZI,U,V,W,P,PSI_ZN,DENF_Z
     &,QVOL_ORI,VOL1)
C ===================================================================== 
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      REAL PSI_ZN(0:M1,0:M2,0:M3)
      REAL DENF_Z(0:M1,0:M2,0:M3)
 
      REAL DENF_XI(0:M1,0:M2,0:M3),DENF_YI(0:M1,0:M2,0:M3)
     &                                         ,DENF_ZI(0:M1,0:M2,0:M3) 

      REAL VOL1
      REAL QVOL_ORI
      
      INTEGER I,J,K,JM1,KM1
      REAL QVOLH,PHCAP
      REAL FUNCBODY

      IF(ICH.EQ.1) THEN 
       CALL QVOLCALC(QVOLH,2,W,PSI_ZN,DENF_Z,VOL1)
       PHCAP=(QVOLH-QVOL_ORI)/VOL1*DTCONSTI
      ELSE
       PHCAP=0.
      ENDIF

      !FRACTIONAL STEP METHOD OF CHOI ET AL(1994)
!$OMP PARALLEL DO private(I,J)
      DO 22 K=1,N3M
      DO 22 J=1,N2M
      DO 22 I=IBG,N1M
         U(I,J,K)=U(I,J,K)
     &         +DTCONST*(P(I,J,K)-P(IMV(I),J,K))*VVDX(I)*DENF_XI(I,J,K)
   22 CONTINUE

!$OMP PARALLEL DO private(I,J,JM1)
      DO 32 K=1,N3M
      DO 32 J=JBG,N2M
        JM1=JMV(J)
      DO 32 I=1,N1M
         V(I,J,K)=V(I,J,K)
     &         +DTCONST*(P(I,J,K)-P(I,JM1,K))*VVDY(J)*DENF_YI(I,J,K)
   32 CONTINUE

       IF ( N3M .NE. 1 ) THEN
!$OMP PARALLEL DO private(I,J,KM1)
      DO 42 K=KBG,N3M
         KM1=KMV(K)
      DO 42 J=1,N2M
      DO 42 I=1,N1M
       IF ( FUNCBODY(XP(I),YP(J),Z(K)) .GT. 0. ) THEN
         W(I,J,K)=W(I,J,K)
     &   +DTCONST*( (P(I,J,K)-P(I,J,KM1))*VVDZ(K)-PHCAP )*DENF_ZI(I,J,K)
       ELSE
         W(I,J,K)=W(I,J,K)
     &       +DTCONST*(P(I,J,K)-P(I,J,KM1))*VVDZ(K)*DENF_ZI(I,J,K)
       ENDIF
   42  CONTINUE
        ENDIF

       RETURN
       END

C =====================================================================      
      SUBROUTINE PRESSURE_CORRECTION(IPOISS,U,V,W,P,DENF_Y,DENF_C
     &,DENF_XI,DENF_YI,DENF_ZI)
C =====================================================================        
      USE PARAM_VAR
      
      USE FLOW_GEOM_VAR
      USE FLOW_VAR
      
      IMPLICIT NONE
      
      REAL DENF_XI(0:M1,0:M2,0:M3),DENF_YI(0:M1,0:M2,0:M3)
     &                                         ,DENF_ZI(0:M1,0:M2,0:M3)

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      REAL DIVGSUM(M1M,M2M,M3M)

      REAL DENF_Y(0:M1,0:M2,0:M3)
      REAL DENF_C(0:M1,0:M2,0:M3)
      
      INTEGER IPOISS
      REAL FUNCBODY

      INTEGER I,J,K
      REAL PHIREF,DVOL,VOL
 
       CALL DIVGS(U(0,0,0),V(0,0,0),W(0,0,0),DIVGSUM(1,1,1))
 
      IF (IPOISS .EQ. 0) THEN  !1D FFT, 1D MG, 1D TDMA
         write(*,*) 'compile with "poiss_ft_ffte_omp_novec.f"'
         write(*,*) 'CALL_POISINIT IS NEEDED."'
         write(*,*) 'z: periodic. x,y: non-periodic is only possible"'
         stop
!       CALL POISSON(P,DIVGSUM)
      ELSE IF (IPOISS .EQ. 1) THEN !PRECONDITIONED CONJUGATE GRADIENT
       CALL POISSON_PCG(P,DIVGSUM,DENF_XI,DENF_YI,DENF_ZI,DENF_Y,DENF_C)
      ENDIF

       CALL UCALC(DENF_XI,DENF_YI,DENF_ZI,U,V,W,P)

C     SET THE AVERAGE PHI AT THE UPPER WALL TO BE ZERO.
      PHIREF=0.
      VOL=0.
!$OMP PARALLEL DO private(I,J,DVOL)
!$OMP&reduction(+:PHIREF,VOL)
      DO 80 K=1,N3M
      DO 80 J=1,N2M
      DO 80 I=1,N1M
       IF ( FUNCBODY(XP(I),YP(J),ZP(K)) .GT. 0. ) THEN
        DVOL=SDX(I)*SDY(J)*SDZ(K)
        PHIREF=PHIREF+P(I,J,K)*DVOL
        VOL=VOL+DVOL
       ENDIF
   80 CONTINUE
      PHIREF=PHIREF/VOL

!$OMP PARALLEL DO private(I,J)
      DO 93 K=1,N3M
      DO 93 J=1,N2M
      DO 93 I=1,N1M
        P(I,J,K)=P(I,J,K)-PHIREF
   93 CONTINUE

      IF (IPX .NE. 1) THEN
!$OMP PARALLEL DO private(J)
      DO K=1,N3M
      DO J=1,N2M
         P(0,J,K)=P(1,J,K)
      ENDDO
      ENDDO
      ENDIF
      IF (IPY .NE. 1) THEN
!$OMP PARALLEL DO private(I)
      DO K=1,N3M
      DO I=1,N1M
         P(I,0,K)=P(I,1,K)
      ENDDO
      ENDDO
      ENDIF
      IF (IPZ .NE. 1) THEN
!$OMP PARALLEL DO private(I)
      DO J=1,N2M
      DO I=1,N1M
         P(I,J,0) =P(I,J,1)
      ENDDO
      ENDDO
      ENDIF

       RETURN
       END

C*******************************************************************
      SUBROUTINE UCALC(DENF_XI,DENF_YI,DENF_ZI,U,V,W,P)
C*******************************************************************
C     CALCULATING U FROM UHAT
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      REAL DENF_XI(0:M1,0:M2,0:M3),DENF_YI(0:M1,0:M2,0:M3)
     &                                         ,DENF_ZI(0:M1,0:M2,0:M3)
     
      INTEGER I,J,K,JM1,KM1
      REAL FUNCBODY

C-----CORRECTION
!$OMP PARALLEL DO private(I,J)
      DO 22 K=1,N3M
      DO 22 J=1,N2M
      DO 22 I=IBG,N1M
         U(I,J,K)=U(I,J,K)
     &      -DTCONST*(P(I,J,K)-P(IMV(I),J,K))*VVDX(I)*DENF_XI(I,J,K)
   22 CONTINUE

!$OMP PARALLEL DO private(I,J,JM1)
      DO 32 K=1,N3M
      DO 32 J=JBG,N2M
        JM1=JMV(J)
      DO 32 I=1,N1M
         V(I,J,K)=V(I,J,K)
     &         -DTCONST*(P(I,J,K)-P(I,JM1,K))*VVDY(J)*DENF_YI(I,J,K)
   32 CONTINUE

       IF ( N3M .NE. 1 ) THEN
!$OMP PARALLEL DO private(I,J,KM1)
      DO 42 K=KBG,N3M
         KM1=KMV(K)
      DO 42 J=1,N2M
      DO 42 I=1,N1M
!       IF ( FUNCBODY(XP(I),YP(J),Z(K)) .GT. 0. ) THEN
!         W(I,J,K)=W(I,J,K)
!     &           -( DTCONST*( (P(I,J,K)-P(I,J,KM1))*VVDZ(K) )
!     &           +PHCAP )*DENF_ZI(I,J,K)
!       ELSE
         !!THIS CAN GENERATE LOCAL RESIDE ERROR, CAUSE PHCAP(CONSTANT) IS NOT ADDED ALL DOMAIN
         W(I,J,K)=W(I,J,K)
     &        -DTCONST*(P(I,J,K)-P(I,J,KM1))*VVDZ(K)*DENF_ZI(I,J,K)
!       ENDIF
   42  CONTINUE

        ENDIF

      !When Neumann condition is used, p(1,j,k)-p(0,j,k)=0
      !boundary condition like U(1,J,K)=U(2,J,K) is not needed!!
      !U(1,J,K)=U(2,J,K) IS UADATED WITH INTERMEDIATE VELOCITY.
!BOUNDARY_CONDITION
        IF (IPX .EQ. 0) THEN
!$OMP PARALLEL DO private(J)
          DO K=1,N3M
          DO J=1,N2M
          V(0,J,K)=V(1,J,K)
          W(0,J,K)=W(1,J,K)
          V(N1,J,K)=V(N1M,J,K)
          W(N1,J,K)=W(N1M,J,K)
         ENDDO
         ENDDO
        ENDIF
        IF (IPY .EQ. 0) THEN
!$OMP PARALLEL DO private(I)
          DO K=1,N3M
          DO I=1,N1M
          U(I,0,K)=U(I,1,K)
          W(I,0,K)=W(I,1,K)
          U(I,N2,K)=U(I,N2M,K)
          W(I,N2,K)=W(I,N2M,K) 
         ENDDO
         ENDDO
        ENDIF
        IF (IPZ .EQ. 0) THEN
!$OMP PARALLEL DO private(I)
          DO J=1,N2M
          DO I=1,N1M
          U(I,J,0)=U(I,J,1)
          V(I,J,0)=V(I,J,1)
          U(I,J,N3)=U(I,J,N3M)
          V(I,J,N3)=V(I,J,N3M)
         ENDDO
         ENDDO
        ENDIF

      RETURN
      END


C*******************************************************************
      SUBROUTINE IBMINIT
C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE IBM_VAR
      
      IMPLICIT NONE
      
      INTEGER I,J,K,N,L,LL

      ALLOCATE(NFC_INTP(3),NFC_INNER(3))
      ALLOCATE(IFC(MBODY,3),JFC(MBODY,3),KFC(MBODY,3))
      ALLOCATE(FCV(MBODY,3),FCVAVG(MBODY,3))
      ALLOCATE(MPI(MBODY0,3,3))
      ALLOCATE(GFI(MBODY0,3,0:2,0:2,0:2))
      ALLOCATE(FORCESUM(10))

       DO L=1,3
        NFC_INTP(L)=0 ;NFC_INNER(L)=0
       ENDDO
      
       DO L=1,3
!$OMP PARALLEL DO
        DO LL=1,MBODY
         IFC(LL,L)=0 ;JFC(LL,L)=0 ;KFC(LL,L)=0
         FCV(LL,L)=0. ;FCVAVG(LL,L)=0.
        ENDDO

!$OMP PARALLEL DO
        DO LL=1,MBODY0
         MPI(LL,L,1)=0 ;MPI(LL,L,2)=0 ;MPI(LL,L,3)=0
        ENDDO

!$OMP PARALLEL DO private(I,J,K)
        DO LL=1,MBODY0
         DO K=0,2 ;DO J=0,2 ;DO I=0,2
          GFI(LL,L,I,J,K)=0.
         ENDDO ;ENDDO ;ENDDO
        ENDDO
       ENDDO
      
       DO L=1,10
        FORCESUM(L)=0.
       ENDDO


      OPEN(11,FILE='ibmpre_intp.bin',FORM='UNFORMATTED')
      READ(11) NFC_INTP(1),NFC_INTP(2),NFC_INTP(3)
      READ(11) NFC_INNER(1),NFC_INNER(2),NFC_INNER(3)
      READ(11) (IFC(N,1),N=1,NFC_INTP(1)+NFC_INNER(1))
      READ(11) (IFC(N,2),N=1,NFC_INTP(2)+NFC_INNER(2))
      READ(11) (IFC(N,3),N=1,NFC_INTP(3)+NFC_INNER(3))
      READ(11) (JFC(N,1),N=1,NFC_INTP(1)+NFC_INNER(1))
      READ(11) (JFC(N,2),N=1,NFC_INTP(2)+NFC_INNER(2))
      READ(11) (JFC(N,3),N=1,NFC_INTP(3)+NFC_INNER(3))
      READ(11) (KFC(N,1),N=1,NFC_INTP(1)+NFC_INNER(1))
      READ(11) (KFC(N,2),N=1,NFC_INTP(2)+NFC_INNER(2))
      READ(11) (KFC(N,3),N=1,NFC_INTP(3)+NFC_INNER(3))
      READ(11) ((MPI(N,1,L),N=1,NFC_INTP(1)),L=1,3)
      READ(11) ((MPI(N,2,L),N=1,NFC_INTP(2)),L=1,3)
      READ(11) ((MPI(N,3,L),N=1,NFC_INTP(3)),L=1,3)
      READ(11) ((((GFI(N,1,I,J,K),N=1,NFC_INTP(1))
     &                                           ,I=0,2),J=0,2),K=0,2)
      READ(11) ((((GFI(N,2,I,J,K),N=1,NFC_INTP(2))
     &                                           ,I=0,2),J=0,2),K=0,2)
      READ(11) ((((GFI(N,3,I,J,K),N=1,NFC_INTP(3))
     &                                           ,I=0,2),J=0,2),K=0,2)
      CLOSE(11)

      WRITE(*,35) NFC_INTP(1),NFC_INTP(2),NFC_INTP(3)
      WRITE(*,36) NFC_INNER(1),NFC_INNER(2),NFC_INNER(3)
 35   FORMAT('FU_intp :',I12,'  FV_intp :',I12,'  FW_intp :',I12)
 36   FORMAT('FU_inner:',I12,'  FV_inner:',I12,'  FW_inner:',I12)

      FCV=0.
      FCVAVG=0.

      RETURN
      END



C******************************************************************
C     SECTION FOR SOLVING UNCOUPLED MOMENTUM EQS. TO GET UHAT
C     - LHS-EI.F
C     - FOR SPECIFIC FEATURES, SEE 'rhsnlhs-ei.f'
C     - THIS IS FASTER THAN 'LHS FOR LES' BECAUSE THIS DOESN'T
C       CALCULATE COEFS. EVERYTIME
C     - CAN HANDLE OPEN( 2 TYPES ) OR CLOSED TOP ACCORDING TO ITOPEN.
C
C                                  LAST REVISED BY SEONGWON KANG
C                                            1998.07.06.   13:00
C******************************************************************

C*******************************************************************
      SUBROUTINE TRDIAG1(A,B,C,R,UU,L1,L2,LL1,LL2)
C*******************************************************************
C     SOLVE THE TRIDIAGONAL MATRIX (X-1 TRIDIAGONAL) WITH
C     DIFFERENT COEFFICIENTS.

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      DIMENSION GAM(M2,M1),A(M2,M1),B(M2,M1),C(M2,M1),
     &          R(M2,M1),UU(M2,M1),BET(M2)

      DO 10 I=LL1,LL2
      BET(I)=1./B(I,L1)
      UU(I,L1)=R(I,L1)*BET(I)
   10 CONTINUE

      DO 20 J=L1+1,L2
      DO 20 I=LL1,LL2
      GAM(I,J)=C(I,J-1)*BET(I)
      BET(I)=1./(B(I,J)-A(I,J)*GAM(I,J))
      UU(I,J)=(R(I,J)-A(I,J)*UU(I,J-1))*BET(I)
   20 CONTINUE
      DO 30 J=L2-1,L1,-1
      DO 30 I=LL1,LL2
      UU(I,J)=UU(I,J)-GAM(I,J+1)*UU(I,J+1)
   30 CONTINUE

      RETURN
      END

C*******************************************************************
      SUBROUTINE TRDIAG2(A,B,C,R,UU,L1,L2,LL1,LL2)
C*******************************************************************
C     SOLVE THE TRIDIAGONAL MATRIX (X-2 TRIDIAGONAL) WITH
C     DIFFERENT COEFFICIENTS.

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      DIMENSION GAM(M1,M2),A(M1,M2),B(M1,M2),C(M1,M2),
     &          R(M1,M2),UU(M1,M2),BET(M1)

      DO 10 I=LL1,LL2
      BET(I)=1./B(I,L1)
      UU(I,L1)=R(I,L1)*BET(I)
   10 CONTINUE

      DO 20 J=L1+1,L2
      DO 20 I=LL1,LL2
      GAM(I,J)=C(I,J-1)*BET(I)
      BET(I)=1./(B(I,J)-A(I,J)*GAM(I,J))
      UU(I,J)=(R(I,J)-A(I,J)*UU(I,J-1))*BET(I)
   20 CONTINUE
      DO 30 J=L2-1,L1,-1
      DO 30 I=LL1,LL2
      UU(I,J)=UU(I,J)-GAM(I,J+1)*UU(I,J+1)
   30 CONTINUE

      RETURN
      END

C*******************************************************************
      SUBROUTINE TRDIAG3(A,B,C,R,UU,L1,L2,LL1,LL2)
C*******************************************************************
C     SOLVE THE TRIDIAGONAL MATRIX (X-3 TRIDIAGONAL) WITH
C     DIFFERENT COEFFICIENTS.

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      DIMENSION GAM(M1,M3),A(M1,M3),B(M1,M3),C(M1,M3),
     &          R(M1,M3),UU(M1,M3),BET(M1)

      DO 10 I=LL1,LL2
      BET(I)=1./B(I,L1)
      UU(I,L1)=R(I,L1)*BET(I)
   10 CONTINUE

      DO 20 J=L1+1,L2
      DO 20 I=LL1,LL2
      GAM(I,J)=C(I,J-1)*BET(I)
      BET(I)=1./(B(I,J)-A(I,J)*GAM(I,J))
      UU(I,J)=(R(I,J)-A(I,J)*UU(I,J-1))*BET(I)
   20 CONTINUE
      DO 30 J=L2-1,L1,-1
      DO 30 I=LL1,LL2
      UU(I,J)=UU(I,J)-GAM(I,J+1)*UU(I,J+1)
   30 CONTINUE

      RETURN
      END

C*******************************************************************
      SUBROUTINE TRDIAG1P(A,B,C,F,J1,J2,L1,L2)
C*******************************************************************
C     INVERT A PERIODIC MATRIX (X-3 TRIDIAGONAL) WITH
C     DIFFERENT COEFFICIENTS.

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      REAL A(M2,M1),B(M2,M1),C(M2,M1),F(M2,M1),Q(M2,M1),
     %     S(M2,M1),QE(M2,M1),FN(M2),PN(M2)

      JA=J1+1
      JJ=J1+J2
      DO 10 K=L1,L2
      BINV=1./B(K,J1)
      Q(K,J1)=-C(K,J1)*BINV
      S(K,J1)=-A(K,J1)*BINV
      FN(K)=F(K,J2)
      F(K,J1)=F(K,J1)*BINV
   10 CONTINUE

C     FORWARD ELIMINATION SWEEP
      DO 20 J=JA,J2
      DO 20 K=L1,L2
      PN(K)=1./(B(K,J)+A(K,J)*Q(K,J-1))
      Q(K,J)=-C(K,J)*PN(K)
      S(K,J)=-A(K,J)*S(K,J-1)*PN(K)
      F(K,J)=(F(K,J)-A(K,J)*F(K,J-1))*PN(K)
   20 CONTINUE

C     BACKWARD PASS
      DO 30 K=L1,L2
      S(K,J2)=1.
      QE(K,J2)=0.
   30 CONTINUE
      DO 40 I=JA,J2
      J=JJ-I
      DO 40 K=L1,L2
      S(K,J)=S(K,J)+Q(K,J)*S(K,J+1)
      QE(K,J)=F(K,J)+Q(K,J)*QE(K,J+1)
   40 CONTINUE
      DO 50 K=L1,L2
      F(K,J2)=(FN(K)-C(K,J2)*QE(K,J1)-A(K,J2)*QE(K,J2-1))
     %       /(C(K,J2)*S(K,J1)+A(K,J2)*S(K,J2-1)+B(K,J2))
   50 CONTINUE

C     BACKWARD ELIMINATION PASS
      DO 60 I=JA,J2
      J=JJ-I
      DO 60 K=L1,L2
      F(K,J)=F(K,J2)*S(K,J)+QE(K,J)
   60 CONTINUE

      RETURN
      END

C*******************************************************************
      SUBROUTINE TRDIAG2P(A,B,C,F,J1,J2,L1,L2)
C*******************************************************************
C     INVERT A PERIODIC MATRIX (X-2 TRIDIAGONAL) WITH
C     DIFFERENT COEFFICIENTS.

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      REAL A(M1,M2),B(M1,M2),C(M1,M2),F(M1,M2),Q(M1,M2),
     %     S(M1,M2),QE(M1,M2),FN(M1),PN(M1)

      JA=J1+1
      JJ=J1+J2
      DO 10 K=L1,L2
      BINV=1./B(K,J1)
      Q(K,J1)=-C(K,J1)*BINV
      S(K,J1)=-A(K,J1)*BINV
      FN(K)=F(K,J2)
      F(K,J1)=F(K,J1)*BINV
   10 CONTINUE

C     FORWARD ELIMINATION SWEEP
      DO 20 J=JA,J2
      DO 20 K=L1,L2
      PN(K)=1./(B(K,J)+A(K,J)*Q(K,J-1))
      Q(K,J)=-C(K,J)*PN(K)
      S(K,J)=-A(K,J)*S(K,J-1)*PN(K)
      F(K,J)=(F(K,J)-A(K,J)*F(K,J-1))*PN(K)
   20 CONTINUE

C     BACKWARD PASS
      DO 30 K=L1,L2
      S(K,J2)=1.
      QE(K,J2)=0.
   30 CONTINUE
      DO 40 I=JA,J2
      J=JJ-I
      DO 40 K=L1,L2
      S(K,J)=S(K,J)+Q(K,J)*S(K,J+1)
      QE(K,J)=F(K,J)+Q(K,J)*QE(K,J+1)
   40 CONTINUE
      DO 50 K=L1,L2
      F(K,J2)=(FN(K)-C(K,J2)*QE(K,J1)-A(K,J2)*QE(K,J2-1))
     %       /(C(K,J2)*S(K,J1)+A(K,J2)*S(K,J2-1)+B(K,J2))
   50 CONTINUE

C     BACKWARD ELIMINATION PASS
      DO 60 I=JA,J2
      J=JJ-I
      DO 60 K=L1,L2
      F(K,J)=F(K,J2)*S(K,J)+QE(K,J)
   60 CONTINUE

      RETURN
      END

C*******************************************************************
      SUBROUTINE TRDIAG3P(A,B,C,F,J1,J2,L1,L2)
C*******************************************************************
C     INVERT A PERIODIC MATRIX (X-3 TRIDIAGONAL) WITH
C     DIFFERENT COEFFICIENTS.

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      REAL A(M1,M3),B(M1,M3),C(M1,M3),F(M1,M3),Q(M1,M3),
     %     S(M1,M3),QE(M1,M3),FN(M1),PN(M1)

      JA=J1+1
      JJ=J1+J2
      DO 10 K=L1,L2
      BINV=1./B(K,J1)
      Q(K,J1)=-C(K,J1)*BINV
      S(K,J1)=-A(K,J1)*BINV
      FN(K)=F(K,J2)
      F(K,J1)=F(K,J1)*BINV
   10 CONTINUE

C     FORWARD ELIMINATION SWEEP
      DO 20 J=JA,J2
      DO 20 K=L1,L2
      PN(K)=1./(B(K,J)+A(K,J)*Q(K,J-1))
      Q(K,J)=-C(K,J)*PN(K)
      S(K,J)=-A(K,J)*S(K,J-1)*PN(K)
      F(K,J)=(F(K,J)-A(K,J)*F(K,J-1))*PN(K)
   20 CONTINUE

C     BACKWARD PASS
      DO 30 K=L1,L2
      S(K,J2)=1.
      QE(K,J2)=0.
   30 CONTINUE
      DO 40 I=JA,J2
      J=JJ-I
      DO 40 K=L1,L2
      S(K,J)=S(K,J)+Q(K,J)*S(K,J+1)
      QE(K,J)=F(K,J)+Q(K,J)*QE(K,J+1)
   40 CONTINUE
      DO 50 K=L1,L2
      F(K,J2)=(FN(K)-C(K,J2)*QE(K,J1)-A(K,J2)*QE(K,J2-1))
     %       /(C(K,J2)*S(K,J1)+A(K,J2)*S(K,J2-1)+B(K,J2))
   50 CONTINUE

C     BACKWARD ELIMINATION PASS
      DO 60 I=JA,J2
      J=JJ-I
      DO 60 K=L1,L2
      F(K,J)=F(K,J2)*S(K,J)+QE(K,J)
   60 CONTINUE

      RETURN
      END

!C*******************************************************************
!      SUBROUTINE TRDIAG1(A,B,C,R,UU,L1,L2,LL1,LL2)
!C*******************************************************************
!C     Tridiagonal Matrix Solver (X-1 direction)
!C     CONSTANT coefficients
!C     L1,L2 : TDMA (J) direction
!C     LL1,LL2 : Vectorization (I) direction
!
!      INCLUDE'param.h'
!      DIMENSION GAM(M1),A(M1),B(M1),C(M1)
!     &         ,R(M2,M1),UU(M2,M1),BET(M1)
!
!      BET(L1)=B(L1)
!      DO 10 I=LL1,LL2
!      UU(I,L1)=R(I,L1)/BET(L1)
!   10 CONTINUE
!      DO 11 J=L1+1,L2
!      GAM(J)=C(J-1)/BET(J-1)
!      BET(J)=B(J)-A(J)*GAM(J)
!   11 CONTINUE
!
!      DO 20 J=L1+1,L2
!      DO 20 I=LL1,LL2
!      UU(I,J)=(R(I,J)-A(J)*UU(I,J-1))/BET(J)
!   20 CONTINUE
!
!      DO 30 J=L2-1,L1,-1
!      DO 30 I=LL1,LL2
!      UU(I,J)=UU(I,J)-GAM(J+1)*UU(I,J+1)
!   30 CONTINUE
!
!      RETURN
!      END
!
!C*******************************************************************
!      SUBROUTINE TRDIAG2(A,B,C,R,UU,L1,L2,LL1,LL2)
!C*******************************************************************
!C     Tridiagonal Matrix Solver (X-2 direction)
!C     CONSTANT coefficients
!C     L1,L2 : TDMA (J) direction
!C     LL1,LL2 : Vectorization (I) direction
!
!      INCLUDE'param.h'
!      DIMENSION GAM(M2),A(M2),B(M2),C(M2)
!     &         ,R(M1,M2),UU(M1,M2),BET(M2)
!
!      BET(L1)=B(L1)
!      DO 10 I=LL1,LL2
!      UU(I,L1)=R(I,L1)/BET(L1)
!   10 CONTINUE
!      DO 11 J=L1+1,L2
!      GAM(J)=C(J-1)/BET(J-1)
!      BET(J)=B(J)-A(J)*GAM(J)
!   11 CONTINUE
!
!      DO 20 J=L1+1,L2
!      DO 20 I=LL1,LL2
!      UU(I,J)=(R(I,J)-A(J)*UU(I,J-1))/BET(J)
!   20 CONTINUE
!
!      DO 30 J=L2-1,L1,-1
!      DO 30 I=LL1,LL2
!      UU(I,J)=UU(I,J)-GAM(J+1)*UU(I,J+1)
!   30 CONTINUE
!
!      RETURN
!      END
!
!C*******************************************************************
!      SUBROUTINE TRDIAG3(A,B,C,R,UU,L1,L2,LL1,LL2)
!C*******************************************************************
!C     Tridiagonal Matrix Solver (X-3 direction)
!C     CONSTANT coefficients
!C     L1,L2 : TDMA (J) direction
!C     LL1,LL2 : Vectorization (I) direction
!
!      INCLUDE'param.h'
!      DIMENSION GAM(M3),A(M3),B(M3),C(M3)
!     &         ,R(M1,M3),UU(M1,M3),BET(M3)
!
!      BET(L1)=B(L1)
!      DO 10 I=LL1,LL2
!      UU(I,L1)=R(I,L1)/BET(L1)
!   10 CONTINUE
!      DO 11 J=L1+1,L2
!      GAM(J)=C(J-1)/BET(J-1)
!      BET(J)=B(J)-A(J)*GAM(J)
!   11 CONTINUE
!
!      DO 20 J=L1+1,L2
!      DO 20 I=LL1,LL2
!      UU(I,J)=(R(I,J)-A(J)*UU(I,J-1))/BET(J)
!   20 CONTINUE
!
!      DO 30 J=L2-1,L1,-1
!      DO 30 I=LL1,LL2
!      UU(I,J)=UU(I,J)-GAM(J+1)*UU(I,J+1)
!   30 CONTINUE
!
!      RETURN
!      END


C----------------------------- PTDIAG1A -------------------------------
C     SOLVE THE PENTADIAGONAL MATRIX (X-1 PENTADIAGONAL)

      SUBROUTINE PTDIAG1A(A,B,C,D,E,F,UU,I1,IN,J1,JN,O,Q,R)

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      DIMENSION A(M2,M1),B(M2,M1),C(M2,M1),D(M2,M1),E(M2,M1)
      DIMENSION O(M2,M1),Q(M2,M1),R(M2,M1),F(M2,M1),UU(M2,M1)
      DIMENSION PDENO2(M2),PDENOJ(M2)

      J2=J1+1
      J3=J2+1
      JNM=JN-1
      DO 1 I=I1,IN
      O(I,J1)=-B(I,J1)/A(I,J1)
      Q(I,J1)=-C(I,J1)/A(I,J1)
      R(I,J1)=F(I,J1)/A(I,J1)
    1 CONTINUE
      DO 2 I=I1,IN
      PDENO2(I)=1./(A(I,J2)+D(I,J2)*O(I,J1))
      O(I,J2)=-(B(I,J2)+D(I,J2)*Q(I,J1))*PDENO2(I)
      Q(I,J2)=-C(I,J2)*PDENO2(I)
      R(I,J2)=(F(I,J2)-D(I,J2)*R(I,J1))*PDENO2(I)
    2 CONTINUE
      DO 10 J=J3,JN
      JM=J-1
      JMM=JM-1
      DO 10 I=I1,IN
      PDENOJ(I)=1./(A(I,J)+E(I,J)*Q(I,JMM)
     1             +(D(I,J)+E(I,J)*O(I,JMM))*O(I,JM))
      O(I,J)=-(B(I,J)+(D(I,J)+E(I,J)*O(I,JMM))*Q(I,JM))*PDENOJ(I)
      Q(I,J)=-C(I,J)*PDENOJ(I)
      R(I,J)=(F(I,J)-E(I,J)*R(I,JMM)-(D(I,J)+E(I,J)*O(I,JMM))*R(I,JM))
     1      *PDENOJ(I)
   10 CONTINUE
      DO 11 I=I1,IN
      UU(I,JN)=R(I,JN)
   11 CONTINUE
      DO 12 I=I1,IN
      UU(I,JNM)=O(I,JNM)*UU(I,JN)+R(I,JNM)
   12 CONTINUE
      DO 20 J=JNM-1,J1,-1
      DO 20 I=I1,IN
      UU(I,J)=O(I,J)*UU(I,J+1)+Q(I,J)*UU(I,J+2)+R(I,J)
   20 CONTINUE

      RETURN
      END

C*******************************************************************
      SUBROUTINE MASSCHECK(QMMAX)
C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE IBM_VAR

      QMMAX=0.
!$OMP PARALLEL DO private(I,J)
!$OMP&reduction(MAX:QMMAX)
      DO 100 K=1,N3M
      DO 100 J=1,N2M
      DO 100 I=1,N1M
      QMMAX=MAX(QMMAX,ABS(QMASS(I,J,K)))  !THIS QMASS IS NOT Q WHEN TWO-PHASE SIMULATION.
                                            !DENSITY SHOULD BE 
 100  CONTINUE

      RETURN
      END

C******************************************************************
      SUBROUTINE CAL_RHSPS(U,V,W,ANPSO,PSI_CN,DF_C)
C******************************************************************
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE IBM_VAR
      USE HEAT_VAR
      
      USE TWO_PHASE_PROPERTY
      
      REAL PSI_CN(0:M1,0:M2,0:M3)
      REAL DF_C(0:M1,0:M2,0:M3,3)

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL ANPSO(M1M,M2M,M3M)

!!!!! FOR MEAN TEMPERATURE
       TGR_TMP=0.
!$OMP PARALLEL DO private(I,J,KP1,DENSC,UU)
!$OMP&reduction(+:TGR_TMP)
       DO K=1,N3M
       	KP1=KPV(K)
       DO J=1,N2M
       DO I=1,N1M
       	IF (FUNCBODY(XP(I),YP(J),ZP(K)) .GT. 0.) THEN
         DENSC=DENSCM+DENSC_DIFF*PSI_CN(I,J,K)
         UU=0.5*(W(I,J,KP1)+W(I,J,K))

         TGR_TMP=TGR_TMP+DENSC*UU*(SDX(I)*SDY(J)*SDZ(K))
        ENDIF
       ENDDO
       ENDDO
       ENDDO

       TGR=-HEAT_SUM/TGR_TMP !0.5*YL IS R (CHARACTERISTIC LENGTH)

      REPRI=1./(RE*PRM)
      
       IF (DEN_DIFF .EQ. 0.) THEN
      	DEN_DIFFI=1.
       ELSE
        DEN_DIFFI=1./DEN_DIFF
       ENDIF
!$OMP PARALLEL DO private(I,J,IP1,IM1,JP1,JM1,KP1,KM1,PSI_TMP)
!$OMP&private(DENSCE,DENSCW,DENSCN,DENSCS,DENSCT,DENSCB)
!$OMP&private(TCE,TCW,TCN,TCS,TCT,TCB)
!$OMP&private(PSE,PSW,PSC,PSF,PSN,PSS,ANPSX,ANPSY,ANPSZ,ANPS)
!$OMP&private(ALPSE,ALPSW,ALPSN,ALPSS,ALPSC,ALPSF,ALPSX,ALPSY,ALPSZ)
!$OMP&private(ALPS,DENSC,UU,AMEANT1,AMEANT2)
      DO K=1,N3M
      	KP1=KPV(K)
      	KM1=KMV(K)
      DO J=1,N2M
      	JP1=JPV(J)
      	JM1=JMV(J)
      DO I=1,N1M
      	IP1=IPV(I)
      	IM1=IMV(I)

       PSI_TMP=(DF_C(I,J,K,1)-DENM)*DEN_DIFFI
       DENSCE=DENSCM+DENSC_DIFF*PSI_TMP
     	 TCE=TCM+TC_DIFF*PSI_TMP
       PSI_TMP=(DF_C(IM1,J,K,1)-DENM)*DEN_DIFFI
       DENSCW=DENSCM+DENSC_DIFF*PSI_TMP
     	 TCW=TCM+TC_DIFF*PSI_TMP
       PSI_TMP=(DF_C(I,J,K,2)-DENM)*DEN_DIFFI
       DENSCN=DENSCM+DENSC_DIFF*PSI_TMP
     	 TCN=TCM+TC_DIFF*PSI_TMP
       PSI_TMP=(DF_C(I,JM1,K,2)-DENM)*DEN_DIFFI
       DENSCS=DENSCM+DENSC_DIFF*PSI_TMP
     	 TCS=TCM+TC_DIFF*PSI_TMP
       PSI_TMP=(DF_C(I,J,K,3)-DENM)*DEN_DIFFI
       DENSCT=DENSCM+DENSC_DIFF*PSI_TMP
     	 TCT=TCM+TC_DIFF*PSI_TMP
       PSI_TMP=(DF_C(I,J,KM1,3)-DENM)*DEN_DIFFI
       DENSCB=DENSCM+DENSC_DIFF*PSI_TMP
     	 TCB=TCM+TC_DIFF*PSI_TMP

      !CONVECTION
      PSE=0.5*(SDX(IP1)*T(I,J,K)+SDX(I)*T(IP1,J,K))*VVDX(IP1)
      PSW=0.5*(SDX(IM1)*T(I,J,K)+SDX(I)*T(IM1,J,K))*VVDX(I)
      PSC=0.5*(SDZ(KP1)*T(I,J,K)+SDZ(K)*T(I,J,KP1))*VVDZ(KP1)
      PSF=0.5*(SDZ(KM1)*T(I,J,K)+SDZ(K)*T(I,J,KM1))*VVDZ(K)
      PSN=0.5*(SDY(JP1)*T(I,J,K)+SDY(J)*T(I,JP1,K))*VVDY(JP1)
     &   *(1.-FIXJU(J))+T(I,N2,K)*FIXJU(J)
      PSS=0.5*(SDY(JM1)*T(I,J,K)+SDY(J)*T(I,JM1,K))*VVDY(J)
     &   *(1.-FIXJL(J))+T(I,0,K)*FIXJL(J)
      ANPSX=(DENSCE*U(IP1,J,K)*PSE-DENSCW*U(I,J,K)*PSW)*SSDX(I)
      ANPSY=(DENSCN*V(I,JP1,K)*PSN-DENSCS*V(I,J,K)*PSS)*SSDY(J)
      ANPSZ=(DENSCT*W(I,J,KP1)*PSC-DENSCB*W(I,J,K)*PSF)*SSDZ(K)
      ANPS=-(ANPSX+ANPSY+ANPSZ)

      !DIFFUSION
      ALPSE=(T(IP1,J,K)-T(I,J,K))*VVDX(IP1)
      ALPSW=(T(I,J,K)-T(IM1,J,K))*VVDX(I)
      ALPSN=(T(I,JP1,K)-T(I,J,K))*VVDY(JP1)
      ALPSS=(T(I,J,K)-T(I,JM1,K))*VVDY(J)
      ALPSC=(T(I,J,KP1)-T(I,J,K))*VVDZ(KP1)
      ALPSF=(T(I,J,K)-T(I,J,KM1))*VVDZ(K)
      ALPSX=(TCE*ALPSE-TCW*ALPSW)*SSDX(I)
      ALPSY=(TCN*ALPSN-TCS*ALPSS)*SSDY(J)
      ALPSZ=(TCT*ALPSC-TCB*ALPSF)*SSDZ(K)
      ALPS=REPRI*(ALPSX+ALPSY+ALPSZ)

!      !LES
!      TALPHE=0.5*(SDX(I+1)*TALPH(I,J,K)+SDX(I)*TALPH(I+1,J,K))*VVDX(I+1)
!      TALPHW=0.5*(SDX(I-1)*TALPH(I,J,K)+SDX(I)*TALPH(I-1,J,K))*VVDX(I)
!      TALPHC=0.5*(SDZ(K+1)*TALPH(I,J,K)+SDZ(K)*TALPH(I,J,K+1))*VVDZ(K+1)
!      TALPHF=0.5*(SDZ(K-1)*TALPH(I,J,K)+SDZ(K)*TALPH(I,J,K-1))*VVDZ(K)
!      TALPHN=0.5*(SDY(J+1)*TALPH(I,J,K)+SDY(J)*TALPH(I,J+1,K))*VVDY(J+1)
!     &   *(1.-FIXJU(J))+TALPH(I,N2,K)*FIXJU(J)
!      TALPHS=0.5*(SDY(J-1)*TALPH(I,J,K)+SDY(J)*TALPH(I,J-1,K))*VVDY(J)
!     &   *(1.-FIXJL(J))+TALPH(I,0,K)*FIXJL(J)
!      TALPSX=(TALPHE*ALPSE-TALPHW*ALPSW)*SSDX(I)
!      TALPSY=(TALPHN*ALPSN-TALPHS*ALPSS)*SSDY(J)
!      TALPSZ=(TALPHC*ALPSC-TALPHF*ALPSF)*SSDZ(K)
!      ALPS=ALPS+FLOAT(ILES)*(TALPSX+TALPSY+TALPSZ)

!        write(*,*) sce,scw,scn,scs,sct,scb
!        write(*,*) tce,tcw,tcn,tcs,tct,tcb

      IF (FUNCBODY(XP(I),YP(J),ZP(K)) .GT. 0.) THEN
       DENSC=DENSCM+DENSC_DIFF*PSI_CN(I,J,K)
       UU=0.5*(W(I,J,KP1)+W(I,J,K))
       AMEANT1=DENSC*UU*TGR

       AMEANT2=REPRI*(TCE-TCW)*SSDX(I)

      RHSPS(I,J,K)=( ( GAMMA(MSUB)*ANPS+RO(MSUB)*ANPSO(I,J,K) )
     &              +2.*ALPHA_RK3*ALPS
     &              +2.*ALPHA_RK3*(AMEANT1-AMEANT2) )*DT
      
      ELSE
      RHSPS(I,J,K)=( ( GAMMA(MSUB)*ANPS+RO(MSUB)*ANPSO(I,J,K) )
     &              +2.*ALPHA_RK3*ALPS )*DT
      ENDIF

      RHSPS_IBM(I,J,K)=RHSPS(I,J,K)

      ANPSO(I,J,K)=ANPS

      ENDDO
      ENDDO
      ENDDO

      RETURN
      END


C*******************************************************************
      SUBROUTINE IBMINIT_PS
C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE IBM_VAR

      ALLOCATE(IFC_PS(MBODY),JFC_PS(MBODY),KFC_PS(MBODY))
      ALLOCATE(FCV_PS(MBODY),FCVAVG_PS(MBODY))
      ALLOCATE(MPI_PS(MBODY0,3))
      ALLOCATE(GFI_PS(MBODY0,0:2,0:2,0:2))

      ALLOCATE(RHSPS_IBM(M1M,M2M,M3M))

!$OMP PARALLEL DO
        DO L=1,MBODY
        IFC_PS(L)=0 ;JFC_PS(L)=0 ;KFC_PS(L)=0
        FCV_PS(L)=0. ;FCVAVG_PS(L)=0.
        ENDDO
!$OMP PARALLEL DO
        DO L=1,MBODY0
         MPI_PS(L,1)=0 ;MPI_PS(L,2)=0 ;MPI_PS(L,3)=0
        ENDDO
!$OMP PARALLEL DO private(I,J,K)
        DO L=1,MBODY0
         DO K=0,2 ;DO J=0,2 ;DO I=0,2
         GFI_PS(L,I,J,K)=0.
         ENDDO ;ENDDO ;ENDDO
        ENDDO

!$OMP PARALLEL DO private(I,J)
        DO K=1,N3M ;DO J=1,N2M ;DO I=1,N1M
         RHSPS_IBM(I,J,K)=0.
        ENDDO ;ENDDO ;ENDDO

      OPEN(11,FILE='ibmpre_intp_ps.bin',FORM='UNFORMATTED')
      READ(11) NFC_INTP_PS
      READ(11) NFC_INNER_PS
      READ(11) (IFC_PS(N),N=1,NFC_INTP_PS+NFC_INNER_PS)
      READ(11) (JFC_PS(N),N=1,NFC_INTP_PS+NFC_INNER_PS)
      READ(11) (KFC_PS(N),N=1,NFC_INTP_PS+NFC_INNER_PS)
      READ(11) ((MPI_PS(N,L),N=1,NFC_INTP_PS),L=1,3)
      READ(11) ((((GFI_PS(N,I,J,K),N=1,NFC_INTP_PS),I=0,2),J=0,2),K=0,2)
      CLOSE(11)

      WRITE(*,35) NFC_INTP_PS
      WRITE(*,36) NFC_INNER_PS
 35   FORMAT('FU_PS_intp :',I12)
 36   FORMAT('FU_PS_inner:',I12)

      RETURN
      END

C*******************************************************************
      SUBROUTINE RHS_IBM_PS(PSI_CN,PSI_C)
C*******************************************************************
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE IBM_VAR

      USE HEAT_VAR
      
      USE TWO_PHASE_PROPERTY

      REAL PSI_CN(0:M1,0:M2,0:M3)
      REAL PSI_C(0:M1,0:M2,0:M3)

C-----compute target velocities & forcing values at forcing points
!      TBODY=1.
      DTI=1./DT

!$OMP PARALLEL DO private(II,JJ,KK,I2,J2,K2,I3,J3,K3)
!$OMP&private(R1,R2,R3,A11,A12,A13,A21,A22,A23,A31,A32,A33)
!$OMP&private(DET1,A1,A2,A3,TTARG,DENSC,DENSCF)
      DO N=1,NFC_INTP_PS
         II =IFC_PS(N)
         JJ =JFC_PS(N)
         KK =KFC_PS(N)

       I2=INT(GFI_PS(N,0,0,1))
       J2=INT(GFI_PS(N,1,0,1))
       K2=INT(GFI_PS(N,2,0,1))
       I3=INT(GFI_PS(N,0,1,1))
       J3=INT(GFI_PS(N,1,1,1))
       K3=INT(GFI_PS(N,2,1,1))

       !RHS PART
       R1=1./(TCM+TC_DIFF*PSI_CN(II,JJ,KK))
       DENSC=DENSCM+DENSC_DIFF*PSI_CN(I2,J2,K2)
       DENSCF=DENSCM+DENSC_DIFF*PSI_C(I2,J2,K2)
       R2=(DENSC*T(I2,J2,K2)+RHSPS_IBM(I2,J2,K2))/DENSCF
       DENSC=DENSCM+DENSC_DIFF*PSI_CN(I3,J3,K3)
       DENSCF=DENSCM+DENSC_DIFF*PSI_C(I3,J3,K3)
       R3=(DENSC*T(I3,J3,K3)+RHSPS_IBM(I3,J3,K3))/DENSCF

       !CONSTRUCT MATRIX
       A11=GFI_PS(N,0,2,1)
       A12=GFI_PS(N,1,2,1)
       A13=GFI_PS(N,2,2,1)
       A21=1.
       A22=XP(I2)
       A23=YP(J2)
       A31=1.
       A32=XP(I3)
       A33=YP(J3)

!http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
       !SOLVE MATRIX
       DETI=1./(A11*A22*A33+A21*A32*A13+A31*A12*A23
     &                           -A11*A32*A23-A31*A22*A13-A21*A12*A33)
       A1=((A22*A33-A23*A32)*R1+(A13*A32-A12*A33)*R2
     &                                    +(A12*A23-A13*A22)*R3)*DETI
       A2=((A23*A31-A21*A33)*R1+(A11*A33-A13*A31)*R2
     &                                    +(A13*A21-A11*A23)*R3)*DETI
       A3=((A21*A32-A22*A31)*R1+(A12*A31-A11*A32)*R2
     &                                    +(A11*A22-A12*A21)*R3)*DETI

       TTARG=A1+A2*XP(II)+A3*YP(JJ)
        
!        WRITE(*,*) II,I2,I3
!        WRITE(*,*) JJ,J2,J3
!        IF (KK.NE.K2 .OR. KK.NE.K3) WRITE(*,*) 'HERE'
        
!        write(*,*) r1,r2,r3
!        write(*,*) a1,a2,a3
!        write(*,*) xp(ii),yp(jj),ttarg

       DENSC=DENSCM+DENSC_DIFF*PSI_CN(II,JJ,KK)
       DENSCF=DENSCM+DENSC_DIFF*PSI_C(II,JJ,KK)
       FCV_PS(N)=(DENSCF*TTARG-(DENSC*T(II,JJ,KK)+RHSPS_IBM(II,JJ,KK)))
     &                                                             *DTI
       RHSPS(II,JJ,KK)=RHSPS(II,JJ,KK)+FCV_PS(N)*DT

      ENDDO
       !FCV_PS IS ACTUALLY 2*ALPHA_RK3*HEAT_SOURCE

!$OMP PARALLEL DO private(II,JJ,KK)
!$OMP&private(DENSC,DENSCF,TTARG)
      DO N=NFC_INTP_PS+1,NFC_INTP_PS+NFC_INNER_PS
         II =IFC_PS(N)
         JJ =JFC_PS(N)
         KK =KFC_PS(N)
         TTARG=0.!TBODY
         
       DENSC=DENSCM+DENSC_DIFF*PSI_CN(II,JJ,KK)
       DENSCF=DENSCM+DENSC_DIFF*PSI_C(II,JJ,KK)

       FCV_PS(N)=(DENSCF*TTARG-(DENSC*T(II,JJ,KK)+RHSPS_IBM(II,JJ,KK)))
     &                                                              *DTI
       RHSPS(II,JJ,KK)=RHSPS(II,JJ,KK)+FCV_PS(N)*DT

      ENDDO

      IF (MSUB .EQ. 1) HEAT_SUM=0.
!$OMP PARALLEL DO reduction(+:HEAT_SUM)
      DO N=1,NFC_INTP_PS!+NFC_INNER_PS
       HEAT_SUM=HEAT_SUM+FCV_PS(N)
     &                    *SDX(IFC_PS(N))*SDY(JFC_PS(N))*SDZ(KFC_PS(N))
      ENDDO
!       IF (ALPHA_RK3 .NE. 0.) HEAT_SUM=HEAT_SUM/(2.*ALPHA_RK3)  !AT INITIAL STEP, ALPHA_RK3=0.      
!       WRITE(*,*) 'Total_HEAT_SUM=',HEAT_SUM

      RETURN
      END

!C*******************************************************************
!      SUBROUTINE RHS_IBM_PS(PSI_CN,PSI_C)
!C*******************************************************************
!      USE FLOW_VAR
!
!      USE PARAM_VAR
!      USE FLOW_GEOM_VAR
!
!      USE IBM_VAR
!
!      USE HEAT_VAR
!      
!      USE TWO_PHASE_PROPERTY
!
!      REAL PSI_CN(0:M1,0:M2,0:M3)
!      REAL PSI_C(0:M1,0:M2,0:M3)
!
!C-----compute target velocities & forcing values at forcing points
!!      TBODY=1.
!      DTI=1./DT
!      REPR=RE*PRM
!
!        IF (DENR .EQ. 1) THEN
!         DEN_DIFFI=1.
!        ELSE
!         DEN_DIFFI=1./DEN_DIFF
!        ENDIF
!
!!$OMP PARALLEL DO private(II,JJ,KK,IP,JP,KP)
!!$OMP&private(DENSC,DENSCF,T_TMP1,T_TMP2,T_TMP3,THETA_TMP,TCC,TTARG)
!      DO N=1,NFC_INTP_PS
!         II =IFC_PS(N)
!         JJ =JFC_PS(N)
!         KK =KFC_PS(N)
!         IP =IFC_PS(N)+MPI_PS(N,1)
!         JP =JFC_PS(N)+MPI_PS(N,2)
!         KP =KFC_PS(N)+MPI_PS(N,3)
!!         IPP=IFC_PS(N)+MPI_PS(N,1)*2
!!         JPP=JFC_PS(N)+MPI_PS(N,2)*2
!!         KPP=KFC_PS(N)+MPI_PS(N,3)*2
!!         TTARG=
!!     &   GFI_PS(N,0,0,0)* TBODY
!!     &  +GFI_PS(N,0,0,1)*(T(II ,JJ ,KP )+RHSPS_IBM(II ,JJ ,KP ))
!!     &  +GFI_PS(N,0,0,2)*(T(II ,JJ ,KPP)+RHSPS_IBM(II ,JJ ,KPP))
!!     &  +GFI_PS(N,0,1,0)*(T(II ,JP ,KK )+RHSPS_IBM(II ,JP ,KK ))
!!     &  +GFI_PS(N,0,1,1)*(T(II ,JP ,KP )+RHSPS_IBM(II ,JP ,KP ))
!!     &  +GFI_PS(N,0,1,2)*(T(II ,JP ,KPP)+RHSPS_IBM(II ,JP ,KPP))
!!     &  +GFI_PS(N,0,2,0)*(T(II ,JPP,KK )+RHSPS_IBM(II ,JPP,KK ))
!!     &  +GFI_PS(N,0,2,1)*(T(II ,JPP,KP )+RHSPS_IBM(II ,JPP,KP ))
!!     &  +GFI_PS(N,0,2,2)*(T(II ,JPP,KPP)+RHSPS_IBM(II ,JPP,KPP))
!!     &  +GFI_PS(N,1,0,0)*(T(IP ,JJ ,KK )+RHSPS_IBM(IP ,JJ ,KK ))
!!     &  +GFI_PS(N,1,0,1)*(T(IP ,JJ ,KP )+RHSPS_IBM(IP ,JJ ,KP ))
!!     &  +GFI_PS(N,1,0,2)*(T(IP ,JJ ,KPP)+RHSPS_IBM(IP ,JJ ,KPP))
!!     &  +GFI_PS(N,1,1,0)*(T(IP ,JP ,KK )+RHSPS_IBM(IP ,JP ,KK ))
!!     &  +GFI_PS(N,1,1,1)*(T(IP ,JP ,KP )+RHSPS_IBM(IP ,JP ,KP ))
!!     &  +GFI_PS(N,1,1,2)*(T(IP ,JP ,KPP)+RHSPS_IBM(IP ,JP ,KPP))
!!     &  +GFI_PS(N,1,2,0)*(T(IP ,JPP,KK )+RHSPS_IBM(IP ,JPP,KK ))
!!     &  +GFI_PS(N,1,2,1)*(T(IP ,JPP,KP )+RHSPS_IBM(IP ,JPP,KP ))
!!     &  +GFI_PS(N,1,2,2)*(T(IP ,JPP,KPP)+RHSPS_IBM(IP ,JPP,KPP))
!!     &  +GFI_PS(N,2,0,0)*(T(IPP,JJ ,KK )+RHSPS_IBM(IPP,JJ ,KK ))
!!     &  +GFI_PS(N,2,0,1)*(T(IPP,JJ ,KP )+RHSPS_IBM(IPP,JJ ,KP ))
!!     &  +GFI_PS(N,2,0,2)*(T(IPP,JJ ,KPP)+RHSPS_IBM(IPP,JJ ,KPP))
!!     &  +GFI_PS(N,2,1,0)*(T(IPP,JP ,KK )+RHSPS_IBM(IPP,JP ,KK ))
!!     &  +GFI_PS(N,2,1,1)*(T(IPP,JP ,KP )+RHSPS_IBM(IPP,JP ,KP ))
!!     &  +GFI_PS(N,2,1,2)*(T(IPP,JP ,KPP)+RHSPS_IBM(IPP,JP ,KPP))
!!     &  +GFI_PS(N,2,2,0)*(T(IPP,JPP,KK )+RHSPS_IBM(IPP,JPP,KK ))
!!     &  +GFI_PS(N,2,2,1)*(T(IPP,JPP,KP )+RHSPS_IBM(IPP,JPP,KP ))
!!     &  +GFI_PS(N,2,2,2)*(T(IPP,JPP,KPP)+RHSPS_IBM(IPP,JPP,KPP))
!
!       DENSC=DENSCM+DENSC_DIFF*PSI_CN(IP,JJ,KK)
!       DENSCF=DENSCM+DENSC_DIFF*PSI_C(IP,JJ,KK)
!       T_TMP1=(DENSC*T(IP,JJ,KK)+RHSPS_IBM(IP,JJ,KK))/DENSCF
!
!       DENSC=DENSCM+DENSC_DIFF*PSI_CN(IP,JP,KK)
!       DENSCF=DENSCM+DENSC_DIFF*PSI_C(IP,JP,KK)
!       T_TMP2=(DENSC*T(IP,JP,KK)+RHSPS_IBM(IP,JP,KK))/DENSCF
!
!       DENSC=DENSCM+DENSC_DIFF*PSI_CN(II,JP,KK)
!       DENSCF=DENSCM+DENSC_DIFF*PSI_C(II,JP,KK)
!       T_TMP3=(DENSC*T(II,JP,KK)+RHSPS_IBM(II,JP,KK))/DENSCF
!
!       THETA_TMP=GFI_PS(N,0,0,0)*T_TMP1
!     &              +GFI_PS(N,0,0,1)*T_TMP2+GFI_PS(N,0,0,2)*T_TMP3
!       TCC=TCM+TC_DIFF*PSI_CN(II,JJ,KK)
!       TTARG=THETA_TMP-GFI_PS(N,0,1,0)/TCC
!
!!       THETA_TMP=T(IP,JP,KK)
!!       TCC=TCM+TC_DIFF*PSI_CN(II,JJ,KK)
!!       TTARG=THETA_TMP-sqrt((xp(ip)-xp(ii))**2+(yp(jp)-yp(jj))**2)/TCC
!
!       DENSC=DENSCM+DENSC_DIFF*PSI_CN(II,JJ,KK)
!       DENSCF=DENSCM+DENSC_DIFF*PSI_C(II,JJ,KK)
!
!       FCV_PS(N)=(DENSCF*TTARG-(DENSC*T(II,JJ,KK)+RHSPS_IBM(II,JJ,KK)))
!     &                                                              *DTI
!       RHSPS(II,JJ,KK)=RHSPS(II,JJ,KK)+FCV_PS(N)*DT
!
!!      if ( kk.eq.1) then
!!      if (msub .eq. 1) write(*,*) ii,jj,GFI_PS(N,0,1,0)!fcv_ps(n)/(2.*alpha_rk3)
!!!      write(*,*) ttarg,t(ii,jj,kk),rhsps_ibm(ii,jj,kk)/(2.*alpha_rk3)
!!!      write(*,*) TTARG-T(II,JJ,KK),TTARG-T(II,JJ,KK)-RHSPS_IBM(II,JJ,KK)
!!!      write(*,*)  (TTARG-T(II,JJ,KK)-RHSPS_IBM(II,JJ,KK))*dti
!!!      write(*,*)  (TTARG-T(II,JJ,KK)-RHSPS_IBM(II,JJ,KK))*dti
!!     & /(2.*alpha_rk3)
!!      endif
!
!      ENDDO
!       !FCV_PS IS ACTUALLY 2*ALPHA_RK3*HEAT_SOURCE
!
!!$OMP PARALLEL DO private(II,JJ,KK)
!!$OMP&private(DENSC,DENSCF,TTARG)
!      DO N=NFC_INTP_PS+1,NFC_INTP_PS+NFC_INNER_PS
!         II =IFC_PS(N)
!         JJ =JFC_PS(N)
!         KK =KFC_PS(N)
!         TTARG=0.!TBODY
!         
!       DENSC=DENSCM+DENSC_DIFF*PSI_CN(II,JJ,KK)
!       DENSCF=DENSCM+DENSC_DIFF*PSI_C(II,JJ,KK)
!
!       FCV_PS(N)=(DENSCF*TTARG-(DENSC*T(II,JJ,KK)+RHSPS_IBM(II,JJ,KK)))
!     &                                                              *DTI
!       RHSPS(II,JJ,KK)=RHSPS(II,JJ,KK)+FCV_PS(N)*DT
!
!      ENDDO
!
!      IF (MSUB .EQ. 1) HEAT_SUM=0.
!!$OMP PARALLEL DO reduction(+:HEAT_SUM)
!      DO N=1,NFC_INTP_PS!+NFC_INNER_PS
!       HEAT_SUM=HEAT_SUM+FCV_PS(N)
!     &                    *SDX(IFC_PS(N))*SDY(JFC_PS(N))*SDZ(KFC_PS(N))
!      ENDDO
!!       IF (ALPHA_RK3 .NE. 0.) HEAT_SUM=HEAT_SUM/(2.*ALPHA_RK3)  !AT INITIAL STEP, ALPHA_RK3=0.      
!!       WRITE(*,*) 'Total_HEAT_SUM=',HEAT_SUM
!
!      RETURN
!      END

C------------------------------- LHSPS --------------------------------
C******************************************************************
      SUBROUTINE LHSPS(U,V,W,PSI_CN,PSI_C,DF_C)
C******************************************************************
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE HEAT_VAR
      
      USE TWO_PHASE_PROPERTY

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL PSI_CN(0:M1,0:M2,0:M3)
      REAL DF_C(0:M1,0:M2,0:M3,3)
      REAL PSI_C(0:M1,0:M2,0:M3)

      real, dimension (:,:), allocatable :: AI,BI,CI,GI
      real, dimension (:,:), allocatable :: AJ,BJ,CJ,GJ
      real, dimension (:,:), allocatable :: AK,BK,CK,GK

       IF (DEN_DIFF .EQ. 0.) THEN
      	DEN_DIFFI=1.
       ELSE
        DEN_DIFFI=1./DEN_DIFF
       ENDIF

!ccccc_explicit_euler
!       DO K=1,N3M
!       DO J=1,N2M
!       DO I=1,N1M
!        DENSC=DENSCM+DENSC_DIFF*PSI_CN(I,J,K)
!        DENSCF=DENSCM+DENSC_DIFF*PSI_C(I,J,K)
!        T(I,J,K)=(DENSC*T(I,J,K)+RHSPS(I,J,K))/DENSCF
!       ENDDO
!       ENDDO
!       ENDDO
!       WRITE(*,*) 'ENERGY EQN. IS SOLVED WITH EXPLICIT EULER'
!ccccc_explicit_euler

      CREPR=0.!FLOAT(ILES)*RE*PRA
      ACOPS=ALPHA_RK3*DT/(RE*PRM)!0.5*DT*REI*PRAI

C     Z-DIRECTION
      IF (N3M.EQ.1) GOTO 1000
!$OMP PARALLEL private(I,K,DK1,DK2,DK3)
!$OMP&private(AK,CK,BK,GK,DENSCF,COEF,PSI_TMP,TCT,TCB,TALPHC,TALPHF)
      allocate(AK(M1,M3),BK(M1,M3),CK(M1,M3),GK(M1,M3))
!$OMP DO
      DO J=1,N2M

      DO K=1,N3M
      DO I=1,N1M
       	DENSCF=DENSCM+DENSC_DIFF*PSI_C(I,J,K)
       	COEF=ACOPS/DENSCF

        PSI_TMP=(DF_C(I,J,K,3)-DENM)*DEN_DIFFI
       	TCT=TCM+TC_DIFF*PSI_TMP
        PSI_TMP=(DF_C(I,J,KMV(K),3)-DENM)*DEN_DIFFI
       	TCB=TCM+TC_DIFF*PSI_TMP

      TALPHC=0.!0.5*(SDZ(K+1)*TALPH(I,J,K)+SDZ(K)*TALPH(I,J,K+1))*VVDZ(K+1)
      TALPHF=0.!0.5*(SDZ(K-1)*TALPH(I,J,K)+SDZ(K)*TALPH(I,J,K-1))*VVDZ(K)
      AK(I,K)=COEF*TCB*AKUV(K)*(1.+CREPR*TALPHF)
      CK(I,K)=COEF*TCT*CKUV(K)*(1.+CREPR*TALPHC)
      BK(I,K)=1.-(AK(I,K)+CK(I,K))
      GK(I,K)=RHSPS(I,J,K)
      ENDDO
      ENDDO

      IF (IPZ .EQ. 1) THEN
       CALL TRDIAG3P(AK,BK,CK,GK,1,N3M,1,N1M)
      ELSE
       CALL TRDIAG3(AK,BK,CK,GK,GK,1,N3M,1,N1M)
      ENDIF

      DO K=1,N3M
      DO I=1,N1M
      RHSPS(I,J,K)=GK(I,K)
      ENDDO
      ENDDO
      ENDDO
!$OMP END DO
      deallocate(AK,BK,CK,GK)
!$OMP END PARALLEL

 1000 CONTINUE

!$OMP PARALLEL private(I,J,DJ1,DJ2,DJ3,DI1,DI2,DI3)
!$OMP&private(AJ,BJ,CJ,GJ,AI,BI,CI,GI,DENSCF,COEF,PSI_TMP,TCN,TCS)
!$OMP&private(TALPHN,TALPHS,TCE,TCW,TALPHE,TALPHW)
      allocate(AI(M2,M1),BI(M2,M1),CI(M2,M1),GI(M2,M1))
      allocate(AJ(M1,M2),BJ(M1,M2),CJ(M1,M2),GJ(M1,M2))
!$OMP DO
      DO K=1,N3M

C     Y-DIRECTION
      DO J=1,N2M
      DO I=1,N1M
       	DENSCF=DENSCM+DENSC_DIFF*PSI_C(I,J,K)
       	COEF=ACOPS/DENSCF

        PSI_TMP=(DF_C(I,J,K,2)-DENM)*DEN_DIFFI
       	TCN=TCM+TC_DIFF*PSI_TMP
        PSI_TMP=(DF_C(I,JMV(J),K,2)-DENM)*DEN_DIFFI
       	TCS=TCM+TC_DIFF*PSI_TMP
      	
      TALPHN=0.!0.5*(SDY(J+1)*TALPH(I,J,K)+SDY(J)*TALPH(I,J+1,K))*VVDY(J+1)
!     &   *(1.-FIXJU(J))+TALPH(I,N2,K)*FIXJU(J)
      TALPHS=0.!0.5*(SDY(J-1)*TALPH(I,J,K)+SDY(J)*TALPH(I,J-1,K))*VVDY(J)
!     &   *(1.-FIXJL(J))+TALPH(I,0,K)*FIXJL(J)
      AJ(I,J)=COEF*TCS*AJUW(J)*(1.+CREPR*TALPHS)
      CJ(I,J)=COEF*TCN*CJUW(J)*(1.+CREPR*TALPHN)     
      BJ(I,J)=1.-(AJ(I,J)+CJ(I,J))
      GJ(I,J)=RHSPS(I,J,K)
      ENDDO
      ENDDO

      IF (IPY .EQ. 1) THEN
       CALL TRDIAG2P(AJ,BJ,CJ,GJ,1,N2M,1,N1M)
      ELSE
       CALL TRDIAG2(AJ,BJ,CJ,GJ,GJ,1,N2M,1,N1M)
      ENDIF
      
C     X-DIRECTION
      DO J=1,N2M
      DO I=1,N1M
       	DENSCF=DENSCM+DENSC_DIFF*PSI_C(I,J,K)
       	COEF=ACOPS/DENSCF

        PSI_TMP=(DF_C(I,J,K,1)-DENM)*DEN_DIFFI
       	TCE=TCM+TC_DIFF*PSI_TMP
        PSI_TMP=(DF_C(IMV(I),J,K,1)-DENM)*DEN_DIFFI
       	TCW=TCM+TC_DIFF*PSI_TMP

      TALPHE=0.!0.5*(SDX(I+1)*TALPH(I,J,K)+SDX(I)*TALPH(I+1,J,K))*VVDX(I+1)
      TALPHW=0.!0.5*(SDX(I-1)*TALPH(I,J,K)+SDX(I)*TALPH(I-1,J,K))*VVDX(I)
      AI(J,I)=COEF*TCW*AIVW(I)*(1.+CREPR*TALPHW)
      CI(J,I)=COEF*TCE*CIVW(I)*(1.+CREPR*TALPHE)
      BI(J,I)=1.-(AI(J,I)+CI(J,I))
      GI(J,I)=GJ(I,J)
      ENDDO
      ENDDO

      IF (IPX .EQ. 1) THEN
       CALL TRDIAG1P(AI,BI,CI,GI,1,N1M,1,N2M)
      ELSE
       CALL TRDIAG1(AI,BI,CI,GI,GI,1,N1M,1,N2M)
      ENDIF

C     UPDATE T
      DO J=1,N2M
      DO I=1,N1M
      T(I,J,K)=GI(J,I)+T(I,J,K)
      ENDDO
      ENDDO

      ENDDO
!$OMP END DO
      deallocate(AI,BI,CI,GI)
      deallocate(AJ,BJ,CJ,GJ)
!$OMP END PARALLEL

!     PERIODICITY
      IF (IPZ .EQ. 1) THEN
!$OMP PARALLEL DO private(I)
      DO J=1,N2M
      DO I=1,N1M
         T(I,J,0) =T(I,J,N3M)
         T(I,J,N3)=T(I,J,1)
      ENDDO
      ENDDO
      ELSE 
!$OMP PARALLEL DO private(I)
      DO J=1,N2M
      DO I=1,N1M
         T(I,J,0) =T(I,J,1)
         T(I,J,N3)=T(I,J,N3M)
      ENDDO
      ENDDO
      ENDIF

      IF (IPY .EQ. 1) THEN
!$OMP PARALLEL DO private(I)
      DO K=1,N3M
      DO I=1,N1M
         T(I,0 ,K)=T(I ,N2M,K)
         T(I,N2,K)=T(I ,1  ,K)
      ENDDO
      ENDDO
      ELSE
!$OMP PARALLEL DO private(I)
      DO K=1,N3M
      DO I=1,N1M
         T(I,0 ,K)=T(I ,1,K)
         T(I,N2,K)=T(I ,N2M,K)
      ENDDO
      ENDDO
      ENDIF
      
      IF (IPX .EQ. 1) THEN
!$OMP PARALLEL DO private(J)
      DO K=1,N3M
      DO J=1,N2M
         T(0 ,J,K)=T(N1M,J,K)
         T(N1,J,K)=T(1  ,J,K)
      ENDDO
      ENDDO
      ELSE
!$OMP PARALLEL DO private(J)
      DO K=1,N3M
      DO J=1,N2M
         T(0 ,J,K)=T(1,J,K)
         T(N1,J,K)=T(N1M,J,K)
      ENDDO
      ENDDO
      ENDIF

      RETURN
      END

C******************************************************************
!                      FIELD AVERAGING ROUTINES
!                                                    JUNGIL LEE
!                                                     AUG. 2006
C******************************************************************
C******************************************************************
      SUBROUTINE FIELD_AVG_INIT_RT(IHIST)
C******************************************************************
C     SUBROUTINE FOR FIELD AVERAGING
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE FLD_AVG

      USE HEAT_VAR

      ALLOCATE (URAVGZ(M1M,M2M),UTAVGZ(M1M,M2M))
      ALLOCATE (WAVGZ(M1M,M2M),UIUJAVGZ(M1M,M2M,6),
     &          PAVGZ(M1M,M2M),P2AVGZ(M1M,M2M),
     &          TNUAVGZ(M1M,M2M),VOLFZ(M1M,M2M))
!     &         SIJAVGZ(M1M,M2M,6),VORAVGZ(M1M,M2M,3),
!     &         SSAVGZ(M1M,M2M),SGSDISSZ(M1M,M2M),
      ALLOCATE (URAVGZ0(M1M,M2M),UTAVGZ0(M1M,M2M))
      ALLOCATE (WAVGZ0(M1M,M2M),UIUJAVGZ0(M1M,M2M,6),
     &          PAVGZ0(M1M,M2M),P2AVGZ0(M1M,M2M),
     &          TIME0(M1M,M2M))
      ALLOCATE (URAVGZ1(M1M,M2M),UTAVGZ1(M1M,M2M))
      ALLOCATE (WAVGZ1(M1M,M2M),UIUJAVGZ1(M1M,M2M,6),
     &          PAVGZ1(M1M,M2M),P2AVGZ1(M1M,M2M),
     &          TIME1(M1M,M2M))

      ALLOCATE( PSAVGZ(M1M,M2M),PSRMSZ(M1M,M2M),
     &          PSUIAVGZ(M1M,M2M,3),TALPHAVGZ(M1M,M2M))
      ALLOCATE( PSAVGZ0(M1M,M2M),PSRMSZ0(M1M,M2M),PSUIAVGZ0(M1M,M2M,3))
      ALLOCATE( PSAVGZ1(M1M,M2M),PSRMSZ1(M1M,M2M),PSUIAVGZ1(M1M,M2M,3))
     
      URAVGZ=0.
      UTAVGZ=0.
      WAVGZ=0.
      UIUJAVGZ=0.
      PAVGZ=0.
      P2AVGZ=0.
      PMIAVG=0.
      TNUAVGZ=0.
      SIJAVGZ=0.
      VORAVGZ=0.
      VOLFZ=0.
      SSAVGZ=0.
      SGSDISSZ=0.

      URAVGZ0=0.
      UTAVGZ0=0.
      WAVGZ0=0.
      UIUJAVGZ0=0.
      PAVGZ0=0.
      P2AVGZ0=0.
      TIME0=0.

      URAVGZ1=0.
      UTAVGZ1=0.
      WAVGZ1=0.
      UIUJAVGZ1=0.
      PAVGZ1=0.
      P2AVGZ1=0.
      TIME1=0.

      IF (IPS.EQ.1) THEN
       PSAVGZ=0.
       PSRMSZ=0.
       PSUIAVGZ=0.
       TALPHAVGZ=0.
       
       PSAVGZ0=0.
       PSRMSZ0=0.
       PSUIAVGZ0=0.

       PSAVGZ1=0.
       PSRMSZ1=0.
       PSUIAVGZ1=0.
      ENDIF

      TIMEINIT=TIME
      IHISTINIT=IHIST

      RETURN
      END

C******************************************************************
      SUBROUTINE FIELD_AVG_RT(IHIST,NAV,U,V,W,P,PSI_CN)
C******************************************************************
C     SUBROUTINE FOR FIELD AVERAGING
C     AVERAGED VARIABLES ARE DEFINED AT CELL CENTER
C     VARIABLES ARE LINEARLY INTERPOLATED
C     FROM STAGGERED GRID STRUCTURE TO DETERMINE CELL CENTER VALUES
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE HEAT_VAR

      USE FLD_AVG
      
      USE LES_VAR

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      REAL PSI_CN(0:M1,0:M2,0:M3)

      CHARACTER*16 tname
      CHARACTER*3  tfn1
      CHARACTER*1  tfnh
      CHARACTER*6  tfn2,tfn3
      REAL SR(6)
      real, dimension (:,:), allocatable :: TMP

      tfnh='-'

!$OMP PARALLEL DO private(I,K,UCC,VCC,WCC,AA,COS1,SIN1,UR,UT)
      DO J=1,N2M
      DO I=1,N1M

      DO K=1,N3M
      UCC=0.5*(U(I,J,K)+U(IPV(I),J,K))
      VCC=0.5*(V(I,J,K)+V(I,JPV(J),K))
      WCC=0.5*(W(I,J,K)+W(I,J,KPV(K)))
        AA=1./SQRT(XP(I)**2+YP(J)**2)    !FOR TE CASE OF RT CENTER IS AT DOMAIN CENTER
          COS1=XP(I)*AA
          SIN1=YP(J)*AA
        UR=UCC*COS1+VCC*SIN1
        UT=-UCC*SIN1+VCC*COS1

         URAVGZ(I,J)    =URAVGZ(I,J)    +UR          *DT
         UTAVGZ(I,J)    =UTAVGZ(I,J)    +UT          *DT
         WAVGZ(I,J)     =WAVGZ(I,J)     +WCC         *DT
         UIUJAVGZ(I,J,1)=UIUJAVGZ(I,J,1)+UR**2.      *DT
         UIUJAVGZ(I,J,2)=UIUJAVGZ(I,J,2)+UR*UT       *DT
         UIUJAVGZ(I,J,3)=UIUJAVGZ(I,J,3)+UR*WCC      *DT
         UIUJAVGZ(I,J,4)=UIUJAVGZ(I,J,4)+UT**2.      *DT
         UIUJAVGZ(I,J,5)=UIUJAVGZ(I,J,5)+UT*WCC      *DT
         UIUJAVGZ(I,J,6)=UIUJAVGZ(I,J,6)+WCC**2.     *DT
         PAVGZ(I,J)     =PAVGZ(I,J)     +P(I,J,K)    *DT
         P2AVGZ(I,J)    =P2AVGZ(I,J)    +P(I,J,K)**2.*DT
         VOLFZ(I,J)     =VOLFZ(I,J)     +PSI_CN(I,J,K) *DT
!         TNUAVGZ(I,J)   =TNUAVGZ(I,J)   +TNU(I,J,K)  *DT

      !PSI=0 FLUID PROPERTIES
      IF (PSI_CN(I,J,K) .EQ. 0.) THEN
         URAVGZ0(I,J)    =URAVGZ0(I,J)    +UR          *DT
         UTAVGZ0(I,J)    =UTAVGZ0(I,J)    +UT          *DT
         WAVGZ0(I,J)     =WAVGZ0(I,J)     +WCC         *DT
         UIUJAVGZ0(I,J,1)=UIUJAVGZ0(I,J,1)+UR**2.      *DT
         UIUJAVGZ0(I,J,2)=UIUJAVGZ0(I,J,2)+UR*UT       *DT
         UIUJAVGZ0(I,J,3)=UIUJAVGZ0(I,J,3)+UR*WCC      *DT
         UIUJAVGZ0(I,J,4)=UIUJAVGZ0(I,J,4)+UT**2.      *DT
         UIUJAVGZ0(I,J,5)=UIUJAVGZ0(I,J,5)+UT*WCC      *DT
         UIUJAVGZ0(I,J,6)=UIUJAVGZ0(I,J,6)+WCC**2.     *DT
         PAVGZ0(I,J)     =PAVGZ0(I,J)     +P(I,J,K)    *DT
         P2AVGZ0(I,J)    =P2AVGZ0(I,J)    +P(I,J,K)**2.*DT
         TIME0(I,J)=TIME0(I,J)+DT
      ENDIF

      !PSI=1 FLUID PROPERTIES
      IF (PSI_CN(I,J,K) .EQ. 1.) THEN
         URAVGZ1(I,J)    =URAVGZ1(I,J)    +UR          *DT
         UTAVGZ1(I,J)    =UTAVGZ1(I,J)    +UT          *DT
         WAVGZ1(I,J)     =WAVGZ1(I,J)     +WCC         *DT
         UIUJAVGZ1(I,J,1)=UIUJAVGZ1(I,J,1)+UR**2.      *DT
         UIUJAVGZ1(I,J,2)=UIUJAVGZ1(I,J,2)+UR*UT       *DT
         UIUJAVGZ1(I,J,3)=UIUJAVGZ1(I,J,3)+UR*WCC      *DT
         UIUJAVGZ1(I,J,4)=UIUJAVGZ1(I,J,4)+UT**2.      *DT
         UIUJAVGZ1(I,J,5)=UIUJAVGZ1(I,J,5)+UT*WCC      *DT
         UIUJAVGZ1(I,J,6)=UIUJAVGZ1(I,J,6)+WCC**2.     *DT
         PAVGZ1(I,J)     =PAVGZ1(I,J)     +P(I,J,K)    *DT
         P2AVGZ1(I,J)    =P2AVGZ1(I,J)    +P(I,J,K)**2.*DT
         TIME1(I,J)=TIME1(I,J)+DT
      ENDIF

      ENDDO
      ENDDO
      ENDDO
         PMIAVG         =PMIAVG         +PMI_DUDY         *DT

!      DO K=1,N3M
!      DO J=1,N2M
!      DO I=1,N1M
!
!      VG11=SSDX(I)*(U(I+1,J,K)-U(I,J,K))
!
!      UP=VVDY(J+1)*0.25
!     &   *(SDY(J+1)*(U(I,J,K)+U(I+1,J,K))
!     &    +SDY(J)*(U(I,J+1,K)+U(I+1,J+1,K)))*(1.-FIXJU(J))
!     &   +0.5*(U(I,N2,K)+U(I+1,N2,K))*FIXJU(J)
!      UM=VVDY(J)*0.25
!     &   *(SDY(J)*(U(I,J-1,K)+U(I+1,J-1,K))
!     &    +SDY(J-1)*(U(I,J,K)+U(I+1,J,K)))*(1.-FIXJL(J))
!     &   +0.5*(U(I,0,K)+U(I+1,0,K))*FIXJL(J)
!      VG12=SSDY(J)*(UP-UM)
!
!      UP=VVDZ(K+1)*0.25
!     &   *(SDZ(K+1)*(U(I,J,K)+U(I+1,J,K))
!     &    +SDZ(K)*(U(I,J,K+1)+U(I+1,J,K+1)))*(1.-FIXKU(K))
!     &   +0.5*(U(I,J,N3)+U(I+1,J,N3))*FIXKU(K)
!      UM=VVDZ(K)*0.25
!     &   *(SDZ(K)*(U(I,J,K-1)+U(I+1,J,K-1))
!     &    +SDZ(K-1)*(U(I,J,K)+U(I+1,J,K)))*(1.-FIXKL(K))
!     &   +0.5*(U(I,J,0)+U(I+1,J,0))*FIXKL(K)
!      VG13=SSDZ(K)*(UP-UM)
!
!      VP=VVDX(I+1)*0.25
!     &   *(SDX(I+1)*(V(I,J,K)+V(I,J+1,K))
!     &    +SDX(I)*(V(I+1,J,K)+V(I+1,J+1,K)))*(1.-FIXIU(I))
!     &   +0.5*(V(N1,J,K)+V(N1,J+1,K))*FIXIU(I)
!      VM=VVDX(I)*0.25
!     &   *(SDX(I)*(V(I-1,J,K)+V(I-1,J+1,K))
!     &    +SDX(I-1)*(V(I,J,K)+V(I,J+1,K)))*(1.-FIXIL(I))
!     &   +0.5*(V(0,J,K)+V(0,J+1,K))*FIXIL(I)
!      VG21=SSDX(I)*(VP-VM)
!
!      VG22=SSDY(J)*(V(I,J+1,K)-V(I,J,K))
!
!      VP=VVDZ(K+1)*0.25
!     &   *(SDZ(K+1)*(V(I,J,K)+V(I,J+1,K))
!     &    +SDZ(K)*(V(I,J,K+1)+V(I,J+1,K+1)))*(1.-FIXKU(K))
!     &   +0.5*(V(I,J,N3)+V(I,J+1,N3))*FIXKU(K)
!      VM=VVDZ(K)*0.25
!     &   *(SDZ(K)*(V(I,J,K-1)+V(I,J+1,K-1))
!     &    +SDZ(K-1)*(V(I,J,K)+V(I,J+1,K)))*(1.-FIXKL(K))
!     &   +0.5*(V(I,J,0)+V(I,J+1,0))*FIXKL(K)
!      VG23=SSDZ(K)*(VP-VM)
!
!      WP=VVDX(I+1)*0.25
!     &   *(SDX(I+1)*(W(I,J,K)+W(I,J,K+1))
!     &    +SDX(I)*(W(I+1,J,K)+W(I+1,J,K+1)))*(1.-FIXIU(I))
!     &   +0.5*(W(N1,J,K)+W(N1,J,K+1))*FIXIU(I)
!      WM=VVDX(I)*0.25
!     &   *(SDX(I)*(W(I-1,J,K)+W(I-1,J,K+1))
!     &    +SDX(I-1)*(W(I,J,K)+W(I,J,K+1)))*(1.-FIXIL(I))
!     &   +0.5*(W(0,J,K)+W(0,J,K+1))*FIXIL(I)
!      VG31=SSDX(I)*(WP-WM)
!
!      WP=VVDY(J+1)*0.25
!     &   *(SDY(J+1)*(W(I,J,K)+W(I,J,K+1))
!     &    +SDY(J)*(W(I,J+1,K)+W(I,J+1,K+1)))*(1.-FIXJU(J))
!     &   +0.5*(W(I,N2,K)+W(I,N2,K+1))*FIXJU(J)
!      WM=VVDY(J)*0.25
!     &   *(SDY(J)*(W(I,J-1,K)+W(I,J-1,K+1))
!     &    +SDY(J-1)*(W(I,J,K)+W(I,J,K+1)))*(1.-FIXJL(J))
!     &   +0.5*(W(I,0,K)+W(I,0,K+1))*FIXJL(J)
!      VG32=SSDY(J)*(WP-WM)
!
!      VG33=SSDZ(K)*(W(I,J,K+1)-W(I,J,K))
!
!      SR(1)=VG11
!      SR(2)=0.5*(VG12+VG21)
!      SR(3)=0.5*(VG13+VG31)
!      SR(4)=VG22
!      SR(5)=0.5*(VG23+VG32)
!      SR(6)=VG33
!
!      SIJAVGZ(I,J,1)=SIJAVGZ(I,J,1)+SR(1)*DT
!      SIJAVGZ(I,J,2)=SIJAVGZ(I,J,2)+SR(2)*DT
!      SIJAVGZ(I,J,3)=SIJAVGZ(I,J,3)+SR(3)*DT
!      SIJAVGZ(I,J,4)=SIJAVGZ(I,J,4)+SR(4)*DT
!      SIJAVGZ(I,J,5)=SIJAVGZ(I,J,5)+SR(5)*DT
!      SIJAVGZ(I,J,6)=SIJAVGZ(I,J,6)+SR(6)*DT
!
!      VORAVGZ(I,J,1)=VORAVGZ(I,J,1)+(VG32-VG23)*DT
!      VORAVGZ(I,J,2)=VORAVGZ(I,J,2)+(VG13-VG31)*DT
!      VORAVGZ(I,J,3)=VORAVGZ(I,J,3)+(VG21-VG12)*DT
!
!      SRSR=2.*SR(2)**2.+SR(1)**2.
!     &    +2.*SR(3)**2.+SR(4)**2.
!     &    +2.*SR(5)**2.+SR(6)**2.
!      SSAVGZ(I,J)=SSAVGZ(I,J)+SRSR*DT
!      SGSDISSZ(I,J)=SGSDISSZ(I,J)+SRSR*TNU(I,J,K)*DT
!
!      ENDDO
!      ENDDO
!      ENDDO
!
      IF (IPS.EQ.1) THEN
!$OMP PARALLEL DO private(I,K,UCC,VCC,WCC,AA,COS1,SIN1,UR,UT)
      DO J=1,N2M
      DO I=1,N1M

      DO K=1,N3M
      UCC=0.5*(U(I,J,K)+U(IPV(I),J,K))
      VCC=0.5*(V(I,J,K)+V(I,JPV(J),K))
      WCC=0.5*(W(I,J,K)+W(I,J,KPV(K)))
        AA=1./SQRT(XP(I)**2+YP(J)**2)    !FOR TE CASE OF RT CENTER IS AT DOMAIN CENTER
          COS1=XP(I)*AA
          SIN1=YP(J)*AA
        UR=UCC*COS1+VCC*SIN1
        UT=-UCC*SIN1+VCC*COS1

      PSAVGZ(I,J)    =PSAVGZ(I,J)    +T(I,J,K)    *DT
      PSRMSZ(I,J)    =PSRMSZ(I,J)    +T(I,J,K)**2 *DT
      PSUIAVGZ(I,J,1)=PSUIAVGZ(I,J,1)+T(I,J,K)*UR *DT
      PSUIAVGZ(I,J,2)=PSUIAVGZ(I,J,2)+T(I,J,K)*UT *DT
      PSUIAVGZ(I,J,3)=PSUIAVGZ(I,J,3)+T(I,J,K)*WCC*DT
      
      !PSI=0 FLUID PROPERTIES
      IF (PSI_CN(I,J,K) .EQ. 0.) THEN
      PSAVGZ0(I,J)    =PSAVGZ0(I,J)    +T(I,J,K)    *DT
      PSRMSZ0(I,J)    =PSRMSZ0(I,J)    +T(I,J,K)**2 *DT
      PSUIAVGZ0(I,J,1)=PSUIAVGZ0(I,J,1)+T(I,J,K)*UR *DT
      PSUIAVGZ0(I,J,2)=PSUIAVGZ0(I,J,2)+T(I,J,K)*UT *DT
      PSUIAVGZ0(I,J,3)=PSUIAVGZ0(I,J,3)+T(I,J,K)*WCC*DT
      !!TIME0(I,J)=TIME0(I,J)+DT !ALREADY DONE
      ENDIF

      !PSI=1 FLUID PROPERTIES
      IF (PSI_CN(I,J,K) .EQ. 1.) THEN
      PSAVGZ1(I,J)    =PSAVGZ1(I,J)    +T(I,J,K)    *DT
      PSRMSZ1(I,J)    =PSRMSZ1(I,J)    +T(I,J,K)**2 *DT
      PSUIAVGZ1(I,J,1)=PSUIAVGZ1(I,J,1)+T(I,J,K)*UR *DT
      PSUIAVGZ1(I,J,2)=PSUIAVGZ1(I,J,2)+T(I,J,K)*UT *DT
      PSUIAVGZ1(I,J,3)=PSUIAVGZ1(I,J,3)+T(I,J,K)*WCC*DT
      !!TIME1(I,J)=TIME1(I,J)+DT !ALREADY DONE
      ENDIF
      
      ENDDO
      ENDDO
      ENDDO
      ENDIF

      IF (MOD(NTIME,NPRIAVG).EQ.0) THEN

      TIMEEND=TIME
      IHISTEND=IHIST

      tfn1='fav'
      idg1=ihistinit/100000
      idg2=(ihistinit-idg1*100000)/10000
      idg3=(ihistinit-idg1*100000-idg2*10000)/1000
      idg4=(ihistinit-idg1*100000-idg2*10000-idg3*1000)/100
      idg5=(ihistinit-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
      idg6=ihistinit-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10
      tfn2=char(idg1+48)//char(idg2+48)//
     &     char(idg3+48)//char(idg4+48)//char(idg5+48)//char(idg6+48)
      idg1=ihist/100000
      idg2=(ihist-idg1*100000)/10000
      idg3=(ihist-idg1*100000-idg2*10000)/1000
      idg4=(ihist-idg1*100000-idg2*10000-idg3*1000)/100
      idg5=(ihist-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
      idg6=ihist-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10
      tfn3=char(idg1+48)//char(idg2+48)//
     &     char(idg3+48)//char(idg4+48)//char(idg5+48)//char(idg6+48)
      tname=tfn1//tfn2//tfnh//tfn3

      OPEN(NAV,FILE=tname)
      WRITE(NAV,101)N1M,N2M,N3M,RE
      WRITE(NAV,102)IPS,PRA
      WRITE(NAV,103)TIMEINIT,TIMEEND,IHISTINIT,IHISTEND
      WRITE(NAV,104)(( URAVGZ(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( UTAVGZ(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( WAVGZ(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(((UIUJAVGZ(I,J,L),I=1,N1M),J=1,N2M),L=1,6)
      WRITE(NAV,104)(( PAVGZ(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( P2AVGZ(I,J)    ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( VOLFZ(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,105) PMIAVG
!      IF (IPS.EQ.1) THEN
!      WRITE(NAV)(( PSAVGZ(I,J)    ,I=1,N1M),J=1,N2M)
!      WRITE(NAV)(( PSRMSZ(I,J)    ,I=1,N1M),J=1,N2M)
!      WRITE(NAV)(((PSUIAVGZ(I,J,L),I=1,N1M),J=1,N2M),L=1,3)
!      WRITE(NAV)(( TALPHAVGZ(I,J) ,I=1,N1M),J=1,N2M)
!      ENDIF

      !PSI=0 FLUID PROPERTIES
      WRITE(NAV,104)(( URAVGZ0(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( UTAVGZ0(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( WAVGZ0(I,J)      ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(((UIUJAVGZ0(I,J,L) ,I=1,N1M),J=1,N2M),L=1,6)
      WRITE(NAV,104)(( PAVGZ0(I,J)      ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( P2AVGZ0(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( TIME0(I,J)       ,I=1,N1M),J=1,N2M)

      !PSI=1 FLUID PROPERTIES
      WRITE(NAV,104)(( URAVGZ1(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( UTAVGZ1(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( WAVGZ1(I,J)      ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(((UIUJAVGZ1(I,J,L) ,I=1,N1M),J=1,N2M),L=1,6)
      WRITE(NAV,104)(( PAVGZ1(I,J)      ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( P2AVGZ1(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( TIME1(I,J)       ,I=1,N1M),J=1,N2M)
      
      !HEAT_TRANSFER
      WRITE(NAV,104)(( PSAVGZ(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( PSRMSZ(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(((PSUIAVGZ(I,J,L),I=1,N1M),J=1,N2M),L=1,3)
      
      WRITE(NAV,104)(( PSAVGZ0(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( PSRMSZ0(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(((PSUIAVGZ0(I,J,L),I=1,N1M),J=1,N2M),L=1,3)
      
      WRITE(NAV,104)(( PSAVGZ1(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( PSRMSZ1(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(((PSUIAVGZ1(I,J,L),I=1,N1M),J=1,N2M),L=1,3)
      
      CLOSE(NAV)

      
      NAV=NAV+1
      URAVGZ=0.
      UTAVGZ=0.
      WAVGZ=0.
      UIUJAVGZ=0.
      PAVGZ=0.
      P2AVGZ=0.
      TNUAVGZ=0.
      SIJAVGZ=0.
      VORAVGZ=0.
      SSAVGZ=0.
      VOLFZ=0.
      SGSDISSZ=0.
      PMIAVG=0.

      URAVGZ0=0.
      UTAVGZ0=0.
      WAVGZ0=0.
      UIUJAVGZ0=0.
      PAVGZ0=0.
      P2AVGZ0=0.
      TIME0=0.

      URAVGZ1=0.
      UTAVGZ1=0.
      WAVGZ1=0.
      UIUJAVGZ1=0.
      PAVGZ1=0.
      P2AVGZ1=0.
      TIME1=0.

      IF (IPS.EQ.1) THEN
       PSAVGZ=0.
       PSRMSZ=0.
       PSUIAVGZ=0.
       TALPHAVGZ=0.
       
       PSAVGZ0=0.
       PSRMSZ0=0.
       PSUIAVGZ0=0.

       PSAVGZ1=0.
       PSRMSZ1=0.
       PSUIAVGZ1=0.
      ENDIF

      TIMEINIT=TIME
      IHISTINIT=IHIST

 101  FORMAT(3I10,1ES20.12)
 102  FORMAT(1I8,1ES20.12)
 103  FORMAT(2ES20.12,2I8)
 104  FORMAT(30ES25.12)
 105  FORMAT(ES20.12)

      ENDIF

      RETURN
      END

C******************************************************************
      SUBROUTINE FIELD_AVG_INIT_XYZ(IHIST)
C******************************************************************
C     SUBROUTINE FOR FIELD AVERAGING
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE FLD_AVG

      USE HEAT_VAR

      ALLOCATE (UAVGZ(M1M,M2M),VAVGZ(M1M,M2M))
      ALLOCATE (WAVGZ(M1M,M2M),UIUJAVGZ(M1M,M2M,6),
     &          PAVGZ(M1M,M2M),P2AVGZ(M1M,M2M),
     &          TNUAVGZ(M1M,M2M),VOLFZ(M1M,M2M))
!     &         SIJAVGZ(M1M,M2M,6),VORAVGZ(M1M,M2M,3),
!     &         SSAVGZ(M1M,M2M),SGSDISSZ(M1M,M2M),
      ALLOCATE (UAVGZ0(M1M,M2M),VAVGZ0(M1M,M2M))
      ALLOCATE (WAVGZ0(M1M,M2M),UIUJAVGZ0(M1M,M2M,6),
     &          PAVGZ0(M1M,M2M),P2AVGZ0(M1M,M2M),
     &          TIME0(M1M,M2M))
      ALLOCATE (UAVGZ1(M1M,M2M),VAVGZ1(M1M,M2M))
      ALLOCATE (WAVGZ1(M1M,M2M),UIUJAVGZ1(M1M,M2M,6),
     &          PAVGZ1(M1M,M2M),P2AVGZ1(M1M,M2M),
     &          TIME1(M1M,M2M))
!      ALLOCATE( PSAVGZ(M1M,M2M),PSRMSZ(M1M,M2M),
!     &          PSUIAVGZ(M1M,M2M,3),TALPHAVGZ(M1M,M2M))

      UAVGZ=0.
      VAVGZ=0.
      WAVGZ=0.
      UIUJAVGZ=0.
      PAVGZ=0.
      P2AVGZ=0.
      PMIAVG=0.
      TNUAVGZ=0.
      SIJAVGZ=0.
      VORAVGZ=0.
      VOLFZ=0.
      SSAVGZ=0.
      SGSDISSZ=0.

      UAVGZ0=0.
      VAVGZ0=0.
      WAVGZ0=0.
      UIUJAVGZ0=0.
      PAVGZ0=0.
      P2AVGZ0=0.
      TIME0=0.

      UAVGZ1=0.
      VAVGZ1=0.
      WAVGZ1=0.
      UIUJAVGZ1=0.
      PAVGZ1=0.
      P2AVGZ1=0.
      TIME1=0.

      IF (IPS.EQ.1) PSAVGZ=0.
      IF (IPS.EQ.1) PSRMSZ=0.
      IF (IPS.EQ.1) PSUIAVGZ=0.
      IF (IPS.EQ.1) TALPHAVGZ=0.
      TIMEINIT=TIME
      IHISTINIT=IHIST

      RETURN
      END

C******************************************************************
      SUBROUTINE FIELD_AVG_XYZ(IHIST,NAV,U,V,W,P,PSI_CN)
C******************************************************************
C     SUBROUTINE FOR FIELD AVERAGING
C     AVERAGED VARIABLES ARE DEFINED AT CELL CENTER
C     VARIABLES ARE LINEARLY INTERPOLATED
C     FROM STAGGERED GRID STRUCTURE TO DETERMINE CELL CENTER VALUES
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE HEAT_VAR

      USE FLD_AVG
      
      USE LES_VAR

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL P(0:M1,0:M2,0:M3)

      REAL PSI_CN(0:M1,0:M2,0:M3)

      CHARACTER*16 tname
      CHARACTER*3  tfn1
      CHARACTER*1  tfnh
      CHARACTER*6  tfn2,tfn3
      REAL SR(6)
      real, dimension (:,:), allocatable :: TMP

      tfnh='-'

!$OMP PARALLEL DO private(I,J,UCC,VCC,WCC)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
      UCC=0.5*(U(I,J,K)+U(IPV(I),J,K))
      VCC=0.5*(V(I,J,K)+V(I,JPV(J),K))
      WCC=0.5*(W(I,J,K)+W(I,J,KPV(K)))
         UAVGZ(I,J)     =UAVGZ(I,J)     +UCC         *DT
         VAVGZ(I,J)     =VAVGZ(I,J)     +VCC         *DT
         WAVGZ(I,J)     =WAVGZ(I,J)     +WCC         *DT
         UIUJAVGZ(I,J,1)=UIUJAVGZ(I,J,1)+UCC**2.     *DT
         UIUJAVGZ(I,J,2)=UIUJAVGZ(I,J,2)+UCC*VCC     *DT
         UIUJAVGZ(I,J,3)=UIUJAVGZ(I,J,3)+UCC*WCC     *DT
         UIUJAVGZ(I,J,4)=UIUJAVGZ(I,J,4)+VCC**2.     *DT
         UIUJAVGZ(I,J,5)=UIUJAVGZ(I,J,5)+VCC*WCC     *DT
         UIUJAVGZ(I,J,6)=UIUJAVGZ(I,J,6)+WCC**2.     *DT
         PAVGZ(I,J)     =PAVGZ(I,J)     +P(I,J,K)    *DT
         P2AVGZ(I,J)    =P2AVGZ(I,J)    +P(I,J,K)**2.*DT
         VOLFZ(I,J)     =VOLFZ(I,J)     +PSI_CN(I,J,K) *DT

      !PSI=0 FLUID PROPERTIES
      IF (PSI_CN(I,J,K) .EQ. 0.) THEN
         UAVGZ0(I,J)     =UAVGZ0(I,J)     +UCC         *DT
         VAVGZ0(I,J)     =VAVGZ0(I,J)     +VCC         *DT
         WAVGZ0(I,J)     =WAVGZ0(I,J)     +WCC         *DT
         UIUJAVGZ0(I,J,1)=UIUJAVGZ0(I,J,1)+UCC**2.     *DT
         UIUJAVGZ0(I,J,2)=UIUJAVGZ0(I,J,2)+UCC*VCC     *DT
         UIUJAVGZ0(I,J,3)=UIUJAVGZ0(I,J,3)+UCC*WCC     *DT
         UIUJAVGZ0(I,J,4)=UIUJAVGZ0(I,J,4)+VCC**2.     *DT
         UIUJAVGZ0(I,J,5)=UIUJAVGZ0(I,J,5)+VCC*WCC     *DT
         UIUJAVGZ0(I,J,6)=UIUJAVGZ0(I,J,6)+WCC**2.     *DT
         PAVGZ0(I,J)     =PAVGZ0(I,J)     +P(I,J,K)    *DT
         P2AVGZ0(I,J)    =P2AVGZ0(I,J)    +P(I,J,K)**2.*DT
         TIME0(I,J)=TIME0(I,J)+DT
      ENDIF

      !PSI=1 FLUID PROPERTIES
      IF (PSI_CN(I,J,K) .EQ. 1.) THEN
         UAVGZ1(I,J)     =UAVGZ1(I,J)     +UCC         *DT
         VAVGZ1(I,J)     =VAVGZ1(I,J)     +VCC         *DT
         WAVGZ1(I,J)     =WAVGZ1(I,J)     +WCC         *DT
         UIUJAVGZ1(I,J,1)=UIUJAVGZ1(I,J,1)+UCC**2.     *DT
         UIUJAVGZ1(I,J,2)=UIUJAVGZ1(I,J,2)+UCC*VCC     *DT
         UIUJAVGZ1(I,J,3)=UIUJAVGZ1(I,J,3)+UCC*WCC     *DT
         UIUJAVGZ1(I,J,4)=UIUJAVGZ1(I,J,4)+VCC**2.     *DT
         UIUJAVGZ1(I,J,5)=UIUJAVGZ1(I,J,5)+VCC*WCC     *DT
         UIUJAVGZ1(I,J,6)=UIUJAVGZ1(I,J,6)+WCC**2.     *DT
         PAVGZ1(I,J)     =PAVGZ1(I,J)     +P(I,J,K)    *DT
         P2AVGZ1(I,J)    =P2AVGZ1(I,J)    +P(I,J,K)**2.*DT
         TIME1(I,J)=TIME1(I,J)+DT
      ENDIF

      ENDDO
      ENDDO
      ENDDO
         PMIAVG         =PMIAVG         +PMI_DUDY         *DT

!      IF (IPS.EQ.1) THEN
!      DO K=1,N3M
!      DO J=1,N2M
!      DO I=1,N1M
!      UCC=0.5*(U(I,J,K)+U(IPV(I),J,K))
!      VCC=0.5*(V(I,J,K)+V(I,JPV(J),K))
!      WCC=0.5*(W(I,J,K)+W(I,J,KPV(K)))
!      PSAVGZ(I,J)    =PSAVGZ(I,J)    +T(I,J,K)    *DT
!      PSRMSZ(I,J)    =PSRMSZ(I,J)    +T(I,J,K)**2 *DT
!      PSUIAVGZ(I,J,1)=PSUIAVGZ(I,J,1)+T(I,J,K)*UCC*DT
!      PSUIAVGZ(I,J,2)=PSUIAVGZ(I,J,2)+T(I,J,K)*VCC*DT
!      PSUIAVGZ(I,J,3)=PSUIAVGZ(I,J,3)+T(I,J,K)*WCC*DT
!      TALPHAVGZ(I,J) =TALPHAVGZ(I,J) +TALPH(I,J,K)*DT
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDIF

      IF (MOD(NTIME,NPRIAVG).EQ.0) THEN

      TIMEEND=TIME
      IHISTEND=IHIST

      tfn1='fav'
      idg1=ihistinit/100000
      idg2=(ihistinit-idg1*100000)/10000
      idg3=(ihistinit-idg1*100000-idg2*10000)/1000
      idg4=(ihistinit-idg1*100000-idg2*10000-idg3*1000)/100
      idg5=(ihistinit-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
      idg6=ihistinit-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10
      tfn2=char(idg1+48)//char(idg2+48)//
     &     char(idg3+48)//char(idg4+48)//char(idg5+48)//char(idg6+48)
      idg1=ihist/100000
      idg2=(ihist-idg1*100000)/10000
      idg3=(ihist-idg1*100000-idg2*10000)/1000
      idg4=(ihist-idg1*100000-idg2*10000-idg3*1000)/100
      idg5=(ihist-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
      idg6=ihist-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10
      tfn3=char(idg1+48)//char(idg2+48)//
     &     char(idg3+48)//char(idg4+48)//char(idg5+48)//char(idg6+48)
      tname=tfn1//tfn2//tfnh//tfn3

      OPEN(NAV,FILE=tname)
      WRITE(NAV,101)N1M,N2M,N3M,RE
      WRITE(NAV,102)IPS,PRA
      WRITE(NAV,103)TIMEINIT,TIMEEND,IHISTINIT,IHISTEND
      WRITE(NAV,104)(( UAVGZ(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( VAVGZ(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( WAVGZ(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(((UIUJAVGZ(I,J,L),I=1,N1M),J=1,N2M),L=1,6)
      WRITE(NAV,104)(( PAVGZ(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( P2AVGZ(I,J)    ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( VOLFZ(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,105) PMIAVG    
!      IF (IPS.EQ.1) THEN
!      WRITE(NAV)(( PSAVGZ(I,J)    ,I=1,N1M),J=1,N2M)
!      WRITE(NAV)(( PSRMSZ(I,J)    ,I=1,N1M),J=1,N2M)
!      WRITE(NAV)(((PSUIAVGZ(I,J,L),I=1,N1M),J=1,N2M),L=1,3)
!      WRITE(NAV)(( TALPHAVGZ(I,J) ,I=1,N1M),J=1,N2M)
!      ENDIF

      !PSI=0 FLUID PROPERTIES
      WRITE(NAV,104)(( UAVGZ0(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( VAVGZ0(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( WAVGZ0(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(((UIUJAVGZ0(I,J,L),I=1,N1M),J=1,N2M),L=1,6)
      WRITE(NAV,104)(( PAVGZ0(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( P2AVGZ0(I,J)    ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( TIME0(I,J)       ,I=1,N1M),J=1,N2M)

      !PSI=1 FLUID PROPERTIES
      WRITE(NAV,104)(( UAVGZ1(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( VAVGZ1(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( WAVGZ1(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(((UIUJAVGZ1(I,J,L),I=1,N1M),J=1,N2M),L=1,6)
      WRITE(NAV,104)(( PAVGZ1(I,J)     ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( P2AVGZ1(I,J)    ,I=1,N1M),J=1,N2M)
      WRITE(NAV,104)(( TIME1(I,J)       ,I=1,N1M),J=1,N2M)
      CLOSE(NAV)
      NAV=NAV+1

      UAVGZ=0.
      VAVGZ=0.
      WAVGZ=0.
      UIUJAVGZ=0.
      PAVGZ=0.
      P2AVGZ=0.
      PMIAVG=0.
      TNUAVGZ=0.
      SIJAVGZ=0.
      VORAVGZ=0.
      VOLFZ=0.
      SSAVGZ=0.
      SGSDISSZ=0.

      UAVGZ0=0.
      VAVGZ0=0.
      WAVGZ0=0.
      UIUJAVGZ0=0.
      PAVGZ0=0.
      P2AVGZ0=0.
      TIME0=0.

      UAVGZ1=0.
      VAVGZ1=0.
      WAVGZ1=0.
      UIUJAVGZ1=0.
      PAVGZ1=0.
      P2AVGZ1=0.
      TIME1=0.

      IF (IPS.EQ.1) PSAVGZ=0.
      IF (IPS.EQ.1) PSRMSZ=0.
      IF (IPS.EQ.1) PSUIAVGZ=0.
      IF (IPS.EQ.1) TALPHAVGZ=0.
      TIMEINIT=TIME
      IHISTINIT=IHIST

 101  FORMAT(3I10,1ES20.12)
 102  FORMAT(1I8,1ES20.12)
 103  FORMAT(2ES20.12,2I8)
 104  FORMAT(30ES25.12)
 105  FORMAT(ES20.12)

      ENDIF

      RETURN
      END

!*****************************************************************
!
!! REAL_TIME returns the real time in seconds.
!
!  Modified:
!
!    07 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) SECONDS, a reading of the real time clock,
!    in seconds.
!
!      subroutine real_time(seconds)
!
!        implicit none
!
!        integer(8) :: clock_count
!        integer(8) :: clock_max
!        integer(8) :: clock_rate
!        real seconds
!
!        call system_clock ( clock_count, clock_rate, clock_max )
!
!        seconds = real (clock_count)
!     &    / real (clock_rate)
!
!!        write(*,*)clock_count, clock_rate
!
!        return
!      end

      subroutine real_time(seconds)

      implicit none
      integer clock_count,clock_max,clock_rate
      real seconds

      call system_clock ( clock_count, clock_rate, clock_max )

      seconds=real(clock_count)/real(clock_rate)

      return
      end


!*******************************************************************************
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
      subroutine timestamp
        implicit none

        character ( len = 8 ) ampm
        integer d
        integer h
        integer m
        integer mm
        character ( len = 9 ), parameter, dimension(12) :: month = (/
     &    'January  ', 'February ', 'March    ', 'April    ',
     &    'May      ', 'June     ', 'July     ', 'August   ',
     &    'September', 'October  ', 'November ', 'December ' /)
        integer n
        integer s
        integer values(8)
        integer y

        call date_and_time ( values = values )

        y = values(1)
        m = values(2)
        d = values(3)
        h = values(5)
        n = values(6)
        s = values(7)
        mm = values(8)

        if ( h < 12 ) then
          ampm = 'AM'
        else if ( h == 12 ) then
          if ( n == 0 .and. s == 0 ) then
            ampm = 'Noon'
          else
            ampm = 'PM'
          end if
        else
          h = h - 12
          if ( h < 12 ) then
            ampm = 'PM'
          else if ( h == 12 ) then
            if ( n == 0 .and. s == 0 ) then
              ampm = 'Midnight'
            else
              ampm = 'AM'
            end if
          end if
        end if

      write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     &d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )
      write(*,*)' '

        return
      end

!C*******************************************************************
!        SUBROUTINE PRTINIT
!C*******************************************************************
!C       read parameters about particle tracking
!        INCLUDE 'particle.h'
!        CHARACTER*30 DUMMY
!
!        PI=ACOS(-1.)
!
!        OPEN(14,FILE='particle.in')
!        READ(14,*) DUMMY
!        READ(14,*) NSEEDST,NWPRTST
!        READ(14,*) DUMMY
!        READ(14,*) NRSEED,NTSEED,NSEEDST
!        READ(14,*) DUMMY
!        READ(14,*) IRPRT,IWPRT,NWPRTST
!
!
!        PRINT*, 'XSEED=',XSEED,' RSEED=',RSEED,' TSEED=',TSEED
!        PRINT*, 'NRSEED=',NRSEED,' NTSEED=',NTSEED,' NSEEDST=',NSEEDST
!        PRINT*, 'IRPRT=',IRPRT,' NWPRT=',IWPRT,' NWPRTST=',NWPRTST
!
!        TSEED=TSEED*PI/180.
!        NPRT=NRSEED*NTSEED
!        NPRTFILE=901  ! initial file name for particle position
!
!        IF (IRPRT.EQ.1) THEN
!           CALL READPRT
!        ELSE
!           NSEEDING=0
!           CALL SEEDPRT
!        ENDIF
!
!        RETURN
!        END
!
!C*******************************************************************
!        SUBROUTINE READPRT
!C*******************************************************************
!C       read particle positions
!        INCLUDE 'particle.h'
!
!        OPEN(50,FILE='part.bin',FORM='UNFORMATTED')
!        READ(50) NPRT,NSEEDING
!        READ(50) (((XPRT(I,J,K),I=1,NPRT),J=1,NSEEDING),K=1,3)
!        READ(50) NOUTPRT,NIBMPRT
!        READ(50) NRSEED,NTSEED
!        CLOSE(50)
!
!        PRINT*, 'NPRT = ',NPRT,' NSEEDING = ',NSEEDING
!        PRINT*, 'NOUTPRT = ',NOUTPRT,' NIBMPRT = ',NIBMPRT
!        PRINT*, 'NTSEED = ',NTSEED,' NRSEED = ',NRSEED
!
!        RETURN
!        END
!
!C*******************************************************************
!        SUBROUTINE SEEDPRT
!C*******************************************************************
!C       set initial positons of particles
!        INCLUDE 'particle.h'
!
!        PI=ACOS(-1.)
!
!        NSEEDING=NSEEDING+1
!
!C       unform distribution in radial & azimuthal directions
!        DRSEED=RSEED/FLOAT(NRSEED)
!        DTSEED=2.*PI/FLOAT(NTSEED)
!
!        II=0
!        DO 200 K=1,NTSEED
!        DO 200 J=1,NRSEED
!        II=II+1
!        XPRT(II,NSEEDING,1)=XSEED
!        XPRT(II,NSEEDING,2)=DRSEED*FLOAT(J)
!c        XPRT(II,NSEEDING,3)=DTSEED*FLOAT(K-1)
!        XPRT(II,NSEEDING,3)=DTSEED*FLOAT(K-1)+TSEED
!C#################
!C       seeding position (K=1) is aligned with lift direction
!        IF (XPRT(II,NSEEDING,3).GT.(2.*PI)) THEN
!           XPRT(II,NSEEDING,3)=XPRT(II,NSEEDING,3)-2.*PI
!        ELSE IF (XPRT(II,NSEEDING,3).LT.0.0) THEN
!           XPRT(II,NSEEDING,3)=XPRT(II,NSEEDING,3)+2.*PI
!        ENDIF
!C#################
! 200    CONTINUE
!
!c        PRINT*, 'NPRT = ',NPRT,' NSEEDING = ',NSEEDING
!c        DO 500 J=1,NSEEDING
!c        DO 500 I=1,NPRT
!c        WRITE(*,511) XPRT(I,J,1),XPRT(I,J,2),XPRT(I,J,3)
!c 500    CONTINUE
!c 511    FORMAT(3E13.5)
!
!
!        RETURN
!        END
!
!C*******************************************************************
!        SUBROUTINE WRTPRT
!C*******************************************************************
!C       write final particle positions
!        INCLUDE 'particle.h'
!
!        OPEN(NPRTFILE,FORM='UNFORMATTED')
!        WRITE(NPRTFILE) NPRT,NSEEDING
!        WRITE(NPRTFILE) (((XPRT(I,J,K),I=1,NPRT),J=1,NSEEDING),K=1,3)
!        WRITE(NPRTFILE) (((VELPRT(I,J,K),I=1,NPRT),J=1,NSEEDING),K=1,3)
!        WRITE(NPRTFILE) NOUTPRT,NIBMPRT
!        WRITE(NPRTFILE) NRSEED,NTSEED
!        CLOSE(NPRTFILE)
!
!        NPRTFILE=NPRTFILE+1
!
!
!
!        RETURN
!        END
!
!C*******************************************************************
!        SUBROUTINE MOVEPRT
!C*******************************************************************
!C       move the particles using the velocity
!        INCLUDE 'common.h'
!        INCLUDE 'particle.h'
!
!        NOUTPRT=0
!        NIBMPRT=0
!
!
!        DO 100 J=1,NSEEDING
!        DO 100 I=1,NPRT
!
!C------ fix the position of the particle getting out of the comput. domain
!        IF (XPRT(I,J,1) .GT. XP(N1M)) THEN
!           NOUTPRT=NOUTPRT+1
!c           XPRT(I,J,1)=XL*0.5
!c           XPRT(I,J,2)=0.
!c           XPRT(I,J,3)=0.
!           GO TO 999
!        ENDIF
!        IF (XPRT(I,J,2) .GT. RP(N2M)) THEN
!           NOUTPRT=NOUTPRT+1
!c           XPRT(I,J,1)=0.
!c           XPRT(I,J,2)=RL
!c           XPRT(I,J,3)=0.
!           GO TO 999
!        ENDIF
!
!C------ fix the position of the particle entering the immersed body
!C       because the velocity is not exactly zero on the immersed boundary
!C       sphere
!        RDIST=SQRT(XPRT(I,J,1)**2+XPRT(I,J,2)**2)
!        IF (RDIST .LE. 0.5) THEN
!           NIBMPRT=NIBMPRT+1
!c           XPRT(I,J,1)=0.
!c           XPRT(I,J,2)=0.
!c           XPRT(I,J,3)=0.
!           GO TO 999
!        ENDIF
!
!C------ compute the geometric factors
!C       X-direction: non-uniform grid assumed
!        DO 110 II=1,N1M-1
!        IF (XP(II+1) .GT. XPRT(I,J,1)) THEN
!           FXGEO=(XPRT(I,J,1)-XP(II))/(XP(II+1)-XP(II))
!           IGEO=II
!           GO TO 111
!        ENDIF
! 110    CONTINUE
! 111    CONTINUE
!
!C       R-direction: non-uniform grid assumed
!        DO 120 JJ=0,N2M-1
!        IF (RP(JJ+1) .GT. XPRT(I,J,2)) THEN
!           FRGEO=(XPRT(I,J,2)-RP(JJ))/(RP(JJ+1)-RP(JJ))
!           JGEO=JJ
!           GO TO 121
!        ENDIF
! 120    CONTINUE
! 121    CONTINUE
!
!C       T-direction: uniform grid assumed
!        KGEO=INT(XPRT(I,J,3)/DTL)+1
!        FTGEO=(XPRT(I,J,3)-TP(KGEO))/DTL
!
!C------ index
!        ICC=IGEO
!        IP1=IP(ICC)
!        IP2=IP(IP1)
!        JCC=JGEO
!        JP1=JP(JCC)
!        JP2=JP(JP1)
!        KCC=KGEO
!        KP1=KP(KCC)
!        KP2=KP(KP1)
!
!C------ UX
!C       UX1~UX8: cell-center streamwise velocity
!        UX1=(UX(ICC,JCC,KCC)+UX(IP1,JCC,KCC))*0.5
!        UX2=(UX(IP1,JCC,KCC)+UX(IP2,JCC,KCC))*0.5
!        UX3=(UX(ICC,JP1,KCC)+UX(IP1,JP1,KCC))*0.5
!        UX4=(UX(IP1,JP1,KCC)+UX(IP2,JP1,KCC))*0.5
!        UX5=(UX(ICC,JCC,KP1)+UX(IP1,JCC,KP1))*0.5
!        UX6=(UX(IP1,JCC,KP1)+UX(IP2,JCC,KP1))*0.5
!        UX7=(UX(ICC,JP1,KP1)+UX(IP1,JP1,KP1))*0.5
!        UX8=(UX(IP1,JP1,KP1)+UX(IP2,JP1,KP1))*0.5
!        UXPRT=(1.-FXGEO)*(1.-FRGEO)*(1.-FTGEO)*UX1
!     &       +(   FXGEO)*(1.-FRGEO)*(1.-FTGEO)*UX2
!     &       +(1.-FXGEO)*(   FRGEO)*(1.-FTGEO)*UX3
!     &       +(   FXGEO)*(   FRGEO)*(1.-FTGEO)*UX4
!     &       +(1.-FXGEO)*(1.-FRGEO)*(   FTGEO)*UX5
!     &       +(   FXGEO)*(1.-FRGEO)*(   FTGEO)*UX6
!     &       +(1.-FXGEO)*(   FRGEO)*(   FTGEO)*UX7
!     &       +(   FXGEO)*(   FRGEO)*(   FTGEO)*UX8
!
!C------ UR
!C       UR1~UR8: cell-center radial velocity
!        IF (JCC.NE.0) THEN
!        UR1=(UR(ICC,JCC,KCC)+UR(ICC,JP1,KCC))*0.5
!        UR2=(UR(IP1,JCC,KCC)+UR(IP1,JP1,KCC))*0.5
!        UR3=(UR(ICC,JP1,KCC)+UR(ICC,JP2,KCC))*0.5
!        UR4=(UR(IP1,JP1,KCC)+UR(IP1,JP2,KCC))*0.5
!        UR5=(UR(ICC,JCC,KP1)+UR(ICC,JP1,KP1))*0.5
!        UR6=(UR(IP1,JCC,KP1)+UR(IP1,JP1,KP1))*0.5
!        UR7=(UR(ICC,JP1,KP1)+UR(ICC,JP2,KP1))*0.5
!        UR8=(UR(IP1,JP1,KP1)+UR(IP1,JP2,KP1))*0.5
!        ELSE
!C##############
!C       CAUTION: UR is not defined when J=0
!        UR1=UR(ICC,1,KCC)
!        UR2=UR(IP1,1,KCC)
!        UR3=(UR(ICC,1,KCC)+UR(ICC,2,KCC))*0.5
!        UR4=(UR(IP1,1,KCC)+UR(IP1,2,KCC))*0.5
!        UR5=UR(ICC,1,KP1)
!        UR6=UR(IP1,1,KP1)
!        UR7=(UR(ICC,1,KP1)+UR(ICC,2,KP1))*0.5
!        UR8=(UR(IP1,1,KP1)+UR(IP1,2,KP1))*0.5
!C##############
!        ENDIF
!        URPRT=(1.-FXGEO)*(1.-FRGEO)*(1.-FTGEO)*UR1
!     &       +(   FXGEO)*(1.-FRGEO)*(1.-FTGEO)*UR2
!     &       +(1.-FXGEO)*(   FRGEO)*(1.-FTGEO)*UR3
!     &       +(   FXGEO)*(   FRGEO)*(1.-FTGEO)*UR4
!     &       +(1.-FXGEO)*(1.-FRGEO)*(   FTGEO)*UR5
!     &       +(   FXGEO)*(1.-FRGEO)*(   FTGEO)*UR6
!     &       +(1.-FXGEO)*(   FRGEO)*(   FTGEO)*UR7
!     &       +(   FXGEO)*(   FRGEO)*(   FTGEO)*UR8
!
!C------ UT
!C       UT1~UT8: cell-center azimuthal velocity
!        UT1=(UT(ICC,JCC,KCC)+UT(ICC,JCC,KP1))*0.5
!        UT2=(UT(IP1,JCC,KCC)+UT(IP1,JCC,KP1))*0.5
!        UT3=(UT(ICC,JP1,KCC)+UT(ICC,JP1,KP1))*0.5
!        UT4=(UT(IP1,JP1,KCC)+UT(IP1,JP1,KP1))*0.5
!        UT5=(UT(ICC,JCC,KP1)+UT(ICC,JCC,KP2))*0.5
!        UT6=(UT(IP1,JCC,KP1)+UT(IP1,JCC,KP2))*0.5
!        UT7=(UT(ICC,JP1,KP1)+UT(ICC,JP1,KP2))*0.5
!        UT8=(UT(IP1,JP1,KP1)+UT(IP1,JP1,KP2))*0.5
!        UTPRT=(1.-FXGEO)*(1.-FRGEO)*(1.-FTGEO)*UT1
!     &       +(   FXGEO)*(1.-FRGEO)*(1.-FTGEO)*UT2
!     &       +(1.-FXGEO)*(   FRGEO)*(1.-FTGEO)*UT3
!     &       +(   FXGEO)*(   FRGEO)*(1.-FTGEO)*UT4
!     &       +(1.-FXGEO)*(1.-FRGEO)*(   FTGEO)*UT5
!     &       +(   FXGEO)*(1.-FRGEO)*(   FTGEO)*UT6
!     &       +(1.-FXGEO)*(   FRGEO)*(   FTGEO)*UT7
!     &       +(   FXGEO)*(   FRGEO)*(   FTGEO)*UT8
!C        UTPRT=UTPRT/XPRT(I,J,2)    ! angular velocity
!
!        X_PRT=XPRT(I,J,1)
!        Y_PRT=XPRT(I,J,2)*COS(XPRT(I,J,3))
!        Z_PRT=XPRT(I,J,2)*SIN(XPRT(I,J,3))
!        UYPRT=URPRT*COS(XPRT(I,J,3))-UTPRT*SIN(XPRT(I,J,3))
!        UZPRT=URPRT*SIN(XPRT(I,J,3))+UTPRT*COS(XPRT(I,J,3))
!
!        X_PRT=X_PRT+GAMMA(ISUB)*DT*UXPRT
!     &                         +RO(ISUB)*DT*VELPRT(I,J,1)
!        Y_PRT=Y_PRT+GAMMA(ISUB)*DT*UYPRT
!     &                         +RO(ISUB)*DT*VELPRT(I,J,2)
!        Z_PRT=Z_PRT+GAMMA(ISUB)*DT*UZPRT
!     &                         +RO(ISUB)*DT*VELPRT(I,J,3)
!        R_PRT=SQRT(Y_PRT**2+Z_PRT**2)
!        T_PRT=ACOS(Y_PRT/R_PRT)
!C############## CAUTION
!        IF (Z_PRT.LT.0.) T_PRT=2.*PI-T_PRT
!C############## CAUTION
!
!        XPRT(I,J,1)=X_PRT
!        XPRT(I,J,2)=R_PRT
!        XPRT(I,J,3)=T_PRT
!
!        VELPRT(I,J,1)=UXPRT
!        VELPRT(I,J,2)=UYPRT
!        VELPRT(I,J,3)=UZPRT
!
!
! 999    CONTINUE
! 100    CONTINUE
!
!
!        RETURN
!        END
