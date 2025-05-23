MODULE CMF_CTRL_LEVEE_MOD
!==========================================================
!* PURPOSE: CaMa-Flood flood stage with levee scheme
!
! (C) Y. Tanaka & D.Yamazaki (U-Tokyo)  Mar 2022
!
!* CONTAINS:
! -- CMF_LEV_NMLIST      : Read setting from namelist
! -- CMF_LEV_INIT        : Initialize levee scheme
! -- CMF_CALC_FLDSTG_LEV : Calculate flood stage considering levee
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
USE PARKIND1,                ONLY: JPIM, JPRB, JPRM, JPRD
USE YOS_CMF_INPUT,           ONLY: LOGNAM, IMIS
!============================
IMPLICIT NONE
SAVE
!*** NAMELIST/NLEVEE/
!*** Levee Map
CHARACTER(LEN=256)             ::  CLEVHGT         !! LEVEE HEIGHT from RIVER
CHARACTER(LEN=256)             ::  CLEVFRC         !! Unprotected fraction. Relative Levee distance from RIVER

NAMELIST/NLEVEE/  CLEVHGT, CLEVFRC

!*** Levee Parameters from map
REAL(KIND=JPRB),ALLOCATABLE    ::  D2LEVHGT(:,:)        !! LEVEE HEIGHT [M] (levee croen elevation above elevtn.bin-river elevation)
REAL(KIND=JPRB),ALLOCATABLE    ::  D2LEVFRC(:,:)        !! Unprotected fraction = RELATIVE DISTANCE between LEVEE and RIVER [0-1].
                                                        !!  0 = just aside channel, 1 = edge of catchment

!*** Levee stage parameter (calculated)
REAL(KIND=JPRB),ALLOCATABLE    ::  D2BASHGT(:,:)        !! LEVEE Base height [M] (levee base elevation above elevtn.bin-river elev)
REAL(KIND=JPRB),ALLOCATABLE    ::  D2LEVDST(:,:)        !! Absolute DISTANCE between LEVEE and RIVER [0-1]. 
                                                        !! 0 = just aside channel, 1 = edge of catchment

REAL(KIND=JPRB),ALLOCATABLE    ::  D2LEVBASSTO(:,:)  !! MAXIMUM STORAGE under LEVEE BASE [M3]
REAL(KIND=JPRB),ALLOCATABLE    ::  D2LEVTOPSTO(:,:)  !! MAXIMUM STORAGE at LEVEE TOP [M3] (only river side)
REAL(KIND=JPRB),ALLOCATABLE    ::  D2LEVFILSTO(:,:) !! MAXIMUM STORAGE at LEVEE TOP [M3] (both river & protected side are filled)

CONTAINS
!####################################################################
!* CONTAINS:
! -- CMF_LEVEE_NMLIST         : Read setting from namelist
! -- CMF_LEVEE_INIT           : Initialize dam data
! -- CMF_LEVEE_FLDSTG         : Calculate inflow and outflow at dam
! -- CMF_LEVEE_OPT_PTHOUT     : Bifurcation scheme with levee consideration
!####################################################################
SUBROUTINE CMF_LEVEE_NMLIST
! reed setting from namelist
! -- Called from CMF_DRV_NMLIST
USE YOS_CMF_INPUT,      ONLY: CSETFILE,NSETFILE
USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID
IMPLICIT NONE
!================================================
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!---------------------!"

!*** 1. open namelist
NSETFILE=INQUIRE_FID()
OPEN(NSETFILE,FILE=CSETFILE,STATUS="OLD")
WRITE(LOGNAM,*) "CMF::LEVEE_NMLIST: namelist OPEN in unit: ", TRIM(CSETFILE), NSETFILE 

!*** 2. default value
CLEVHGT   ="NONE"
CLEVFRC   ="NONE"

!*** 3. read namelist
REWIND(NSETFILE)
READ(NSETFILE,NML=NLEVEE)

WRITE(LOGNAM,*)   "=== NAMELIST, NLEVEE ==="
WRITE(LOGNAM,*)   "CLEVHGT   : ", CLEVHGT   
WRITE(LOGNAM,*)   "CLEVFRC   : ", CLEVFRC   

CLOSE(NSETFILE)

WRITE(LOGNAM,*) "CMF::LEVEE_NMLIST: end" 

END SUBROUTINE CMF_LEVEE_NMLIST
!####################################################################





!####################################################################
SUBROUTINE CMF_LEVEE_INIT
USE YOS_CMF_INPUT,      ONLY: TMPNAM, NX, NY, NLFP
USE YOS_CMF_MAP,        ONLY: NSEQALL
USE YOS_CMF_MAP,        ONLY: D2GRAREA, D2RIVLEN, D2RIVWTH, D2FLDHGT, &
                             & D2FLDGRD, D2RIVSTOMAX, D2FLDSTOMAX, DFRCINC
USE CMF_UTILS_MOD,      ONLY: mapR2vecD, INQUIRE_FID
!
IMPLICIT NONE
!* local variables
REAL(KIND=JPRM)            ::  R2TEMP(NX,NY)
! SAVE for OpenMP
INTEGER(KIND=JPIM),SAVE    ::  ISEQ, I, ILEV
REAL(KIND=JPRB),SAVE       ::  DSTONOW,DSTOPRE,DHGTPRE,DWTHINC,DWTHPRE,DWTHNOW,DHGTNOW,DHGTDIF
!$OMP THREADPRIVATE    (I,ILEV,DSTONOW,DSTOPRE,DHGTPRE,DWTHINC,DWTHPRE,DWTHNOW,DHGTNOW,DHGTDIF)
!####################################################################
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!---------------------!"
WRITE(LOGNAM,*) "CMF::LEVEE_INIT: initialize levee"

!********************
! [1] Read Levee Parameter Map
WRITE(LOGNAM,*) "CMF::LEVEE_INIT: read levee parameter files"

ALLOCATE( D2LEVHGT(NSEQALL,1) )
ALLOCATE( D2LEVFRC(NSEQALL,1) )
D2LEVHGT(:,:)   =0._JPRB
D2LEVFRC(:,:)   =0._JPRB

TMPNAM=INQUIRE_FID()

WRITE(LOGNAM,*)'INIT_LEVEE: levee crown height : ',TRIM(CLEVHGT)
OPEN(TMPNAM,FILE=CLEVHGT,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
READ(TMPNAM,REC=1) R2TEMP(:,:)
CALL mapR2vecD(R2TEMP,D2LEVHGT)
CLOSE(TMPNAM)

WRITE(LOGNAM,*)'INIT_LEVEE: distance from levee to river : ',TRIM(CLEVFRC)
OPEN(TMPNAM,FILE=CLEVFRC,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
READ(TMPNAM,REC=1) R2TEMP(:,:)
CALL mapR2vecD(R2TEMP,D2LEVFRC)
CLOSE(TMPNAM)

!*******************************
! [2] Calculate Levee Stage Parameter
WRITE(LOGNAM,*) "CMF::LEVEE_INIT: flood stage parameters considering levee"

ALLOCATE( D2BASHGT(NSEQALL,1) )
ALLOCATE( D2LEVDST(NSEQALL,1) )

ALLOCATE( D2LEVBASSTO(NSEQALL,1) )
ALLOCATE( D2LEVTOPSTO(NSEQALL,1) )
ALLOCATE( D2LEVFILSTO(NSEQALL,1) )

D2FLDSTOMAX(:,:,:) = 0._JPRD   !! max floodplain  storage  at each layer
D2FLDGRD(:,:,:)    = 0._JPRB   !! floodplain topo gradient of each layer
DFRCINC=REAL(NLFP,KIND=JPRB)**(-1._JPRB)   !! fration of each layer

D2LEVBASSTO(:,:)= 0._JPRB      !! storage at levee base     (levee protection start)
D2LEVTOPSTO(:,:)= 0._JPRB      !! storage at levee top      (levee protection end)
D2LEVFILSTO(:,:)= 0._JPRB      !! storage when levee filled (protected-side depth reach levee top)

!$OMP PARALLEL DO SIMD
DO ISEQ=1, NSEQALL
  IF( D2LEVHGT(ISEQ,1)<=0._JPRB )THEN
    D2LEVHGT(ISEQ,1)=0._JPRB
    D2LEVFRC(ISEQ,1)=1._JPRB   !! If no levee, all area is unprotected/
  ENDIF
  D2LEVFRC(ISEQ,1)=MAX(0._JPRB,MIN(1._JPRB,D2LEVFRC(ISEQ,1)))
END DO
!$OMP END PARALLEL DO SIMD

!$OMP PARALLEL DO
DO ISEQ=1, NSEQALL
! calculate floodplain parameters (without levee, same as SET_FLDSTG)
  DSTOPRE = D2RIVSTOMAX(ISEQ,1)
  DHGTPRE = 0._JPRB
  DWTHINC = D2GRAREA(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.) * DFRCINC  !! width increlment for each layer
  DO I=1, NLFP
    DSTONOW = D2RIVLEN(ISEQ,1) * ( D2RIVWTH(ISEQ,1) + DWTHINC*(REAL(I,KIND=JPRB)-0.5) ) * (D2FLDHGT(ISEQ,1,I)-DHGTPRE)  
    D2FLDSTOMAX(ISEQ,1,I) = DSTOPRE + DSTONOW
    D2FLDGRD(ISEQ,1,I) = (D2FLDHGT(ISEQ,1,I)-DHGTPRE) * DWTHINC**(-1.)
    DSTOPRE = D2FLDSTOMAX(ISEQ,1,I)
    DHGTPRE = D2FLDHGT(ISEQ,1,I)
  END DO

! Levee parameters calculation
  IF( D2LEVHGT(ISEQ,1) == 0._JPRB )THEN ! Grid without levee, treat everything as unprotected
    D2BASHGT(ISEQ,1) = 1.E18
    D2LEVDST(ISEQ,1) = 1.E18
    D2LEVBASSTO(ISEQ,1) = 1.E18
    D2LEVTOPSTO(ISEQ,1) = 1.E18
    D2LEVFILSTO(ISEQ,1) = 1.E18
  ELSE  !! levee exist
    !!*********
    !! [1] levee base storage & levee top storage (water only in river side)

    DSTOPRE = D2RIVSTOMAX(ISEQ,1)
    DHGTPRE = 0._JPRB
    DWTHPRE = 0._JPRB
    D2LEVDST(ISEQ,1) = D2LEVFRC(ISEQ,1) * DWTHINC*NLFP !! distance from channel to levee [m]

    ILEV=INT( D2LEVFRC(ISEQ,1)*NLFP )+1 !! which layer levee exist
    IF( ILEV>=2 )THEN
      DSTOPRE = D2FLDSTOMAX(ISEQ,1,ILEV-1)
      DHGTPRE = D2FLDHGT(ISEQ,1,ILEV-1)
      DWTHPRE = DWTHINC * (ILEV-1)
    ENDIF

    IF( ILEV<=NLFP )THEN
      !! levee in floodplain layer ILEV
      DWTHNOW = D2LEVDST(ISEQ,1) - DWTHPRE
      DHGTNOW = DWTHNOW * D2FLDGRD(ISEQ,1,ILEV) !! levee height above lower floodplain profile point
      D2BASHGT(ISEQ,1) = DHGTNOW + DHGTPRE
      D2LEVHGT(ISEQ,1) = max( D2LEVHGT(ISEQ,1), D2BASHGT(ISEQ,1) ) !! levee height >= base height

      DSTONOW = ( DWTHNOW*0.5 + DWTHPRE + D2RIVWTH(ISEQ,1) ) * DHGTNOW * D2RIVLEN(ISEQ,1) 
      D2LEVBASSTO(ISEQ,1) = DSTOPRE + DSTONOW

      DHGTDIF = D2LEVHGT(ISEQ,1) - D2BASHGT(ISEQ,1)
      D2LEVTOPSTO(ISEQ,1) = D2LEVBASSTO(ISEQ,1) + ( D2LEVDST(ISEQ,1)+D2RIVWTH(ISEQ,1) ) * DHGTDIF * D2RIVLEN(ISEQ,1)
    ELSE
      !! levee on the floodplain edge (ILEV=NLEV+1)
      D2BASHGT(ISEQ,1) = DHGTPRE
      D2LEVHGT(ISEQ,1) = max( D2LEVHGT(ISEQ,1), D2BASHGT(ISEQ,1) ) !! levee height >= base height

      D2LEVBASSTO(ISEQ,1) = DSTOPRE

      DHGTDIF = D2LEVHGT(ISEQ,1) - D2BASHGT(ISEQ,1)
      D2LEVTOPSTO(ISEQ,1) = D2LEVBASSTO(ISEQ,1) + ( D2LEVDST(ISEQ,1)+D2RIVWTH(ISEQ,1) ) * DHGTDIF * D2RIVLEN(ISEQ,1)
    ENDIF

    !!*********
    !! [2] levee fill storage (water in both river side & protected side)
    I=1
    DSTOPRE = D2RIVSTOMAX(ISEQ,1)
    DWTHPRE = D2RIVWTH(ISEQ,1)
    DHGTPRE = 0._JPRB

    !! check which layer levee top belongs
    DO WHILE( D2LEVHGT(ISEQ,1) > D2FLDHGT(ISEQ,1,I) .AND. I<=NLFP )
      DSTOPRE = D2FLDSTOMAX(ISEQ,1,I)
      DWTHPRE = DWTHPRE + DWTHINC
      DHGTPRE = D2FLDHGT(ISEQ,1,I)
      I=I+1
      IF( I>NLFP ) EXIT
    END DO

    !! calculate levee fill volume
    IF( I<=NLFP )THEN 
      !! levee top height collesponds to layer I
      DHGTNOW = D2LEVHGT(ISEQ,1) - DHGTPRE
      DWTHNOW = DHGTNOW * D2FLDGRD(ISEQ,1,I)**(-1.)

      DSTONOW = ( DWTHNOW*0.5 + DWTHPRE ) * DHGTNOW * D2RIVLEN(ISEQ,1)
      D2LEVFILSTO(ISEQ,1) = DSTOPRE + DSTONOW
    ELSE
      !! levee higher than catchment boundary height
      DHGTNOW = D2LEVHGT(ISEQ,1) - DHGTPRE
      DSTONOW = DWTHPRE * DHGTNOW * D2RIVLEN(ISEQ,1)
      D2LEVFILSTO(ISEQ,1) = DSTOPRE + DSTONOW
    ENDIF
  ENDIF

END DO
!$OMP END PARALLEL DO

END SUBROUTINE CMF_LEVEE_INIT
!####################################################################




!####################################################################
SUBROUTINE CMF_LEVEE_FLDSTG
! ================================================
! calculate river and floodplain staging considering levee
! ================================================
USE YOS_CMF_INPUT  ,ONLY: NLFP
USE YOS_CMF_MAP    ,ONLY: NSEQALL
USE YOS_CMF_MAP    ,ONLY: D2GRAREA, D2RIVLEN, D2RIVWTH, D2RIVELV, D2RIVSTOMAX, D2FLDSTOMAX, D2FLDGRD, DFRCINC, D2FLDHGT
USE YOS_CMF_PROG   ,ONLY: P2RIVSTO, P2FLDSTO
USE YOS_CMF_DIAG   ,ONLY: D2RIVDPH, D2FLDDPH, D2FLDFRC, D2FLDARE, D2SFCELV, D2STORGE
USE YOS_CMF_DIAG   ,ONLY: P0GLBSTOPRE2, P0GLBSTONEW2, P0GLBRIVSTO, P0GLBFLDSTO, P0GLBLEVSTO, P0GLBFLDARE

!! levee specific data
USE YOS_CMF_PROG   ,ONLY: P2LEVSTO  !! flood storage in protected side (P2FLDSTO for storage betwen river & levee)
USE YOS_CMF_DIAG   ,ONLY: D2LEVDPH  !! flood depth in protected side   (D2FLDDPH for water depth betwen river & levee)
IMPLICIT NONE

!*** LOCAL
! Save for OpenMP
INTEGER(KIND=JPIM),SAVE ::  ISEQ, I, ILEV
REAL(KIND=JPRD),SAVE    ::  PSTOALL
REAL(KIND=JPRB),SAVE    ::  DSTOALL, DSTONOW, DSTOPRE, DWTHNOW, DWTHPRE, DDPHPRE, DDPHNOW, DWTHINC, DSTOADD
!$OMP THREADPRIVATE (I,ILEV,DSTOALL, DSTONOW, DSTOPRE, DWTHNOW, DWTHPRE, DDPHPRE, DDPHNOW, DWTHINC, DSTOADD, PSTOALL)
!!==============================
P0GLBSTOPRE2=0._JPRD
P0GLBSTONEW2=0._JPRD
P0GLBRIVSTO =0._JPRD
P0GLBFLDSTO =0._JPRD
P0GLBLEVSTO =0._JPRD
P0GLBFLDARE =0._JPRD

!$OMP PARALLEL DO REDUCTION(+:P0GLBSTOPRE2,P0GLBSTONEW2,P0GLBRIVSTO,P0GLBFLDSTO,P0GLBLEVSTO,P0GLBFLDARE)
DO ISEQ=1, NSEQALL
!
  PSTOALL = P2RIVSTO(ISEQ,1) + P2FLDSTO(ISEQ,1) + P2LEVSTO(ISEQ,1)
  DSTOALL = REAL( PSTOALL, KIND=JPRB)
  DWTHINC = D2GRAREA(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.) * DFRCINC    !! width of each layer [m]
  IF( DSTOALL > D2RIVSTOMAX(ISEQ,1) )THEN
    !**********
    ! [Case-1] Water surface is under levee base (all water is between river-levee)
    IF( DSTOALL < D2LEVBASSTO(ISEQ,1) )THEN 
      I=1
      DSTOPRE = D2RIVSTOMAX(ISEQ,1)
      DWTHPRE = D2RIVWTH(ISEQ,1)
      DDPHPRE = 0._JPRB

      ! which layer current water level is
      DO WHILE( DSTOALL > D2FLDSTOMAX(ISEQ,1,I) .AND. I<=NLFP )
        DSTOPRE = D2FLDSTOMAX(ISEQ,1,I)
        DWTHPRE = DWTHPRE + DWTHINC
        DDPHPRE = DDPHPRE + D2FLDGRD(ISEQ,1,I) * DWTHINC
        I=I+1
        IF( I>NLFP ) EXIT
      END DO

      ! water depth at unprotected area
      IF( I<=NLFP )THEN
        DSTONOW =  DSTOALL - DSTOPRE
        DWTHNOW = -DWTHPRE + ( DWTHPRE**2. + 2. * DSTONOW * D2RIVLEN(ISEQ,1)**(-1.) * D2FLDGRD(ISEQ,1,I)**(-1.) )**0.5
        D2FLDDPH(ISEQ,1) = DDPHPRE + D2FLDGRD(ISEQ,1,I) * DWTHNOW
      ELSE
        DSTONOW = DSTOALL - DSTOPRE
        DWTHNOW = 0._JPRB
        D2FLDDPH(ISEQ,1) = DDPHPRE + DSTONOW * DWTHPRE**(-1.) * D2RIVLEN(ISEQ,1)**(-1.)
      ENDIF

      P2RIVSTO(ISEQ,1) = D2RIVSTOMAX(ISEQ,1) + D2RIVLEN(ISEQ,1) * D2RIVWTH(ISEQ,1) * D2FLDDPH(ISEQ,1)
      D2RIVDPH(ISEQ,1) = REAL(P2RIVSTO(ISEQ,1),KIND=JPRB) * D2RIVLEN(ISEQ,1)**(-1.) * D2RIVWTH(ISEQ,1)**(-1.)
!
      P2FLDSTO(ISEQ,1) = PSTOALL - P2RIVSTO(ISEQ,1)
      P2FLDSTO(ISEQ,1) = MAX( P2FLDSTO(ISEQ,1), 0._JPRD )
      D2FLDFRC(ISEQ,1) = (-D2RIVWTH(ISEQ,1) + DWTHPRE + DWTHNOW ) * (DWTHINC*NLFP)**(-1.)
      D2FLDFRC(ISEQ,1) = MAX( D2FLDFRC(ISEQ,1),0._JPRB )
      D2FLDFRC(ISEQ,1) = MIN( D2FLDFRC(ISEQ,1),1._JPRB )
      D2FLDARE(ISEQ,1) = D2GRAREA(ISEQ,1)*D2FLDFRC(ISEQ,1)
!
      P2LEVSTO(ISEQ,1) = 0._JPRD  !! no flooding in protected area
      D2LEVDPH(ISEQ,1) = 0._JPRB

    !**********
    ! [Case-2]  River-side water surface is under levee crown (water only in river side)
    ELSEIF( DSTOALL < D2LEVTOPSTO(ISEQ,1) )THEN 

      DSTONOW = DSTOALL - D2LEVBASSTO(ISEQ,1)
      DWTHNOW = D2LEVDST(ISEQ,1) + D2RIVWTH(ISEQ,1)
      D2FLDDPH(ISEQ,1) = D2BASHGT(ISEQ,1) + DSTONOW * DWTHNOW**(-1.) * D2RIVLEN(ISEQ,1)**(-1.)

      P2RIVSTO(ISEQ,1) = D2RIVSTOMAX(ISEQ,1) + D2RIVLEN(ISEQ,1) * D2RIVWTH(ISEQ,1) * D2FLDDPH(ISEQ,1)
      D2RIVDPH(ISEQ,1) = REAL(P2RIVSTO(ISEQ,1),KIND=JPRB) * D2RIVLEN(ISEQ,1)**(-1.) * D2RIVWTH(ISEQ,1)**(-1.)
  !
      P2FLDSTO(ISEQ,1) = PSTOALL - P2RIVSTO(ISEQ,1)
      P2FLDSTO(ISEQ,1) = MAX( P2FLDSTO(ISEQ,1), 0._JPRD )
      D2FLDFRC(ISEQ,1) = D2LEVFRC(ISEQ,1)
      D2FLDARE(ISEQ,1) = D2GRAREA(ISEQ,1)*D2FLDFRC(ISEQ,1)
  ! 
      P2LEVSTO(ISEQ,1) = 0._JPRD  !! no flooding in protected area
      D2LEVDPH(ISEQ,1) = 0._JPRB

    !**********
    ! [Case-3] River side is full, protected side is under levee crown height (Water both in river side & protected side)
    ELSEIF( DSTOALL < D2LEVFILSTO(ISEQ,1) )THEN 
      ! river side stage = levee height
      D2FLDDPH(ISEQ,1) = D2LEVHGT(ISEQ,1)
      P2RIVSTO(ISEQ,1) = D2RIVSTOMAX(ISEQ,1) + D2RIVLEN(ISEQ,1) * D2RIVWTH(ISEQ,1) * D2FLDDPH(ISEQ,1)
      D2RIVDPH(ISEQ,1) = REAL(P2RIVSTO(ISEQ,1),KIND=JPRB) * D2RIVLEN(ISEQ,1)**(-1.) * D2RIVWTH(ISEQ,1)**(-1.)

      P2FLDSTO(ISEQ,1) = D2LEVTOPSTO(ISEQ,1) - P2RIVSTO(ISEQ,1)
      P2FLDSTO(ISEQ,1) = MAX( P2FLDSTO(ISEQ,1), 0._JPRD )

      !! protected side storate calculation
      P2LEVSTO(ISEQ,1) = PSTOALL - P2RIVSTO(ISEQ,1) - P2FLDSTO(ISEQ,1)
      P2LEVSTO(ISEQ,1) = MAX( P2LEVSTO(ISEQ,1), 0._JPRD )

      !!****
      !! protected side stage calculation
      ILEV=INT( D2LEVFRC(ISEQ,1)*NLFP )+1 !! levee relative distance -> floodplain layer with levee
      DSTOPRE = D2LEVTOPSTO(ISEQ,1)
      DWTHPRE = 0._JPRB
      DDPHPRE = 0._JPRB
      !! which layer current water level is
      I=ILEV
      DO WHILE( I<=NLFP )
        DSTOADD = ( D2LEVDST(ISEQ,1)+D2RIVWTH(ISEQ,1) ) * ( D2LEVHGT(ISEQ,1)-D2FLDHGT(ISEQ,1,I) ) * D2RIVLEN(ISEQ,1) 
        IF( DSTOALL < D2FLDSTOMAX(ISEQ,1,I) + DSTOADD ) EXIT
        DSTOPRE = D2FLDSTOMAX(ISEQ,1,I) + DSTOADD
        DWTHPRE = DWTHINC*I - D2LEVDST(ISEQ,1)
        DDPHPRE = D2FLDHGT(ISEQ,1,I) - D2BASHGT(ISEQ,1)
        I=I+1
        IF( I>NLFP ) EXIT
      END DO

      IF( I<=NLFP )THEN
        DSTONOW = DSTOALL - DSTOPRE
        DWTHNOW = -DWTHPRE + ( DWTHPRE**2. + 2. * DSTONOW*D2RIVLEN(ISEQ,1)**(-1.) * D2FLDGRD(ISEQ,1,I)**(-1.) )**0.5
        DDPHNOW = DWTHNOW * D2FLDGRD(ISEQ,1,I)
        D2LEVDPH(ISEQ,1) = D2BASHGT(ISEQ,1) + DDPHPRE + DDPHNOW

        D2FLDFRC(ISEQ,1) = ( DWTHPRE + D2LEVDST(ISEQ,1) ) * (DWTHINC*NLFP)**(-1.)
        D2FLDFRC(ISEQ,1) = MAX( D2FLDFRC(ISEQ,1),0._JPRB)
        D2FLDFRC(ISEQ,1) = MIN( D2FLDFRC(ISEQ,1),1._JPRB)
        D2FLDARE(ISEQ,1) = D2GRAREA(ISEQ,1)*D2FLDFRC(ISEQ,1)
      ELSE
        DSTONOW = DSTOALL - DSTOPRE
        DDPHNOW = DSTONOW * DWTHPRE**(-1.) * D2RIVLEN(ISEQ,1)**(-1.)
        D2LEVDPH(ISEQ,1) = D2BASHGT(ISEQ,1) + DDPHPRE + DDPHNOW

        D2FLDFRC(ISEQ,1) = 1._JPRB
        D2FLDARE(ISEQ,1) = D2GRAREA(ISEQ,1)*D2FLDFRC(ISEQ,1)
      ENDIF

    !**********
    ! [Case-4] Water level above levee crown (Both river side and protected side exceed levee crown height)
    ELSE 
      I=1
      DSTOPRE = D2RIVSTOMAX(ISEQ,1)
      DWTHPRE = D2RIVWTH(ISEQ,1)
      DDPHPRE = 0._JPRB
      DO WHILE( DSTOALL > D2FLDSTOMAX(ISEQ,1,I) .AND. I<=NLFP)
        DSTOPRE = D2FLDSTOMAX(ISEQ,1,I)
        DWTHPRE = DWTHPRE + DWTHINC
        DDPHPRE = DDPHPRE + D2FLDGRD(ISEQ,1,I) * DWTHINC
        I=I+1
        IF( I>NLFP ) EXIT
      END DO

      IF( I<=NLFP )THEN
        DSTONOW =  DSTOALL - DSTOPRE
        DWTHNOW = -DWTHPRE + ( DWTHPRE**2. + 2. * DSTONOW * D2RIVLEN(ISEQ,1)**(-1.) * D2FLDGRD(ISEQ,1,I)**(-1.) )**0.5
        D2FLDDPH(ISEQ,1) = DDPHPRE + D2FLDGRD(ISEQ,1,I) * DWTHNOW
      ELSE
        DSTONOW = DSTOALL - DSTOPRE
        DWTHNOW = 0._JPRB
        D2FLDDPH(ISEQ,1) = DDPHPRE + DSTONOW * DWTHPRE**(-1.) * D2RIVLEN(ISEQ,1)**(-1.)
      ENDIF

      D2FLDFRC(ISEQ,1) = (-D2RIVWTH(ISEQ,1) + DWTHPRE + DWTHNOW ) * (DWTHINC*NLFP)**(-1.)
      D2FLDARE(ISEQ,1) = D2GRAREA(ISEQ,1)*D2FLDFRC(ISEQ,1)

      !! river channel storage
      P2RIVSTO(ISEQ,1) = D2RIVSTOMAX(ISEQ,1) + D2RIVLEN(ISEQ,1) * D2RIVWTH(ISEQ,1) * D2FLDDPH(ISEQ,1)
      D2RIVDPH(ISEQ,1) = REAL(P2RIVSTO(ISEQ,1),KIND=JPRB) * D2RIVLEN(ISEQ,1)**(-1.) * D2RIVWTH(ISEQ,1)**(-1.)
!
      DSTOADD = ( D2FLDDPH(ISEQ,1)-D2LEVHGT(ISEQ,1) ) * (D2LEVDST(ISEQ,1)+D2RIVWTH(ISEQ,1)) * D2RIVLEN(ISEQ,1)
      P2FLDSTO(ISEQ,1) = D2LEVTOPSTO(ISEQ,1) + DSTOADD - P2RIVSTO(ISEQ,1)
      P2FLDSTO(ISEQ,1) = MAX( P2FLDSTO(ISEQ,1), 0._JPRD )

      P2LEVSTO(ISEQ,1) = PSTOALL - P2RIVSTO(ISEQ,1) - P2FLDSTO(ISEQ,1)
      P2LEVSTO(ISEQ,1) = MAX( P2LEVSTO(ISEQ,1), 0._JPRD )
      D2LEVDPH(ISEQ,1) = D2FLDDPH(ISEQ,1)
    ENDIF

  ! [Case-0] Water only in river channel
  ELSE
    P2RIVSTO(ISEQ,1) = PSTOALL
    D2RIVDPH(ISEQ,1) = DSTOALL * D2RIVLEN(ISEQ,1)**(-1.) * D2RIVWTH(ISEQ,1)**(-1.)
    D2RIVDPH(ISEQ,1) = MAX( D2RIVDPH(ISEQ,1), 0._JPRB )
    P2FLDSTO(ISEQ,1) = 0._JPRD
    D2FLDDPH(ISEQ,1) = 0._JPRB
    D2FLDFRC(ISEQ,1) = 0._JPRB
    D2FLDARE(ISEQ,1) = 0._JPRB
    P2LEVSTO(ISEQ,1) = 0._JPRD
    D2LEVDPH(ISEQ,1) = 0._JPRB
  ENDIF
  D2SFCELV(ISEQ,1)     = D2RIVELV(ISEQ,1) + D2RIVDPH(ISEQ,1)
  D2STORGE(ISEQ,1) =REAL(P2RIVSTO(ISEQ,1) + P2FLDSTO(ISEQ,1)+ P2LEVSTO(ISEQ,1) ,KIND=JPRB)

  P0GLBSTOPRE2     = P0GLBSTOPRE2 + PSTOALL
  P0GLBSTONEW2     = P0GLBSTONEW2 + P2RIVSTO(ISEQ,1) + P2FLDSTO(ISEQ,1) + P2LEVSTO(ISEQ,1)
  P0GLBRIVSTO      = P0GLBRIVSTO  + P2RIVSTO(ISEQ,1)
  P0GLBFLDSTO      = P0GLBFLDSTO  + P2FLDSTO(ISEQ,1)
  P0GLBLEVSTO      = P0GLBLEVSTO  + P2LEVSTO(ISEQ,1)
  P0GLBFLDARE      = P0GLBFLDARE  + D2FLDARE(ISEQ,1)
END DO
!$OMP END PARALLEL DO

END SUBROUTINE CMF_LEVEE_FLDSTG
!####################################################################




!####################################################################
SUBROUTINE CMF_LEVEE_OPT_PTHOUT
! realistic bifurcation considering levee
USE PARKIND1,           ONLY: JPIM, JPRB, JPRD
USE YOS_CMF_INPUT,      ONLY: DT, PGRV, DMIS
USE YOS_CMF_MAP,        ONLY: NSEQALL, NPTHOUT, NPTHLEV, PTH_UPST, PTH_DOWN, PTH_DST, &
                            & PTH_ELV, PTH_WTH, PTH_MAN, I2MASK
USE YOS_CMF_MAP,        ONLY: D2ELEVTN, D2RIVELV
USE YOS_CMF_PROG,       ONLY: D1PTHFLW, D1PTHFLW_PRE, D2RIVDPH_PRE
USE YOS_CMF_DIAG,       ONLY: D2LEVDPH, D2SFCELV
IMPLICIT NONE
!*** Local
REAL(KIND=JPRB)    ::  D2SFCELV_LEV(NSEQALL,1)                  !! water surface elev protected [m]
REAL(KIND=JPRB)    ::  D2SFCELV_PRE(NSEQALL,1)                  !! water surface elev (t-1) [m] (for stable calculation)

! SAVE for OpenMP
INTEGER(KIND=JPIM),SAVE ::  IPTH, ILEV, ISEQ, ISEQP, JSEQP
REAL(KIND=JPRB),SAVE    ::  DSLOPE, DFLW, DOUT_PRE, DFLW_PRE, DFLW_IMP
!$OMP THREADPRIVATE        (DSLOPE, DFLW, DOUT_PRE, DFLW_PRE, DFLW_IMP, ILEV, ISEQP, JSEQP)
!================================================
!$OMP PARALLEL DO SIMD
DO ISEQ=1, NSEQALL
  IF( D2LEVFRC(ISEQ,1)<1.0 )THEN
    D2SFCELV_LEV(ISEQ,1) = D2ELEVTN(ISEQ,1)+D2LEVDPH(ISEQ,1)  !! levee exist, calculate pthout based on levee protected depth
  ELSE
    D2SFCELV_LEV(ISEQ,1) = D2SFCELV(ISEQ,1)
  ENDIF

  D2SFCELV_PRE(ISEQ,1) = D2RIVELV(ISEQ,1)+D2RIVDPH_PRE(ISEQ,1)
END DO
!$OMP END PARALLEL DO SIMD

D1PTHFLW(:,:) = DMIS
!$OMP PARALLEL DO
DO IPTH=1, NPTHOUT  
  ISEQP=PTH_UPST(IPTH)
  JSEQP=PTH_DOWN(IPTH)
  !! Avoid calculation outside of domain
  IF (ISEQP<=0 .OR. JSEQP<=0 ) CYCLE
  IF (I2MASK(ISEQP,1) == 1 .OR. I2MASK(JSEQP,1) == 1 ) CYCLE  !! I2MASK is for kinematic-inertial mixed flow scheme. 

!! [1] for channel bifurcation, use river surface elevation  
  DSLOPE  = (D2SFCELV(ISEQP,1)-D2SFCELV(JSEQP,1)) * PTH_DST(IPTH)**(-1.)
  DSLOPE = max(-0.005_JPRB,min(0.005_JPRB,DSLOPE))                                    !! v390 stabilization

  ILEV=1 !! for river channek
    DFLW = MAX(D2SFCELV(ISEQP,1),D2SFCELV(JSEQP,1)) - PTH_ELV(IPTH,ILEV) 
    DFLW = MAX(DFLW,0._JPRB)

    DFLW_PRE = MAX(D2SFCELV_PRE(ISEQP,1),D2SFCELV_PRE(JSEQP,1)) - PTH_ELV(IPTH,ILEV)
    DFLW_PRE = MAX(DFLW_PRE,0._JPRB)

    DFLW_IMP = (DFLW*DFLW_PRE)**0.5                                       !! semi implicit flow depth
    IF( DFLW_IMP<=0._JPRB ) DFLW_IMP=DFLW

    IF( DFLW_IMP>1.E-5 )THEN                         !! local inertial equation, see [Bates et al., 2010, J.Hydrol.]
      DOUT_PRE = D1PTHFLW_PRE(IPTH,ILEV) * PTH_WTH(IPTH,ILEV)**(-1.)           !! outflow (t-1) [m2/s] (unit width)
      D1PTHFLW(IPTH,ILEV) = PTH_WTH(IPTH,ILEV) * ( DOUT_PRE + PGRV*DT*DFLW_IMP*DSLOPE ) &
                         * ( 1. + PGRV*DT*PTH_MAN(ILEV)**2. * abs(DOUT_PRE)*DFLW_IMP**(-7./3.) )**(-1.)
    ELSE
      D1PTHFLW(IPTH,ILEV) = 0._JPRB
    ENDIF

!! [1] for overland bifurcation, use levee protected surface elevation
  IF( NPTHLEV<=1 ) CYCLE

  DSLOPE  = (D2SFCELV_LEV(ISEQP,1)-D2SFCELV_LEV(JSEQP,1)) * PTH_DST(IPTH)**(-1.)
  DSLOPE = max(-0.005_JPRB,min(0.005_JPRB,DSLOPE))      

  DO ILEV=2, NPTHLEV
    DFLW = MAX(D2SFCELV_LEV(ISEQP,1),D2SFCELV_LEV(JSEQP,1)) - PTH_ELV(IPTH,ILEV) 
    DFLW = MAX(DFLW,0._JPRB)

    DFLW_IMP=DFLW  !! do not consider implicit flow depth for overland bifurcation
    IF( DFLW_IMP>1.E-5 )THEN                         !! local inertial equation, see [Bates et al., 2010, J.Hydrol.]
      DOUT_PRE = D1PTHFLW_PRE(IPTH,ILEV) * PTH_WTH(IPTH,ILEV)**(-1.)            !! outflow (t-1) [m2/s] (unit width)
      D1PTHFLW(IPTH,ILEV) = PTH_WTH(IPTH,ILEV) * ( DOUT_PRE + PGRV*DT*DFLW_IMP*DSLOPE ) &
                         * ( 1. + PGRV*DT*PTH_MAN(ILEV)**2. * abs(DOUT_PRE)*DFLW_IMP**(-7./3.) )**(-1.)
    ELSE
      D1PTHFLW(IPTH,ILEV) = 0._JPRB
    ENDIF
  END DO
END DO
!$OMP END PARALLEL DO

END SUBROUTINE CMF_LEVEE_OPT_PTHOUT
!################################################################

END MODULE CMF_CTRL_LEVEE_MOD
