MODULE CMF_CTRL_MAPS_MOD
!==========================================================
!* PURPOSE: Control CaMa-Flood map/topography data
!
!* CONTAINS:
! -- CMF_MAPS_NMLIST   : configuration from namelist
! -- CMF_RIVMAP_INIT  : read & set river network map 
! -- CMF_TOPO_INIT    : read & set topography
!
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  Aug 2019
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
! shared variables in module
USE PARKIND1,                ONLY: JPIM, JPRB, JPRD, JPRM
USE YOS_CMF_INPUT,           ONLY: LOGNAM
IMPLICIT NONE
SAVE
!*** NAMELIST/NMAP/ from inputnam
CHARACTER(LEN=256)              :: CNEXTXY         !! river network nextxy
CHARACTER(LEN=256)              :: CGRAREA         !! catchment area
CHARACTER(LEN=256)              :: CELEVTN         !! bank top elevation
CHARACTER(LEN=256)              :: CNXTDST         !! distance to next outlet
CHARACTER(LEN=256)              :: CRIVLEN         !! river channel length
CHARACTER(LEN=256)              :: CFLDHGT         !! floodplain elevation profile
!* river channel parameters
CHARACTER(LEN=256)              :: CRIVWTH         !! channel width
CHARACTER(LEN=256)              :: CRIVHGT         !! channel depth
CHARACTER(LEN=256)              :: CRIVMAN         !! river manning coefficient
!* optional maps
CHARACTER(LEN=256)              :: CPTHOUT         !! bifurcation channel table
CHARACTER(LEN=256)              :: CGDWDLY         !! Groundwater Delay Parameter
CHARACTER(LEN=256)              :: CMEANSL         !! mean sea level
!* MPI parallelization
CHARACTER(LEN=256)              :: CMPIREG         !! MPI region map
!* netCDF map
LOGICAL                         :: LMAPCDF         !! true for netCDF map input
CHARACTER(LEN=256)              :: CRIVCLINC       !! river map netcdf
CHARACTER(LEN=256)              :: CRIVPARNC       !! river parameter netcdf (WIDTH,HEIGHT, Manning, ground wateer delay)
CHARACTER(LEN=256)              :: CMEANSLNC       !! mean sea level netCDF
CHARACTER(LEN=256)              :: CMPIREGNC       !! MPI region map in netcdf

NAMELIST/NMAP/     CNEXTXY,  CGRAREA,  CELEVTN,  CNXTDST, CRIVLEN, CFLDHGT, &
                   CRIVWTH,  CRIVHGT,  CRIVMAN,  CPTHOUT, CGDWDLY, CMEANSL, &
                   CMPIREG,  LMAPCDF,  CRIVCLINC,CRIVPARNC,CMEANSLNC,CMPIREGNC


CONTAINS
!####################################################################
! -- CMF_MAP_NMLIST   : configuration from namelist
! -- CMF_RIVMAP_INIT  : read & set river network map 
! -- CMF_TOPO_INIT    : read & set topography
!
!
!####################################################################
SUBROUTINE CMF_MAPS_NMLIST
! reed setting from namelist
! -- Called from CMF_DRV_NMLIST
USE YOS_CMF_INPUT,      ONLY: CSETFILE,NSETFILE,LMEANSL,LGDWDLY
USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID
IMPLICIT NONE
!================================================
!*** 1. open namelist
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!---------------------!"

NSETFILE=INQUIRE_FID()
OPEN(NSETFILE,FILE=CSETFILE,STATUS="OLD")
WRITE(LOGNAM,*) "CMF::MAP_NMLIST: namelist OPEN in unit: ", TRIM(CSETFILE), NSETFILE 

!*** 2. default value
CNEXTXY="./nextxy.bin"
CGRAREA="./ctmare.bin"
CELEVTN="./elevtn.bin"
CNXTDST="./nxtdst.bin"
CRIVLEN="./rivlen.bin"
CFLDHGT="./fldhgt.bin"

CRIVWTH="./rivwth.bin"
CRIVHGT="./rivhgt.bin"
CRIVMAN="./rivman.bin"

CPTHOUT="./bifprm.txt"
CGDWDLY="NONE"
CMEANSL="NONE"

CMPIREG="NONE"

LMAPCDF=.FALSE.
CRIVCLINC="NONE"
CRIVPARNC="NONE"
CMEANSLNC="NONE"
CMPIREGNC="NONE"

!*** 3. read namelist
REWIND(NSETFILE)
READ(NSETFILE,NML=NMAP)

WRITE(LOGNAM,*)     "=== NAMELIST, NMAP ==="
WRITE(LOGNAM,*)     "LMAPCDF:   ", LMAPCDF
IF( LMAPCDF )THEN
  WRITE(LOGNAM,*)   "CRIVCLINC: ", TRIM(CRIVCLINC)
  WRITE(LOGNAM,*)   "CRIVPARNC: ", TRIM(CRIVPARNC)
  IF( LMEANSL ) THEN
    WRITE(LOGNAM,*) "CMEANSLNC: ", TRIM(CMEANSLNC)
  ENDIF
#ifdef UseMPI_CMF
    WRITE(LOGNAM,*) "CMPIREGNC:   ", TRIM(CMPIREGNC)
#endif
ELSE
  WRITE(LOGNAM,*)   "CNEXTXY:   ", TRIM(CNEXTXY)
  WRITE(LOGNAM,*)   "CGRAREA:   ", TRIM(CGRAREA)
  WRITE(LOGNAM,*)   "CELEVTN:   ", TRIM(CELEVTN)
  WRITE(LOGNAM,*)   "CNXTDST:   ", TRIM(CNXTDST)
  WRITE(LOGNAM,*)   "CRIVLEN:   ", TRIM(CRIVLEN)
  WRITE(LOGNAM,*)   "CFLDHGT:   ", TRIM(CFLDHGT)

  WRITE(LOGNAM,*)   "CRIVWTH:   ", TRIM(CRIVWTH)
  WRITE(LOGNAM,*)   "CRIVHGT:   ", TRIM(CRIVHGT)
  WRITE(LOGNAM,*)   "CRIVMAN:   ", TRIM(CRIVMAN)

  WRITE(LOGNAM,*)   "CPTHOUT:   ", TRIM(CPTHOUT)
  IF( LGDWDLY )THEN
    WRITE(LOGNAM,*) "CGDWDLY:    ",TRIM(CGDWDLY)
  ENDIF
  IF( LMEANSL )THEN
    WRITE(LOGNAM,*) "CMEANSL:   ", TRIM(CMEANSL)
  ENDIF
#ifdef UseMPI_CMF
    WRITE(LOGNAM,*) "CMPIREG:   ", TRIM(CMPIREG)
#endif
ENDIF

CLOSE(NSETFILE)

WRITE(LOGNAM,*) "CMF::MAP_NMLIST: end"

END SUBROUTINE CMF_MAPS_NMLIST
!####################################################################





!####################################################################
SUBROUTINE CMF_RIVMAP_INIT
! read & set river network map 
! -- call from CMF_DRV_INIT
USE YOS_CMF_INPUT,      ONLY: TMPNAM, NX,NY,NLFP, LPTHOUT
USE YOS_CMF_MAP,        ONLY: I2NEXTX,I2NEXTY, I2REGION, REGIONALL,REGIONTHIS, &
                            & I1SEQX, I1SEQY,  I1NEXT,  I2VECTOR, D1LON,    D1LAT,      &
                            & NSEQRIV,NSEQALL, NSEQMAX, NPTHOUT,  NPTHLEV
USE YOS_CMF_MAP,        ONLY: I1UPST, I1UPN, I1P_OUT, I1P_OUTN, I1P_INF, I1P_INFN
USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID
IMPLICIT NONE
!================================================
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!---------------------!"
WRITE(LOGNAM,*) 'CMF::RIVMAP_INIT: river network initialization'

! *** 1. ALLOCATE ARRAYS
ALLOCATE( I2NEXTX(NX,NY) )
ALLOCATE( I2NEXTY(NX,NY) )
ALLOCATE( I2REGION(NX,NY) )
ALLOCATE( D1LON(NX) )
ALLOCATE( D1LAT(NY) )

!============================
!*** 2a. read river network map
WRITE(LOGNAM,*) 'CMF::RIVMAP_INIT: read nextXY & set lat lon'
IF( LMAPCDF )THEN
#ifdef UseCDF_CMF
  CALL READ_MAP_CDF
#endif
ELSE
  CALL READ_MAP_BIN
ENDIF

!*** 2b. calculate river sequence & regions
WRITE(LOGNAM,*) 'CMF::RIVMAP_INIT: calc region'
CALL CALC_REGION

!============================
!*** 3. conversion 2D map -> 1D vector
WRITE(LOGNAM,*) 'CMF::RIVMAP_INIT: calculate 1d river sequence'

CALL CALC_1D_SEQ                                  !! 2D map to 1D vector conversion. for faster calculation

WRITE(LOGNAM,*) '  NSEQRIV=',NSEQRIV
WRITE(LOGNAM,*) '  NSEQALL=',NSEQALL
WRITE(LOGNAM,*) '  NSEQMAX=',NSEQMAX

!*** 3c. Write Map Data                                       !! used for combining mpi distributed output into one map
IF( REGIONTHIS==1 )THEN
  TMPNAM=INQUIRE_FID()
  OPEN(TMPNAM,FILE='./mapdata.txt',FORM='FORMATTED')
  WRITE(TMPNAM,*) 'NX',        NX
  WRITE(TMPNAM,*) 'NY',        NY
  WRITE(TMPNAM,*) 'NLFP',      NLFP
  WRITE(TMPNAM,*) 'REGIONALL', REGIONALL
  WRITE(TMPNAM,*) 'NSEQMAX',   NSEQMAX
  CLOSE(TMPNAM)
ENDIF

!============================
!*** 4.  bifurcation channel parameters
NPTHOUT=1  !! tentative number
NPTHLEV=1
IF( LPTHOUT )THEN
  WRITE(LOGNAM,*) 'CMF::RIVMAP_INIT: read bifurcation channel setting'
  CALL READ_BIFPARAM
ENDIF

DEALLOCATE( I2NEXTX,I2NEXTY,I2REGION )

WRITE(LOGNAM,*) 'CMF::RIVMAP_INIT: end'

CONTAINS
!==========================================================
!+ READ_MAP_BIN
!+ READ_MAP_CDF
!+ CALC_REGION
!+ CALC_1D_SEQ
!+ READ_BIFPRM
!==========================================================
SUBROUTINE READ_MAP_BIN
USE YOS_CMF_INPUT,      ONLY: TMPNAM, LMAPEND
USE YOS_CMF_INPUT,      ONLY: WEST,EAST,NORTH,SOUTH
USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID, CONV_ENDI
IMPLICIT NONE
!* local variables
INTEGER(KIND=JPIM),SAVE    :: IX,IY
!==========================================================
!*** read river map
WRITE(LOGNAM,*)'RIVMAP_INIT: nextxy binary: ',TRIM(CNEXTXY)
TMPNAM=INQUIRE_FID()
OPEN(TMPNAM,FILE=CNEXTXY,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
READ(TMPNAM,REC=1) I2NEXTX
READ(TMPNAM,REC=2) I2NEXTY
CLOSE(TMPNAM)

IF ( LMAPEND )THEN
  CALL CONV_ENDI(I2NEXTX,NX,NY)
  CALL CONV_ENDI(I2NEXTY,NX,NY)
ENDIF

!*** calculate lat, lon
IF( WEST>=-180._JPRB .and. EAST<=360._JPRB .and. SOUTH>=-180._JPRB .and. NORTH<=180._JPRB )THEN  !! bugfix_v396a
  !$OMP PARALLEL DO
  DO IX=1,NX
    D1LON(IX)=WEST +(REAL(IX,KIND=JPRB)-0.5_JPRB)*(EAST-WEST)  /REAL(NX,KIND=JPRB)
  ENDDO
  !$OMP END PARALLEL DO
  !$OMP PARALLEL DO
  DO IY=1,NY
    D1LAT(IY)=NORTH-(REAL(IY,KIND=JPRB)-0.5_JPRB)*(NORTH-SOUTH)/REAL(NY,KIND=JPRB)
  ENDDO
  !$OMP END PARALLEL DO
ENDIF

END SUBROUTINE READ_MAP_BIN
!==========================================================
!+
!+
!+
!==========================================================
#ifdef UseCDF_CMF
SUBROUTINE READ_MAP_CDF
USE CMF_UTILS_MOD  ,ONLY: NCERROR
USE NETCDF
IMPLICIT NONE
!* local variables
INTEGER(KIND=JPIM)              :: NCID,VARID
!================================================
WRITE(LOGNAM,*)'RIVMAP_INIT: nextxy netCDF: ', TRIM(CRIVCLINC)

CALL NCERROR (NF90_OPEN(CRIVCLINC,NF90_NOWRITE,NCID),'opening '//TRIM(CRIVCLINC) )

!*** next xy
CALL NCERROR ( NF90_INQ_VARID(NCID,'nextx',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,I2NEXTX),'reading data' ) 

CALL NCERROR ( NF90_INQ_VARID(NCID,'nexty',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,I2NEXTY),'reading data' )

!*** lat, lon
CALL NCERROR ( NF90_INQ_VARID(NCID,'lat',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,D1LAT),'reading data' )

CALL NCERROR ( NF90_INQ_VARID(NCID,'lon',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,D1LON),'reading data' )

CALL NCERROR( NF90_CLOSE(NCID))

END SUBROUTINE READ_MAP_CDF
#endif
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE CALC_REGION    !! evenly allocate pixels to mpi nodes (updated in v4.03. MPI region given from file)
USE YOS_CMF_INPUT,           ONLY: IMIS
#ifdef UseCDF_CMF
USE CMF_UTILS_MOD,           ONLY: NCERROR
USE NETCDF
#endif
IMPLICIT NONE
!* local variables
INTEGER(KIND=JPIM),ALLOCATABLE  :: REGIONGRID(:)
!
INTEGER(KIND=JPIM),SAVE         :: IX,IY
INTEGER(KIND=JPIM),SAVE         :: IREGION
INTEGER(KIND=JPIM),SAVE         :: NSEQMAX_MPI
#ifdef UseMPI_CMF
#ifdef UseCDF_CMF
INTEGER(KIND=JPIM)              :: NCID
INTEGER(KIND=JPIM)              :: VARID
#endif
#endif
!$OMP THREADPRIVATE               (IX)
!================================================
WRITE(LOGNAM,*) 'RIVMAP_INIT: region code'

!*** when MPI is not used. Treat all grid as region=1
REGIONALL=1
I2REGION(:,:)=IMIS
!$OMP PARALLEL DO
DO IY=1, NY
  DO IX=1, NX
    IF( I2NEXTX(IX,IY)/=IMIS ) THEN
      I2REGION(IX,IY)=1
    ENDIF
  END DO
END DO
!$OMP END PARALLEL DO

!*** read MPI region map
!! Use MPI: read MPI region map, allocate regions to MPI nodes
#ifdef UseMPI_CMF
  IF ( LMAPCDF ) THEN
#ifdef UseCDF_CMF
    CALL NCERROR (NF90_OPEN(CMPIREGNC,NF90_NOWRITE,NCID),'opening '//TRIM(CMPIREGNC) )
    CALL NCERROR (NF90_INQ_VARID(NCID, 'mpireg',VARID),'getting id' )
    CALL NCERROR (NF90_GET_VAR(NCID,VARID,I2REGION),'reading data' )
    CALL NCERROR (NF90_CLOSE(NCID))
#endif
  ELSE
    WRITE(LOGNAM,*)'RIVMAP_INIT: read MPI region: ',TRIM(CNEXTXY)
    TMPNAM=INQUIRE_FID()
    OPEN(TMPNAM,FILE=CMPIREG,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
    READ(TMPNAM,REC=1) I2REGION
    CLOSE(TMPNAM)
  ENDIF

  REGIONALL=1
!$OMP PARALLEL DO REDUCTION(max:REGIONALL)
  DO IY=1, NY
    DO IX=1, NX
      REGIONALL=MAX( REGIONALL, I2REGION(IX,IY) )
    END DO
  END DO
!$OMP END PARALLEL DO
#endif


WRITE(LOGNAM,*)'RIVMAP_INIT: count number of grid in each region: '
ALLOCATE(REGIONGRID(REGIONALL))
REGIONGRID(:)=0
!! OMP reduction operation for array might not be available in some environment
DO IY=1, NY
  DO IX=1, NX
    IF( I2REGION(IX,IY)>0 ) THEN
      IREGION=I2REGION(IX,IY)
      REGIONGRID(IREGION)=REGIONGRID(IREGION)+1
    ENDIF
  END DO
END DO

NSEQMAX=0
NSEQMAX_MPI=0
DO IREGION=1, REGIONALL
  IF( IREGION==REGIONTHIS ) NSEQMAX=REGIONGRID(IREGION)
  NSEQMAX_MPI=MAX(NSEQMAX,REGIONGRID(IREGION))  !! maximum nseqmax among all MPI region
END DO

WRITE(LOGNAM,*) 'CALC_REGION: REGIONALL= ', REGIONALL
WRITE(LOGNAM,*) 'CALC_REGION: NSEQMAX='   , NSEQMAX
WRITE(LOGNAM,*) 'CALC_REGION: NSEQMAX='   , NSEQMAX_MPI

END SUBROUTINE CALC_REGION
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE CALC_1D_SEQ
! OpenMP is not used, because results of this subroutine highly depents on calculation order
USE YOS_CMF_INPUT,           ONLY: IMIS
IMPLICIT NONE
!* local variables
INTEGER(KIND=JPIM)              :: IX,IY,JX,JY,ISEQ,JSEQ,ISEQ1,ISEQ2,AGAIN
INTEGER(KIND=JPIM)              :: UPNMAX

INTEGER(KIND=JPIM),ALLOCATABLE  :: NUPST(:,:), UPNOW(:,:)
!================================================
WRITE(LOGNAM,*) 'RIVMAP_INIT: convert 2D map to 1D sequence'

ALLOCATE( NUPST(NX,NY) )
ALLOCATE( UPNOW(NX,NY) )

ALLOCATE( I1SEQX(NSEQMAX) )
ALLOCATE( I1SEQY(NSEQMAX) )
ALLOCATE( I1NEXT(NSEQMAX) )
ALLOCATE( I2VECTOR(NX,NY) )
I1SEQX(:)=0
I1SEQY(:)=0
I1NEXT(:)=0
I2VECTOR(:,:)=0

! count number of upstream 
NUPST(:,:)=0
UPNOW(:,:)=0
UPNMAX=0
DO IY=1, NY
  DO IX=1, NX
    IF( I2NEXTX(IX,IY).GT.0 .and. I2REGION(IX,IY)==REGIONTHIS )THEN
      JX=I2NEXTX(IX,IY)
      JY=I2NEXTY(IX,IY)
      NUPST(JX,JY)=NUPST(JX,JY)+1
      UPNMAX=max(UPNMAX,NUPST(JX,JY))
    ENDIF
  END DO
END DO

! register upmost grid in 1d sequence
ISEQ=0
DO IY=1, NY
  DO IX=1, NX
    IF( I2NEXTX(IX,IY).GT.0 .and. I2REGION(IX,IY)==REGIONTHIS )THEN
      IF( NUPST(IX,IY)==UPNOW(IX,IY) )THEN
        ISEQ=ISEQ+1
        I1SEQX(ISEQ)=IX
        I1SEQY(ISEQ)=IY
        I2VECTOR(IX,IY)=ISEQ
      ENDIF
    ENDIF
  END DO
END DO
ISEQ1=1
ISEQ2=ISEQ

AGAIN=1
DO WHILE( AGAIN==1 )
  AGAIN=0
  JSEQ=ISEQ2
  DO ISEQ=ISEQ1, ISEQ2
    IX=I1SEQX(ISEQ)
    IY=I1SEQY(ISEQ)
    JX=I2NEXTX(IX,IY)
    JY=I2NEXTY(IX,IY)
    UPNOW(JX,JY)=UPNOW(JX,JY)+1
    IF( UPNOW(JX,JY)==NUPST(JX,JY) .and. I2NEXTX(JX,JY)>0 )THEN !! if all upstream calculated, register to 1D sequence
      JSEQ=JSEQ+1
      I1SEQX(JSEQ)=JX
      I1SEQY(JSEQ)=JY
      I2VECTOR(JX,JY)=JSEQ
      AGAIN=1
    ENDIF
  END DO
  ISEQ1=ISEQ2+1
  ISEQ2=JSEQ
END DO
NSEQRIV=JSEQ    !! END OF RIVER-LINK GRID

ISEQ=NSEQRIV
DO IY=1, NY
  DO IX=1, NX
    IF( I2NEXTX(IX,IY).LT.0 .AND. I2NEXTX(IX,IY).NE.IMIS .AND. I2REGION(IX,IY)==REGIONTHIS )THEN
      ISEQ=ISEQ+1
      I1SEQX(ISEQ)=IX
      I1SEQY(ISEQ)=IY
      I2VECTOR(IX,IY)=ISEQ
    ENDIF
  END DO
END DO
NSEQALL=ISEQ    !! END OF RIVER-MOUTH GRID

DO ISEQ=1, NSEQALL
  IX=I1SEQX(ISEQ)
  IY=I1SEQY(ISEQ)
  IF( I2NEXTX(IX,IY)>0 )THEN
    JX=I2NEXTX(IX,IY)
    JY=I2NEXTY(IX,IY)
    I1NEXT(ISEQ)=I2VECTOR(JX,JY)
  ELSE
    I1NEXT(ISEQ)=I2NEXTX(IX,IY)
  ENDIF
END DO

DEALLOCATE(NUPST,UPNOW)

WRITE(LOGNAM,*) 'CALC_1DSEQ: NSEQRIV='   , NSEQRIV
WRITE(LOGNAM,*) 'CALC_1DSEQ: NSEQALL='   , NSEQALL

!====
!!IF( LSPAMAT )THEN
  WRITE(LOGNAM,*) 'RIVMAP_INIT: Upstream Matrix size', UPNMAX
  ALLOCATE(I1UPST(NSEQMAX,UPNMAX),I1UPN(NSEQMAX))
  I1UPST(:,:)=-9999
  I1UPN(:)=0
  DO ISEQ=1, NSEQRIV
    JSEQ=I1NEXT(ISEQ)
    I1UPN(JSEQ)=I1UPN(JSEQ)+1
    I1UPST(JSEQ,I1UPN(JSEQ))=ISEQ
  END DO
!!ENDIF

END SUBROUTINE CALC_1D_SEQ
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE READ_BIFPARAM    !! evenly allocate pixels to mpi nodes (not used in vcurrent version)
USE YOS_CMF_INPUT,      ONLY: PMANRIV, PMANFLD
USE YOS_CMF_MAP,        ONLY: NPTHOUT, NPTHLEV, PTH_UPST, PTH_DOWN,&
                            & PTH_DST, PTH_ELV, PTH_WTH,  PTH_MAN, NSEQMAX
USE CMF_UTILS_MOD,      ONLY: INQUIRE_FID
IMPLICIT NONE
!* local variables
INTEGER(KIND=JPIM)         :: IX,IY, JX,JY, ISEQ, JSEQ
INTEGER(KIND=JPIM)         :: IPTH,  ILEV,  NPTHOUT1
REAL(KIND=JPRB)            :: PELV,  PWTH,  PDPH
INTEGER(KIND=JPIM)         :: ONMAX,  INMAX
INTEGER(KIND=JPIM)         :: OCOUNT(NSEQMAX),  ICOUNT(NSEQMAX)
!================================================
WRITE(LOGNAM,*)"RIVMAP_INIT: Bifuraction channel:", TRIM(CPTHOUT)

TMPNAM=INQUIRE_FID()
OPEN(TMPNAM,FILE=CPTHOUT,FORM='FORMATTED')
READ(TMPNAM,*) NPTHOUT,NPTHLEV

WRITE(LOGNAM,*) "Bifurcation channel dimantion", NPTHOUT, NPTHLEV

ALLOCATE( PTH_UPST(NPTHOUT) )
ALLOCATE( PTH_DOWN(NPTHOUT) )
ALLOCATE( PTH_DST(NPTHOUT)  )
ALLOCATE( PTH_ELV(NPTHOUT,NPTHLEV) )
ALLOCATE( PTH_WTH(NPTHOUT,NPTHLEV) )
ALLOCATE( PTH_MAN(NPTHLEV)  )

OCOUNT(:)=0
ICOUNT(:)=0

NPTHOUT1=0
DO IPTH=1, NPTHOUT
  READ(TMPNAM,*) IX, IY, JX, JY, PTH_DST(IPTH), PELV, PDPH, (PTH_WTH(IPTH,ILEV),ILEV=1,NPTHLEV)
  PTH_UPST(IPTH)=I2VECTOR(IX,IY)
  PTH_DOWN(IPTH)=I2VECTOR(JX,JY)
  IF (PTH_UPST(IPTH) > 0 .AND. PTH_DOWN(IPTH) > 0) THEN  ! both upst and down grids are within domain
    NPTHOUT1=NPTHOUT1+1
    ISEQ=PTH_UPST(IPTH)
    JSEQ=PTH_DOWN(IPTH)
    OCOUNT(ISEQ)=OCOUNT(ISEQ)+1
    ICOUNT(JSEQ)=ICOUNT(JSEQ)+1
  ENDIF
  DO ILEV=1, NPTHLEV
    IF( ILEV==1 )THEN            !!ILEV=1: water channel bifurcation. consider bifurcation channel depth
      PWTH=PTH_WTH(IPTH,ILEV)
      IF( PWTH>0 )then
        PTH_ELV(IPTH,ILEV)=PELV - PDPH
      ELSE
        PTH_ELV(IPTH,ILEV)=1.E20
      ENDIF
    ELSE
      PWTH=PTH_WTH(IPTH,ILEV)
      IF( PWTH>0 )then
        PTH_ELV(IPTH,ILEV)=PELV + ILEV - 2.0    !! ILEV=2: bank top level 
      ELSE
        PTH_ELV(IPTH,ILEV)=1.E20
      ENDIF
    ENDIF
  END DO
END DO
CLOSE(TMPNAM)

DO ILEV=1, NPTHLEV
  IF( ILEV==1 )THEN
    PTH_MAN(ILEV)=PMANRIV
  ELSE
    PTH_MAN(ILEV)=PMANFLD
  ENDIF
END DO

IF (NPTHOUT /= NPTHOUT1) THEN
  WRITE(LOGNAM,*)"Bifuraction channel outside of domain. Only valid:", NPTHOUT1
ENDIF

!=====
!!IF( LSPAMAT )THEN
  ONMAX=0
  INMAX=0
  DO ISEQ=1, NSEQMAX
    ONMAX=max(ONMAX,OCOUNT(ISEQ))
    INMAX=max(INMAX,ICOUNT(ISEQ))
  END DO
  WRITE(LOGNAM,*) "Max # of Bifuraction channel in each grid (Out, Inf):", ONMAX, INMAX
  
  ALLOCATE(I1P_OUT(NSEQMAX,ONMAX))
  ALLOCATE(I1P_INF(NSEQMAX,INMAX))
  ALLOCATE(I1P_OUTN(NSEQMAX))
  ALLOCATE(I1P_INFN(NSEQMAX))
  I1P_OUT(:,:)=0
  I1P_INF(:,:)=0
  I1P_OUTN(:)=0
  I1P_INFN(:)=0

  DO IPTH=1, NPTHOUT
    IF (PTH_UPST(IPTH) > 0 .AND. PTH_DOWN(IPTH) > 0) THEN  ! both upst and down grids are within domain
      ISEQ=PTH_UPST(IPTH)
      JSEQ=PTH_DOWN(IPTH)
  
      I1P_OUTN(ISEQ)=I1P_OUTN(ISEQ)+1
      I1P_OUT(ISEQ,I1P_OUTN(ISEQ))=IPTH
  
      I1P_INFN(JSEQ)=I1P_INFN(JSEQ)+1
      I1P_INF(JSEQ,I1P_INFN(JSEQ))=IPTH
    ENDIF
  END DO
!!ENDIF

END SUBROUTINE READ_BIFPARAM
!==========================================================


END SUBROUTINE CMF_RIVMAP_INIT
!####################################################################





!####################################################################
SUBROUTINE CMF_TOPO_INIT
! read & set topography map 
! -- call from CMF_DRV_INIT
USE YOS_CMF_INPUT,  ONLY: TMPNAM,   NX,NY,NLFP, LMAPEND,  &
                        & LFPLAIN,  LMEANSL,  LGDWDLY,  LSLPMIX, LSLOPEMOUTH
USE YOS_CMF_MAP,    ONLY: D2NXTDST, D2GRAREA, D2ELEVTN, D2RIVLEN, &
                        & D2RIVWTH, D2RIVHGT, D2FLDHGT, D2RIVELV, &
                        & D2FLDGRD, D2RIVMAN, D2RIVSTOMAX, D2FLDSTOMAX,  &
                        & DFRCINC,  NSEQMAX,  D2MEANSL, D2DWNELV, &
                        & D2GDWDLY, I2MASK
IMPLICIT NONE
!================================================
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!---------------------!"

WRITE(LOGNAM,*) 'CMF::TOPO_INIT: topography map initialization'

! *** 1. ALLOCATE ARRAYS
ALLOCATE( D2GRAREA(NSEQMAX,1) )
ALLOCATE( D2ELEVTN(NSEQMAX,1) )
ALLOCATE( D2NXTDST(NSEQMAX,1) )
ALLOCATE( D2RIVLEN(NSEQMAX,1) )
ALLOCATE( D2RIVWTH(NSEQMAX,1) )
ALLOCATE( D2RIVHGT(NSEQMAX,1) )
ALLOCATE( D2FLDHGT(NSEQMAX,1,NLFP) )
ALLOCATE( D2RIVMAN(NSEQMAX,1) )
ALLOCATE( D2MEANSL(NSEQMAX,1) )
ALLOCATE( D2DWNELV(NSEQMAX,1) )
ALLOCATE( D2GDWDLY(NSEQMAX,1) )
ALLOCATE( I2MASK(NSEQMAX,1) )

D2GRAREA(:,:)  =0._JPRB
D2ELEVTN(:,:)  =0._JPRB
D2NXTDST(:,:)  =0._JPRB
D2RIVLEN(:,:)  =0._JPRB
D2RIVWTH(:,:)  =0._JPRB
D2RIVHGT(:,:)  =0._JPRB
D2FLDHGT(:,:,:)=0._JPRB
D2RIVMAN(:,:)  =0._JPRB
D2MEANSL(:,:)  =0._JPRB
D2DWNELV(:,:)  =0._JPRB
D2GDWDLY(:,:)  =0._JPRB
I2MASK(:,:)    =0._JPIM     !! mask for calculation (IFS slopemix: Kinemacti Wave for Mask=1; Reservoir: dam=2, dam upstream=1)

!============================
! *** 2. Read topo map
WRITE(LOGNAM,*) 'CMF::TOPO_INIT: read topography maps'
IF ( .not. LMAPCDF ) THEN
  CALL READ_TOPO_BIN
ELSE
  CALL READ_TOPO_CDF
ENDIF

I2MASK(:,:)=0_JPIM
IF ( LSLPMIX ) THEN
  CALL SET_SLOPEMIX
ENDIF

!============================
! *** 3a. Calc Channel Parameters
WRITE(LOGNAM,*) 'TOPO_INIT: calc river channel parameters'

ALLOCATE(D2RIVSTOMAX(NSEQMAX,1))
ALLOCATE(D2RIVELV(NSEQMAX,1))

IF ( LFPLAIN ) THEN
  D2RIVSTOMAX(:,:) = D2RIVLEN(:,:) * D2RIVWTH(:,:) * D2RIVHGT(:,:)
ELSE
  WRITE(LOGNAM,*) 'TOPO_INIT: no floodplain (rivstomax=1.D18)'
  D2RIVSTOMAX(:,:) = 1.E18
ENDIF
D2RIVELV(:,:) = D2ELEVTN(:,:) - D2RIVHGT(:,:)

!*** 3b. Calc Channel Parameters
WRITE(LOGNAM,*) 'TOPO_INIT: calc floodplain parameters'

ALLOCATE(D2FLDSTOMAX(NSEQMAX,1,NLFP))
ALLOCATE(D2FLDGRD(NSEQMAX,1,NLFP))
CALL SET_FLDSTG

!*** 3c. Calc downstream boundary
WRITE(LOGNAM,*) 'TOPO_INIT: calc downstream boundary elevation'
D2DWNELV(:,:)=D2ELEVTN(:,:)
IF( LMEANSL ) THEN
  D2DWNELV(:,:)=D2ELEVTN(:,:)+D2MEANSL(:,:)
ENDIF

CONTAINS
!==========================================================
!+ READ_TOPO_BIN
!+ READ_TOPO_CDF
!+ SET_FLDSTG
!+ SET_SLOPEMIX
!==========================================================
SUBROUTINE READ_TOPO_BIN
USE CMF_UTILS_MOD,       ONLY: mapR2vecD, CONV_END,  INQUIRE_FID
IMPLICIT NONE
!* local variables
INTEGER(KIND=JPIM)          :: ILFP
REAL(KIND=JPRM),ALLOCATABLE :: R2TEMP(:,:)
REAL(KIND=JPRB),ALLOCATABLE :: D2TEMP(:,:)
!================================================
ALLOCATE(R2TEMP(NX,NY))
ALLOCATE(D2TEMP(NSEQMAX,1))

TMPNAM=INQUIRE_FID()

WRITE(LOGNAM,*)'TOPO_INIT: unit-catchment area : ',TRIM(CGRAREA) 
OPEN(TMPNAM,FILE=CGRAREA,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
READ(TMPNAM,REC=1) R2TEMP(:,:)
  IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
CALL mapR2vecD(R2TEMP,D2GRAREA)
CLOSE(TMPNAM)

WRITE(LOGNAM,*)'TOPO_INIT: ground elevation : ',TRIM(CELEVTN)
OPEN(TMPNAM,FILE=CELEVTN,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
READ(TMPNAM,REC=1) R2TEMP(:,:)
  IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
CALL mapR2vecD(R2TEMP,D2ELEVTN)
CLOSE(TMPNAM)

WRITE(LOGNAM,*)'TOPO_INIT: downstream distance : ',TRIM(CNXTDST)
OPEN(TMPNAM,FILE=CNXTDST,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
READ(TMPNAM,REC=1) R2TEMP(:,:)
  IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
CALL mapR2vecD(R2TEMP,D2NXTDST)
CLOSE(TMPNAM)

WRITE(LOGNAM,*)'TOPO_INIT: river channel length : ',TRIM(CRIVLEN)
OPEN(TMPNAM,FILE=CRIVLEN,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
READ(TMPNAM,REC=1) R2TEMP(:,:)
  IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
CALL mapR2vecD(R2TEMP,D2RIVLEN)
CLOSE(TMPNAM)

WRITE(LOGNAM,*)'TOPO_INIT: floodplain elevation profile : ',TRIM(CFLDHGT)
OPEN(TMPNAM,FILE=TRIM(CFLDHGT),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
DO ILFP=1,NLFP
  READ(TMPNAM,REC=ILFP) R2TEMP
  IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
  CALL mapR2vecD(R2TEMP,D2TEMP)
  D2FLDHGT(:,:,ILFP)= D2TEMP(:,:)
ENDDO
CLOSE(TMPNAM)

!*** river channel / groundwater parameters)

WRITE(LOGNAM,*)'TOPO_INIT: river channel depth : ',TRIM(CRIVHGT)
OPEN(TMPNAM,FILE=CRIVHGT,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
READ(TMPNAM,REC=1) R2TEMP(:,:)
  IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
CALL mapR2vecD(R2TEMP,D2RIVHGT)
CLOSE(TMPNAM)

WRITE(LOGNAM,*)'TOPO_INIT: river channel width : ',TRIM(CRIVWTH)
OPEN(TMPNAM,FILE=CRIVWTH,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
READ(TMPNAM,REC=1) R2TEMP(:,:)
  IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
CALL mapR2vecD(R2TEMP,D2RIVWTH)
CLOSE(TMPNAM)

WRITE(LOGNAM,*)'TOPO_INIT: manning coefficient river: ',TRIM(CRIVMAN)
OPEN(TMPNAM,FILE=TRIM(CRIVMAN),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
READ(TMPNAM,REC=1) R2TEMP(:,:)
  IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
CALL mapR2vecD(R2TEMP,D2RIVMAN)
CLOSE(TMPNAM)

IF( LGDWDLY )THEN
  WRITE(LOGNAM,*)'TOPO_INIT: groundwater delay parameter: ',TRIM(CGDWDLY)
  OPEN(TMPNAM,FILE=TRIM(CGDWDLY),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=4*NX*NY)
  READ(TMPNAM,REC=1) R2TEMP(:,:)
    IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
  CALL mapR2vecD(R2TEMP,D2GDWDLY)
  CLOSE(TMPNAM)
ENDIF

IF( LSLPMIX )THEN
  WRITE(LOGNAM,*)'TOPO_INIT: LSLPMIX only used in IFS, not availabke with binary map'
ENDIF
IF( LSLOPEMOUTH )THEN
  WRITE(LOGNAM,*)'TOPO_INIT: LSLOPEMOUTH only used in IFS, not availabke with binary map'
ENDIF

! ==========

IF( LMEANSL ) THEN
  WRITE(LOGNAM, *)'TOPO_INIT: mean sea level: ', TRIM(CMEANSL)
  OPEN(TMPNAM, FILE=CMEANSL, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=4*NX*NY)
  READ(TMPNAM, REC=1) R2TEMP(:,:)
    IF( LMAPEND ) CALL CONV_END(R2TEMP,NX,NY)
  CALL mapR2vecD(R2TEMP, D2MEANSL)
  CLOSE(TMPNAM)
ENDIF

DEALLOCATE(R2TEMP)
DEALLOCATE(D2TEMP)

END SUBROUTINE READ_TOPO_BIN
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE READ_TOPO_CDF
#ifdef UseCDF_CMF
USE NETCDF 
USE CMF_UTILS_MOD,            ONLY: NCERROR,mapR2vecD
USE YOS_CMF_MAP,              ONLY: D2ELEVSLOPE     !! only used in ECMWF
IMPLICIT NONE
!* local variables
INTEGER(KIND=JPIM)               :: NCID,VARID,STATUS
INTEGER(KIND=JPIM)               :: ILEV
REAL(KIND=JPRM),ALLOCATABLE      :: R2TEMP(:,:)
REAL(KIND=JPRB),ALLOCATABLE      :: D2TEMP(:,:)
!================================================
ALLOCATE(R2TEMP(NX,NY))
ALLOCATE(D2TEMP(NSEQMAX,1))

!! CLIM FILE
CALL NCERROR (NF90_OPEN(CRIVCLINC,NF90_NOWRITE,NCID),'opening '//TRIM(CRIVCLINC) )

WRITE(LOGNAM,*)'TOPO_INIT: ctmare:',TRIM(CRIVCLINC)
CALL NCERROR ( NF90_INQ_VARID(NCID,'ctmare',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' )
CALL mapR2vecD(R2TEMP,D2GRAREA)

WRITE(LOGNAM,*)'TOPO_INIT: elevtn:',TRIM(CRIVCLINC)
CALL NCERROR ( NF90_INQ_VARID(NCID,'elevtn',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
CALL mapR2vecD(R2TEMP,D2ELEVTN)

WRITE(LOGNAM,*)'TOPO_INIT: nxtdst:',TRIM(CRIVCLINC)
CALL NCERROR ( NF90_INQ_VARID(NCID,'nxtdst',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
CALL mapR2vecD(R2TEMP,D2NXTDST)

WRITE(LOGNAM,*)'TOPO_INIT: rivlen:',TRIM(CRIVCLINC)
CALL NCERROR ( NF90_INQ_VARID(NCID,'rivlen',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
CALL mapR2vecD(R2TEMP,D2RIVLEN)

WRITE(LOGNAM,*)'TOPO_INIT: fldhgt:',TRIM(CRIVCLINC)
CALL NCERROR ( NF90_INQ_VARID(NCID,'fldhgt',VARID),'getting id' )
DO ILEV=1,NLFP
  CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP,(/1,1,ILEV/),(/NX,NY,1/)),'reading data' ) 
  CALL mapR2vecD(R2TEMP,D2TEMP)
  D2FLDHGT(:,:,ILEV)=D2TEMP(:,:)
ENDDO

CALL NCERROR( NF90_CLOSE(NCID))

IF ( LSLOPEMOUTH ) THEN
  ALLOCATE( D2ELEVSLOPE(NSEQMAX,1) )
  WRITE(LOGNAM,*)'TOPO_INIT: elevslope:',TRIM(CRIVPARNC)
  STATUS = NF90_INQ_VARID(NCID,'elevslope',VARID)
  IF (STATUS /= 0 ) THEN
    WRITE(LOGNAM,*)'TOPO_INIT: elevslope: not present, aborting'
    STOP 9 
  ELSE
    CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
  ENDIF 
  CALL mapR2vecD(R2TEMP,D2ELEVSLOPE)
ENDIF

!!========== 
!! PAR FILE (river channel / groundwater parameters)
CALL NCERROR (NF90_OPEN(CRIVPARNC,NF90_NOWRITE,NCID),'opening '//TRIM(CRIVPARNC) )

WRITE(LOGNAM,*)'TOPO_INIT: rivwth:',TRIM(CRIVPARNC)
CALL NCERROR ( NF90_INQ_VARID(NCID,'rivwth',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
CALL mapR2vecD(R2TEMP,D2RIVWTH)

WRITE(LOGNAM,*)'TOPO_INIT: rivhgt:',TRIM(CRIVPARNC)
CALL NCERROR ( NF90_INQ_VARID(NCID,'rivhgt',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
CALL mapR2vecD(R2TEMP,D2RIVHGT)

WRITE(LOGNAM,*)'TOPO_INIT: rivman:',TRIM(CRIVPARNC)
CALL NCERROR ( NF90_INQ_VARID(NCID,'rivman',VARID),'getting id' )
CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
CALL mapR2vecD(R2TEMP,D2RIVMAN)

IF ( LGDWDLY ) THEN
  WRITE(LOGNAM,*)'TOPO_INIT: GDWDLY:',TRIM(CRIVPARNC)
  STATUS = NF90_INQ_VARID(NCID,'gdwdly',VARID)
  IF (STATUS /= 0 ) THEN
    WRITE(LOGNAM,*)'TOPO_INIT: GDWDLY: not present, setting to zero'
    R2TEMP(:,:) = 0._JPRB
  ELSE
    CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
  ENDIF 
  CALL mapR2vecD(R2TEMP,D2GDWDLY)
ENDIF

CALL NCERROR( NF90_CLOSE(NCID))

!!========== 
!! MEAN SEA LEVEL FILE
IF( LMEANSL ) THEN
  CALL NCERROR (NF90_OPEN(CMEANSLNC,NF90_NOWRITE,NCID),'opening '//TRIM(CMEANSLNC) )
  WRITE(LOGNAM,*)'TOPO_INIT: rivhgt:',TRIM(CMEANSLNC)
  CALL NCERROR ( NF90_INQ_VARID(NCID,'meansl',VARID),'getting id' )
  CALL NCERROR ( NF90_GET_VAR(NCID,VARID,R2TEMP),'reading data' ) 
  CALL mapR2vecD ( R2TEMP,D2MEANSL  )
  CALL NCERROR ( NF90_CLOSE(NCID) )
ENDIF 

DEALLOCATE(R2TEMP)
DEALLOCATE(D2TEMP)
#endif
END SUBROUTINE READ_TOPO_CDF
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE SET_FLDSTG
IMPLICIT NONE
!* local variables
INTEGER(KIND=JPIM),SAVE  ::  ISEQ, I
REAL(KIND=JPRB),SAVE     ::  DSTONOW
REAL(KIND=JPRB),SAVE     ::  DSTOPRE
REAL(KIND=JPRB),SAVE     ::  DHGTPRE
REAL(KIND=JPRB),SAVE     ::  DWTHINC
!$OMP THREADPRIVATE               (I,DSTONOW,DSTOPRE,DHGTPRE,DWTHINC)
!================================================
D2FLDSTOMAX(:,:,:) = 0._JPRD
D2FLDGRD(:,:,:)    = 0._JPRB
DFRCINC=REAL(NLFP,KIND=JPRB)**(-1.)
!
!$OMP PARALLEL DO
DO ISEQ=1, NSEQMAX
  DSTOPRE = D2RIVSTOMAX(ISEQ,1)
  DHGTPRE = 0._JPRB
  DWTHINC = D2GRAREA(ISEQ,1) * D2RIVLEN(ISEQ,1)**(-1.) * DFRCINC
  DO I=1, NLFP
    DSTONOW = D2RIVLEN(ISEQ,1) * ( D2RIVWTH(ISEQ,1) + DWTHINC*(REAL(I,KIND=JPRB)-0.5) ) * (D2FLDHGT(ISEQ,1,I)-DHGTPRE)
    D2FLDSTOMAX(ISEQ,1,I) = DSTOPRE + DSTONOW
    D2FLDGRD(ISEQ,1,I) = (D2FLDHGT(ISEQ,1,I)-DHGTPRE) * DWTHINC**(-1.)
    DSTOPRE = D2FLDSTOMAX(ISEQ,1,I)
    DHGTPRE = D2FLDHGT(ISEQ,1,I)
  END DO
END DO
!$OMP END PARALLEL DO

!
END SUBROUTINE SET_FLDSTG
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE SET_SLOPEMIX    !! only used in IFS0
#ifdef UseCDF_CMF
USE NETCDF 
USE CMF_UTILS_MOD,           ONLY: NCERROR,mapI2vecI

IMPLICIT NONE
INTEGER(KIND=JPIM),ALLOCATABLE  :: I2TEMP(:,:)
INTEGER(KIND=JPIM)              :: ISEQ, I0, I1
INTEGER(KIND=JPIM)              :: NCID,VARID,STATUS

  ALLOCATE(I2TEMP(NX,NY))
  WRITE(LOGNAM,*)'TOPO_INIT: mask_slope:',TRIM(CRIVPARNC)
  STATUS =  NF90_INQ_VARID(NCID,'mask_slope',VARID)
  IF (STATUS /= 0 ) THEN
    WRITE(LOGNAM,*)'TOPO_INIT: mask_slope: LSLPMIX should be set to FALSE: ABORTING!'
    STOP 9
  ENDIF 
  CALL NCERROR ( NF90_GET_VAR(NCID,VARID,I2TEMP),'reading data' ) 
  CALL mapI2vecI(I2TEMP,I2MASK)
  I0=0
  I1=0
  DO ISEQ=1,NSEQMAX
    IF (I2MASK(ISEQ,1) == 1 ) THEN  !! kinematic wave applied
      I1=I1+1
    ENDIF
    IF (I2MASK(ISEQ,1) == 0 ) THEN 
      I0=I0+1
    ENDIF
  ENDDO
  WRITE(LOGNAM,*)'TOPO_INIT: sum(mask==0), sum(mask==1)',I0,I1
  IF ( I0+I1 .NE. NSEQMAX ) THEN 
     WRITE(LOGNAM,*)'TOPO_INIT: mask==0 + mask == 1 does not match NSEQMAX.. something wrong, aborting'
     STOP 9
  ENDIF 

  DEALLOCATE(I2TEMP)
#endif
END SUBROUTINE SET_SLOPEMIX

END SUBROUTINE CMF_TOPO_INIT
!####################################################################


END MODULE CMF_CTRL_MAPS_MOD
