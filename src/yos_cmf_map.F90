MODULE YOS_CMF_MAP
!==========================================================
!* PURPOSE: Shared map/topography variables for CaMa-Flood
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
USE PARKIND1, ONLY: JPIM, JPRM, JPRB,JPRD
IMPLICIT NONE
SAVE 
!================================================
!*** river network
INTEGER(KIND=JPIM),ALLOCATABLE           ::  I2NEXTX(:,:)       !! POINT DOWNSTREAM HORIZONTAL
INTEGER(KIND=JPIM),ALLOCATABLE           ::  I2NEXTY(:,:)       !! POINT DOWNSTREAM VERTICAL

INTEGER(KIND=JPIM),ALLOCATABLE           ::  I1SEQX(:)          !! 1D SEQUENCE HORIZONTAL
INTEGER(KIND=JPIM),ALLOCATABLE           ::  I1SEQY(:)          !! 1D SEQUENCE VERTICAL
INTEGER(KIND=JPIM),ALLOCATABLE           ::  I1NEXT(:)          !! 1D DOWNSTREAM
INTEGER(KIND=JPIM)                       ::  NSEQRIV            !! END OF RIVER-LINK  SEQUENCE (1 ~ NSEQRIV)
INTEGER(KIND=JPIM)                       ::  NSEQALL            !! END OF RIVER-MOUTH SEQUENCE (NSEQRUV+1 ~ NSEQALL)
INTEGER(KIND=JPIM)                       ::  NSEQMAX            !! TOTAL GRIDS ON RIVER MAP (lakes added in future dev)

INTEGER(KIND=JPIM),ALLOCATABLE           ::  I2VECTOR(:,:)      !! VECTOR INDEX
INTEGER(KIND=JPIM),ALLOCATABLE           ::  I2REGION(:,:)      !! REGION INDEX
INTEGER(KIND=JPIM)                       ::  REGIONALL          !! REGION TOTAL
INTEGER(KIND=JPIM)                       ::  REGIONTHIS         !! REGION THIS CPU
INTEGER(KIND=JPIM)                       ::  MPI_COMM_CAMA      !! MPI COMMUNICATOR

INTEGER(KIND=JPIM),ALLOCATABLE           ::  I1UPST(:,:)        !! UPSTREAM SEQ
INTEGER(KIND=JPIM),ALLOCATABLE           ::  I1UPN(:)           !! MAX UPSTREAM NUMBER

INTEGER(KIND=JPIM),ALLOCATABLE           ::  I1P_OUT(:,:)        !! BIF PTH OUT
INTEGER(KIND=JPIM),ALLOCATABLE           ::  I1P_OUTN(:)           !! MAX UPSTREAM NUMBER

INTEGER(KIND=JPIM),ALLOCATABLE           ::  I1P_INF(:,:)        !! UPSTREAM SEQ
INTEGER(KIND=JPIM),ALLOCATABLE           ::  I1P_INFN(:)           !! MAX UPSTREAM NUMBER


!================================================
!*** lat, lon
REAL(KIND=JPRB),ALLOCATABLE              ::  D1LON(:)           !! longitude [degree_east]
REAL(KIND=JPRB),ALLOCATABLE              ::  D1LAT(:)           !! latitude  [degree_north]

!================================================
!*** River + Floodplain topography (map)
REAL(KIND=JPRB),ALLOCATABLE              ::  D2GRAREA(:,:)      !! GRID AREA [M2]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2ELEVTN(:,:)      !! ELEVATION [M]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2NXTDST(:,:)      !! DISTANCE TO THE NEXT GRID [M]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2RIVLEN(:,:)      !! RIVER LENGTH [M]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2RIVWTH(:,:)      !! RIVER WIDTH [M]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2RIVMAN(:,:)      !! RIVER MANNING COEFFICIENT
REAL(KIND=JPRB),ALLOCATABLE              ::  D2RIVHGT(:,:)      !! RIVER HEIGHT [M]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2FLDHGT(:,:,:)    !! FLOODPLAIN HEIGHT [M]

REAL(KIND=JPRB),ALLOCATABLE              ::  D2GDWDLY(:,:)      !! Ground water delay
REAL(KIND=JPRB),ALLOCATABLE              ::  D2ELEVSLOPE(:,:)   !! River bed slope
INTEGER(KIND=JPIM),ALLOCATABLE           ::  I2MASK(:,:)        !! Mask 

!================================================
!*** Floodplain Topography (diagnosed)
REAL(KIND=JPRB),ALLOCATABLE              ::  D2RIVSTOMAX(:,:)   !! maximum river storage [m3]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2FLDSTOMAX(:,:,:) !! MAXIMUM FLOODPLAIN STORAGE [M3]

REAL(KIND=JPRB),ALLOCATABLE              ::  D2RIVELV(:,:)      !! elevation of river bed [m3]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2FLDGRD(:,:,:)    !! FLOODPLAIN GRADIENT
REAL(KIND=JPRB)                          ::  DFRCINC            !! FLOODPLAIN FRACTION INCREMENT [-] (1/NLFP)

!================================================
!*** Downstream boundary
REAL(KIND=JPRB),ALLOCATABLE              ::  D2MEANSL(:,:)      !! MEAN SEA LEVEL [M]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2SEALEV(:,:)        !! sea level variation [m]
REAL(KIND=JPRB),ALLOCATABLE              ::  D2DWNELV(:,:)        !! downstream boundary elevation [m]

!================================================
!*** bifurcation channel
INTEGER(KIND=JPIM)                       ::  NPTHOUT            !! NUMBER OF FLOODPLAIN PATH
INTEGER(KIND=JPIM)                       ::  NPTHLEV            !! NUMBER OF FLOODPLAIN PATH LAYER
INTEGER(KIND=JPIM),ALLOCATABLE           ::  PTH_UPST(:)        !! FLOOD PATHWAY UPSTREAM   ISEQ
INTEGER(KIND=JPIM),ALLOCATABLE           ::  PTH_DOWN(:)        !! FLOOD PATHWAY DOWNSTREAM JSEQ
REAL(KIND=JPRB),ALLOCATABLE              ::  PTH_DST(:)         !! FLOOD PATHWAY DISTANCE [m]
REAL(KIND=JPRB),ALLOCATABLE              ::  PTH_ELV(:,:)         !! FLOOD PATHWAY ELEVATION [m]
REAL(KIND=JPRB),ALLOCATABLE              ::  PTH_WTH(:,:)         !! FLOOD PATHWAY WIDTH [m]
REAL(KIND=JPRB),ALLOCATABLE              ::  PTH_MAN(:)         !! FLOOD PATHWAY Manning

!================================================
! input matrix (converted from NX:NY*INPN to NSEQMAX*INPN)
INTEGER(KIND=JPIM),ALLOCATABLE           :: INPX(:,:)        !! INPUT GRID XIN
INTEGER(KIND=JPIM),ALLOCATABLE           :: INPY(:,:)        !! INPUT GRID YIN
REAL(KIND=JPRB),ALLOCATABLE              :: INPA(:,:)        !! INPUT AREA

! input matrix Inverse
INTEGER(KIND=JPIM),ALLOCATABLE           :: INPXI(:,:,:)        !! OUTPUT GRID XOUT
INTEGER(KIND=JPIM),ALLOCATABLE           :: INPYI(:,:,:)        !! OUTPUT GRID YOUT
REAL(KIND=JPRB),ALLOCATABLE              :: INPAI(:,:,:)        !! OUTPUT AREA
INTEGER(KIND=JPIM)                       :: INPNI               !! MAX INPUT NUMBER for inverse interpolation

DATA REGIONALL  /1/
DATA REGIONTHIS /1/

END MODULE YOS_CMF_MAP
