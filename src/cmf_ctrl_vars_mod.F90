MODULE CMF_CTRL_VARS_MOD
!==========================================================
!* PURPOSE: Manage prognostic/diagnostic variables in CaMa-Flood
!
!* CONTAINS:
! -- CMF_PROG_INIT      : Initialize Prognostic variables (include restart data handling)
! -- CMF_DIAG_INIT      : Initialize Diagnostic variables
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
USE PARKIND1,                ONLY: JPIM, JPRM, JPRB, JPRD
USE YOS_CMF_INPUT,           ONLY: LOGNAM, LPTHOUT, LDAMOUT, LLEVEE, LWEVAP, LOUTINS
IMPLICIT NONE
CONTAINS 
!####################################################################
! -- CMF_PROG_INIT      : Initialize Prognostic variables (include restart data handling)
! -- CMF_DIAG_INIT      : Initialize Diagnostic variables
!
!####################################################################
SUBROUTINE CMF_PROG_INIT
USE YOS_CMF_MAP,             ONLY: NSEQALL, NPTHOUT, NPTHLEV
USE YOS_CMF_PROG,            ONLY: D2RUNOFF,     D2ROFSUB,     &
                                 & P2RIVSTO,     P2FLDSTO,     D2RIVOUT,     D2FLDOUT,     &
                                 & D2RIVOUT_PRE, D2FLDOUT_PRE, D2RIVDPH_PRE, D2FLDSTO_PRE, &
                                 & D1PTHFLW,     D1PTHFLW_PRE, P2GDWSTO,     D2GDWRTN,     &
                                 & P2DAMSTO,     P2DAMINF,     P2LEVSTO ,    D2WEVAP,      &   !! optional
                                 & D2DAMMY,      D2COPY
IMPLICIT NONE
!================================================
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!---------------------!"

WRITE(LOGNAM,*) "CMF::PROG_INIT: prognostic variable initialization"

!*** 1. ALLOCATE 
! runoff input
ALLOCATE( D2RUNOFF(NSEQALL,1)     )
ALLOCATE( D2ROFSUB(NSEQALL,1)     )

! river+floodplain storage
ALLOCATE( P2RIVSTO(NSEQALL,1)     )
ALLOCATE( P2FLDSTO(NSEQALL,1)     )

! discharge calculation
ALLOCATE( D2RIVOUT(NSEQALL,1)     )
ALLOCATE( D2FLDOUT(NSEQALL,1)     )
ALLOCATE( D2RIVOUT_PRE(NSEQALL,1)     )
ALLOCATE( D2FLDOUT_PRE(NSEQALL,1)     )
ALLOCATE( D2RIVDPH_PRE(NSEQALL,1)     )
ALLOCATE( D2FLDSTO_PRE(NSEQALL,1)     )

D2RUNOFF(:,:)=0._JPRB
D2ROFSUB(:,:)=0._JPRB

P2RIVSTO(:,:)=0._JPRD
P2FLDSTO(:,:)=0._JPRD

D2RIVOUT(:,:)=0._JPRB
D2FLDOUT(:,:)=0._JPRB
D2RIVOUT_PRE(:,:)=0._JPRB
D2FLDOUT_PRE(:,:)=0._JPRB
D2RIVDPH_PRE(:,:)=0._JPRB
D2FLDSTO_PRE(:,:)=0._JPRB

!! prognostics for bifurcation (if no bifurcation, NPTHOUT=1, NPTHLEV=1)
ALLOCATE( D1PTHFLW(NPTHOUT,NPTHLEV)     )
ALLOCATE( D1PTHFLW_PRE(NPTHOUT,NPTHLEV) )
D1PTHFLW(:,:)=0._JPRB
D1PTHFLW_PRE(:,:)=0._JPRB

! reservoir operation
IF( LDAMOUT )THEN
  ALLOCATE( P2DAMSTO(NSEQALL,1)     )
  ALLOCATE( P2DAMINF(NSEQALL,1)     )
  P2DAMSTO(:,:)=0._JPRD
  P2DAMINF(:,:)=0._JPRD
ENDIF

IF( LLEVEE ) THEN  !! additional prognostics for LLEVEE
  ALLOCATE( P2LEVSTO(NSEQALL,1)     )
  P2LEVSTO(:,:)=0._JPRD
ENDIF

!! Used in ECMWF
IF( LWEVAP ) THEN  !! additional prognostics for LLEVEE
  ALLOCATE( D2WEVAP(NSEQALL,1)     )
  D2WEVAP(:,:)=0._JPRB
ENDIF

!! keep these variables even when LGDWDLY is not used. (for simplifying runoff calculation)
ALLOCATE( P2GDWSTO(NSEQALL,1)     )
ALLOCATE( D2GDWRTN(NSEQALL,1)     )
P2GDWSTO(:,:)=0._JPRD
D2GDWRTN(:,:)=0._JPRB

!! dammy variable for data handling
ALLOCATE( D2DAMMY(NSEQALL,1)) !! Float64/32 switch (Dammy for unused var)
ALLOCATE( D2COPY(NSEQALL,1))  !! Float64/32 switch (Dammy for output)
D2DAMMY(:,:)=0._JPRB
D2COPY(:,:) =0._JPRB

!============================
!***  2. set initial water surface elevation to sea surface level
WRITE(LOGNAM,*) 'PROG_INIT: fill channels below downstream boundary'
CALL STORAGE_SEA_SURFACE


WRITE(LOGNAM,*) "CMF::PROG_INIT: end"

CONTAINS
!==========================================================
!+ STORAGE_SEA_SURFACE: set initial storage, assuming water surface not lower than downstream sea surface elevation
!+
!+
! ==================================================
SUBROUTINE STORAGE_SEA_SURFACE
! set initial storage, assuming water surface not lower than downstream sea surface elevation
USE YOS_CMF_MAP,  ONLY: NSEQRIV,  NSEQALL,  I1NEXT
USE YOS_CMF_MAP,  ONLY: D2DWNELV, D2RIVELV,D2RIVHGT,D2RIVWTH,D2RIVLEN,D2RIVSTOMAX
IMPLICIT NONE
! local variables
INTEGER(KIND=JPIM)   :: ISEQ, JSEQ
!
REAL(KIND=JPRB),SAVE :: DSEAELV, DDPH
!$OMP THREADPRIVATE    (DSEAELV, DDPH)
!!=================
! For River Mouth Grid
!$OMP PARALLEL DO SIMD
DO ISEQ=NSEQRIV+1,NSEQALL
  DSEAELV=D2DWNELV(ISEQ,1) !! downstream boundary elevation

  !! set initial water level to sea level if river bed is lower than sea level
  DDPH=MAX( DSEAELV-D2RIVELV(ISEQ,1),0._JPRB )
  DDPH=MIN( DDPH,D2RIVHGT(ISEQ,1) )
  P2RIVSTO(ISEQ,1)=DDPH*D2RIVLEN(ISEQ,1)*D2RIVWTH(ISEQ,1)
  P2RIVSTO(ISEQ,1)=MIN( P2RIVSTO(ISEQ,1),D2RIVSTOMAX(ISEQ,1)*1._JPRD )
  D2RIVDPH_PRE(ISEQ,1)=DDPH
END DO
!$OMP END PARALLEL DO SIMD

!! For Usual River Grid (from downstream to upstream). OMP cannot be applied
DO ISEQ=NSEQRIV,1, -1
  JSEQ=I1NEXT(ISEQ)
  DSEAELV=D2RIVELV(JSEQ,1)+D2RIVDPH_PRE(JSEQ,1)

  !! set initial water level to sea level if river bed is lower than sea level
  DDPH=MAX( DSEAELV-D2RIVELV(ISEQ,1),0._JPRB )
  DDPH=MIN( DDPH,D2RIVHGT(ISEQ,1) )

  P2RIVSTO(ISEQ,1)=DDPH*D2RIVLEN(ISEQ,1)*D2RIVWTH(ISEQ,1)
  P2RIVSTO(ISEQ,1)=MIN(  P2RIVSTO(ISEQ,1),D2RIVSTOMAX(ISEQ,1)*1._JPRD )
  D2RIVDPH_PRE(ISEQ,1)=DDPH
END DO

  
! old version before v4.02 (too slow)
!DO ISEQ=1, NSEQALL
!  JSEQ=ISEQ
!  DO WHILE( I1NEXT(JSEQ)>0 )
!    KSEQ=JSEQ
!    JSEQ=I1NEXT(KSEQ)
!  END DO
!
!  DSEAELV=D2DWNELV(JSEQ,1) !! downstream boundary elevation
!  !! set initial water level to sea level if river bed is lower than sea level
!  DDPH=MAX( DSEAELV-D2RIVELV(ISEQ,1),0._JPRB )
!  DDPH=MIN( DDPH,D2RIVHGT(ISEQ,1) )
!  P2RIVSTO(ISEQ,1)=DDPH*D2RIVLEN(ISEQ,1)*D2RIVWTH(ISEQ,1)
!END DO
    
END SUBROUTINE STORAGE_SEA_SURFACE
! ==================================================

END SUBROUTINE CMF_PROG_INIT
!####################################################################






!####################################################################
SUBROUTINE CMF_DIAG_INIT

USE YOS_CMF_MAP,        ONLY: NSEQALL,NPTHOUT,NPTHLEV
USE YOS_CMF_DIAG,       ONLY: D2RIVINF, D2RIVDPH, D2RIVVEL, D2FLDINF, D2FLDDPH, D2FLDFRC, D2FLDARE, &
                            & D2PTHOUT, D2PTHINF, D2SFCELV, D2OUTFLW, D2STORGE, D2OUTINS, D2LEVDPH, &
                            & D1PTHFLWSUM,   D2WEVAPEX, P2STOOUT, P2RIVINF, P2FLDINF, P2PTHOUT, D2RATE,  &
                            & D2SFCELV_PRE,  D2DWNELV_PRE,  D2FLDDPH_PRE

USE YOS_CMF_DIAG,       ONLY: D2RIVOUT_oAVG, D2FLDOUT_oAVG, D2OUTFLW_oAVG, D2RIVVEL_oAVG, D2PTHOUT_oAVG, &
                            & D2GDWRTN_oAVG, D2RUNOFF_oAVG, D2ROFSUB_oAVG, D1PTHFLW_oAVG, D2WEVAPEX_oAVG,&
                            & D2DAMINF_oAVG, D2STORGE_oMAX, D2OUTFLW_oMAX, D2RIVDPH_oMAX, NADD_out
USE YOS_CMF_DIAG,       ONLY: D2RIVOUT_aAVG, D2FLDOUT_aAVG, D2OUTFLW_aAVG, D2RIVVEL_aAVG, D2PTHOUT_aAVG, &
                            & D2GDWRTN_aAVG, D2RUNOFF_aAVG, D2ROFSUB_aAVG, D1PTHFLW_aAVG, D2WEVAPEX_aAVG,&
                            & D2DAMINF_aAVG, D2STORGE_aMAX, D2OUTFLW_aMAX, D2RIVDPH_aMAX, NADD_adp,      &
                            & D1PTHFLWSUM_aAVG
IMPLICIT NONE
!================================================
WRITE(LOGNAM,*) ""
WRITE(LOGNAM,*) "!---------------------!"

WRITE(LOGNAM,*) "CMF::DIAG_INIT: initialize diagnostic variables"

!*** 1. snapshot 2D diagnostics
ALLOCATE(D2RIVINF(NSEQALL,1))
ALLOCATE(D2RIVDPH(NSEQALL,1))
ALLOCATE(D2RIVVEL(NSEQALL,1))
ALLOCATE(D2FLDINF(NSEQALL,1))
ALLOCATE(D2FLDDPH(NSEQALL,1))
ALLOCATE(D2FLDFRC(NSEQALL,1))
ALLOCATE(D2FLDARE(NSEQALL,1))
ALLOCATE(D2PTHOUT(NSEQALL,1))
ALLOCATE(D2PTHINF(NSEQALL,1))
ALLOCATE(D2SFCELV(NSEQALL,1))
ALLOCATE(D2OUTFLW(NSEQALL,1))
ALLOCATE(D2STORGE(NSEQALL,1))
D2RIVINF(:,:)=0._JPRB
D2RIVDPH(:,:)=0._JPRB
D2RIVVEL(:,:)=0._JPRB
D2FLDINF(:,:)=0._JPRB
D2FLDDPH(:,:)=0._JPRB
D2FLDFRC(:,:)=0._JPRB
D2FLDARE(:,:)=0._JPRB
D2PTHOUT(:,:)=0._JPRB
D2PTHINF(:,:)=0._JPRB
D2SFCELV(:,:)=0._JPRB
D2OUTFLW(:,:)=0._JPRB
D2STORGE(:,:)=0._JPRB

ALLOCATE(D1PTHFLWSUM(NPTHOUT))
D1PTHFLWSUM(:)=0._JPRB

IF ( LLEVEE  )THEN
  ALLOCATE(D2LEVDPH(NSEQALL,1))
  D2LEVDPH(:,:)=0._JPRB
ENDIF
IF ( LWEVAP  )THEN
  ALLOCATE(D2WEVAPEX(NSEQALL,1))
  D2WEVAPEX(:,:)=0._JPRB
ENDIF
IF ( LOUTINS  )THEN
  ALLOCATE(D2OUTINS (NSEQALL,1))
  D2OUTINS (:,:)=0._JPRB
ENDIF


!! temporally variables in subroutines
ALLOCATE(D2SFCELV_PRE(NSEQALL,1))
ALLOCATE(D2DWNELV_PRE(NSEQALL,1))
ALLOCATE(D2FLDDPH_PRE(NSEQALL,1))
D2SFCELV_PRE(:,:)=0._JPRB
D2DWNELV_PRE(:,:)=0._JPRB
D2FLDDPH_PRE(:,:)=0._JPRB

ALLOCATE(P2STOOUT(NSEQALL,1))
ALLOCATE(P2RIVINF(NSEQALL,1))
ALLOCATE(P2FLDINF(NSEQALL,1))
ALLOCATE(P2PTHOUT(NSEQALL,1))
ALLOCATE(D2RATE(NSEQALL,1))
P2STOOUT(:,:)=0._JPRD
P2RIVINF(:,:)=0._JPRD
P2FLDINF(:,:)=0._JPRD
P2PTHOUT(:,:)=0._JPRD
D2RATE(:,:)  =1._JPRB


!============================
!*** 2a. time-average 2D diagnostics for adaptive time step
NADD_adp=0

ALLOCATE(D2RIVOUT_aAVG(NSEQALL,1))
ALLOCATE(D2FLDOUT_aAVG(NSEQALL,1))
ALLOCATE(D2OUTFLW_aAVG(NSEQALL,1))
ALLOCATE(D2RIVVEL_aAVG(NSEQALL,1))
ALLOCATE(D2PTHOUT_aAVG(NSEQALL,1))
ALLOCATE(D2GDWRTN_aAVG(NSEQALL,1))
ALLOCATE(D2RUNOFF_aAVG(NSEQALL,1))
ALLOCATE(D2ROFSUB_aAVG(NSEQALL,1))
D2RIVOUT_aAVG(:,:)=0._JPRB
D2FLDOUT_aAVG(:,:)=0._JPRB
D2OUTFLW_aAVG(:,:)=0._JPRB
D2RIVVEL_aAVG(:,:)=0._JPRB
D2PTHOUT_aAVG(:,:)=0._JPRB
D2GDWRTN_aAVG(:,:)=0._JPRB
D2RUNOFF_aAVG(:,:)=0._JPRB
D2ROFSUB_aAVG(:,:)=0._JPRB

IF ( LDAMOUT ) THEN
  ALLOCATE(D2DAMINF_aAVG(NSEQALL,1))
  D2DAMINF_aAVG(:,:)=0._JPRB
ENDIF
IF ( LWEVAP ) THEN
  ALLOCATE(D2WEVAPEX_aAVG(NSEQALL,1))
  D2WEVAPEX_aAVG(:,:)=0._JPRB
ENDIF

!*** 2b time-average 1D Diagnostics (bifurcation channel) for adaptive time step
ALLOCATE(D1PTHFLW_aAVG(NPTHOUT,NPTHLEV))
ALLOCATE(D1PTHFLWSUM_aAVG(NPTHOUT))
D1PTHFLW_aAVG(:,:)  = 0._JPRB 
D1PTHFLWSUM_aAVG(:) = 0._JPRB 

!*** 2c. Maximum 2D Diagnostics 

ALLOCATE(D2STORGE_aMAX(NSEQALL,1))
ALLOCATE(D2OUTFLW_aMAX(NSEQALL,1))
ALLOCATE(D2RIVDPH_aMAX(NSEQALL,1))
D2STORGE_aMAX(:,:)=0._JPRB
D2OUTFLW_aMAX(:,:)=0._JPRB
D2RIVDPH_aMAX(:,:)=0._JPRB

!============
!*** 3a. time-average 2D diagnostics for output
NADD_out=0

ALLOCATE(D2RIVOUT_oAVG(NSEQALL,1))
ALLOCATE(D2FLDOUT_oAVG(NSEQALL,1))
ALLOCATE(D2OUTFLW_oAVG(NSEQALL,1))
ALLOCATE(D2RIVVEL_oAVG(NSEQALL,1))
ALLOCATE(D2PTHOUT_oAVG(NSEQALL,1))
ALLOCATE(D2GDWRTN_oAVG(NSEQALL,1))
ALLOCATE(D2RUNOFF_oAVG(NSEQALL,1))
ALLOCATE(D2ROFSUB_oAVG(NSEQALL,1))
D2RIVOUT_oAVG(:,:)=0._JPRB
D2FLDOUT_oAVG(:,:)=0._JPRB
D2OUTFLW_oAVG(:,:)=0._JPRB
D2RIVVEL_oAVG(:,:)=0._JPRB
D2PTHOUT_oAVG(:,:)=0._JPRB
D2GDWRTN_oAVG(:,:)=0._JPRB
D2RUNOFF_oAVG(:,:)=0._JPRB
D2ROFSUB_oAVG(:,:)=0._JPRB

IF ( LDAMOUT ) THEN
  ALLOCATE(D2DAMINF_oAVG(NSEQALL,1))
  D2DAMINF_oAVG(:,:)=0._JPRB
ENDIF
IF ( LWEVAP ) THEN
  ALLOCATE(D2WEVAPEX_oAVG(NSEQALL,1))
  D2WEVAPEX_oAVG(:,:)=0._JPRB
ENDIF

!*** ab time-average 1D Diagnostics (bifurcation channel)
ALLOCATE(D1PTHFLW_oAVG(NPTHOUT,NPTHLEV))
D1PTHFLW_oAVG(:,:) = 0._JPRB 

!*** 3c. Maximum 2D Diagnostics 

ALLOCATE(D2STORGE_oMAX(NSEQALL,1))
ALLOCATE(D2OUTFLW_oMAX(NSEQALL,1))
ALLOCATE(D2RIVDPH_oMAX(NSEQALL,1))
D2STORGE_oMAX(:,:)=0._JPRB
D2OUTFLW_oMAX(:,:)=0._JPRB
D2RIVDPH_oMAX(:,:)=0._JPRB

WRITE(LOGNAM,*) "CMF::DIAG_INIT: end"

END SUBROUTINE CMF_DIAG_INIT
!####################################################################

END MODULE CMF_CTRL_VARS_MOD
