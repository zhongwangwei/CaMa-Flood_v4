MODULE CMF_CTRL_MPI_MOD
!! contains nothing is UseMPI is not defined
#ifdef UseMPI
!==========================================================
!* PURPOSE: modules related to MPI usage 
!
!* CONTAINS:
! -- CMF_MPI_INIT     : MPI Initialization
! -- CMF_MPI_END      : MPI Finalization
!
! (C) D.Yamazaki (U-Tokyo)  Oct 2021
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
!** shared variables in module
USE PARKIND1,                ONLY: JPIM, JPRB, JPRM
USE YOS_CMF_INPUT,           ONLY: LOGNAM
USE YOS_CMF_MAP,             ONLY: REGIONALL, REGIONTHIS
!$ USE OMP_LIB
IMPLICIT NONE
!** local variables
SAVE
!** MPI setting
INTEGER(KIND=JPIM)              :: ierr, Nproc, Nid
INTEGER(KIND=JPIM)              :: iOMP, nOMP
!==========================================================
CONTAINS
!####################################################################
! -- CMF_DRV_INPUT    : Set namelist & logfile
! -- CMF_DRV_INIT     : Initialize        CaMa-Flood
! -- CMF_DRV_END      : Finalize          CaMa-Flood
!
!####################################################################
SUBROUTINE CMF_MPI_INIT
USE MPI
IMPLICIT NONE
!================================================
!*** 0. MPI specific setting
REGIONTHIS=1
CALL MPI_Init(ierr)

CALL MPI_Comm_size(MPI_COMM_WORLD, Nproc, ierr)
CALL MPI_Comm_rank(MPI_COMM_WORLD, Nid, ierr)

REGIONALL =Nproc
REGIONTHIS=Nid+1

! For BUGFIX: Check MPI  / OpenMPI is working or not.
! Write to standard output (log file is not opened yet)
#ifdef _OPENMP
nOMP = omp_get_max_threads();
!$OMP PARALLEL DO
DO iOMP=1, nOMP
  print *, 'MPI: ', REGIONTHIS, REGIONALL, ' OMP: ', omp_get_thread_num(), nOMP
END DO
!$OMP END PARALLEL DO
#endif

END SUBROUTINE CMF_MPI_INIT
!####################################################################



!####################################################################
SUBROUTINE CMF_MPI_END
USE MPI
IMPLICIT NONE
INTEGER(KIND=JPIM)              :: ierr
!================================================
CALL MPI_Finalize(ierr)
END SUBROUTINE CMF_MPI_END
!####################################################################


!####################################################################
SUBROUTINE MPI_REDUCE_R2MAP(R2MAP)
USE MPI
USE YOS_CMF_INPUT,           ONLY: RMIS, NX,NY
IMPLICIT NONE
!* input/output
REAL(KIND=JPRM),INTENT(INOUT)   :: R2MAP(NX,NY)
!* local variable
REAL(KIND=JPRM)                 :: R2TMP(NX,NY)
!================================================
! gather to master node
  R2TMP(:,:)=RMIS
  CALL MPI_Reduce(R2MAP,R2TMP,NX*NY,MPI_REAL4,MPI_MIN,0,mpi_comm_world,ierr)
  R2MAP(:,:)=R2TMP(:,:)
END SUBROUTINE MPI_REDUCE_R2MAP
!####################################################################




!####################################################################
SUBROUTINE MPI_REDUCE_R1PTH(R1PTH)
USE MPI
USE YOS_CMF_INPUT,           ONLY: RMIS
USE YOS_CMF_MAP,             ONLY: NPTHOUT, NPTHLEV
IMPLICIT NONE
!* input/output
REAL(KIND=JPRM),INTENT(INOUT)   :: R1PTH(NPTHOUT,NPTHLEV)
!* local variable
REAL(KIND=JPRM)                 :: R1PTMP(NPTHOUT,NPTHLEV)
!================================================
! gather to master node
  R1PTMP(:,:)=RMIS
  CALL MPI_Reduce(R1PTH,R1PTMP,NPTHOUT*NPTHLEV,MPI_REAL4,MPI_MIN,0,mpi_comm_world,ierr)
  R1PTH(:,:)=R1PTMP(:,:)
END SUBROUTINE MPI_REDUCE_R1PTH
!####################################################################


!####################################################################
SUBROUTINE MPI_ADPSTP(DT_MIN)
USE MPI
USE YOS_CMF_INPUT,           ONLY: LOGNAM
IMPLICIT NONE
!* input/output
REAL(KIND=JPRB),INTENT(INOUT)   :: DT_MIN
!* local variable
REAL(KIND=JPRB)                 :: DT_LOC
!================================================
!*** MPI: use same DT in all node
DT_LOC=DT_MIN

CALL MPI_AllReduce(DT_LOC, DT_MIN, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD,ierr)
WRITE(LOGNAM,'(A,2F10.2)') "ADPSTP (MPI_AllReduce): DT_LOC->DTMIN", DT_LOC, DT_MIN

END SUBROUTINE MPI_ADPSTP
!####################################################################

#endif
END MODULE CMF_CTRL_MPI_MOD