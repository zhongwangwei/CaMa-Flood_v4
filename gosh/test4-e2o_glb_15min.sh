#!/bin/sh
#==========================================================
# CaMa-Flood sample go script (4) global 15min simulation with "external runoff data"
# Please follow the doc/Guideline.Md about detail, You need to download external runoff data not included in the package to execute this script.
# -- Long-term  simulations (1980 spinup -> 1980 ~ 2014)
# -- Daily runoff forcing (netCDF) at 15min resolution from "E2O WRR2 WCMEF" is used
# -- input matrix is same as the sample ruoff data (test-15min_nc)
#
# (C) X. Zhou, D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  Feb 2021
#
# Licensed under the Apache License, Version 2.0 (the "License");
#   You may not use this file except in compliance with the License.
#   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is 
#  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the License for the specific language governing permissions and limitations under the License.
#==========================================================

#*** PBS setting when needed
#PBS -q F20
#PBS -l select=1:ncpus=20:mem=10gb
#PBS -j oe
#PBS -m ea
#PBS -V

#================================================
# (0) Basic Setting (for workstation)

#*** 0a. Set CaMa-Flood base directory
BASE=`pwd`/..
# BASE="/home/yamadai/work/CaMa_v424_simd/cmf_v424_pkg"  # setting for PBS in cluster

echo $BASE

#*** 0b. Set dynamic library if needed
export IFORTLIB="/opt/intel/lib:/opt/intel/mkl/lib"
export HDF5LIB="/opt/local/hdf5/lib"
export DYLD_LIBRARY_PATH="${HDF5LIB}:${IFORTLIB}:${DYLD_LIBRARY_PATH}"

#*** 0c. OpenMP thread number
export OMP_NUM_THREADS=16                  # OpenMP cpu num
#export OMP_NUM_THREADS=20                  # OpenMP cpu num


#================================================
# (1) Experiment setting
# -- some non-default options can be modified in NAMELIST section 

#============================
#*** 1a. Experiment directory setting
EXP="test4-e2o_ecmwf-glb_15min"                       # experiment name (output directory name)
RDIR=${BASE}/out/${EXP}                     # directory to run CaMa-Flood
EXE="MAIN_cmf"                              # Execute file name
PROG=${BASE}/src/${EXE}                     # location of Fortran main program
NMLIST="./input_cmf.nam"                    # standard namelist
LOGOUT="./log_CaMa.txt"                     # standard log output


#============================
#*** 1b. Model physics option
DT=3600                                     # base DT (modified in physics loop by LADPSTP)
LADPSTP=".TRUE."                            # .TRUE. for adaptive time step
LPTHOUT=".TRUE."                            # .TRUE. to activate bifurcation flow, mainly for delta simulation

#============================
#*** 1c. simulation time
YSTA=1980                                   # start year ( from YSTA / Jan  1st _ 00:00)
YEND=2014                                   # end   year (until YEND / Dec 31st _ 24:00)
SPINUP=0                                    # [0]: zero-storage start, [1]: from restart file
NSP=3                                       # spinup repeat time


#============================
#*** 1d. spinup setting

#* input restart file
LRESTART="" # see (3) set each year         # TRUE. to use restart initial condition
CRESTSTO="" # see (3) set each year         # input restart FIle
LSTOONLY=".FALSE."                          # .TRUE. for storage only restart (for assimilation)

#* output restart file
CRESTDIR="./"                               # output restart file directory
CVNREST="restart"                           # output restart file prefix
LRESTCDF=".FALSE."                          # .TRUE. to use netCDF restart file
IFRQ_RST="0"                                # output restat frequency.
                                            # [0]: only at last time, [1,2,3,...,24] hourly restart, [30]: monthly restart

#============================
#*** 1e. river map & topography
FMAP="${BASE}/map/glb_15min"                # map directory

#----- for plain binary map input
#** basic topography
LMAPCDF=".FALSE."                           # .TRUE. for netCDF map
CNEXTXY="${FMAP}/nextxy.bin"                # downstream xy (river network map)
CGRAREA="${FMAP}/ctmare.bin"                # unit-catchment area   [m2]
CELEVTN="${FMAP}/elevtn.bin"                # channel top elevation [m]
CNXTDST="${FMAP}/nxtdst.bin"                # downstream distance   [m]
CRIVLEN="${FMAP}/rivlen.bin"                # channel length        [m]
CFLDHGT="${FMAP}/fldhgt.bin"                # floodplain elevation profile (height above 'elevtn') [m]

#** channel parameter
###CRIVWTH=${FMAP}/rivwth.bin"              # channel width [m] (empirical power-low)
CRIVWTH="${FMAP}/rivwth_gwdlr.bin"          # channel width [m] (GWD-LR + filled with empirical)
CRIVHGT="${FMAP}/rivhgt.bin"                # channel depth [m] (empirical power-low)
CRIVMAN="${FMAP}/rivman.bin"                # manning coefficient river (The one in flood plain is a global parameter; set $PMANFLD below.)

#** bifurcation channel info
CPTHOUT="${FMAP}/bifprm.txt"                #   bifurcation channel list

#============================
#*** 1f. forcing setting
#CDIMINFO="${FMAP}/diminfo_test-1deg.txt"    # dimention information file
#CINPMAT=${FMAP}/inpmat_test-1deg.bin        # runoff input matrix for interporlation
CDIMINFO="${FMAP}/diminfo_test-15min_nc.txt" # dimention information file
CINPMAT=${FMAP}/inpmat_test-15min_nc.bin     # runoff input matrix for interporlation

IFRQ_INP="24"                               # input forcing frequency: [1,2,3,...,24] hour
DROFUNIT="86400000"   # [mm/day->m/s]       # runoff unit conversion

#----- for netCDF runoff forcing ###
LINPCDF=".TRUE."                              # true for netCDF runoff
LINTERP=".TRUE."                              # .TRUE. to interporlate with input matrix
LINTERPCDF=".FALSE."                          # .TRUE. to use netCDF input matrix
CROFDIR="${BASE}/inp/E2O/"          # runoff directory
CROFPRE="e2o_ecmwf_wrr2_glob15_day_Runoff_"   # runoff prefix/suffix  
CROFCDF=""     # see (3) set each year        # netCDF runoff file
CVNROF="Runoff"                               # netCDF runoff    variable name
###SYEARIN=""     # see (3) set each year        #   netCDF runoff file, start date
###SMONIN=""      # see (3) set each year
###SDAYIN=""      # see (3) set each year
###SHOURIN=""     # see (3) set each year


#============================
#*** 1h. Output Settings 
IFRQ_OUT=24                                 # output frequency: [1,2,3,...,24] hour

LOUTCDF=".FALSE."                           # .TRUE. netCDF output, .FALSE. plain binary output
COUTDIR="./"                                # output directory 
#CVARSOUT="outflw,storge,fldfrc,maxdph,flddph" # list output variable (comma separated)
CVARSOUT="rivout,rivsto,rivdph,rivvel,fldout,fldsto,flddph,fldfrc,fldare,sfcelv,outflw,storge,pthflw,pthout,maxsto,maxflw,maxdph" # list output variable (comma separated)
COUTTAG=""  # see (3) set each year         #   output tag $(COUTDIR)/$(VARNAME)$(OUTTAG).bin

##### Model Parameters ################
PMANRIV="0.03D0"                            # manning coefficient river
PMANFLD="0.10D0"                            # manning coefficient floodplain
PCADP="0.7"                                 # satety coefficient for CFL condition
PDSTMTH="10000.D0"                          # downstream distance at river mouth [m]



#================================================
# (2) Initial setting

#*** 2a. create running dir 
mkdir -p ${RDIR}
cd ${RDIR}

#*** 2b. for new simulation, remove old files in running directory

if [ ${SPINUP} -eq 0 ]; then
  rm -rf ${RDIR}/????-sp*
  rm -rf ${RDIR}/*.bin
  rm -rf ${RDIR}/*.pth
  rm -rf ${RDIR}/*.vec
  rm -rf ${RDIR}/*.nc
  rm -rf ${RDIR}/*.log
  rm -rf ${RDIR}/*.txt
  rm -rf ${RDIR}/restart*
else
  NSP=0  # restart, no spinup
fi



#================================================
# (3) For each simulation year, modify setting
#--  loop 1-year simulation from $YSTART to $YEND

ISP=1           ## spinup count
IYR=${YSTA}     ## curent year
while [ ${IYR} -le ${YEND} ];
do 
  CYR=`printf %04d ${IYR}`   ## update file name, bugfix in v3.96a

  #*** 3a. modify restart setting
  if [ ${SPINUP} -eq 0 ];then
    LRESTART=".FALSE."                  ## from zero storage
    CRESTSTO=""
  else
    LRESTART=".TRUE."
    CRESTSTO="${CVNREST}${CYR}010100.bin"    ## from restart file
    if [ ${LRESTCDF} = ".TRUE." ]; then
      CRESTSTO="${CVNREST}${CYR}010100.nc"    ## from restart file
    fi
  fi

  #*** 3b. update start-end year
  SYEAR=$IYR
  SMON=1
  SDAY=1
  SHOUR=0

  EYEAR=`expr $SYEAR + 1`
  EMON=1
  EDAY=1
  EHOUR=0

  ln -sf $PROG $EXE

  #*** 3c. update input / output file data
  CSYEAR=`printf %04d ${SYEAR}`
  COUTTAG=${CSYEAR}                  # output file tag

  CROFCDF="${CROFDIR}/${CROFPRE}${CSYEAR}.nc"  # input netCDF runoff file
  SYEARIN=$IYR
  SMONIN=1
  SDAYIN=1
  SHOURIN=0



#================================================
# (4) Create NAMELIST for simulation year
# it is OK to remove optional variables (set to default in CaMa-Flood)

rm -f ${NMLIST}

#*** 0. config
cat >> ${NMLIST} << EOF
&NRUNVER
LADPSTP  = ${LADPSTP}                  ! true: use adaptive time step
LPTHOUT  = ${LPTHOUT}                  ! true: activate bifurcation scheme
LRESTART = ${LRESTART}                 ! true: initial condition from restart file
/
&NDIMTIME
CDIMINFO = "${CDIMINFO}"               ! text file for dimention information
DT       = ${DT}                       ! time step length (sec)
IFRQ_INP = ${IFRQ_INP}                 ! input forcing update frequency (hour)
/
&NPARAM
PMANRIV  = ${PMANRIV}                  ! manning coefficient river
PMANFLD  = ${PMANFLD}                  ! manning coefficient floodplain
PDSTMTH  = ${PDSTMTH}                  ! downstream distance at river mouth [m]
PCADP    = ${PCADP}                    ! CFL coefficient
/
EOF

#*** 1. time
cat >> ${NMLIST} << EOF
&NSIMTIME
SYEAR   = ${SYEAR}                     ! start year
SMON    = ${SMON}                      !  month 
SDAY    = ${SDAY}                      !  day 
SHOUR   = ${SHOUR}                     !  houe
EYEAR   = ${EYEAR}                     ! end year
EMON    = ${EMON}                      !  month 
EDAY    = ${EDAY}                      !  day 
EHOUR   = ${EHOUR}                     !  hour
/
EOF

#*** 2. map
cat >> ${NMLIST} << EOF
&NMAP
LMAPCDF    = ${LMAPCDF}                ! * true for netCDF map input
CNEXTXY    = "${CNEXTXY}"              ! river network nextxy
CGRAREA    = "${CGRAREA}"              ! catchment area
CELEVTN    = "${CELEVTN}"              ! bank top elevation
CNXTDST    = "${CNXTDST}"              ! distance to next outlet
CRIVLEN    = "${CRIVLEN}"              ! river channel length
CFLDHGT    = "${CFLDHGT}"              ! floodplain elevation profile
CRIVWTH    = "${CRIVWTH}"              ! channel width
CRIVHGT    = "${CRIVHGT}"              ! channel depth
CRIVMAN    = "${CRIVMAN}"              ! river manning coefficient
CPTHOUT    = "${CPTHOUT}"              ! bifurcation channel table
/
EOF

#*** 3. restart
cat >> ${NMLIST} << EOF
&NRESTART
CRESTSTO = "${CRESTSTO}"               ! restart file
CRESTDIR = "${CRESTDIR}"               ! restart directory
CVNREST  = "${CVNREST}"                ! restart variable name
LRESTCDF = ${LRESTCDF}                 ! * true for netCDF restart file (double precision)
IFRQ_RST = ${IFRQ_RST}                 ! restart write frequency (1-24: hour, 0:end of run)
/
EOF

#*** 4. forcing
if [ ${LINPCDF} = ".FALSE." ]; then
cat >> ${NMLIST} << EOF
&NFORCE
LINPCDF  = ${LINPCDF}                  ! true for netCDF runoff
LINTERP  = ${LINTERP}                  ! true for runoff interpolation using input matrix
CINPMAT  = "${CINPMAT}"                ! input matrix file name
DROFUNIT = ${DROFUNIT}                 ! runoff unit conversion
CROFDIR  = "${CROFDIR}"                ! runoff             input directory
CROFPRE  = "${CROFPRE}"                ! runoff             input prefix
CROFSUF  = "${CROFSUF}"                ! runoff             input suffix
/
EOF

elif [ ${LINPCDF} = ".TRUE." ]; then
cat >> ${NMLIST} << EOF
&NFORCE
LINPCDF  = ${LINPCDF}                  ! true for netCDF runoff
LINTERP  = ${LINTERP}                  ! true for runoff interpolation using input matrix
CINPMAT  = "${CINPMAT}"                ! input matrix file name
DROFUNIT = ${DROFUNIT}                 ! runoff unit conversion
CROFCDF  = "${CROFCDF}"                ! * netCDF input runoff file name
CVNROF   = "${CVNROF}"                 ! * netCDF input runoff variable name
SYEARIN  = ${SYEARIN}                  ! * netCDF input start year
SMONIN   = ${SMONIN}                   ! * netCDF input start year
SDAYIN   = ${SDAYIN}                   ! * netCDF input start year
SHOURIN  = ${SHOURIN}                  ! * netCDF input start year
/
EOF

fi # (if LINPCDF)

#*** 5. outputs
cat >> ${NMLIST} << EOF
&NOUTPUT
COUTDIR  = "${COUTDIR}"                ! OUTPUT DIRECTORY
CVARSOUT = "${CVARSOUT}"               ! Comma-separated list of output variables to save 
COUTTAG  = "${COUTTAG}"                ! Output Tag Name for each experiment
LOUTCDF  = ${LOUTCDF}                  ! * true for netcdf outptu false for binary
NDLEVEL  = 0                           ! * NETCDF DEFLATION LEVEL 
IFRQ_OUT = ${IFRQ_OUT}                 ! output data write frequency (hour)
/
EOF

#### 6. sea level (optional) 
#cat >> ${NMLIST} << EOF
#&NBOUND
#LSEALEVCDF =  ${LSEALEVCDF}            ! * true : netCDF sea level boundary
#CSEALEVDIR = "${CSEALEVDIR}"           ! Sea level boundary DIRECTORY
#CSEALEVPRE = "${CSEALEVPRE}"           ! Sea level boundary PREFIX
#CSEALEVSUF = "${CSEALEVSUF}"           ! Sea level boundary SUFFIX
#CSEALEVCDF = "${CSEALEVCDF}"           ! * Sea level netCDF file name
#CVNSEALEV  = "${CVNSEALEV}"            ! * Sea Level netCDF variable name
#SYEARSL    = ${SYEARSL}                ! * netCDF sea level start year
#SMONSL     = ${SMONSL}                 ! * netCDF sea level start year
#SDAYSL     = ${SDAYSL}                 ! * netCDF sea level start year
#SHOURSL    = ${SHOURSL}                ! * netCDF sea level start year
#NSTATIONS  = ${NSTATIONS}              ! sea level data points
#CSLMAP     = "${CSLMAP}                ! station to XY conversion table
#IFRQ_SL    = ${IFRQ_SL}                ! sea level boundary update frequency (min)
#/
#EOF

#================================================
# (5) Execute main program

echo "start: ${SYEAR}" `date`  >> log.txt
time ./${EXE}                  >> log.txt 
echo "end:   ${SYEAR}" `date`  >> log.txt

mv ${LOGOUT} log_CaMa-${CYR}.txt


#================================================
# (6) manage spin up

# if curent spinup time $ISP < required spinup time $NSP
#   copy the restart file restart$(IYR+1) to restart${IYR}
#   copy the outputs to directory "${IYR}-sp1"

SPINUP=1
if [ ${IYR} -eq ${YSTA} ];
then
  if [ ${ISP} -le ${NSP} ];
    then
    IYR1=`expr ${IYR} + 1`
    CYR1=`printf %04d ${IYR1}`
    cp -f ${CVNREST}${CYR1}010100.bin ${CVNREST}${CYR}010100.bin-sp${ISP}         2> /dev/null
    mv -f ${CVNREST}${CYR1}010100.bin ${CVNREST}${CYR}010100.bin                  2> /dev/null

    cp -f ${CVNREST}${CYR1}010100.bin.pth ${CVNREST}${CYR}010100.bin.pth-sp${ISP} 2> /dev/null
    mv -f ${CVNREST}${CYR1}010100.bin.pth ${CVNREST}${CYR}010100.bin.pth          2> /dev/null

    cp -f ${CVNREST}${CYR1}010100.nc ${CVNREST}${CYR}010100.nc-sp${ISP}           2> /dev/null
    mv -f ${CVNREST}${CYR1}010100.nc ${CVNREST}${CYR}010100.nc                    2> /dev/null

    mkdir -p ${CYR}-sp${ISP}
    mv -f ./${CVNREST}${CYR}010100.bin-sp${ISP}      ${CYR}-sp${ISP}  2> /dev/null
    mv -f ./${CVNREST}${CYR}010100.bin.pth-sp${ISP}  ${CYR}-sp${ISP}  2> /dev/null
    mv -f ./${CVNREST}${CYR}010100.nc-sp${ISP}       ${CYR}-sp${ISP}  2> /dev/null
    mv -f ./*${CYR}.bin                              ${CYR}-sp${ISP}  2> /dev/null
    mv -f ./*${CYR}.pth                              ${CYR}-sp${ISP}  2> /dev/null
    mv -f ./o_*${CYR}.nc                             ${CYR}-sp${ISP}  2> /dev/null
    mv -f ./*${CYR}.log                              ${CYR}-sp${ISP}  2> /dev/null
    mv -f ./log_CaMa-${CYR}.txt                      ${CYR}-sp${ISP}  2> /dev/null

    ISP=`expr ${ISP} + 1`
  else
    ISP=0
    IYR=`expr ${IYR} + 1`
  fi
else
  IYR=`expr ${IYR} + 1`
fi

#================================================
# (7) End of each year loop. Back to (3)

done # loop to next year simulation

exit 0
