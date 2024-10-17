import numpy as np
from netCDF4 import Dataset

class CMF_CTRL_OUTPUT_MOD:
    def __init__(self):
        self.COUTDIR = "./"
        self.CVARSOUT = "outflw,storge,rivdph"
        self.COUTTAG = "_cmf"
        self.LOUTCDF = False
        self.NDLEVEL = 0
        self.LOUTVEC = False
        self.IFRQ_OUT = 24
        self.LOUTTXT = False
        self.CGAUTXT = "None"
        self.NVARS = 100
        self.NVARSOUT = 0
        self.IRECOUT = 0
        self.VAROUT = []

    def CMF_OUTPUT_NMLIST(self, CSETFILE):
        print("")
        print("!---------------------!")

        with open(CSETFILE, 'r') as file:
            print("CMF::OUTPUT_NMLIST: namelist OPEN in unit:", CSETFILE)

            self.COUTDIR = "./"
            self.CVARSOUT = "outflw,storge,rivdph"
            self.COUTTAG = "_cmf"
            self.LOUTCDF = False
            self.NDLEVEL = 0
            self.LOUTVEC = False
            self.IFRQ_OUT = 24
            self.LOUTTXT = False
            self.CGAUTXT = "None"

            file.seek(0)
            for line in file:
                if "NOUTPUT" in line:
                    break

            print("=== NAMELIST, NOUTPUT ===")
            print("COUTDIR:  ", self.COUTDIR)
            print("CVARSOUT: ", self.CVARSOUT)
            print("COUTTAG:  ", self.COUTTAG)
            print("LOUTCDF:  ", self.LOUTCDF)
            if self.LOUTCDF:
                print("NDLEVEL:  ", self.NDLEVEL)
            if self.LOUTVEC:
                print("LOUTVEC:  ", self.LOUTVEC)
            print("IFRQ_OUT  ", self.IFRQ_OUT)
            print("IFRQ_OUT  ", self.LOUTTXT)
            print("CGAUTXRT  ", self.CGAUTXT)

        print("CMF::OUTPUT_NMLIST: end")

    def CMF_OUTPUT_INIT(self, ISYYYY, ISMM, ISDD, ISHOUR, ISMIN, NSEQMAX, NPTHOUT, NPTHLEV, REGIONTHIS, NX, NY):
        print("")
        print("!---------------------!")

        print("CMF::OUTPUT_INIT: check output variables")
        self.NVARSOUT = 0
        J0 = 0
        CVNAMES = []

        for J in range(len(self.CVARSOUT)):
            if J > J0 and self.CVARSOUT[J] == ',':
                CTMP = self.CVARSOUT[J0:J].strip()
                if len(CTMP) > 0:
                    self.NVARSOUT += 1
                    CVNAMES.append(CTMP)
                J0 = J + 1

        if J0 <= len(self.CVARSOUT):
            J = len(self.CVARSOUT)
            CTMP = self.CVARSOUT[J0:J].strip()
            if len(CTMP) > 0:
                self.NVARSOUT += 1
                CVNAMES.append(CTMP)

        if self.NVARSOUT == 0:
            print("CMF::OUTPUT_INIT: No output files will be produced!")
            return

        self.VAROUT = [self.TVAROUT() for _ in range(self.NVARSOUT)]
        CTIME = f'seconds since {ISYYYY}-{ISMM:02d}-{ISDD:02d} {ISHOUR:02d}:{ISMIN:02d}'

        for JF in range(self.NVARSOUT):
            print("Creating output for variable:", CVNAMES[JF])
            if CVNAMES[JF] == 'rivout':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'river discharge'
                self.VAROUT[JF].CVUNITS = 'm3/s'
            elif CVNAMES[JF] == 'rivsto':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'river storage'
                self.VAROUT[JF].CVUNITS = 'm3'
            elif CVNAMES[JF] == 'rivdph':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'river depth'
                self.VAROUT[JF].CVUNITS = 'm'
            elif CVNAMES[JF] == 'rivvel':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'river velocity'
                self.VAROUT[JF].CVUNITS = 'm/s'
            elif CVNAMES[JF] == 'fldout':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'floodplain discharge'
                self.VAROUT[JF].CVUNITS = 'm3/s'
            elif CVNAMES[JF] == 'fldsto':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'floodplain storage'
                self.VAROUT[JF].CVUNITS = 'm3'
            elif CVNAMES[JF] == 'flddph':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'floodplain depth'
                self.VAROUT[JF].CVUNITS = 'm'
            elif CVNAMES[JF] == 'fldfrc':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'flooded fraction'
                self.VAROUT[JF].CVUNITS = '0-1'
            elif CVNAMES[JF] == 'fldare':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'flooded area'
                self.VAROUT[JF].CVUNITS = 'm2'
            elif CVNAMES[JF] == 'sfcelv':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'water surface elevation'
                self.VAROUT[JF].CVUNITS = 'm'
            elif CVNAMES[JF] == 'totout':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'discharge (river+floodplain)'
                self.VAROUT[JF].CVUNITS = 'm3/s'
            elif CVNAMES[JF] == 'outflw':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'discharge (river+floodplain)'
                self.VAROUT[JF].CVUNITS = 'm3/s'
            elif CVNAMES[JF] == 'totsto':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'total storage (river+floodplain)'
                self.VAROUT[JF].CVUNITS = 'm3'
            elif CVNAMES[JF] == 'storge':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'total storage (river+floodplain)'
                self.VAROUT[JF].CVUNITS = 'm3'
            elif CVNAMES[JF] == 'pthflw':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'bifurcation channel discharge'
                self.VAROUT[JF].CVUNITS = 'm3/s'
            elif CVNAMES[JF] == 'pthout':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'net bifurcation discharge'
                self.VAROUT[JF].CVUNITS = 'm3/s'
            elif CVNAMES[JF] == 'maxsto':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'daily maximum storage'
                self.VAROUT[JF].CVUNITS = 'm3'
            elif CVNAMES[JF] == 'maxflw':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'daily maximum discharge'
                self.VAROUT[JF].CVUNITS = 'm3/s'
            elif CVNAMES[JF] == 'maxdph':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'daily maximum river depth'
                self.VAROUT[JF].CVUNITS = 'm'
            elif CVNAMES[JF] == 'runoff':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'Surface runoff'
                self.VAROUT[JF].CVUNITS = 'm3/s'
            elif CVNAMES[JF] == 'runoffsub':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'sub-surface runoff'
                self.VAROUT[JF].CVUNITS = 'm3/s'
            elif CVNAMES[JF] == 'damsto':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'reservoir storage'
                self.VAROUT[JF].CVUNITS = 'm3'
            elif CVNAMES[JF] == 'daminf':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'reservoir inflow'
                self.VAROUT[JF].CVUNITS = 'm3/s'
            elif CVNAMES[JF] == 'levsto':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'protected area storage'
                self.VAROUT[JF].CVUNITS = 'm3'
            elif CVNAMES[JF] == 'levdph':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'protected area depth'
                self.VAROUT[JF].CVUNITS = 'm'
            elif CVNAMES[JF] == 'gdwsto':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'ground water storage'
                self.VAROUT[JF].CVUNITS = 'm3'
            elif CVNAMES[JF] == 'gwsto':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'ground water storage'
                self.VAROUT[JF].CVUNITS = 'm3'
            elif CVNAMES[JF] == 'gwout':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'ground water discharge'
                self.VAROUT[JF].CVUNITS = 'm3/s'
            elif CVNAMES[JF] == 'wevap':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'water evaporation'
                self.VAROUT[JF].CVUNITS = 'm3/s'
            elif CVNAMES[JF] == 'outins':
                self.VAROUT[JF].CVNAME = CVNAMES[JF]
                self.VAROUT[JF].CVLNAME = 'instantaneous discharge'
                self.VAROUT[JF].CVUNITS = 'm3/s'
            else:
                print(CVNAMES[JF], ' Not defined in CMF_CREATE_OUTCDF_MOD')

            self.VAROUT[JF].BINID = self.INQUIRE_FID()

            if self.LOUTCDF:
                if REGIONTHIS == 1:
                    self.CREATE_OUTCDF(JF, CTIME, NX, NY)
            else:
                self.CREATE_OUTBIN(JF, NSEQMAX, NPTHOUT, NPTHLEV, REGIONTHIS, NX, NY)

        self.IRECOUT = 0

    def CREATE_OUTBIN(self, JF, NSEQMAX, NPTHOUT, NPTHLEV, REGIONTHIS, NX, NY):
        if self.VAROUT[JF].CVNAME == 'pthflw':
            if REGIONTHIS == 1:
                self.VAROUT[JF].CFILE = f"{self.COUTDIR}{self.VAROUT[JF].CVNAME}{self.COUTTAG}.pth"
                self.VAROUT[JF].BINID = open(self.VAROUT[JF].CFILE, 'wb')
                print("output file opened in unit: ", self.VAROUT[JF].CFILE)
        elif self.LOUTVEC:
            self.VAROUT[JF].CFILE = f"{self.COUTDIR}{self.VAROUT[JF].CVNAME}{self.COUTTAG}.vec"
            self.VAROUT[JF].BINID = open(self.VAROUT[JF].CFILE, 'wb')
            print("output file opened in unit: ", self.VAROUT[JF].CFILE)
        else:
            if REGIONTHIS == 1:
                self.VAROUT[JF].CFILE = f"{self.COUTDIR}{self.VAROUT[JF].CVNAME}{self.COUTTAG}.bin"
                self.VAROUT[JF].BINID = open(self.VAROUT[JF].CFILE, 'wb')
                print("output file opened in unit: ", self.VAROUT[JF].CFILE)

    def CREATE_OUTCDF(self, JF, CTIME, NX, NY):
        self.VAROUT[JF].IRECNC = 1
        self.VAROUT[JF].CFILE = f"{self.COUTDIR}o_{self.VAROUT[JF].CVNAME}{self.COUTTAG}.nc"
        self.VAROUT[JF].NCID = Dataset(self.VAROUT[JF].CFILE, 'w', format='NETCDF4')
        self.VAROUT[JF].NCID.createDimension('time', None)
        self.VAROUT[JF].NCID.createDimension('lat', NY)
        self.VAROUT[JF].NCID.createDimension('lon', NX)
        lat = self.VAROUT[JF].NCID.createVariable('lat', 'f4', ('lat',))
        lon = self.VAROUT[JF].NCID.createVariable('lon', 'f4', ('lon',))
        time = self.VAROUT[JF].NCID.createVariable('time', 'f8', ('time',))
        var = self.VAROUT[JF].NCID.createVariable(self.VAROUT[JF].CVNAME, 'f4', ('lon', 'lat', 'time',), zlib=True, complevel=self.NDLEVEL)
        lat.long_name = 'latitude'
        lat.units = 'degrees_north'
        lon.long_name = 'longitude'
        lon.units = 'degrees_east'
        time.long_name = 'time'
        time.units = CTIME
        var.long_name = self.VAROUT[JF].CVLNAME
        var.units = self.VAROUT[JF].CVUNITS
        var._FillValue = 1.0e20
        print('CFILE: ', self.VAROUT[JF].CFILE, ' CVAR:', self.VAROUT[JF].CVNAME, ' CLNAME: ', self.VAROUT[JF].CVLNAME, ' CUNITS: ', self.VAROUT[JF].CVUNITS)
        print('OPEN IN UNIT: ', self.VAROUT[JF].NCID)

    def CMF_OUTPUT_WRITE(self, JYYYYMMDD, JHHMM, JHOUR, JMIN, KSTEP, NX, NY, LOUTINI, NSEQMAX, NPTHOUT, NPTHLEV, REGIONTHIS, D2VEC, R2OUT, R1POUT, D2OUTFLW_oAVG, D2RIVDPH, D2RIVVEL_oAVG, D2FLDOUT_oAVG, D2FLDDPH, D2FLDFRC, D2FLDARE, D2SFCELV, D2STORGE, D2PTHOUT_oAVG, D1PTHFLW_oAVG, D2OUTFLW_oMAX, D2RIVDPH_oMAX, D2STORGE_oMAX, D2OUTINS, P2GDWSTO, D2GDWRTN_oAVG, D2RUNOFF_oAVG, D2ROFSUB_oAVG, D2WEVAPEX_oAVG, P2DAMSTO, D2DAMINF_oAVG, P2LEVSTO, D2LEVDPH):
        print("")
        print("!---------------------!")

        if JHOUR % self.IFRQ_OUT == 0 and JMIN == 0:
            self.IRECOUT += 1
            print('CMF::OUTPUT_WRITE: write at time: ', JYYYYMMDD, JHHMM, self.IRECOUT)

            for JF in range(self.NVARSOUT):
                if self.VAROUT[JF].CVNAME == 'rivsto':
                    D2VEC = P2RIVSTO
                elif self.VAROUT[JF].CVNAME == 'fldsto':
                    D2VEC = P2FLDSTO
                elif self.VAROUT[JF].CVNAME == 'rivout':
                    D2VEC = D2OUTFLW_oAVG
                elif self.VAROUT[JF].CVNAME == 'rivdph':
                    D2VEC = D2RIVDPH
                elif self.VAROUT[JF].CVNAME == 'rivvel':
                    D2VEC = D2RIVVEL_oAVG
                elif self.VAROUT[JF].CVNAME == 'fldout':
                    D2VEC = D2FLDOUT_oAVG
                elif self.VAROUT[JF].CVNAME == 'flddph':
                    D2VEC = D2FLDDPH
                elif self.VAROUT[JF].CVNAME == 'fldfrc':
                    D2VEC = D2FLDFRC
                elif self.VAROUT[JF].CVNAME == 'fldare':
                    D2VEC = D2FLDARE
                elif self.VAROUT[JF].CVNAME == 'sfcelv':
                    D2VEC = D2SFCELV
                elif self.VAROUT[JF].CVNAME == 'totout':
                    D2VEC = D2OUTFLW_oAVG
                elif self.VAROUT[JF].CVNAME == 'outflw':
                    D2VEC = D2OUTFLW_oAVG
                elif self.VAROUT[JF].CVNAME == 'totsto':
                    D2VEC = D2STORGE
                elif self.VAROUT[JF].CVNAME == 'storge':
                    D2VEC = D2STORGE
                elif self.VAROUT[JF].CVNAME == 'pthout':
                    if not self.LPTHOUT:
                        continue
                    D2VEC = D2PTHOUT_oAVG
                elif self.VAROUT[JF].CVNAME == 'pthflw':
                    if not self.LPTHOUT:
                        continue
                elif self.VAROUT[JF].CVNAME == 'maxflw':
                    D2VEC = D2OUTFLW_oMAX
                elif self.VAROUT[JF].CVNAME == 'maxdph':
                    D2VEC = D2RIVDPH_oMAX
                elif self.VAROUT[JF].CVNAME == 'maxsto':
                    D2VEC = D2STORGE_oMAX
                elif self.VAROUT[JF].CVNAME == 'outins':
                    if not self.LOUTINS:
                        continue
                    D2VEC = D2OUTINS
                elif self.VAROUT[JF].CVNAME == 'gwsto':
                    if not self.LGDWDLY:
                        continue
                    D2VEC = P2GDWSTO
                elif self.VAROUT[JF].CVNAME == 'gdwsto':
                    if not self.LGDWDLY:
                        continue
                    D2VEC = P2GDWSTO
                elif self.VAROUT[JF].CVNAME == 'gwout':
                    if not self.LGDWDLY:
                        continue
                    D2VEC = D2GDWRTN_oAVG
                elif self.VAROUT[JF].CVNAME == 'gdwrtn':
                    if not self.LGDWDLY:
                        continue
                    D2VEC = D2GDWRTN_oAVG
                elif self.VAROUT[JF].CVNAME == 'runoff':
                    D2VEC = D2RUNOFF_oAVG
                elif self.VAROUT[JF].CVNAME == 'runoffsub':
                    if not self.LROSPLIT:
                        continue
                    D2VEC = D2ROFSUB_oAVG
                elif self.VAROUT[JF].CVNAME == 'rofsfc':
                    D2VEC = D2RUNOFF_oAVG
                elif self.VAROUT[JF].CVNAME == 'rofsub':
                    D2VEC = D2ROFSUB_oAVG
                elif self.VAROUT[JF].CVNAME == 'wevap':
                    if not self.LWEVAP:
                        continue
                    D2VEC = D2WEVAPEX_oAVG
                elif self.VAROUT[JF].CVNAME == 'damsto':
                    if not self.LDAMOUT:
                        continue
                    D2VEC = P2DAMSTO
                elif self.VAROUT[JF].CVNAME == 'daminf':
                    if not self.LDAMOUT:
                        continue
                    D2VEC = D2DAMINF_oAVG
                elif self.VAROUT[JF].CVNAME == 'levsto':
                    if not self.LLEVEE:
                        continue
                    D2VEC = P2LEVSTO
                elif self.VAROUT[JF].CVNAME == 'levdph':
                    if not self.LLEVEE:
                        continue
                    D2VEC = D2LEVDPH
                else:
                    continue

                if KSTEP == 0 and LOUTINI:
                    if not self.LOUTCDF:
                        continue
                    if self.VAROUT[JF].CVNAME not in ['rivsto', 'fldsto', 'gwsto']:
                        continue

                if self.VAROUT[JF].CVNAME != 'pthflw':
                    R2OUT = D2VEC
                else:
                    if not self.LPTHOUT:
                        continue
                    R1POUT = D1PTHFLW_oAVG

                if self.LOUTCDF:
                    if REGIONTHIS == 1:
                        self.WRTE_OUTCDF(JF, R2OUT, NX, NY)
                else:
                    if self.VAROUT[JF].CVNAME == 'pthflw':
                        if REGIONTHIS == 1:
                            self.WRTE_OUTPTH(self.VAROUT[JF].BINID, self.IRECOUT, R1POUT)
                    else:
                        if self.LOUTVEC:
                            self.WRTE_OUTVEC(self.VAROUT[JF].BINID, self.IRECOUT, D2VEC, NSEQMAX)
                        else:
                            if REGIONTHIS == 1:
                                self.WRTE_OUTBIN(self.VAROUT[JF].BINID, self.IRECOUT, R2OUT, NX, NY)

            print('CMF::OUTPUT_WRITE: end')

    def WRTE_OUTBIN(self, IFN, IREC, R2OUTDAT, NX, NY):
        IFN.write(R2OUTDAT.tobytes())

    def WRTE_OUTPTH(self, IFN, IREC, R2OUTDAT):
        IFN.write(R2OUTDAT.tobytes())

    def WRTE_OUTVEC(self, IFN, IREC, D2OUTDAT, NSEQMAX):
        R2OUTDAT = D2OUTDAT
        IFN.write(R2OUTDAT.tobytes())

    def WRTE_OUTCDF(self, JF, R2OUT, NX, NY):
        XTIME = (self.KMINNEXT - self.KMINSTART) * 60.0
        self.VAROUT[JF].NCID.variables['time'][self.VAROUT[JF].IRECNC] = XTIME
        self.VAROUT[JF].NCID.variables[self.VAROUT[JF].CVNAME][:, :, self.VAROUT[JF].IRECNC] = R2OUT
        self.VAROUT[JF].IRECNC += 1

    def CMF_OUTPUT_END(self, REGIONTHIS):
        print("")
        print("!---------------------!")
        print("CMF::OUTPUT_END: finalize output module")

        if REGIONTHIS == 1:
            if self.LOUTCDF:
                for JF in range(self.NVARSOUT):
                    self.VAROUT[JF].NCID.close()
                    print("Output netcdf output unit closed:", self.VAROUT[JF].NCID)
            else:
                for JF in range(self.NVARSOUT):
                    self.VAROUT[JF].BINID.close()
                    print("Output binary output unit closed:", self.VAROUT[JF].BINID)
                if self.LOUTVEC:
                    self.WRTE_mapR2vecD()

        print("CMF::OUTPUT_END: end")

    def WRTE_mapR2vecD(self):
        CFILE1 = './ind_xy.vec'
        print("LOUTVEC: write mapR2vecD conversion table", CFILE1)
        TMPNAM = self.INQUIRE_FID()
        with open(CFILE1, 'wb') as file:
            file.write(self.I1SEQX.tobytes())
            file.write(self.I1SEQY.tobytes())

    def CMF_OUTTXT_WRTE(self, D2OUTFLW, IYYYYMMDD, ISYYYY, I2VECTOR):
        if self.LOUTTXT:
            if not hasattr(self, 'IsOpen'):
                self.IsOpen = False

            if not self.IsOpen:
                self.IsOpen = True
                self.NGAUGEX = 0
                self.LOGOUTTXT = self.INQUIRE_FID()
                with open(self.CGAUTXT, 'r') as file:
                    self.NGAUGE = int(file.readline().strip())
                    for _ in range(self.NGAUGE):
                        GID, GNAME, GIX, GIY = file.readline().strip().split()
                        GID, GIX, GIY = int(GID), int(GIX), int(GIY)
                        if I2VECTOR[GIX, GIY] > 0:
                            self.NGAUGEX += 1

                self.WriteID = np.zeros(self.NGAUGEX, dtype=int)
                self.WriteISEQ = np.zeros(self.NGAUGEX, dtype=int)
                self.WriteOut = np.zeros(self.NGAUGEX, dtype=float)
                self.WriteName = np.zeros(self.NGAUGEX, dtype='U9')

                self.NGAUGEX = 0
                with open(self.CGAUTXT, 'r') as file:
                    self.NGAUGE = int(file.readline().strip())
                    for _ in range(self.NGAUGE):
                        GID, GNAME, GIX, GIY = file.readline().strip().split()
                        GID, GIX, GIY = int(GID), int(GIX), int(GIY)
                        if I2VECTOR[GIX, GIY] > 0:
                            self.NGAUGEX += 1
                            self.WriteID[self.NGAUGEX - 1] = GID
                            self.WriteName[self.NGAUGEX - 1] = GNAME
                            self.WriteISEQ[self.NGAUGEX - 1] = I2VECTOR[GIX, GIY]

                self.cYYYY = f"{ISYYYY:04d}"
                self.COUTTXT = f"./outtxt-{self.cYYYY}.txt"
                self.LOGOUTTXT = self.INQUIRE_FID()
                with open(self.COUTTXT, 'w') as file:
                    file.write(f"{self.NGAUGEX} " + " ".join(map(str, self.WriteID)) + "\n")
                    file.write(f"{self.NGAUGEX} " + " ".join(self.WriteName) + "\n")

            for IGAUGE in range(self.NGAUGEX):
                GISEQ = self.WriteISEQ[IGAUGE]
                self.WriteOut[IGAUGE] = D2OUTFLW[GISEQ, 0]

            with open(self.COUTTXT, 'a') as file:
                file.write(f"{IYYYYMMDD} " + " ".join(f"{x:.2f}" for x in self.WriteOut) + "\n")

    class TVAROUT:
        def __init__(self):
            self.CVNAME = ""
            self.CVLNAME = ""
            self.CVUNITS = ""
            self.CFILE = ""
            self.BINID = None
            self.NCID = None
            self.VARID = None
            self.TIMID = None
            self.IRECNC = 0

    def INQUIRE_FID(self):
        return 10
