import numpy as np

class CMF_CTRL_FORCING_MOD:
    def __init__(self):
        self.LINPCDF = False
        self.LINPEND = False
        self.LINTERP = False
        self.LITRPCDF = False
        self.CINPMAT = "NONE"
        self.DROFUNIT = 86400 * 1000.0
        self.CROFDIR = "./runoff/"
        self.CROFPRE = "Roff____"
        self.CROFSUF = ".one"
        self.CSUBDIR = "./runoff/"
        self.CSUBPRE = "Rsub____"
        self.CSUBSUF = ".one"
        self.CROFCDF = "NONE"
        self.CVNROF = "runoff"
        self.CVNSUB = "NONE"
        self.SYEARIN = 0
        self.SMONIN = 0
        self.SDAYIN = 0
        self.SHOURIN = 0
        self.ROFCDF = None
        self.INPX = None
        self.INPY = None
        self.INPA = None
        self.INPXI = None
        self.INPYI = None
        self.INPAI = None

    def CMF_FORCING_NMLIST(self):
        print("CMF::FORCING_NMLIST: namelist OPEN in unit: ", self.CSETFILE, self.NSETFILE)
        self.LINPCDF = False
        self.LINPEND = False
        self.LINTERP = False
        self.LITRPCDF = False
        self.CINPMAT = "NONE"
        self.DROFUNIT = 86400 * 1000.0
        self.CROFDIR = "./runoff/"
        self.CROFPRE = "Roff____"
        self.CROFSUF = ".one"
        self.CSUBDIR = "./runoff/"
        self.CSUBPRE = "Rsub____"
        self.CSUBSUF = ".one"
        self.CROFCDF = "NONE"
        self.CVNROF = "runoff"
        self.CVNSUB = "NONE"
        self.SYEARIN = 0
        self.SMONIN = 0
        self.SDAYIN = 0
        self.SHOURIN = 0
        print("=== NAMELIST, NFORCE ===")
        print("LINPCDF:   ", self.LINPCDF)
        print("LINTERP:   ", self.LINTERP)
        print("LITRPCDF:  ", self.LITRPCDF)
        print("CINPMAT:   ", self.CINPMAT)
        print("LROSPLIT:  ", self.LROSPLIT)
        if not self.LINPCDF:
            print("CROFDIR:   ", self.CROFDIR)
            print("CROFPRE:   ", self.CROFPRE)
            print("CROFSUF:   ", self.CROFSUF)
            if self.LROSPLIT:
                print("CSUBDIR:   ", self.CSUBDIR)
                print("CSUBPRE:   ", self.CSUBPRE)
                print("CSUBSUF:   ", self.CSUBSUF)
        else:
            print("CROFCDF:   ", self.CROFCDF)
            print("CVNROF:    ", self.CVNROF)
            if self.LROSPLIT:
                print("CVNSUB:    ", self.CVNSUB)
            print("SYEARIN,SMONIN,SDAYIN,SHOURIN ", self.SYEARIN, self.SMONIN, self.SDAYIN, self.SHOURIN)
        if self.LINPEND:
            print("LINPEND:   ", self.LINPEND)
        print("CMF::FORCING_NMLIST: end")

    def CMF_FORCING_INIT(self, LECMF2LAKEC=None):
        print("CMF::FORCING_INIT: Initialize runoff forcing file (only for netCDF)")
        if self.LINPCDF:
            self.CMF_FORCING_INIT_CDF()
        if self.LINTERP:
            if self.LITRPCDF:
                self.CMF_INPMAT_INIT_CDF(LECMF2LAKEC)
            else:
                self.CMF_INPMAT_INIT_BIN()
        print("CMF::FORCING_INIT: end")

    def CMF_FORCING_INIT_CDF(self):
        self.KMINSTAIN = self.DATE2MIN(self.SYEARIN * 10000 + self.SMONIN * 100 + self.SDAYIN, self.SHOURIN * 100)
        self.ROFCDF = {
            "CNAME": self.CROFCDF,
            "CVAR": [self.CVNROF, self.CVNSUB, "NONE"],
            "NSTART": self.KMINSTAIN
        }
        if not self.LROSPLIT:
            self.ROFCDF["CVAR"][1] = "NONE"
        if not self.LWEVAP:
            self.ROFCDF["CVAR"][2] = "NONE"
        print("CMF::FORCING_INIT_CDF:", self.ROFCDF["CNAME"], self.ROFCDF["CVAR"][0])
        self.ROFCDF["NCID"] = self.NF90_OPEN(self.ROFCDF["CNAME"], self.NF90_NOWRITE)
        self.ROFCDF["NVARID"] = [self.NF90_INQ_VARID(self.ROFCDF["NCID"], self.ROFCDF["CVAR"][0])]
        if self.LROSPLIT:
            self.ROFCDF["NVARID"].append(self.NF90_INQ_VARID(self.ROFCDF["NCID"], self.ROFCDF["CVAR"][1]))
        if self.LWEVAP:
            self.ROFCDF["NVARID"].append(self.NF90_INQ_VARID(self.ROFCDF["NCID"], self.ROFCDF["CVAR"][2]))
        self.NTIMEID = self.NF90_INQ_DIMID(self.ROFCDF["NCID"], "time")
        self.NCDFSTP = self.NF90_INQUIRE_DIMENSION(self.ROFCDF["NCID"], self.NTIMEID)
        print("CMF::FORCING_INIT_CDF: CNAME,NCID,VARID", self.ROFCDF["CNAME"], self.ROFCDF["NCID"], self.ROFCDF["NVARID"][0])
        if self.KMINSTART < self.KMINSTAIN:
            print("Run start earlier than forcing data", self.ROFCDF["CNAME"], self.KMINSTART, self.KMINSTAIN)
            raise Exception("Run start earlier than forcing data")
        self.KMINENDIN = self.KMINSTAIN + self.NCDFSTP * int(self.DTIN / 60)
        if self.KMINEND > self.KMINENDIN:
            print("Run end later than forcing data", self.ROFCDF["CNAME"], self.KMINEND, self.KMINENDIN)
            raise Exception("Run end later than forcing data")

    def CMF_INPMAT_INIT_CDF(self, LECMF2LAKEC=None):
        print('INPUT MATRIX netCDF', self.CINPMAT)
        self.NCID = self.NF90_OPEN(self.CINPMAT, self.NF90_NOWRITE)
        self.INPX = np.zeros((self.NSEQMAX, self.INPN), dtype=int)
        self.INPY = np.zeros((self.NSEQMAX, self.INPN), dtype=int)
        self.INPA = np.zeros((self.NSEQMAX, self.INPN))
        self.D2TMP = self.NF90_GET_VAR(self.NCID, self.NF90_INQ_VARID(self.NCID, 'inpa'), (1, 1, 1), (self.NX, self.NY, self.INPN))
        for INPI in range(self.INPN):
            self.INPA[:, INPI] = self.mapD2vecD(self.D2TMP[:, :, INPI])
        self.I2TMP = self.NF90_GET_VAR(self.NCID, self.NF90_INQ_VARID(self.NCID, 'inpx'), (1, 1, 1), (self.NX, self.NY, self.INPN))
        for INPI in range(self.INPN):
            self.INPX[:, INPI] = self.mapI2vecI(self.I2TMP[:, :, INPI])
        self.I2TMP = self.NF90_GET_VAR(self.NCID, self.NF90_INQ_VARID(self.NCID, 'inpy'), (1, 1, 1), (self.NX, self.NY, self.INPN))
        for INPI in range(self.INPN):
            self.INPY[:, INPI] = self.mapI2vecI(self.I2TMP[:, :, INPI])
        if LECMF2LAKEC is not None and LECMF2LAKEC != 0:
            self.INPNI = self.NF90_INQUIRE_DIMENSION(self.NCID, self.NF90_INQ_DIMID(self.NCID, 'levI'))
            if self.INPNI == -1:
                print("Could not find levI variable in inpmat.nc: inverse interpolation not available")
                return
            self.INPXI = np.zeros((self.NXIN, self.NYIN, self.INPNI), dtype=int)
            self.INPYI = np.zeros((self.NXIN, self.NYIN, self.INPNI), dtype=int)
            self.INPAI = np.zeros((self.NXIN, self.NYIN, self.INPNI))
            self.INPAI = self.NF90_GET_VAR(self.NCID, self.NF90_INQ_VARID(self.NCID, 'inpaI'), (1, 1, 1), (self.NXIN, self.NYIN, self.INPNI))
            self.INPXI = self.NF90_GET_VAR(self.NCID, self.NF90_INQ_VARID(self.NCID, 'inpxI'), (1, 1, 1), (self.NXIN, self.NYIN, self.INPNI))
            self.INPYI = self.NF90_GET_VAR(self.NCID, self.NF90_INQ_VARID(self.NCID, 'inpyI'), (1, 1, 1), (self.NXIN, self.NYIN, self.INPNI))
            for IX in range(self.NXIN):
                for IY in range(self.NYIN):
                    ZTMP = np.sum(self.INPAI[IX, IY, :])
                    if ZTMP > 0.0:
                        self.INPAI[IX, IY, :] /= ZTMP

    def CMF_INPMAT_INIT_BIN(self):
        print('INPUT MATRIX binary', self.CINPMAT)
        self.INPX = np.zeros((self.NSEQMAX, self.INPN), dtype=int)
        self.INPY = np.zeros((self.NSEQMAX, self.INPN), dtype=int)
        self.INPA = np.zeros((self.NSEQMAX, self.INPN))
        self.I2TMP = np.zeros((self.NX, self.NY), dtype=int)
        self.R2TMP = np.zeros((self.NX, self.NY))
        self.TMPNAM = self.INQUIRE_FID()
        for INPI in range(self.INPN):
            self.I2TMP = self.read_binary_file(self.CINPMAT, INPI)
            self.INPX[:, INPI] = self.mapI2vecI(self.I2TMP)
            self.I2TMP = self.read_binary_file(self.CINPMAT, self.INPN + INPI)
            self.INPY[:, INPI] = self.mapI2vecI(self.I2TMP)
            self.R2TMP = self.read_binary_file(self.CINPMAT, 2 * self.INPN + INPI)
            self.INPA[:, INPI] = self.mapR2vecD(self.R2TMP)

    def CMF_FORCING_GET(self, PBUFF):
        if self.LINPCDF:
            self.CMF_FORCING_GET_CDF(PBUFF)
        else:
            self.CMF_FORCING_GET_BIN(PBUFF)

    def CMF_FORCING_GET_BIN(self, PBUFF):
        ISEC = self.IHOUR * 60 * 60 + self.IMIN * 60
        IRECINP = int(ISEC / self.DTIN) + 1
        CDATE = f"{self.IYYYY:04d}{self.IMM:02d}{self.IDD:02d}"
        CIFNAME = f"{self.CROFDIR}/{self.CROFPRE}{CDATE}{self.CROFSUF}"
        print("CMF::FORCING_GET_BIN:", CIFNAME)
        self.TMPNAM = self.INQUIRE_FID()
        R2TMP = self.read_binary_file(CIFNAME, IRECINP)
        if self.LINPEND:
            R2TMP = self.CONV_END(R2TMP)
        PBUFF[:, :, 0] = R2TMP
        PBUFF[:, :, 1] = 0.0
        if self.LROSPLIT:
            CIFNAME = f"{self.CSUBDIR}/{self.CSUBPRE}{CDATE}{self.CSUBSUF}"
            print("CMF::FORCING_GET_BIN: (sub-surface)", CIFNAME)
            R2TMP = self.read_binary_file(CIFNAME, IRECINP)
            if self.LINPEND:
                R2TMP = self.CONV_END(R2TMP)
            PBUFF[:, :, 1] = R2TMP

    def CMF_FORCING_GET_CDF(self, PBUFF):
        IRECINP = (self.KMIN - self.ROFCDF["NSTART"]) * 60 // int(self.DTIN) + 1
        self.NF90_GET_VAR(self.ROFCDF["NCID"], self.ROFCDF["NVARID"][0], PBUFF[:, :, 0], (1, 1, IRECINP), (self.NXIN, self.NYIN, 1))
        if self.ROFCDF["NVARID"][1] != -1:
            self.NF90_GET_VAR(self.ROFCDF["NCID"], self.ROFCDF["NVARID"][1], PBUFF[:, :, 1], (1, 1, IRECINP), (self.NXIN, self.NYIN, 1))
        print("CMF::FORCING_GET_CDF: read runoff:", self.IYYYYMMDD, self.IHHMM, IRECINP)

    def CMF_FORCING_COM(self, PBUFF):
        D2MAPTMP = self.vecD2mapD(self.D2FLDFRC)
        self.CMF_MPI_AllReduce_D2MAP(D2MAPTMP)
        self.INTERPI(D2MAPTMP, PBUFF[:, :, 0])

    def INTERPI(self, PBUFFIN, PBUFFOUT):
        if self.INPNI == -1:
            print("INPNI==-1, no inverse interpolation possible")
            raise Exception("INPNI==-1, no inverse interpolation possible")
        PBUFFOUT[:, :] = 1.0
        for IXIN in range(self.NXIN):
            for IYIN in range(self.NYIN):
                PBUFFOUT[IXIN, IYIN] = 0.0
                for INP in range(self.INPNI):
                    IX = self.INPXI[IXIN, IYIN, INP]
                    IY = self.INPYI[IXIN, IYIN, INP]
                    if 0 < IX <= self.NX and 0 < IY <= self.NY:
                        PBUFFOUT[IXIN, IYIN] += PBUFFIN[IX, IY] * self.INPAI[IXIN, IYIN, INP]

    def CMF_FORCING_PUT(self, PBUFF):
        if self.LINTERP:
            self.ROFF_INTERP(PBUFF[:, :, 0], self.D2RUNOFF)
            if self.LROSPLIT:
                self.ROFF_INTERP(PBUFF[:, :, 1], self.D2ROFSUB)
            else:
                self.D2ROFSUB[:, :] = 0.0
        else:
            self.CONV_RESOL(PBUFF[:, :, 0], self.D2RUNOFF)
            if self.LROSPLIT:
                self.CONV_RESOL(PBUFF[:, :, 1], self.D2ROFSUB)
            else:
                self.D2ROFSUB[:, :] = 0.0
        if self.LWEVAP:
            if PBUFF.shape[2] == 3:
                self.ROFF_INTERP(PBUFF[:, :, 2], self.D2WEVAP)
            else:
                print("LWEVAP is true but evaporation not provide in input array for interpolation")
                print("CMF_FORCING_PUT(PBUFF), PBUFF should have 3 fields for interpolation ")
                raise Exception("LWEVAP is true but evaporation not provide in input array for interpolation")

    def ROFF_INTERP(self, PBUFFIN, PBUFFOUT):
        for ISEQ in range(self.NSEQALL):
            PBUFFOUT[ISEQ, 0] = 0.0
            for INPI in range(self.INPN):
                IXIN = self.INPX[ISEQ, INPI]
                IYIN = self.INPY[ISEQ, INPI]
                if IXIN > 0:
                    if IXIN > self.NXIN or IYIN > self.NYIN:
                        print("error")
                        print('XXX', ISEQ, INPI, IXIN, IYIN)
                        continue
                    if PBUFFIN[IXIN, IYIN] != self.RMIS and not self.CMF_CheckNanB(PBUFFIN[IXIN, IYIN], 0.0):
                        PBUFFOUT[ISEQ, 0] += PBUFFIN[IXIN, IYIN] * self.INPA[ISEQ, INPI] / self.DROFUNIT
            PBUFFOUT[ISEQ, 0] = max(PBUFFOUT[ISEQ, 0], 0.0)

    def CONV_RESOL(self, PBUFFIN, PBUFFOUT):
        D2TEMP = self.mapD2vecD(PBUFFIN)
        for ISEQ in range(self.NSEQALL):
            if D2TEMP[ISEQ, 0] != self.RMIS:
                PBUFFOUT[ISEQ, 0] = D2TEMP[ISEQ, 0] * self.D2GRAREA[ISEQ, 0] / self.DROFUNIT
                PBUFFOUT[ISEQ, 0] = max(PBUFFOUT[ISEQ, 0], 0.0)
            else:
                PBUFFOUT[ISEQ, 0] = 0.0

    def CMF_FORCING_END(self):
        print("CMF::FORCING_END: Finalize forcing module")
        if self.LINPCDF:
            self.NF90_CLOSE(self.ROFCDF["NCID"])
            print("input netCDF runoff closed:", self.ROFCDF["NCID"])
        print("CMF::FORCING_END: end")

    def DATE2MIN(self, date, hour):
        # Implement the DATE2MIN function here
        pass

    def NF90_OPEN(self, cname, mode):
        # Implement the NF90_OPEN function here
        pass

    def NF90_INQ_VARID(self, ncid, varname):
        # Implement the NF90_INQ_VARID function here
        pass

    def NF90_INQ_DIMID(self, ncid, dimname):
        # Implement the NF90_INQ_DIMID function here
        pass

    def NF90_INQUIRE_DIMENSION(self, ncid, dimid):
        # Implement the NF90_INQUIRE_DIMENSION function here
        pass

    def NF90_GET_VAR(self, ncid, varid, rec):
        # Implement the NF90_GET_VAR function here
        pass

    def NF90_CLOSE(self, ncid):
        # Implement the NF90_CLOSE function here
        pass

    def INQUIRE_FID(self):
        # Implement the INQUIRE_FID function here
        pass

    def read_binary_file(self, filename, rec):
        # Implement the read_binary_file function here
        pass

    def CONV_END(self, data):
        # Implement the CONV_END function here
        pass

    def mapD2vecD(self, data):
        # Implement the mapD2vecD function here
        pass

    def mapI2vecI(self, data):
        # Implement the mapI2vecI function here
        pass

    def vecD2mapD(self, data):
        # Implement the vecD2mapD function here
        pass

    def CMF_MPI_AllReduce_D2MAP(self, data):
        # Implement the CMF_MPI_AllReduce_D2MAP function here
        pass

    def CMF_CheckNanB(self, value, default):
        # Implement the CMF_CheckNanB function here
        pass
