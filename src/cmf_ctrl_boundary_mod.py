import numpy as np

class CMF_CTRL_BOUNDARY_MOD:
    def __init__(self):
        self.LSEALEVCDF = False
        self.CSEALEVDIR = "./sealev/"
        self.CSEALEVPRE = "sealev"
        self.CSEALEVSUF = ".bin"
        self.CSEALEVCDF = "./sealev/"
        self.CVNSEALEV = "variable"
        self.SYEARSL = 0
        self.SMONSL = 0
        self.SDAYSL = 0
        self.SHOURSL = 0
        self.NLINKS = 0
        self.NCDFSTAT = 0
        self.CSLMAP = "./sealev/"
        self.SLCDF = None
        self.R1SLIN = None
        self.I2SLMAP = None

    def CMF_BOUNDARY_NMLIST(self):
        print("CMF::BOUNDARY_NMLIST: namelist OPEN in unit: ", self.CSETFILE, self.NSETFILE)
        self.LSEALEVCDF = False
        self.CSEALEVDIR = "./sealev/"
        self.CSEALEVPRE = "sealev"
        self.CSEALEVSUF = ".bin"
        self.CSEALEVCDF = "./sealev/"
        self.CVNSEALEV = "variable"
        self.CSLMAP = "./sealev/"
        self.SYEARSL = 0
        self.SMONSL = 0
        self.SDAYSL = 0
        self.SHOURSL = 0
        self.IFRQ_SL = 9999
        if self.LSEALEV:
            print("=== NAMELIST, NBOUNDARY ===")
            print("LSEALEVCDF: ", self.LSEALEVCDF)
            if self.LSEALEVCDF:
                print("CSEALEVCDF: ", self.CSEALEVCDF)
                print("CVNSEALEV:  ", self.CVNSEALEV)
                print("SYEARSL:    ", self.SYEARSL)
                print("SMONSL:     ", self.SMONSL)
                print("SDAYSL:     ", self.SDAYSL)
                print("SHOURSL:    ", self.SHOURSL)
                print("CSLMAP:     ", self.CSLMAP)
            else:
                print("CSEALEVDIR: ", self.CSEALEVDIR)
                print("CSEALEVPRE: ", self.CSEALEVPRE)
                print("CSEALEVSUF: ", self.CSEALEVSUF)
            print("IFRQ_SL   ", self.IFRQ_SL)
        self.DTSL = self.IFRQ_SL * 60
        print("CMF::BOUNDARY_NMLIST: end")

    def CMF_BOUNDARY_INIT(self):
        print("CMF::BOUNDARY_INIT: initialize boundary")
        self.D2SEALEV = np.zeros((self.NSEQMAX, 1))
        if self.LSEALEVCDF:
            self.CMF_BOUNDARY_INIT_CDF()
        print("CMF::BOUNDARY_INIT: end")

    def CMF_BOUNDARY_INIT_CDF(self):
        self.KMINSTASL = self.DATE2MIN(self.SYEARSL * 10000 + self.SMONSL * 100 + self.SDAYSL, self.SHOURSL * 100)
        self.SLCDF = {
            "CNAME": self.CSEALEVCDF,
            "CVAR": self.CVNSEALEV,
            "NSTART": self.KMINSTASL
        }
        print("CMF::BOUNRARY_INIT_CDF:", self.SLCDF["CNAME"], self.SLCDF["NSTART"])
        self.SLCDF["NCID"] = self.NF90_OPEN(self.SLCDF["CNAME"], self.NF90_NOWRITE)
        self.SLCDF["NVARID"] = self.NF90_INQ_VARID(self.SLCDF["NCID"], self.SLCDF["CVAR"])
        self.NTIMEID = self.NF90_INQ_DIMID(self.SLCDF["NCID"], "time")
        self.NCDFSTP = self.NF90_INQUIRE_DIMENSION(self.SLCDF["NCID"], self.NTIMEID)
        self.SLCDF["NSTAID"] = self.NF90_INQ_DIMID(self.SLCDF["NCID"], "stations")
        self.NCDFSTAT = self.NF90_INQUIRE_DIMENSION(self.SLCDF["NCID"], self.SLCDF["NSTAID"])
        self.R1SLIN = np.zeros(self.NCDFSTAT)
        print("CMF::BOUNDARY_INIT_CDF: CNAME,NCID,VARID", self.SLCDF["CNAME"], self.SLCDF["NCID"], self.SLCDF["NVARID"])
        if self.KMINSTART < self.KMINSTASL:
            print("Run start earlier than boundary data", self.SLCDF["CNAME"], self.KMINSTART, self.KMINSTASL)
            raise Exception("Run start earlier than boundary data")
        self.KMINENDSL = self.KMINSTASL + self.NCDFSTP * int(self.DTSL / 60)
        if self.KMINEND > self.KMINENDSL:
            print("Run end later than sealev data", self.SLCDF["CNAME"], self.KMINEND, self.KMINENDSL)
            raise Exception("Run end later than sealev data")
        self.TMPNAM = self.INQUIRE_FID()
        self.NLINKS = self.read_nlinks(self.CSLMAP)
        print("Dynamic sea level links", self.NLINKS)
        self.I2SLMAP = np.zeros((3, self.NLINKS), dtype=int)
        for ILNK in range(self.NLINKS):
            IX, IY, IS = self.read_link(self.TMPNAM)
            ISEQ = self.I2VECTOR(IX, IY)
            if ISEQ > 0:
                if self.I1NEXT[ISEQ] != -9:
                    print("Sealev link not at river outlet cell", IX, IY)
                    raise Exception("Sealev link not at river outlet cell")
                elif IS < 1 or IS > self.NCDFSTAT:
                    print("Sealev link outside netcdf index", IS)
                    raise Exception("Sealev link outside netcdf index")
            else:
                print("Sealev link outside land grids", IX, IY)
                raise Exception("Sealev link outside land grids")
            self.I2SLMAP[0, ILNK] = IX
            self.I2SLMAP[1, ILNK] = IY
            self.I2SLMAP[2, ILNK] = IS

    def CMF_BOUNDARY_UPDATE(self):
        IUPDATE = 0
        if int(self.IMIN) % self.IFRQ_SL == 0:
            IUPDATE = 1
        if self.LSEALEV and IUPDATE == 1:
            print("CMF::BOUNDARY_UPDATE: update at time: ", self.IYYYYMMDD, self.IHHMM)
            if self.LSEALEVCDF:
                self.CMF_BOUNDARY_GET_CDF()
            else:
                self.CMF_BOUNDARY_GET_BIN()
        if self.LMEANSL:
            self.D2DWNELV = self.D2ELEVTN + self.D2MEANSL
        else:
            self.D2DWNELV = self.D2ELEVTN
        if self.LSEALEV:
            self.D2DWNELV += self.D2SEALEV

    def CMF_BOUNDARY_GET_BIN(self):
        CDATE = f"{self.IYYYY:04d}{self.IMM:02d}{self.IDD:02d}{self.IHHMM:04d}"
        CIFNAME = f"{self.CSEALEVDIR}/{self.CSEALEVPRE}{CDATE}{self.CSEALEVSUF}"
        print("CMF::BOUNDARY_GET_BIN: read sealev:", CIFNAME)
        self.TMPNAM = self.INQUIRE_FID()
        R2TMP = self.read_binary_file(CIFNAME)
        self.D2SEALEV = self.mapR2vecD(R2TMP)

    def CMF_BOUNDARY_GET_CDF(self):
        IRECSL = (self.KMIN - self.SLCDF["NSTART"]) * 60 // int(self.DTSL) + 1
        print("CMF::BOUNDARY_GET_CDF:", self.SLCDF["CNAME"], IRECSL)
        self.R1SLIN = self.NF90_GET_VAR(self.SLCDF["NCID"], self.SLCDF["NVARID"], IRECSL)
        R2TMP = np.zeros((self.NX, self.NY))
        for ILNK in range(self.NLINKS):
            IX = self.I2SLMAP[0, ILNK]
            IY = self.I2SLMAP[1, ILNK]
            IS = self.I2SLMAP[2, ILNK]
            R2TMP[IX, IY] = self.R1SLIN[IS]
        self.D2SEALEV = self.mapR2vecD(R2TMP)

    def CMF_BOUNDARY_END(self):
        print("CMF::BOUNDARY_END: Finalize boundary module")
        if self.LSEALEVCDF:
            self.NF90_CLOSE(self.SLCDF["NCID"])
            print("Input netcdf sealev closed:", self.SLCDF["NCID"])
        print("CMF::BOUNDARY_END: end")

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

    def read_nlinks(self, cslmap):
        # Implement the read_nlinks function here
        pass

    def read_link(self, tmpnam):
        # Implement the read_link function here
        pass

    def I2VECTOR(self, ix, iy):
        # Implement the I2VECTOR function here
        pass

    def read_binary_file(self, cifname):
        # Implement the read_binary_file function here
        pass

    def mapR2vecD(self, r2tmp):
        # Implement the mapR2vecD function here
        pass
