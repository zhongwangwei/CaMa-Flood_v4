import numpy as np

class CMF_CTRL_MAPS_MOD:
    def __init__(self, NX, NY, NLFP, LPTHOUT, LMAPCDF, LMAPEND, LMEANSL, LGDWDLY, LSLPMIX, LSLOPEMOUTH, WEST, EAST, NORTH, SOUTH, PMANRIV, PMANFLD, IMIS, REGIONTHIS, REGIONALL, NSEQMAX, NSEQALL, NSEQRIV, DFRCINC, CNEXTXY, CGRAREA, CELEVTN, CNXTDST, CRIVLEN, CFLDHGT, CRIVWTH, CRIVHGT, CRIVMAN, CPTHOUT, CGDWDLY, CMEANSL, CMPIREG, CRIVCLINC, CRIVPARNC, CMEANSLNC, CMPIREGNC):
        self.NX = NX
        self.NY = NY
        self.NLFP = NLFP
        self.LPTHOUT = LPTHOUT
        self.LMAPCDF = LMAPCDF
        self.LMAPEND = LMAPEND
        self.LMEANSL = LMEANSL
        self.LGDWDLY = LGDWDLY
        self.LSLPMIX = LSLPMIX
        self.LSLOPEMOUTH = LSLOPEMOUTH
        self.WEST = WEST
        self.EAST = EAST
        self.NORTH = NORTH
        self.SOUTH = SOUTH
        self.PMANRIV = PMANRIV
        self.PMANFLD = PMANFLD
        self.IMIS = IMIS
        self.REGIONTHIS = REGIONTHIS
        self.REGIONALL = REGIONALL
        self.NSEQMAX = NSEQMAX
        self.NSEQALL = NSEQALL
        self.NSEQRIV = NSEQRIV
        self.DFRCINC = DFRCINC
        self.CNEXTXY = CNEXTXY
        self.CGRAREA = CGRAREA
        self.CELEVTN = CELEVTN
        self.CNXTDST = CNXTDST
        self.CRIVLEN = CRIVLEN
        self.CFLDHGT = CFLDHGT
        self.CRIVWTH = CRIVWTH
        self.CRIVHGT = CRIVHGT
        self.CRIVMAN = CRIVMAN
        self.CPTHOUT = CPTHOUT
        self.CGDWDLY = CGDWDLY
        self.CMEANSL = CMEANSL
        self.CMPIREG = CMPIREG
        self.CRIVCLINC = CRIVCLINC
        self.CRIVPARNC = CRIVPARNC
        self.CMEANSLNC = CMEANSLNC
        self.CMPIREGNC = CMPIREGNC

        self.I2NEXTX = np.zeros((NX, NY), dtype=int)
        self.I2NEXTY = np.zeros((NX, NY), dtype=int)
        self.I2REGION = np.zeros((NX, NY), dtype=int)
        self.D1LON = np.zeros(NX)
        self.D1LAT = np.zeros(NY)
        self.I1SEQX = np.zeros(NSEQMAX, dtype=int)
        self.I1SEQY = np.zeros(NSEQMAX, dtype=int)
        self.I1NEXT = np.zeros(NSEQMAX, dtype=int)
        self.I2VECTOR = np.zeros((NX, NY), dtype=int)
        self.D2GRAREA = np.zeros((NSEQMAX, 1))
        self.D2ELEVTN = np.zeros((NSEQMAX, 1))
        self.D2NXTDST = np.zeros((NSEQMAX, 1))
        self.D2RIVLEN = np.zeros((NSEQMAX, 1))
        self.D2RIVWTH = np.zeros((NSEQMAX, 1))
        self.D2RIVHGT = np.zeros((NSEQMAX, 1))
        self.D2FLDHGT = np.zeros((NSEQMAX, 1, NLFP))
        self.D2RIVMAN = np.zeros((NSEQMAX, 1))
        self.D2MEANSL = np.zeros((NSEQMAX, 1))
        self.D2DWNELV = np.zeros((NSEQMAX, 1))
        self.D2GDWDLY = np.zeros((NSEQMAX, 1))
        self.I2MASK = np.zeros((NSEQMAX, 1), dtype=int)
        self.P2RIVSTOMAX = np.zeros((NSEQMAX, 1))
        self.D2RIVELV = np.zeros((NSEQMAX, 1))
        self.P2FLDSTOMAX = np.zeros((NSEQMAX, 1, NLFP))
        self.D2FLDGRD = np.zeros((NSEQMAX, 1, NLFP))

    def CMF_MAPS_NMLIST(self, CSETFILE, NSETFILE, LMEANSL, LGDWDLY):
        print("")
        print("!---------------------!")

        NSETFILE = self.INQUIRE_FID()
        with open(CSETFILE, 'r') as file:
            print("CMF::MAP_NMLIST: namelist OPEN in unit: ", CSETFILE, NSETFILE)

            self.CNEXTXY = "./nextxy.bin"
            self.CGRAREA = "./ctmare.bin"
            self.CELEVTN = "./elevtn.bin"
            self.CNXTDST = "./nxtdst.bin"
            self.CRIVLEN = "./rivlen.bin"
            self.CFLDHGT = "./fldhgt.bin"
            self.CRIVWTH = "./rivwth.bin"
            self.CRIVHGT = "./rivhgt.bin"
            self.CRIVMAN = "./rivman.bin"
            self.CPTHOUT = "./bifprm.txt"
            self.CGDWDLY = "NONE"
            self.CMEANSL = "NONE"
            self.CMPIREG = "NONE"
            self.LMAPCDF = False
            self.CRIVCLINC = "NONE"
            self.CRIVPARNC = "NONE"
            self.CMEANSLNC = "NONE"
            self.CMPIREGNC = "NONE"

            file.seek(0)
            for line in file:
                if "NMAP" in line:
                    break

            print("=== NAMELIST, NMAP ===")
            print("LMAPCDF:   ", self.LMAPCDF)
            if self.LMAPCDF:
                print("CRIVCLINC: ", self.CRIVCLINC)
                print("CRIVPARNC: ", self.CRIVPARNC)
                if LMEANSL:
                    print("CMEANSLNC: ", self.CMEANSLNC)
                print("CMPIREGNC:   ", self.CMPIREGNC)
            else:
                print("CNEXTXY:   ", self.CNEXTXY)
                print("CGRAREA:   ", self.CGRAREA)
                print("CELEVTN:   ", self.CELEVTN)
                print("CNXTDST:   ", self.CNXTDST)
                print("CRIVLEN:   ", self.CRIVLEN)
                print("CFLDHGT:   ", self.CFLDHGT)
                print("CRIVWTH:   ", self.CRIVWTH)
                print("CRIVHGT:   ", self.CRIVHGT)
                print("CRIVMAN:   ", self.CRIVMAN)
                print("CPTHOUT:   ", self.CPTHOUT)
                if LGDWDLY:
                    print("CGDWDLY:    ", self.CGDWDLY)
                if LMEANSL:
                    print("CMEANSL:   ", self.CMEANSL)
                print("CMPIREG:   ", self.CMPIREG)

        print("CMF::MAP_NMLIST: end")

    def CMF_RIVMAP_INIT(self, TMPNAM, NX, NY, NLFP, LPTHOUT, LMAPCDF, LMAPEND, LMEANSL, LGDWDLY, LSLPMIX, LSLOPEMOUTH, WEST, EAST, NORTH, SOUTH, PMANRIV, PMANFLD, IMIS, REGIONTHIS, REGIONALL, NSEQMAX, NSEQALL, NSEQRIV, DFRCINC, CNEXTXY, CGRAREA, CELEVTN, CNXTDST, CRIVLEN, CFLDHGT, CRIVWTH, CRIVHGT, CRIVMAN, CPTHOUT, CGDWDLY, CMEANSL, CMPIREG, CRIVCLINC, CRIVPARNC, CMEANSLNC, CMPIREGNC):
        print("")
        print("!---------------------!")
        print('CMF::RIVMAP_INIT: river network initialization')

        self.I2NEXTX = np.zeros((NX, NY), dtype=int)
        self.I2NEXTY = np.zeros((NX, NY), dtype=int)
        self.I2REGION = np.zeros((NX, NY), dtype=int)
        self.D1LON = np.zeros(NX)
        self.D1LAT = np.zeros(NY)

        print('CMF::RIVMAP_INIT: read nextXY & set lat lon')
        if LMAPCDF:
            self.READ_MAP_CDF()
        else:
            self.READ_MAP_BIN()

        print('CMF::RIVMAP_INIT: calc region')
        self.CALC_REGION()

        print('CMF::RIVMAP_INIT: calculate 1d river sequence')
        self.CALC_1D_SEQ()

        print('  NSEQRIV=', self.NSEQRIV)
        print('  NSEQALL=', self.NSEQALL)

        if self.REGIONTHIS == 1:
            TMPNAM = self.INQUIRE_FID()
            with open('./mapdata.txt', 'w') as file:
                file.write('NX {}\n'.format(self.NX))
                file.write('NY {}\n'.format(self.NY))
                file.write('NLFP {}\n'.format(self.NLFP))
                file.write('REGIONALL {}\n'.format(self.REGIONALL))
                file.write('NSEQMAX {}\n'.format(self.NSEQMAX))

        if self.LPTHOUT:
            print('CMF::RIVMAP_INIT: read bifurcation channel setting')
            self.READ_BIFPARAM()

        print('CMF::RIVMAP_INIT: end')

    def CMF_TOPO_INIT(self, TMPNAM, NX, NY, NLFP, LMAPEND, LFPLAIN, LMEANSL, LGDWDLY, LSLPMIX, LSLOPEMOUTH, D2NXTDST, D2GRAREA, D2ELEVTN, D2RIVLEN, D2RIVWTH, D2RIVHGT, D2FLDHGT, D2RIVELV, D2FLDGRD, D2RIVMAN, P2RIVSTOMAX, P2FLDSTOMAX, DFRCINC, NSEQALL, NSEQMAX, D2MEANSL, D2DWNELV, D2GDWDLY, I2MASK):
        print("")
        print("!---------------------!")
        print('CMF::TOPO_INIT: topography map initialization')

        self.D2GRAREA = np.zeros((NSEQMAX, 1))
        self.D2ELEVTN = np.zeros((NSEQMAX, 1))
        self.D2NXTDST = np.zeros((NSEQMAX, 1))
        self.D2RIVLEN = np.zeros((NSEQMAX, 1))
        self.D2RIVWTH = np.zeros((NSEQMAX, 1))
        self.D2RIVHGT = np.zeros((NSEQMAX, 1))
        self.D2FLDHGT = np.zeros((NSEQMAX, 1, NLFP))
        self.D2RIVMAN = np.zeros((NSEQMAX, 1))
        self.D2MEANSL = np.zeros((NSEQMAX, 1))
        self.D2DWNELV = np.zeros((NSEQMAX, 1))
        self.D2GDWDLY = np.zeros((NSEQMAX, 1))
        self.I2MASK = np.zeros((NSEQMAX, 1), dtype=int)

        print('CMF::TOPO_INIT: read topography maps')
        if not self.LMAPCDF:
            self.READ_TOPO_BIN()
        else:
            self.READ_TOPO_CDF()

        print('TOPO_INIT: calc river channel parameters')
        self.P2RIVSTOMAX = np.zeros((NSEQMAX, 1))
        self.D2RIVELV = np.zeros((NSEQMAX, 1))

        if LFPLAIN:
            self.P2RIVSTOMAX[:, :] = self.D2RIVLEN[:, :] * self.D2RIVWTH[:, :] * self.D2RIVHGT[:, :]
        else:
            print('TOPO_INIT: no floodplain (rivstomax=1.D18)')
            self.P2RIVSTOMAX[:, :] = 1.E18

        self.D2RIVELV[:, :] = self.D2ELEVTN[:, :] - self.D2RIVHGT[:, :]

        print('TOPO_INIT: calc floodplain parameters')
        self.P2FLDSTOMAX = np.zeros((NSEQMAX, 1, NLFP))
        self.D2FLDGRD = np.zeros((NSEQMAX, 1, NLFP))
        self.SET_FLDSTG()

        print('TOPO_INIT: calc downstream boundary elevation')
        self.D2DWNELV[:, :] = self.D2ELEVTN[:, :]
        if self.LMEANSL:
            self.D2DWNELV[:, :] = self.D2ELEVTN[:, :] + self.D2MEANSL[:, :]

    def INQUIRE_FID(self):
        return 10

    def READ_MAP_BIN(self):
        print('RIVMAP_INIT: nextxy binary: ', self.CNEXTXY)
        TMPNAM = self.INQUIRE_FID()
        with open(self.CNEXTXY, 'rb') as file:
            self.I2NEXTX = np.fromfile(file, dtype=np.int32, count=self.NX * self.NY).reshape((self.NX, self.NY))
            self.I2NEXTY = np.fromfile(file, dtype=np.int32, count=self.NX * self.NY).reshape((self.NX, self.NY))

        if self.LMAPEND:
            self.CONV_ENDI(self.I2NEXTX, self.NX, self.NY)
            self.CONV_ENDI(self.I2NEXTY, self.NX, self.NY)

        if self.WEST >= -180.0 and self.EAST <= 360.0 and self.SOUTH >= -180.0 and self.NORTH <= 180.0:
            for IX in range(self.NX):
                self.D1LON[IX] = self.WEST + (IX + 0.5) * (self.EAST - self.WEST) / self.NX
            for IY in range(self.NY):
                self.D1LAT[IY] = self.NORTH - (IY + 0.5) * (self.NORTH - self.SOUTH) / self.NY

    def READ_MAP_CDF(self):
        print('RIVMAP_INIT: nextxy netCDF: ', self.CRIVCLINC)
        import netCDF4 as nc

        with nc.Dataset(self.CRIVCLINC, 'r') as ds:
            self.I2NEXTX = ds.variables['nextx'][:]
            self.I2NEXTY = ds.variables['nexty'][:]
            self.D1LAT = ds.variables['lat'][:]
            self.D1LON = ds.variables['lon'][:]

    def CALC_REGION(self):
        print('RIVMAP_INIT: region code')

        self.REGIONALL = 1
        self.I2REGION[:, :] = self.IMIS
        for IY in range(self.NY):
            for IX in range(self.NX):
                if self.I2NEXTX[IX, IY] != self.IMIS:
                    self.I2REGION[IX, IY] = 1

        print('RIVMAP_INIT: count number of grid in each region: ')
        REGIONGRID = np.zeros(self.REGIONALL, dtype=int)
        for IY in range(self.NY):
            for IX in range(self.NX):
                if self.I2REGION[IX, IY] > 0:
                    IREGION = self.I2REGION[IX, IY]
                    REGIONGRID[IREGION - 1] += 1

        self.NSEQMAX = np.max(REGIONGRID)

        print('CALC_REGION: REGIONALL= ', self.REGIONALL)
        print('CALC_REGION: NSEQMAX=', self.NSEQMAX)
        print('CALC_REGION: NSEQALL=', self.NSEQALL)

    def CALC_1D_SEQ(self):
        print('RIVMAP_INIT: convert 2D map to 1D sequence')

        NUPST = np.zeros((self.NX, self.NY), dtype=int)
        UPNOW = np.zeros((self.NX, self.NY), dtype=int)

        self.I1SEQX = np.zeros(self.NSEQMAX, dtype=int)
        self.I1SEQY = np.zeros(self.NSEQMAX, dtype=int)
        self.I1NEXT = np.zeros(self.NSEQMAX, dtype=int)
        self.I2VECTOR = np.zeros((self.NX, self.NY), dtype=int)

        for IY in range(self.NY):
            for IX in range(self.NX):
                if self.I2NEXTX[IX, IY] > 0 and self.I2REGION[IX, IY] == self.REGIONTHIS:
                    JX = self.I2NEXTX[IX, IY]
                    JY = self.I2NEXTY[IX, IY]
                    NUPST[JX, JY] += 1

        ISEQ = 0
        for IY in range(self.NY):
            for IX in range(self.NX):
                if self.I2NEXTX[IX, IY] > 0 and self.I2REGION[IX, IY] == self.REGIONTHIS:
                    if NUPST[IX, IY] == UPNOW[IX, IY]:
                        ISEQ += 1
                        self.I1SEQX[ISEQ - 1] = IX
                        self.I1SEQY[ISEQ - 1] = IY
                        self.I2VECTOR[IX, IY] = ISEQ

        ISEQ1 = 1
        ISEQ2 = ISEQ
        AGAIN = 1
        while AGAIN == 1:
            AGAIN = 0
            JSEQ = ISEQ2
            for ISEQ in range(ISEQ1 - 1, ISEQ2):
                IX = self.I1SEQX[ISEQ]
                IY = self.I1SEQY[ISEQ]
                JX = self.I2NEXTX[IX, IY]
                JY = self.I2NEXTY[IX, IY]
                UPNOW[JX, JY] += 1
                if UPNOW[JX, JY] == NUPST[JX, JY] and self.I2NEXTX[JX, JY] > 0:
                    JSEQ += 1
                    self.I1SEQX[JSEQ - 1] = JX
                    self.I1SEQY[JSEQ - 1] = JY
                    self.I2VECTOR[JX, JY] = JSEQ
                    AGAIN = 1
            ISEQ1 = ISEQ2 + 1
            ISEQ2 = JSEQ

        self.NSEQRIV = JSEQ

        ISEQ = self.NSEQRIV
        for IY in range(self.NY):
            for IX in range(self.NX):
                if self.I2NEXTX[IX, IY] < 0 and self.I2NEXTX[IX, IY] != self.IMIS and self.I2REGION[IX, IY] == self.REGIONTHIS:
                    ISEQ += 1
                    self.I1SEQX[ISEQ - 1] = IX
                    self.I1SEQY[ISEQ - 1] = IY
                    self.I2VECTOR[IX, IY] = ISEQ

        self.NSEQALL = ISEQ

        for ISEQ in range(self.NSEQALL):
            IX = self.I1SEQX[ISEQ]
            IY = self.I1SEQY[ISEQ]
            if self.I2NEXTX[IX, IY] > 0:
                JX = self.I2NEXTX[IX, IY]
                JY = self.I2NEXTY[IX, IY]
                self.I1NEXT[ISEQ] = self.I2VECTOR[JX, JY]
            else:
                self.I1NEXT[ISEQ] = self.I2NEXTX[IX, IY]

    def READ_BIFPARAM(self):
        print("RIVMAP_INIT: Bifuraction channel:", self.CPTHOUT)

        with open(self.CPTHOUT, 'r') as file:
            self.NPTHOUT, self.NPTHLEV = map(int, file.readline().split())

            print("Bifurcation channel dimantion", self.NPTHOUT, self.NPTHLEV)

            self.PTH_UPST = np.zeros(self.NPTHOUT, dtype=int)
            self.PTH_DOWN = np.zeros(self.NPTHOUT, dtype=int)
            self.PTH_DST = np.zeros(self.NPTHOUT)
            self.PTH_ELV = np.zeros((self.NPTHOUT, self.NPTHLEV))
            self.PTH_WTH = np.zeros((self.NPTHOUT, self.NPTHLEV))
            self.PTH_MAN = np.zeros(self.NPTHLEV)

            NPTHOUT1 = 0
            for IPTH in range(self.NPTHOUT):
                data = list(map(float, file.readline().split()))
                IX, IY, JX, JY, self.PTH_DST[IPTH], PELV, PDPH = data[:7]
                self.PTH_WTH[IPTH, :] = data[7:]
                self.PTH_UPST[IPTH] = self.I2VECTOR[int(IX), int(IY)]
                self.PTH_DOWN[IPTH] = self.I2VECTOR[int(JX), int(JY)]
                if self.PTH_UPST[IPTH] > 0 and self.PTH_DOWN[IPTH] > 0:
                    NPTHOUT1 += 1
                for ILEV in range(self.NPTHLEV):
                    if ILEV == 0:
                        PWTH = self.PTH_WTH[IPTH, ILEV]
                        if PWTH > 0:
                            self.PTH_ELV[IPTH, ILEV] = PELV - PDPH
                        else:
                            self.PTH_ELV[IPTH, ILEV] = 1.E20
                    else:
                        PWTH = self.PTH_WTH[IPTH, ILEV]
                        if PWTH > 0:
                            self.PTH_ELV[IPTH, ILEV] = PELV + ILEV - 2.0
                        else:
                            self.PTH_ELV[IPTH, ILEV] = 1.E20

            for ILEV in range(self.NPTHLEV):
                if ILEV == 0:
                    self.PTH_MAN[ILEV] = self.PMANRIV
                else:
                    self.PTH_MAN[ILEV] = self.PMANFLD

            if self.NPTHOUT != NPTHOUT1:
                print("Bifuraction channel outside of domain. Only valid:", NPTHOUT1)

    def READ_TOPO_BIN(self):
        import struct

        def read_binary_file(filename, shape, dtype):
            with open(filename, 'rb') as file:
                data = file.read()
                return np.array(struct.unpack(dtype * np.prod(shape), data)).reshape(shape)

        print('TOPO_INIT: unit-catchment area : ', self.CGRAREA)
        self.D2GRAREA = read_binary_file(self.CGRAREA, (self.NX, self.NY), 'f')

        print('TOPO_INIT: ground elevation : ', self.CELEVTN)
        self.D2ELEVTN = read_binary_file(self.CELEVTN, (self.NX, self.NY), 'f')

        print('TOPO_INIT: downstream distance : ', self.CNXTDST)
        self.D2NXTDST = read_binary_file(self.CNXTDST, (self.NX, self.NY), 'f')

        print('TOPO_INIT: river channel length : ', self.CRIVLEN)
        self.D2RIVLEN = read_binary_file(self.CRIVLEN, (self.NX, self.NY), 'f')

        print('TOPO_INIT: floodplain elevation profile : ', self.CFLDHGT)
        self.D2FLDHGT = np.zeros((self.NX, self.NY, self.NLFP))
        for ILFP in range(self.NLFP):
            self.D2FLDHGT[:, :, ILFP] = read_binary_file(self.CFLDHGT, (self.NX, self.NY), 'f')

        print('TOPO_INIT: river channel depth : ', self.CRIVHGT)
        self.D2RIVHGT = read_binary_file(self.CRIVHGT, (self.NX, self.NY), 'f')

        print('TOPO_INIT: river channel width : ', self.CRIVWTH)
        self.D2RIVWTH = read_binary_file(self.CRIVWTH, (self.NX, self.NY), 'f')

        print('TOPO_INIT: manning coefficient river: ', self.CRIVMAN)
        self.D2RIVMAN = read_binary_file(self.CRIVMAN, (self.NX, self.NY), 'f')

        if self.LGDWDLY:
            print('TOPO_INIT: groundwater delay parameter: ', self.CGDWDLY)
            self.D2GDWDLY = read_binary_file(self.CGDWDLY, (self.NX, self.NY), 'f')

        if self.LMEANSL:
            print('TOPO_INIT: mean sea level: ', self.CMEANSL)
            self.D2MEANSL = read_binary_file(self.CMEANSL, (self.NX, self.NY), 'f')

    def READ_TOPO_CDF(self):
        import netCDF4 as nc

        def read_netcdf_variable(filename, varname):
            with nc.Dataset(filename, 'r') as ds:
                return ds.variables[varname][:]

        print('TOPO_INIT: ctmare:', self.CRIVCLINC)
        self.D2GRAREA = read_netcdf_variable(self.CRIVCLINC, 'ctmare')

        print('TOPO_INIT: elevtn:', self.CRIVCLINC)
        self.D2ELEVTN = read_netcdf_variable(self.CRIVCLINC, 'elevtn')

        print('TOPO_INIT: nxtdst:', self.CRIVCLINC)
        self.D2NXTDST = read_netcdf_variable(self.CRIVCLINC, 'nxtdst')

        print('TOPO_INIT: rivlen:', self.CRIVCLINC)
        self.D2RIVLEN = read_netcdf_variable(self.CRIVCLINC, 'rivlen')

        print('TOPO_INIT: fldhgt:', self.CRIVCLINC)
        self.D2FLDHGT = np.zeros((self.NX, self.NY, self.NLFP))
        for ILEV in range(self.NLFP):
            self.D2FLDHGT[:, :, ILEV] = read_netcdf_variable(self.CRIVCLINC, 'fldhgt')[:, :, ILEV]

        if self.LSLOPEMOUTH:
            self.D2ELEVSLOPE = np.zeros((self.NSEQMAX, 1))
            print('TOPO_INIT: elevslope:', self.CRIVPARNC)
            self.D2ELEVSLOPE = read_netcdf_variable(self.CRIVPARNC, 'elevslope')

        print('TOPO_INIT: rivwth:', self.CRIVPARNC)
        self.D2RIVWTH = read_netcdf_variable(self.CRIVPARNC, 'rivwth')

        print('TOPO_INIT: rivhgt:', self.CRIVPARNC)
        self.D2RIVHGT = read_netcdf_variable(self.CRIVPARNC, 'rivhgt')

        print('TOPO_INIT: rivman:', self.CRIVPARNC)
        self.D2RIVMAN = read_netcdf_variable(self.CRIVPARNC, 'rivman')

        if self.LGDWDLY:
            print('TOPO_INIT: GDWDLY:', self.CRIVPARNC)
            self.D2GDWDLY = read_netcdf_variable(self.CRIVPARNC, 'gdwdly')

        self.I2MASK = np.zeros((self.NSEQMAX, 1), dtype=int)
        if self.LSLPMIX:
            self.SET_SLOPEMIX()

        if self.LMEANSL:
            print('TOPO_INIT: rivhgt:', self.CMEANSLNC)
            self.D2MEANSL = read_netcdf_variable(self.CMEANSLNC, 'meansl')

    def SET_FLDSTG(self):
        self.P2FLDSTOMAX[:, :, :] = 0.0
        self.D2FLDGRD[:, :, :] = 0.0
        self.DFRCINC = 1.0 / self.NLFP

        for ISEQ in range(self.NSEQALL):
            DSTOPRE = self.P2RIVSTOMAX[ISEQ, 0]
            DHGTPRE = 0.0
            DWTHINC = self.D2GRAREA[ISEQ, 0] * self.D2RIVLEN[ISEQ, 0] * self.DFRCINC
            for I in range(self.NLFP):
                DSTONOW = self.D2RIVLEN[ISEQ, 0] * (self.D2RIVWTH[ISEQ, 0] + DWTHINC * (I + 0.5)) * (self.D2FLDHGT[ISEQ, 0, I] - DHGTPRE)
                self.P2FLDSTOMAX[ISEQ, 0, I] = DSTOPRE + DSTONOW
                self.D2FLDGRD[ISEQ, 0, I] = (self.D2FLDHGT[ISEQ, 0, I] - DHGTPRE) * (1.0 / DWTHINC)
                DSTOPRE = self.P2FLDSTOMAX[ISEQ, 0, I]
                DHGTPRE = self.D2FLDHGT[ISEQ, 0, I]

    def SET_SLOPEMIX(self):
        import netCDF4 as nc

        def read_netcdf_variable(filename, varname):
            with nc.Dataset(filename, 'r') as ds:
                return ds.variables[varname][:]

        self.I2MASK = read_netcdf_variable(self.CRIVPARNC, 'mask_slope')
        I0 = np.sum(self.I2MASK == 0)
        I1 = np.sum(self.I2MASK == 1)
        print('TOPO_INIT: sum(mask==0), sum(mask==1)', I0, I1)
        if I0 + I1 != self.NSEQALL:
            print('TOPO_INIT: mask==0 + mask == 1 does not match NSEQALL.. something wrong, aborting')
            raise ValueError('mask==0 + mask == 1 does not match NSEQALL')

    def CONV_ENDI(self, array, NX, NY):
        array.byteswap(inplace=True)
