import numpy as np

class CMF_CTRL_LEVEE_MOD:
    def __init__(self):
        self.CLEVHGT = "NONE"
        self.CLEVFRC = "NONE"
        self.D2LEVHGT = None
        self.D2LEVFRC = None
        self.D2BASHGT = None
        self.D2LEVDST = None
        self.D2LEVBASSTO = None
        self.D2LEVTOPSTO = None
        self.D2LEVFILSTO = None

    def CMF_LEVEE_NMLIST(self):
        print("CMF::LEVEE_NMLIST: namelist OPEN in unit:", self.CSETFILE, self.NSETFILE)
        self.CLEVHGT = "NONE"
        self.CLEVFRC = "NONE"
        print("=== NAMELIST, NLEVEE ===")
        print("CLEVHGT   :", self.CLEVHGT)
        print("CLEVFRC   :", self.CLEVFRC)
        print("CMF::LEVEE_NMLIST: end")

    def CMF_LEVEE_INIT(self):
        print("CMF::LEVEE_INIT: initialize levee")
        self.D2LEVHGT = np.zeros((self.NSEQMAX, 1))
        self.D2LEVFRC = np.zeros((self.NSEQMAX, 1))
        self.D2BASHGT = np.zeros((self.NSEQMAX, 1))
        self.D2LEVDST = np.zeros((self.NSEQMAX, 1))
        self.D2LEVBASSTO = np.zeros((self.NSEQMAX, 1))
        self.D2LEVTOPSTO = np.zeros((self.NSEQMAX, 1))
        self.D2LEVFILSTO = np.zeros((self.NSEQMAX, 1))

        print('INIT_LEVEE: levee crown height:', self.CLEVHGT)
        R2TEMP = self.read_binary_file(self.CLEVHGT)
        self.D2LEVHGT = self.mapR2vecD(R2TEMP)

        print('INIT_LEVEE: distance from levee to river:', self.CLEVFRC)
        R2TEMP = self.read_binary_file(self.CLEVFRC)
        self.D2LEVFRC = self.mapR2vecD(R2TEMP)

        self.calculate_levee_parameters()

    def calculate_levee_parameters(self):
        for ISEQ in range(self.NSEQALL):
            if self.D2LEVHGT[ISEQ, 0] <= 0.0:
                self.D2LEVHGT[ISEQ, 0] = 0.0
                self.D2LEVFRC[ISEQ, 0] = 1.0
            self.D2LEVFRC[ISEQ, 0] = max(0.0, min(1.0, self.D2LEVFRC[ISEQ, 0]))

        for ISEQ in range(self.NSEQALL):
            DSTOPRE = self.P2RIVSTOMAX[ISEQ, 0]
            DHGTPRE = 0.0
            DWTHINC = self.D2GRAREA[ISEQ, 0] * self.D2RIVLEN[ISEQ, 0]**(-1.0) * self.DFRCINC
            for I in range(self.NLFP):
                DSTONOW = self.D2RIVLEN[ISEQ, 0] * (self.D2RIVWTH[ISEQ, 0] + DWTHINC * (I + 0.5)) * (self.D2FLDHGT[ISEQ, I] - DHGTPRE)
                self.P2FLDSTOMAX[ISEQ, 0, I] = DSTOPRE + DSTONOW
                self.D2FLDGRD[ISEQ, 0, I] = (self.D2FLDHGT[ISEQ, I] - DHGTPRE) * DWTHINC**(-1.0)
                DSTOPRE = self.P2FLDSTOMAX[ISEQ, 0, I]
                DHGTPRE = self.D2FLDHGT[ISEQ, I]

            if self.D2LEVHGT[ISEQ, 0] == 0.0:
                self.D2BASHGT[ISEQ, 0] = 1.0e18
                self.D2LEVDST[ISEQ, 0] = 1.0e18
                self.D2LEVBASSTO[ISEQ, 0] = 1.0e18
                self.D2LEVTOPSTO[ISEQ, 0] = 1.0e18
                self.D2LEVFILSTO[ISEQ, 0] = 1.0e18
            else:
                self.calculate_levee_storage(ISEQ)

    def calculate_levee_storage(self, ISEQ):
        DSTOPRE = self.P2RIVSTOMAX[ISEQ, 0]
        DHGTPRE = 0.0
        DWTHPRE = 0.0
        self.D2LEVDST[ISEQ, 0] = self.D2LEVFRC[ISEQ, 0] * self.DWTHINC * self.NLFP

        ILEV = int(self.D2LEVFRC[ISEQ, 0] * self.NLFP) + 1
        if ILEV >= 2:
            DSTOPRE = self.P2FLDSTOMAX[ISEQ, 0, ILEV - 1]
            DHGTPRE = self.D2FLDHGT[ISEQ, ILEV - 1]
            DWTHPRE = self.DWTHINC * (ILEV - 1)

        if ILEV <= self.NLFP:
            DWTHNOW = self.D2LEVDST[ISEQ, 0] - DWTHPRE
            DHGTNOW = DWTHNOW * self.D2FLDGRD[ISEQ, 0, ILEV]
            self.D2BASHGT[ISEQ, 0] = DHGTNOW + DHGTPRE
            self.D2LEVHGT[ISEQ, 0] = max(self.D2LEVHGT[ISEQ, 0], self.D2BASHGT[ISEQ, 0])

            DSTONOW = (DWTHNOW * 0.5 + DWTHPRE + self.D2RIVWTH[ISEQ, 0]) * DHGTNOW * self.D2RIVLEN[ISEQ, 0]
            self.D2LEVBASSTO[ISEQ, 0] = DSTOPRE + DSTONOW

            DHGTDIF = self.D2LEVHGT[ISEQ, 0] - self.D2BASHGT[ISEQ, 0]
            self.D2LEVTOPSTO[ISEQ, 0] = self.D2LEVBASSTO[ISEQ, 0] + (self.D2LEVDST[ISEQ, 0] + self.D2RIVWTH[ISEQ, 0]) * DHGTDIF * self.D2RIVLEN[ISEQ, 0]
        else:
            self.D2BASHGT[ISEQ, 0] = DHGTPRE
            self.D2LEVHGT[ISEQ, 0] = max(self.D2LEVHGT[ISEQ, 0], self.D2BASHGT[ISEQ, 0])
            self.D2LEVBASSTO[ISEQ, 0] = DSTOPRE
            DHGTDIF = self.D2LEVHGT[ISEQ, 0] - self.D2BASHGT[ISEQ, 0]
            self.D2LEVTOPSTO[ISEQ, 0] = self.D2LEVBASSTO[ISEQ, 0] + (self.D2LEVDST[ISEQ, 0] + self.D2RIVWTH[ISEQ, 0]) * DHGTDIF * self.D2RIVLEN[ISEQ, 0]

        self.calculate_levee_fill_storage(ISEQ)

    def calculate_levee_fill_storage(self, ISEQ):
        I = 1
        DSTOPRE = self.P2RIVSTOMAX[ISEQ, 0]
        DWTHPRE = self.D2RIVWTH[ISEQ, 0]
        DHGTPRE = 0.0

        while self.D2LEVHGT[ISEQ, 0] > self.D2FLDHGT[ISEQ, I] and I <= self.NLFP:
            DSTOPRE = self.P2FLDSTOMAX[ISEQ, 0, I]
            DWTHPRE += self.DWTHINC
            DHGTPRE = self.D2FLDHGT[ISEQ, I]
            I += 1

        if I <= self.NLFP:
            DHGTNOW = self.D2LEVHGT[ISEQ, 0] - DHGTPRE
            DWTHNOW = DHGTNOW * self.D2FLDGRD[ISEQ, 0, I]**(-1.0)
            DSTONOW = (DWTHNOW * 0.5 + DWTHPRE) * DHGTNOW * self.D2RIVLEN[ISEQ, 0]
            self.D2LEVFILSTO[ISEQ, 0] = DSTOPRE + DSTONOW
        else:
            DHGTNOW = self.D2LEVHGT[ISEQ, 0] - DHGTPRE
            DSTONOW = DWTHPRE * DHGTNOW * self.D2RIVLEN[ISEQ, 0]
            self.D2LEVFILSTO[ISEQ, 0] = DSTOPRE + DSTONOW

    def CMF_LEVEE_FLDSTG(self):
        P0GLBSTOPRE2 = 0.0
        P0GLBSTONEW2 = 0.0
        P0GLBRIVSTO = 0.0
        P0GLBFLDSTO = 0.0
        P0GLBLEVSTO = 0.0
        P0GLBFLDARE = 0.0

        for ISEQ in range(self.NSEQALL):
            DSTOALL = self.P2RIVSTO[ISEQ, 0] + self.P2FLDSTO[ISEQ, 0] + self.P2LEVSTO[ISEQ, 0]
            DWTHINC = self.D2GRAREA[ISEQ, 0] * self.D2RIVLEN[ISEQ, 0]**(-1.0) * self.DFRCINC

            if DSTOALL > self.P2RIVSTOMAX[ISEQ, 0]:
                if DSTOALL < self.D2LEVBASSTO[ISEQ, 0]:
                    self.calculate_river_stage(ISEQ, DSTOALL, DWTHINC)
                elif DSTOALL < self.D2LEVTOPSTO[ISEQ, 0]:
                    self.calculate_river_levee_stage(ISEQ, DSTOALL)
                elif DSTOALL < self.D2LEVFILSTO[ISEQ, 0]:
                    self.calculate_river_levee_fill_stage(ISEQ, DSTOALL)
                else:
                    self.calculate_above_levee_stage(ISEQ, DSTOALL, DWTHINC)
            else:
                self.P2RIVSTO[ISEQ, 0] = DSTOALL
                self.D2RIVDPH[ISEQ, 0] = DSTOALL * self.D2RIVLEN[ISEQ, 0]**(-1.0) * self.D2RIVWTH[ISEQ, 0]**(-1.0)
                self.D2RIVDPH[ISEQ, 0] = max(self.D2RIVDPH[ISEQ, 0], 0.0)
                self.P2FLDSTO[ISEQ, 0] = 0.0
                self.D2FLDDPH[ISEQ, 0] = 0.0
                self.D2FLDFRC[ISEQ, 0] = 0.0
                self.D2FLDARE[ISEQ, 0] = 0.0
                self.P2LEVSTO[ISEQ, 0] = 0.0
                self.D2LEVDPH[ISEQ, 0] = 0.0

            self.D2SFCELV[ISEQ, 0] = self.D2RIVELV[ISEQ, 0] + self.D2RIVDPH[ISEQ, 0]
            P0GLBSTOPRE2 += DSTOALL
            P0GLBSTONEW2 += self.P2RIVSTO[ISEQ, 0] + self.P2FLDSTO[ISEQ, 0] + self.P2LEVSTO[ISEQ, 0]
            P0GLBRIVSTO += self.P2RIVSTO[ISEQ, 0]
            P0GLBFLDSTO += self.P2FLDSTO[ISEQ, 0]
            P0GLBLEVSTO += self.P2LEVSTO[ISEQ, 0]
            P0GLBFLDARE += self.D2FLDARE[ISEQ, 0]

    def calculate_river_stage(self, ISEQ, DSTOALL, DWTHINC):
        I = 1
        DSTOPRE = self.P2RIVSTOMAX[ISEQ, 0]
        DWTHPRE = self.D2RIVWTH[ISEQ, 0]
        DDPHPRE = 0.0

        while DSTOALL > self.P2FLDSTOMAX[ISEQ, 0, I] and I <= self.NLFP:
            DSTOPRE = self.P2FLDSTOMAX[ISEQ, 0, I]
            DWTHPRE += DWTHINC
            DDPHPRE += self.D2FLDGRD[ISEQ, 0, I] * DWTHINC
            I += 1

        if I <= self.NLFP:
            DSTONOW = DSTOALL - DSTOPRE
            DWTHNOW = -DWTHPRE + (DWTHPRE**2.0 + 2.0 * DSTONOW * self.D2RIVLEN[ISEQ, 0]**(-1.0) * self.D2FLDGRD[ISEQ, 0, I]**(-1.0))**0.5
            self.D2FLDDPH[ISEQ, 0] = DDPHPRE + self.D2FLDGRD[ISEQ, 0, I] * DWTHNOW
        else:
            DSTONOW = DSTOALL - DSTOPRE
            DWTHNOW = 0.0
            self.D2FLDDPH[ISEQ, 0] = DDPHPRE + DSTONOW * DWTHPRE**(-1.0) * self.D2RIVLEN[ISEQ, 0]**(-1.0)

        self.P2RIVSTO[ISEQ, 0] = self.P2RIVSTOMAX[ISEQ, 0] + self.D2RIVLEN[ISEQ, 0] * self.D2RIVWTH[ISEQ, 0] * self.D2FLDDPH[ISEQ, 0]
        self.D2RIVDPH[ISEQ, 0] = self.P2RIVSTO[ISEQ, 0] * self.D2RIVLEN[ISEQ, 0]**(-1.0) * self.D2RIVWTH[ISEQ, 0]**(-1.0)
        self.P2FLDSTO[ISEQ, 0] = DSTOALL - self.P2RIVSTO[ISEQ, 0]
        self.P2FLDSTO[ISEQ, 0] = max(self.P2FLDSTO[ISEQ, 0], 0.0)
        self.D2FLDFRC[ISEQ, 0] = (-self.D2RIVWTH[ISEQ, 0] + DWTHPRE + DWTHNOW) * (DWTHINC * self.NLFP)**(-1.0)
        self.D2FLDFRC[ISEQ, 0] = max(self.D2FLDFRC[ISEQ, 0], 0.0)
        self.D2FLDFRC[ISEQ, 0] = min(self.D2FLDFRC[ISEQ, 0], 1.0)
        self.D2FLDARE[ISEQ, 0] = self.D2GRAREA[ISEQ, 0] * self.D2FLDFRC[ISEQ, 0]
        self.P2LEVSTO[ISEQ, 0] = 0.0
        self.D2LEVDPH[ISEQ, 0] = 0.0

    def calculate_river_levee_stage(self, ISEQ, DSTOALL):
        DSTONOW = DSTOALL - self.D2LEVBASSTO[ISEQ, 0]
        DWTHNOW = self.D2LEVDST[ISEQ, 0] + self.D2RIVWTH[ISEQ, 0]
        self.D2FLDDPH[ISEQ, 0] = self.D2BASHGT[ISEQ, 0] + DSTONOW * DWTHNOW**(-1.0) * self.D2RIVLEN[ISEQ, 0]**(-1.0)
        self.P2RIVSTO[ISEQ, 0] = self.P2RIVSTOMAX[ISEQ, 0] + self.D2RIVLEN[ISEQ, 0] * self.D2RIVWTH[ISEQ, 0] * self.D2FLDDPH[ISEQ, 0]
        self.D2RIVDPH[ISEQ, 0] = self.P2RIVSTO[ISEQ, 0] * self.D2RIVLEN[ISEQ, 0]**(-1.0) * self.D2RIVWTH[ISEQ, 0]**(-1.0)
        self.P2FLDSTO[ISEQ, 0] = DSTOALL - self.P2RIVSTO[ISEQ, 0]
        self.P2FLDSTO[ISEQ, 0] = max(self.P2FLDSTO[ISEQ, 0], 0.0)
        self.D2FLDFRC[ISEQ, 0] = self.D2LEVFRC[ISEQ, 0]
        self.D2FLDARE[ISEQ, 0] = self.D2GRAREA[ISEQ, 0] * self.D2FLDFRC[ISEQ, 0]
        self.P2LEVSTO[ISEQ, 0] = 0.0
        self.D2LEVDPH[ISEQ, 0] = 0.0

    def calculate_river_levee_fill_stage(self, ISEQ, DSTOALL):
        self.D2FLDDPH[ISEQ, 0] = self.D2LEVHGT[ISEQ, 0]
        self.P2RIVSTO[ISEQ, 0] = self.P2RIVSTOMAX[ISEQ, 0] + self.D2RIVLEN[ISEQ, 0] * self.D2RIVWTH[ISEQ, 0] * self.D2FLDDPH[ISEQ, 0]
        self.D2RIVDPH[ISEQ, 0] = self.P2RIVSTO[ISEQ, 0] * self.D2RIVLEN[ISEQ, 0]**(-1.0) * self.D2RIVWTH[ISEQ, 0]**(-1.0)
        self.P2FLDSTO[ISEQ, 0] = self.D2LEVTOPSTO[ISEQ, 0] - self.P2RIVSTO[ISEQ, 0]
        self.P2FLDSTO[ISEQ, 0] = max(self.P2FLDSTO[ISEQ, 0], 0.0)
        self.P2LEVSTO[ISEQ, 0] = DSTOALL - self.P2RIVSTO[ISEQ, 0] - self.P2FLDSTO[ISEQ, 0]
        self.P2LEVSTO[ISEQ, 0] = max(self.P2LEVSTO[ISEQ, 0], 0.0)

        ILEV = int(self.D2LEVFRC[ISEQ, 0] * self.NLFP) + 1
        DSTOPRE = self.D2LEVTOPSTO[ISEQ, 0]
        DWTHPRE = 0.0
        DDPHPRE = 0.0

        I = ILEV
        while I <= self.NLFP:
            DSTOADD = (self.D2LEVDST[ISEQ, 0] + self.D2RIVWTH[ISEQ, 0]) * (self.D2LEVHGT[ISEQ, 0] - self.D2FLDHGT[ISEQ, I]) * self.D2RIVLEN[ISEQ, 0]
            if DSTOALL < self.P2FLDSTOMAX[ISEQ, 0, I] + DSTOADD:
                break
            DSTOPRE = self.P2FLDSTOMAX[ISEQ, 0, I] + DSTOADD
            DWTHPRE = self.DWTHINC * I - self.D2LEVDST[ISEQ, 0]
            DDPHPRE = self.D2FLDHGT[ISEQ, I] - self.D2BASHGT[ISEQ, 0]
            I += 1

        if I <= self.NLFP:
            DSTONOW = DSTOALL - DSTOPRE
            DWTHNOW = -DWTHPRE + (DWTHPRE**2.0 + 2.0 * DSTONOW * self.D2RIVLEN[ISEQ, 0]**(-1.0) * self.D2FLDGRD[ISEQ, 0, I]**(-1.0))**0.5
            DDPHNOW = DWTHNOW * self.D2FLDGRD[ISEQ, 0, I]
            self.D2LEVDPH[ISEQ, 0] = self.D2BASHGT[ISEQ, 0] + DDPHPRE + DDPHNOW
            self.D2FLDFRC[ISEQ, 0] = (DWTHPRE + self.D2LEVDST[ISEQ, 0]) * (self.DWTHINC * self.NLFP)**(-1.0)
            self.D2FLDFRC[ISEQ, 0] = max(self.D2FLDFRC[ISEQ, 0], 0.0)
            self.D2FLDFRC[ISEQ, 0] = min(self.D2FLDFRC[ISEQ, 0], 1.0)
            self.D2FLDARE[ISEQ, 0] = self.D2GRAREA[ISEQ, 0] * self.D2FLDFRC[ISEQ, 0]
        else:
            DSTONOW = DSTOALL - DSTOPRE
            DDPHNOW = DSTONOW * DWTHPRE**(-1.0) * self.D2RIVLEN[ISEQ, 0]**(-1.0)
            self.D2LEVDPH[ISEQ, 0] = self.D2BASHGT[ISEQ, 0] + DDPHPRE + DDPHNOW
            self.D2FLDFRC[ISEQ, 0] = 1.0
            self.D2FLDARE[ISEQ, 0] = self.D2GRAREA[ISEQ, 0] * self.D2FLDFRC[ISEQ, 0]

    def calculate_above_levee_stage(self, ISEQ, DSTOALL, DWTHINC):
        I = 1
        DSTOPRE = self.P2RIVSTOMAX[ISEQ, 0]
        DWTHPRE = self.D2RIVWTH[ISEQ, 0]
        DDPHPRE = 0.0

        while DSTOALL > self.P2FLDSTOMAX[ISEQ, 0, I] and I <= self.NLFP:
            DSTOPRE = self.P2FLDSTOMAX[ISEQ, 0, I]
            DWTHPRE += DWTHINC
            DDPHPRE += self.D2FLDGRD[ISEQ, 0, I] * DWTHINC
            I += 1

        if I <= self.NLFP:
            DSTONOW = DSTOALL - DSTOPRE
            DWTHNOW = -DWTHPRE + (DWTHPRE**2.0 + 2.0 * DSTONOW * self.D2RIVLEN[ISEQ, 0]**(-1.0) * self.D2FLDGRD[ISEQ, 0, I]**(-1.0))**0.5
            self.D2FLDDPH[ISEQ, 0] = DDPHPRE + self.D2FLDGRD[ISEQ, 0, I] * DWTHNOW
        else:
            DSTONOW = DSTOALL - DSTOPRE
            DWTHNOW = 0.0
            self.D2FLDDPH[ISEQ, 0] = DDPHPRE + DSTONOW * DWTHPRE**(-1.0) * self.D2RIVLEN[ISEQ, 0]**(-1.0)

        self.D2FLDFRC[ISEQ, 0] = (-self.D2RIVWTH[ISEQ, 0] + DWTHPRE + DWTHNOW) * (DWTHINC * self.NLFP)**(-1.0)
        self.D2FLDARE[ISEQ, 0] = self.D2GRAREA[ISEQ, 0] * self.D2FLDFRC[ISEQ, 0]
        self.P2RIVSTO[ISEQ, 0] = self.P2RIVSTOMAX[ISEQ, 0] + self.D2RIVLEN[ISEQ, 0] * self.D2RIVWTH[ISEQ, 0] * self.D2FLDDPH[ISEQ, 0]
        self.D2RIVDPH[ISEQ, 0] = self.P2RIVSTO[ISEQ, 0] * self.D2RIVLEN[ISEQ, 0]**(-1.0) * self.D2RIVWTH[ISEQ, 0]**(-1.0)
        DSTOADD = (self.D2FLDDPH[ISEQ, 0] - self.D2LEVHGT[ISEQ, 0]) * (self.D2LEVDST[ISEQ, 0] + self.D2RIVWTH[ISEQ, 0]) * self.D2RIVLEN[ISEQ, 0]
        self.P2FLDSTO[ISEQ, 0] = self.D2LEVTOPSTO[ISEQ, 0] + DSTOADD - self.P2RIVSTO[ISEQ, 0]
        self.P2FLDSTO[ISEQ, 0] = max(self.P2FLDSTO[ISEQ, 0], 0.0)
        self.P2LEVSTO[ISEQ, 0] = DSTOALL - self.P2RIVSTO[ISEQ, 0] - self.P2FLDSTO[ISEQ, 0]
        self.P2LEVSTO[ISEQ, 0] = max(self.P2LEVSTO[ISEQ, 0], 0.0)
        self.D2LEVDPH[ISEQ, 0] = self.D2FLDDPH[ISEQ, 0]

    def CMF_LEVEE_OPT_PTHOUT(self):
        D2SFCELV_LEV = np.zeros((self.NSEQMAX, 1))
        D2SFCELV_PRE = np.zeros((self.NSEQMAX, 1))

        for ISEQ in range(self.NSEQALL):
            if self.D2LEVFRC[ISEQ, 0] < 1.0:
                D2SFCELV_LEV[ISEQ, 0] = self.D2ELEVTN[ISEQ, 0] + self.D2LEVDPH[ISEQ, 0]
            else:
                D2SFCELV_LEV[ISEQ, 0] = self.D2SFCELV[ISEQ, 0]

            D2SFCELV_PRE[ISEQ, 0] = self.D2RIVELV[ISEQ, 0] + self.D2RIVDPH_PRE[ISEQ, 0]

        self.D1PTHFLW.fill(self.DMIS)

        for IPTH in range(self.NPTHOUT):
            ISEQP = self.PTH_UPST[IPTH]
            JSEQP = self.PTH_DOWN[IPTH]

            if ISEQP <= 0 or JSEQP <= 0:
                continue

            if self.I2MASK[ISEQP, 0] == 1 or self.I2MASK[JSEQP, 0] == 1:
                continue

            DSLOPE = (self.D2SFCELV[ISEQP, 0] - self.D2SFCELV[JSEQP, 0]) * self.PTH_DST[IPTH]**(-1.0)
            DSLOPE = max(-0.005, min(0.005, DSLOPE))

            ILEV = 1
            DFLW = max(self.D2SFCELV[ISEQP, 0], self.D2SFCELV[JSEQP, 0]) - self.PTH_ELV[IPTH, ILEV]
            DFLW = max(DFLW, 0.0)
            DFLW_PRE = max(self.D2SFCELV_PRE[ISEQP, 0], self.D2SFCELV_PRE[JSEQP, 0]) - self.PTH_ELV[IPTH, ILEV]
            DFLW_PRE = max(DFLW_PRE, 0.0)
            DFLW_IMP = (DFLW * DFLW_PRE)**0.5
            if DFLW_IMP <= 0.0:
                DFLW_IMP = DFLW

            if DFLW_IMP > 1.0e-5:
                DOUT_PRE = self.D1PTHFLW_PRE[IPTH, ILEV] * self.PTH_WTH[IPTH, ILEV]**(-1.0)
                self.D1PTHFLW[IPTH, ILEV] = self.PTH_WTH[IPTH, ILEV] * (DOUT_PRE + self.PGRV * self.DT * DFLW_IMP * DSLOPE) * (1.0 + self.PGRV * self.DT * self.PTH_MAN[ILEV]**2.0 * abs(DOUT_PRE) * DFLW_IMP**(-7.0 / 3.0))**(-1.0)
            else:
                self.D1PTHFLW[IPTH, ILEV] = 0.0

            if self.NPTHLEV <= 1:
                continue

            DSLOPE = (D2SFCELV_LEV[ISEQP, 0] - D2SFCELV_LEV[JSEQP, 0]) * self.PTH_DST[IPTH]**(-1.0)
            DSLOPE = max(-0.005, min(0.005, DSLOPE))

            for ILEV in range(2, self.NPTHLEV + 1):
                DFLW = max(D2SFCELV_LEV[ISEQP, 0], D2SFCELV_LEV[JSEQP, 0]) - self.PTH_ELV[IPTH, ILEV]
                DFLW = max(DFLW, 0.0)
                DFLW_IMP = DFLW
                if DFLW_IMP > 1.0e-5:
                    DOUT_PRE = self.D1PTHFLW_PRE[IPTH, ILEV] * self.PTH_WTH[IPTH, ILEV]**(-1.0)
                    self.D1PTHFLW[IPTH, ILEV] = self.PTH_WTH[IPTH, ILEV] * (DOUT_PRE + self.PGRV * self.DT * DFLW_IMP * DSLOPE) * (1.0 + self.PGRV * self.DT * self.PTH_MAN[ILEV]**2.0 * abs(DOUT_PRE) * DFLW_IMP**(-7.0 / 3.0))**(-1.0)
                else:
                    self.D1PTHFLW[IPTH, ILEV] = 0.0

    def read_binary_file(self, filename):
        # Implement the read_binary_file function here
        pass

    def mapR2vecD(self, r2tmp):
        # Implement the mapR2vecD function here
        pass
