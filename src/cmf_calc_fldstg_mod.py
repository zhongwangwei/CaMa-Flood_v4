import numpy as np

class CMF_CALC_FLDSTG_MOD:
    def __init__(self, NLFP, NSEQALL, D2GRAREA, D2RIVLEN, D2RIVWTH, D2RIVELV, P2RIVSTOMAX, P2FLDSTOMAX, D2FLDGRD, DFRCINC, P2RIVSTO, P2FLDSTO, D2RIVDPH, D2FLDDPH, D2FLDFRC, D2FLDARE, D2SFCELV, P0GLBSTOPRE2, P0GLBSTONEW2, P0GLBRIVSTO, P0GLBFLDSTO, P0GLBFLDARE):
        self.NLFP = NLFP
        self.NSEQALL = NSEQALL
        self.D2GRAREA = D2GRAREA
        self.D2RIVLEN = D2RIVLEN
        self.D2RIVWTH = D2RIVWTH
        self.D2RIVELV = D2RIVELV
        self.P2RIVSTOMAX = P2RIVSTOMAX
        self.P2FLDSTOMAX = P2FLDSTOMAX
        self.D2FLDGRD = D2FLDGRD
        self.DFRCINC = DFRCINC
        self.P2RIVSTO = P2RIVSTO
        self.P2FLDSTO = P2FLDSTO
        self.D2RIVDPH = D2RIVDPH
        self.D2FLDDPH = D2FLDDPH
        self.D2FLDFRC = D2FLDFRC
        self.D2FLDARE = D2FLDARE
        self.D2SFCELV = D2SFCELV
        self.P0GLBSTOPRE2 = P0GLBSTOPRE2
        self.P0GLBSTONEW2 = P0GLBSTONEW2
        self.P0GLBRIVSTO = P0GLBRIVSTO
        self.P0GLBFLDSTO = P0GLBFLDSTO
        self.P0GLBFLDARE = P0GLBFLDARE

    def CMF_CALC_FLDSTG_DEF(self):
        self.P0GLBSTOPRE2 = 0.0
        self.P0GLBSTONEW2 = 0.0
        self.P0GLBRIVSTO = 0.0
        self.P0GLBFLDSTO = 0.0
        self.P0GLBFLDARE = 0.0

        for ISEQ in range(self.NSEQALL):
            DSTOALL = self.P2RIVSTO[ISEQ, 0] + self.P2FLDSTO[ISEQ, 0]

            if DSTOALL > self.P2RIVSTOMAX[ISEQ, 0]:
                I = 0
                DSTOPRE = self.P2RIVSTOMAX[ISEQ, 0]
                DWTHPRE = self.D2RIVWTH[ISEQ, 0]
                DDPHPRE = 0.0
                DWTHINC = self.D2GRAREA[ISEQ, 0] * self.D2RIVLEN[ISEQ, 0]**(-1.0) * self.DFRCINC
                while DSTOALL > self.P2FLDSTOMAX[ISEQ, 0, I] and I < self.NLFP:
                    DSTOPRE = self.P2FLDSTOMAX[ISEQ, 0, I]
                    DWTHPRE += DWTHINC
                    DDPHPRE += self.D2FLDGRD[ISEQ, 0, I] * DWTHINC
                    I += 1
                    if I >= self.NLFP:
                        break
                if I >= self.NLFP:
                    DSTONOW = DSTOALL - DSTOPRE
                    DWTHNOW = 0.0
                    self.D2FLDDPH[ISEQ, 0] = DDPHPRE + DSTONOW * DWTHPRE**(-1.0) * self.D2RIVLEN[ISEQ, 0]**(-1.0)
                else:
                    DSTONOW = DSTOALL - DSTOPRE
                    DWTHNOW = -DWTHPRE + (DWTHPRE**2.0 + 2.0 * DSTONOW * self.D2RIVLEN[ISEQ, 0]**(-1.0) * self.D2FLDGRD[ISEQ, 0, I]**(-1.0))**0.5
                    self.D2FLDDPH[ISEQ, 0] = DDPHPRE + self.D2FLDGRD[ISEQ, 0, I] * DWTHNOW
                self.P2RIVSTO[ISEQ, 0] = self.P2RIVSTOMAX[ISEQ, 0] + self.D2RIVLEN[ISEQ, 0] * self.D2RIVWTH[ISEQ, 0] * self.D2FLDDPH[ISEQ, 0]
                self.P2RIVSTO[ISEQ, 0] = min(self.P2RIVSTO[ISEQ, 0], DSTOALL)
                self.D2RIVDPH[ISEQ, 0] = self.P2RIVSTO[ISEQ, 0] * self.D2RIVLEN[ISEQ, 0]**(-1.0) * self.D2RIVWTH[ISEQ, 0]**(-1.0)
                self.P2FLDSTO[ISEQ, 0] = DSTOALL - self.P2RIVSTO[ISEQ, 0]
                self.P2FLDSTO[ISEQ, 0] = max(self.P2FLDSTO[ISEQ, 0], 0.0)
                self.D2FLDFRC[ISEQ, 0] = (-self.D2RIVWTH[ISEQ, 0] + DWTHPRE + DWTHNOW) * (DWTHINC * self.NLFP)**(-1.0)
                self.D2FLDFRC[ISEQ, 0] = max(self.D2FLDFRC[ISEQ, 0], 0.0)
                self.D2FLDFRC[ISEQ, 0] = min(self.D2FLDFRC[ISEQ, 0], 1.0)
                self.D2FLDARE[ISEQ, 0] = self.D2GRAREA[ISEQ, 0] * self.D2FLDFRC[ISEQ, 0]
            else:
                self.P2RIVSTO[ISEQ, 0] = DSTOALL
                self.D2RIVDPH[ISEQ, 0] = DSTOALL * self.D2RIVLEN[ISEQ, 0]**(-1.0) * self.D2RIVWTH[ISEQ, 0]**(-1.0)
                self.D2RIVDPH[ISEQ, 0] = max(self.D2RIVDPH[ISEQ, 0], 0.0)
                self.P2FLDSTO[ISEQ, 0] = 0.0
                self.D2FLDDPH[ISEQ, 0] = 0.0
                self.D2FLDFRC[ISEQ, 0] = 0.0
                self.D2FLDARE[ISEQ, 0] = 0.0
            self.D2SFCELV[ISEQ, 0] = self.D2RIVELV[ISEQ, 0] + self.D2RIVDPH[ISEQ, 0]
            self.P0GLBSTOPRE2 += DSTOALL
            self.P0GLBSTONEW2 += self.P2RIVSTO[ISEQ, 0] + self.P2FLDSTO[ISEQ, 0]
            self.P0GLBRIVSTO += self.P2RIVSTO[ISEQ, 0]
            self.P0GLBFLDSTO += self.P2FLDSTO[ISEQ, 0]
            self.P0GLBFLDARE += self.D2FLDARE[ISEQ, 0]

    def CMF_OPT_FLDSTG_ES(self):
        D2STODWN = np.zeros((self.NSEQALL, 1))
        D2WTHPRE = np.zeros((self.NSEQALL, 1))
        D2WTHINC = np.zeros((self.NSEQALL, 1))

        self.P0GLBRIVSTO = 0.0
        self.P0GLBFLDSTO = 0.0
        self.P0GLBFLDARE = 0.0
        self.P0GLBSTOPRE2 = 0.0
        self.P0GLBSTONEW2 = 0.0

        for ISEQ in range(self.NSEQALL):
            DSTOALL = self.P2RIVSTO[ISEQ, 0] + self.P2FLDSTO[ISEQ, 0]
            self.P2RIVSTO[ISEQ, 0] = DSTOALL
            self.D2RIVDPH[ISEQ, 0] = DSTOALL * self.D2RIVLEN[ISEQ, 0]**(-1.0) * self.D2RIVWTH[ISEQ, 0]**(-1.0)
            self.D2RIVDPH[ISEQ, 0] = max(self.D2RIVDPH[ISEQ, 0], 0.0)
            self.P2FLDSTO[ISEQ, 0] = 0.0
            self.D2FLDDPH[ISEQ, 0] = 0.0
            self.D2FLDFRC[ISEQ, 0] = 0.0
            self.D2FLDARE[ISEQ, 0] = 0.0
            D2STODWN[ISEQ, 0] = self.P2RIVSTOMAX[ISEQ, 0]
            D2WTHPRE[ISEQ, 0] = self.D2RIVWTH[ISEQ, 0]
            D2WTHINC[ISEQ, 0] = self.D2GRAREA[ISEQ, 0] * self.D2RIVLEN[ISEQ, 0]**(-1.0) * self.DFRCINC
            self.P0GLBSTOPRE2 += DSTOALL

        for I in range(self.NLFP):
            for ISEQ in range(self.NSEQALL):
                DSTOALL = self.P2RIVSTO[ISEQ, 0] + self.P2FLDSTO[ISEQ, 0]
                DSTONOW = DSTOALL - D2STODWN[ISEQ, 0]
                DSTONOW = max(DSTONOW, 0.0)
                DWTHNOW = -D2WTHPRE[ISEQ, 0] + (D2WTHPRE[ISEQ, 0]**2.0 + 2.0 * DSTONOW * self.D2RIVLEN[ISEQ, 0]**(-1.0) * self.D2FLDGRD[ISEQ, 0, I]**(-1.0))**0.5
                DWTHNOW = min(DWTHNOW, D2WTHINC[ISEQ, 0])
                DWTHNOW = max(DWTHNOW, 0.0)
                self.D2FLDDPH[ISEQ, 0] += self.D2FLDGRD[ISEQ, 0, I] * DWTHNOW
                self.D2FLDFRC[ISEQ, 0] += DWTHNOW / D2WTHINC[ISEQ, 0] * self.NLFP**(-1.0)
                D2STODWN[ISEQ, 0] = self.P2FLDSTOMAX[ISEQ, 0, I]
                D2WTHPRE[ISEQ, 0] += D2WTHINC[ISEQ, 0]

        for ISEQ in range(self.NSEQALL):
            DSTOALL = self.P2RIVSTO[ISEQ, 0] + self.P2FLDSTO[ISEQ, 0]
            DSTONOW = DSTOALL - D2STODWN[ISEQ, 0]
            DSTONOW = max(DSTONOW, 0.0)
            self.D2FLDDPH[ISEQ, 0] += DSTONOW * D2WTHPRE[ISEQ, 0]**(-1.0) * self.D2RIVLEN[ISEQ, 0]**(-1.0)

        for ISEQ in range(self.NSEQALL):
            if self.D2FLDDPH[ISEQ, 0] > 1.0e-5:
                DSTOALL = self.P2RIVSTO[ISEQ, 0] + self.P2FLDSTO[ISEQ, 0]
                self.D2RIVDPH[ISEQ, 0] = self.D2RIVHGT[ISEQ, 0] + self.D2FLDDPH[ISEQ, 0]
                self.P2RIVSTO[ISEQ, 0] = self.D2RIVLEN[ISEQ, 0] * self.D2RIVWTH[ISEQ, 0] * self.D2RIVDPH[ISEQ, 0]
                self.P2RIVSTO[ISEQ, 0] = min(self.P2RIVSTO[ISEQ, 0], DSTOALL)
                self.P2FLDSTO[ISEQ, 0] = DSTOALL - self.P2RIVSTO[ISEQ, 0]
                self.P2FLDSTO[ISEQ, 0] = max(self.P2FLDSTO[ISEQ, 0], 0.0)
                self.D2FLDFRC[ISEQ, 0] = max(self.D2FLDFRC[ISEQ, 0], 0.0)
                self.D2FLDFRC[ISEQ, 0] = min(self.D2FLDFRC[ISEQ, 0], 1.0)
                self.D2FLDARE[ISEQ, 0] = self.D2GRAREA[ISEQ, 0] * self.D2FLDFRC[ISEQ, 0]

        for ISEQ in range(self.NSEQALL):
            self.D2SFCELV[ISEQ, 0] = self.D2RIVELV[ISEQ, 0] + self.D2RIVDPH[ISEQ, 0]
            self.P0GLBSTONEW2 += self.P2RIVSTO[ISEQ, 0] + self.P2FLDSTO[ISEQ, 0]
            self.P0GLBRIVSTO += self.P2RIVSTO[ISEQ, 0]
            self.P0GLBFLDSTO += self.P2FLDSTO[ISEQ, 0]
            self.P0GLBFLDARE += self.D2FLDARE[ISEQ, 0]
