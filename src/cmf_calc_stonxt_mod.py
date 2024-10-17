import numpy as np

class CMF_CALC_STONXT_MOD:
    def __init__(self, LGDWDLY, DT, LWEVAP, NSEQALL, D2RIVOUT, D2FLDOUT, P2RIVSTO, P2FLDSTO, D2RUNOFF, P2GDWSTO, D2GDWRTN, D2ROFSUB, D2WEVAP, D2RIVINF, D2FLDINF, D2PTHOUT, D2FLDFRC, D2OUTFLW, D2STORGE, D2WEVAPEX, P0GLBSTOPRE, P0GLBSTONXT, P0GLBSTONEW, P0GLBRIVINF, P0GLBRIVOUT):
        self.LGDWDLY = LGDWDLY
        self.DT = DT
        self.LWEVAP = LWEVAP
        self.NSEQALL = NSEQALL
        self.D2RIVOUT = D2RIVOUT
        self.D2FLDOUT = D2FLDOUT
        self.P2RIVSTO = P2RIVSTO
        self.P2FLDSTO = P2FLDSTO
        self.D2RUNOFF = D2RUNOFF
        self.P2GDWSTO = P2GDWSTO
        self.D2GDWRTN = D2GDWRTN
        self.D2ROFSUB = D2ROFSUB
        self.D2WEVAP = D2WEVAP
        self.D2RIVINF = D2RIVINF
        self.D2FLDINF = D2FLDINF
        self.D2PTHOUT = D2PTHOUT
        self.D2FLDFRC = D2FLDFRC
        self.D2OUTFLW = D2OUTFLW
        self.D2STORGE = D2STORGE
        self.D2WEVAPEX = D2WEVAPEX
        self.P0GLBSTOPRE = P0GLBSTOPRE
        self.P0GLBSTONXT = P0GLBSTONXT
        self.P0GLBSTONEW = P0GLBSTONEW
        self.P0GLBRIVINF = P0GLBRIVINF
        self.P0GLBRIVOUT = P0GLBRIVOUT

    def CMF_CALC_STONXT(self):
        if self.LGDWDLY:
            self.CALC_GDWDLY()
        else:
            self.P2GDWSTO.fill(0.0)
            self.D2GDWRTN = self.D2ROFSUB.copy()

        self.P0GLBSTOPRE = 0.0
        self.P0GLBSTONXT = 0.0
        self.P0GLBSTONEW = 0.0
        self.P0GLBRIVINF = 0.0
        self.P0GLBRIVOUT = 0.0

        for ISEQ in range(self.NSEQALL):
            self.P0GLBSTOPRE += self.P2RIVSTO[ISEQ, 0] + self.P2FLDSTO[ISEQ, 0]
            self.P0GLBRIVINF += self.D2RIVINF[ISEQ, 0] * self.DT + self.D2FLDINF[ISEQ, 0] * self.DT
            self.P0GLBRIVOUT += self.D2RIVOUT[ISEQ, 0] * self.DT + self.D2FLDOUT[ISEQ, 0] * self.DT + self.D2PTHOUT[ISEQ, 0] * self.DT

            self.P2RIVSTO[ISEQ, 0] += self.D2RIVINF[ISEQ, 0] * self.DT - self.D2RIVOUT[ISEQ, 0] * self.DT
            if self.P2RIVSTO[ISEQ, 0] < 0.0:
                self.P2FLDSTO[ISEQ, 0] += self.P2RIVSTO[ISEQ, 0]
                self.P2RIVSTO[ISEQ, 0] = 0.0

            self.P2FLDSTO[ISEQ, 0] += self.D2FLDINF[ISEQ, 0] * self.DT - self.D2FLDOUT[ISEQ, 0] * self.DT - self.D2PTHOUT[ISEQ, 0] * self.DT
            if self.P2FLDSTO[ISEQ, 0] < 0.0:
                self.P2RIVSTO[ISEQ, 0] = max(self.P2RIVSTO[ISEQ, 0] + self.P2FLDSTO[ISEQ, 0], 0.0)
                self.P2FLDSTO[ISEQ, 0] = 0.0

            self.P0GLBSTONXT += self.P2RIVSTO[ISEQ, 0] + self.P2FLDSTO[ISEQ, 0]
            self.D2OUTFLW[ISEQ, 0] = self.D2RIVOUT[ISEQ, 0] + self.D2FLDOUT[ISEQ, 0]

            DRIVROF = (self.D2RUNOFF[ISEQ, 0] + self.D2GDWRTN[ISEQ, 0]) * (1.0 - self.D2FLDFRC[ISEQ, 0]) * self.DT
            DFLDROF = (self.D2RUNOFF[ISEQ, 0] + self.D2GDWRTN[ISEQ, 0]) * self.D2FLDFRC[ISEQ, 0] * self.DT
            self.P2RIVSTO[ISEQ, 0] += DRIVROF
            self.P2FLDSTO[ISEQ, 0] += DFLDROF

            if self.LWEVAP:
                DWEVAPEX = min(self.P2FLDSTO[ISEQ, 0], self.D2FLDFRC[ISEQ, 0] * self.DT * self.D2WEVAP[ISEQ, 0])
                self.P2FLDSTO[ISEQ, 0] -= DWEVAPEX
                self.D2WEVAPEX[ISEQ, 0] = DWEVAPEX / self.DT

            self.D2STORGE[ISEQ, 0] = self.P2RIVSTO[ISEQ, 0] + self.P2FLDSTO[ISEQ, 0]
            self.P0GLBSTONEW += self.D2STORGE[ISEQ, 0]

    def CALC_GDWDLY(self):
        ZDTI = 1.0 / self.DT
        for ISEQ in range(self.NSEQALL):
            if self.D2GDWDLY[ISEQ, 0] > 0.0:
                ZMULGW = 1.0 / (ZDTI + 1.0 / self.D2GDWDLY[ISEQ, 0])
                self.P2GDWSTO[ISEQ, 0] = (self.D2ROFSUB[ISEQ, 0] + self.P2GDWSTO[ISEQ, 0] * ZDTI) * ZMULGW
                self.D2GDWRTN[ISEQ, 0] = self.P2GDWSTO[ISEQ, 0] / self.D2GDWDLY[ISEQ, 0]
            else:
                self.P2GDWSTO[ISEQ, 0] = 0.0
                self.D2GDWRTN[ISEQ, 0] = self.D2ROFSUB[ISEQ, 0]
