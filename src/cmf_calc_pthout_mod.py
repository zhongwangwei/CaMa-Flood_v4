import numpy as np

class CMF_CALC_PTHOUT_MOD:
    def __init__(self, DT, PGRV, DMIS, NSEQALL, NSEQMAX, NPTHOUT, NPTHLEV, PTH_UPST, PTH_DOWN, PTH_DST, PTH_ELV, PTH_WTH, PTH_MAN, I2MASK, D2RIVELV, D1PTHFLW, D1PTHFLW_PRE, D2RIVDPH_PRE, D2SFCELV, D1PTHFLWSUM):
        self.DT = DT
        self.PGRV = PGRV
        self.DMIS = DMIS
        self.NSEQALL = NSEQALL
        self.NSEQMAX = NSEQMAX
        self.NPTHOUT = NPTHOUT
        self.NPTHLEV = NPTHLEV
        self.PTH_UPST = PTH_UPST
        self.PTH_DOWN = PTH_DOWN
        self.PTH_DST = PTH_DST
        self.PTH_ELV = PTH_ELV
        self.PTH_WTH = PTH_WTH
        self.PTH_MAN = PTH_MAN
        self.I2MASK = I2MASK
        self.D2RIVELV = D2RIVELV
        self.D1PTHFLW = D1PTHFLW
        self.D1PTHFLW_PRE = D1PTHFLW_PRE
        self.D2RIVDPH_PRE = D2RIVDPH_PRE
        self.D2SFCELV = D2SFCELV
        self.D1PTHFLWSUM = D1PTHFLWSUM

    def CMF_CALC_PTHOUT(self):
        D2SFCELV_PRE = np.zeros((self.NSEQMAX, 1))

        for ISEQ in range(self.NSEQALL):
            D2SFCELV_PRE[ISEQ, 0] = self.D2RIVELV[ISEQ, 0] + self.D2RIVDPH_PRE[ISEQ, 0]

        self.D1PTHFLW.fill(0.0)

        for IPTH in range(self.NPTHOUT):
            ISEQP = self.PTH_UPST[IPTH]
            JSEQP = self.PTH_DOWN[IPTH]

            if ISEQP <= 0 or JSEQP <= 0:
                continue

            if self.I2MASK[ISEQP, 0] > 0 or self.I2MASK[JSEQP, 0] > 0:
                continue

            DSLOPE = (self.D2SFCELV[ISEQP, 0] - self.D2SFCELV[JSEQP, 0]) * self.PTH_DST[IPTH]**(-1.0)
            DSLOPE = max(-0.005, min(0.005, DSLOPE))

            for ILEV in range(self.NPTHLEV):
                DFLW = max(max(self.D2SFCELV[ISEQP, 0], self.D2SFCELV[JSEQP, 0]) - self.PTH_ELV[IPTH, ILEV], 0.0)
                DFLW_PRE = max(max(D2SFCELV_PRE[ISEQP, 0], D2SFCELV_PRE[JSEQP, 0]) - self.PTH_ELV[IPTH, ILEV], 0.0)
                DFLW_IMP = max((DFLW * DFLW_PRE)**0.5, (DFLW * 0.01)**0.5)

                if DFLW_IMP > 1.0e-5:
                    DOUT_PRE = self.D1PTHFLW_PRE[IPTH, ILEV] * self.PTH_WTH[IPTH, ILEV]**(-1.0)
                    self.D1PTHFLW[IPTH, ILEV] = self.PTH_WTH[IPTH, ILEV] * (DOUT_PRE + self.PGRV * self.DT * DFLW_IMP * DSLOPE) * (1.0 + self.PGRV * self.DT * self.PTH_MAN[ILEV]**2.0 * abs(DOUT_PRE) * DFLW_IMP**(-7.0 / 3.0))**(-1.0)
                else:
                    self.D1PTHFLW[IPTH, ILEV] = 0.0

        self.D1PTHFLWSUM.fill(0.0)

        for ILEV in range(self.NPTHLEV):
            self.D1PTHFLWSUM += self.D1PTHFLW[:, ILEV]
