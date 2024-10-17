import numpy as np

class CMF_CALC_OUTFLW_MOD:
    def __init__(self, DT, PDSTMTH, PMANFLD, PGRV, LFLDOUT, LPTHOUT, LSLOPEMOUTH, I1NEXT, NSEQALL, NSEQRIV, NSEQMAX, D2RIVELV, D2ELEVTN, D2NXTDST, D2RIVWTH, D2RIVHGT, D2RIVLEN, D2RIVMAN, D2DWNELV, D2ELEVSLOPE, P2RIVSTO, D2RIVOUT, P2FLDSTO, D2FLDOUT, D2RIVOUT_PRE, D2RIVDPH_PRE, D2FLDOUT_PRE, D2FLDSTO_PRE, D2RIVDPH, D2RIVVEL, D2RIVINF, D2FLDDPH, D2FLDINF, D2SFCELV):
        self.DT = DT
        self.PDSTMTH = PDSTMTH
        self.PMANFLD = PMANFLD
        self.PGRV = PGRV
        self.LFLDOUT = LFLDOUT
        self.LPTHOUT = LPTHOUT
        self.LSLOPEMOUTH = LSLOPEMOUTH
        self.I1NEXT = I1NEXT
        self.NSEQALL = NSEQALL
        self.NSEQRIV = NSEQRIV
        self.NSEQMAX = NSEQMAX
        self.D2RIVELV = D2RIVELV
        self.D2ELEVTN = D2ELEVTN
        self.D2NXTDST = D2NXTDST
        self.D2RIVWTH = D2RIVWTH
        self.D2RIVHGT = D2RIVHGT
        self.D2RIVLEN = D2RIVLEN
        self.D2RIVMAN = D2RIVMAN
        self.D2DWNELV = D2DWNELV
        self.D2ELEVSLOPE = D2ELEVSLOPE
        self.P2RIVSTO = P2RIVSTO
        self.D2RIVOUT = D2RIVOUT
        self.P2FLDSTO = P2FLDSTO
        self.D2FLDOUT = D2FLDOUT
        self.D2RIVOUT_PRE = D2RIVOUT_PRE
        self.D2RIVDPH_PRE = D2RIVDPH_PRE
        self.D2FLDOUT_PRE = D2FLDOUT_PRE
        self.D2FLDSTO_PRE = D2FLDSTO_PRE
        self.D2RIVDPH = D2RIVDPH
        self.D2RIVVEL = D2RIVVEL
        self.D2RIVINF = D2RIVINF
        self.D2FLDDPH = D2FLDDPH
        self.D2FLDINF = D2FLDINF
        self.D2SFCELV = D2SFCELV

    def CMF_CALC_OUTFLW(self):
        D2SFCELV_PRE = np.zeros((self.NSEQMAX, 1))
        D2FLDDPH_PRE = np.zeros((self.NSEQMAX, 1))

        for ISEQ in range(self.NSEQALL):
            self.D2SFCELV[ISEQ, 0] = self.D2RIVELV[ISEQ, 0] + self.D2RIVDPH[ISEQ, 0]
            D2SFCELV_PRE[ISEQ, 0] = self.D2RIVELV[ISEQ, 0] + self.D2RIVDPH_PRE[ISEQ, 0]
            D2FLDDPH_PRE[ISEQ, 0] = max(self.D2RIVDPH_PRE[ISEQ, 0] - self.D2RIVHGT[ISEQ, 0], 0.0)

        for ISEQ in range(self.NSEQRIV):
            JSEQ = self.I1NEXT[ISEQ]
            DSFCMAX = max(self.D2SFCELV[ISEQ, 0], self.D2SFCELV[JSEQ, 0])
            DSFCMAX_PRE = max(D2SFCELV_PRE[ISEQ, 0], D2SFCELV_PRE[JSEQ, 0])
            DSLOPE = (self.D2SFCELV[ISEQ, 0] - self.D2SFCELV[JSEQ, 0]) * self.D2NXTDST[ISEQ, 0]**(-1.0)
            DSLOPE_F = max(-0.005, min(0.005, DSLOPE))

            DFLW = DSFCMAX - self.D2RIVELV[ISEQ, 0]
            DAREA = self.D2RIVWTH[ISEQ, 0] * DFLW

            DFLW_PRE = DSFCMAX_PRE - self.D2RIVELV[ISEQ, 0]
            DFLW_IMP = max((DFLW * DFLW_PRE)**0.5, 1.0e-6)

            if DFLW_IMP > 1.0e-5 and DAREA > 1.0e-5:
                DOUT_PRE = self.D2RIVOUT_PRE[ISEQ, 0] * self.D2RIVWTH[ISEQ, 0]**(-1.0)
                self.D2RIVOUT[ISEQ, 0] = self.D2RIVWTH[ISEQ, 0] * (DOUT_PRE + self.PGRV * self.DT * DFLW_IMP * DSLOPE) * (1.0 + self.PGRV * self.DT * self.D2RIVMAN[ISEQ, 0]**2.0 * abs(DOUT_PRE) * DFLW_IMP**(-7.0 / 3.0))**(-1.0)
                self.D2RIVVEL[ISEQ, 0] = self.D2RIVOUT[ISEQ, 0] * DAREA**(-1.0)
            else:
                self.D2RIVOUT[ISEQ, 0] = 0.0
                self.D2RIVVEL[ISEQ, 0] = 0.0

            if self.LFLDOUT:
                DFLW_F = max(DSFCMAX - self.D2ELEVTN[ISEQ, 0], 0.0)
                DARE_F = self.P2FLDSTO[ISEQ, 0] * self.D2RIVLEN[ISEQ, 0]**(-1.0)
                DARE_F = max(DARE_F - self.D2FLDDPH[ISEQ, 0] * self.D2RIVWTH[ISEQ, 0], 0.0)

                DFLW_PRE_F = DSFCMAX_PRE - self.D2ELEVTN[ISEQ, 0]
                DFLW_IMP_F = max((max(DFLW_F * DFLW_PRE_F, 0.0))**0.5, 1.0e-6)

                DARE_PRE_F = self.D2FLDSTO_PRE[ISEQ, 0] * self.D2RIVLEN[ISEQ, 0]**(-1.0)
                DARE_PRE_F = max(DARE_PRE_F - D2FLDDPH_PRE[ISEQ, 0] * self.D2RIVWTH[ISEQ, 0], 1.0e-6)
                DARE_IMP_F = max((DARE_F * DARE_PRE_F)**0.5, 1.0e-6)

                if DFLW_IMP_F > 1.0e-5 and DARE_IMP_F > 1.0e-5:
                    DOUT_PRE_F = self.D2FLDOUT_PRE[ISEQ, 0]
                    self.D2FLDOUT[ISEQ, 0] = (DOUT_PRE_F + self.PGRV * self.DT * DARE_IMP_F * DSLOPE_F) * (1.0 + self.PGRV * self.DT * self.PMANFLD**2.0 * abs(DOUT_PRE_F) * DFLW_IMP_F**(-4.0 / 3.0) * DARE_IMP_F**(-1.0))**(-1.0)
                else:
                    self.D2FLDOUT[ISEQ, 0] = 0.0

                if self.D2FLDOUT[ISEQ, 0] * self.D2RIVOUT[ISEQ, 0] < 0.0:
                    self.D2FLDOUT[ISEQ, 0] = 0.0

        for ISEQ in range(self.NSEQRIV, self.NSEQALL):
            if self.LSLOPEMOUTH:
                DSLOPE = self.D2ELEVSLOPE[ISEQ, 0]
            else:
                DSLOPE = (self.D2SFCELV[ISEQ, 0] - self.D2DWNELV[ISEQ, 0]) * self.PDSTMTH**(-1.0)
            DSLOPE_F = max(-0.005, min(0.005, DSLOPE))

            DFLW = self.D2RIVDPH[ISEQ, 0]
            DAREA = self.D2RIVWTH[ISEQ, 0] * DFLW

            DFLW_PRE = self.D2RIVDPH_PRE[ISEQ, 0]
            DFLW_IMP = max((DFLW * DFLW_PRE)**0.5, 1.0e-6)

            if DFLW_IMP > 1.0e-5 and DAREA > 1.0e-5:
                DOUT_PRE = self.D2RIVOUT_PRE[ISEQ, 0] * self.D2RIVWTH[ISEQ, 0]**(-1.0)
                self.D2RIVOUT[ISEQ, 0] = self.D2RIVWTH[ISEQ, 0] * (DOUT_PRE + self.PGRV * self.DT * DFLW_IMP * DSLOPE) * (1.0 + self.PGRV * self.DT * self.D2RIVMAN[ISEQ, 0]**2.0 * abs(DOUT_PRE) * DFLW_IMP**(-7.0 / 3.0))**(-1.0)
                self.D2RIVVEL[ISEQ, 0] = self.D2RIVOUT[ISEQ, 0] * DAREA**(-1.0)
            else:
                self.D2RIVOUT[ISEQ, 0] = 0.0
                self.D2RIVVEL[ISEQ, 0] = 0.0

            if self.LFLDOUT:
                DFLW_F = self.D2SFCELV[ISEQ, 0] - self.D2ELEVTN[ISEQ, 0]

                DARE_F = self.P2FLDSTO[ISEQ, 0] * self.D2RIVLEN[ISEQ, 0]**(-1.0)
                DARE_F = max(DARE_F - self.D2FLDDPH[ISEQ, 0] * self.D2RIVWTH[ISEQ, 0], 0.0)

                DFLW_PRE_F = D2SFCELV_PRE[ISEQ, 0] - self.D2ELEVTN[ISEQ, 0]
                DFLW_IMP_F = max((max(DFLW_F * DFLW_PRE_F, 0.0))**0.5, 1.0e-6)

                DARE_PRE_F = self.D2FLDSTO_PRE[ISEQ, 0] * self.D2RIVLEN[ISEQ, 0]**(-1.0)
                DARE_PRE_F = max(DARE_PRE_F - D2FLDDPH_PRE[ISEQ, 0] * self.D2RIVWTH[ISEQ, 0], 1.0e-6)
                DARE_IMP_F = max((DARE_F * DARE_PRE_F)**0.5, 1.0e-6)

                if DFLW_IMP_F > 1.0e-5 and DARE_IMP_F > 1.0e-5:
                    DOUT_PRE_F = self.D2FLDOUT_PRE[ISEQ, 0]
                    self.D2FLDOUT[ISEQ, 0] = (DOUT_PRE_F + self.PGRV * self.DT * DARE_IMP_F * DSLOPE_F) * (1.0 + self.PGRV * self.DT * self.PMANFLD**2.0 * abs(DOUT_PRE_F) * DFLW_IMP_F**(-4.0 / 3.0) * DARE_IMP_F**(-1.0))**(-1.0)
                else:
                    self.D2FLDOUT[ISEQ, 0] = 0.0

                if self.D2FLDOUT[ISEQ, 0] * self.D2RIVOUT[ISEQ, 0] < 0.0:
                    self.D2FLDOUT[ISEQ, 0] = 0.0

    def CMF_CALC_INFLOW(self, NPTHOUT, NPTHLEV, I2MASK, PTH_UPST, PTH_DOWN, D1PTHFLW, D2PTHOUT, D1PTHFLWSUM):
        P2STOOUT = np.zeros((self.NSEQMAX, 1))
        P2RIVINF = np.zeros((self.NSEQMAX, 1))
        P2FLDINF = np.zeros((self.NSEQMAX, 1))
        P2PTHOUT = np.zeros((self.NSEQMAX, 1))
        D2RATE = np.ones((self.NSEQMAX, 1))

        for ISEQ in range(self.NSEQALL):
            P2RIVINF[ISEQ, 0] = 0.0
            P2FLDINF[ISEQ, 0] = 0.0
            P2PTHOUT[ISEQ, 0] = 0.0
            P2STOOUT[ISEQ, 0] = 0.0
            D2RATE[ISEQ, 0] = 1.0

        for ISEQ in range(self.NSEQRIV):
            JSEQ = self.I1NEXT[ISEQ]
            OUT_R1 = max(self.D2RIVOUT[ISEQ, 0], 0.0)
            OUT_R2 = max(-self.D2RIVOUT[ISEQ, 0], 0.0)
            OUT_F1 = max(self.D2FLDOUT[ISEQ, 0], 0.0)
            OUT_F2 = max(-self.D2FLDOUT[ISEQ, 0], 0.0)
            DIUP = (OUT_R1 + OUT_F1) * self.DT
            DIDW = (OUT_R2 + OUT_F2) * self.DT
            P2STOOUT[ISEQ, 0] += DIUP
            P2STOOUT[JSEQ, 0] += DIDW

        for ISEQ in range(self.NSEQRIV, self.NSEQALL):
            OUT_R1 = max(self.D2RIVOUT[ISEQ, 0], 0.0)
            OUT_F1 = max(self.D2FLDOUT[ISEQ, 0], 0.0)
            P2STOOUT[ISEQ, 0] += OUT_R1 * self.DT + OUT_F1 * self.DT

        if self.LPTHOUT:
            for IPTH in range(NPTHOUT):
                ISEQP = PTH_UPST[IPTH]
                JSEQP = PTH_DOWN[IPTH]
                if ISEQP <= 0 or JSEQP <= 0:
                    continue
                if I2MASK[ISEQP, 0] > 0 or I2MASK[JSEQP, 0] > 0:
                    continue
                OUT_R1 = max(D1PTHFLWSUM[IPTH], 0.0)
                OUT_R2 = max(-D1PTHFLWSUM[IPTH], 0.0)
                DIUP = OUT_R1 * self.DT
                DIDW = OUT_R2 * self.DT
                P2STOOUT[ISEQP, 0] += DIUP
                P2STOOUT[JSEQP, 0] += DIDW

        for ISEQ in range(self.NSEQALL):
            if P2STOOUT[ISEQ, 0] > 1.0e-8:
                D2RATE[ISEQ, 0] = min((self.P2RIVSTO[ISEQ, 0] + self.P2FLDSTO[ISEQ, 0]) * P2STOOUT[ISEQ, 0]**(-1.0), 1.0)

        for ISEQ in range(self.NSEQRIV):
            JSEQ = self.I1NEXT[ISEQ]
            if self.D2RIVOUT[ISEQ, 0] >= 0.0:
                self.D2RIVOUT[ISEQ, 0] *= D2RATE[ISEQ, 0]
                self.D2FLDOUT[ISEQ, 0] *= D2RATE[ISEQ, 0]
            else:
                self.D2RIVOUT[ISEQ, 0] *= D2RATE[JSEQ, 0]
                self.D2FLDOUT[ISEQ, 0] *= D2RATE[JSEQ, 0]
            P2RIVINF[JSEQ, 0] += self.D2RIVOUT[ISEQ, 0]
            P2FLDINF[JSEQ, 0] += self.D2FLDOUT[ISEQ, 0]

        for ISEQ in range(self.NSEQRIV, self.NSEQALL):
            self.D2RIVOUT[ISEQ, 0] *= D2RATE[ISEQ, 0]
            self.D2FLDOUT[ISEQ, 0] *= D2RATE[ISEQ, 0]

        if self.LPTHOUT:
            for IPTH in range(NPTHOUT):
                ISEQP = PTH_UPST[IPTH]
                JSEQP = PTH_DOWN[IPTH]
                if ISEQP <= 0 or JSEQP <= 0:
                    continue
                if I2MASK[ISEQP, 0] > 0 or I2MASK[JSEQP, 0] > 0:
                    continue
                for ILEV in range(NPTHLEV):
                    if D1PTHFLW[IPTH, ILEV] >= 0.0:
                        D1PTHFLW[IPTH, ILEV] *= D2RATE[ISEQP, 0]
                    else:
                        D1PTHFLW[IPTH, ILEV] *= D2RATE[JSEQP, 0]
                if D1PTHFLWSUM[IPTH] >= 0.0:
                    D1PTHFLWSUM[IPTH] *= D2RATE[ISEQP, 0]
                else:
                    D1PTHFLWSUM[IPTH] *= D2RATE[JSEQP, 0]
                P2PTHOUT[ISEQP, 0] += D1PTHFLWSUM[IPTH]
                P2PTHOUT[JSEQP, 0] -= D1PTHFLWSUM[IPTH]

        self.D2RIVINF[:, :] = P2RIVINF[:, :]
        self.D2FLDINF[:, :] = P2FLDINF[:, :]
        D2PTHOUT[:, :] = P2PTHOUT[:, :]
