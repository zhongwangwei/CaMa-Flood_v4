import numpy as np

class CMF_CALC_DIAG_MOD:
    def __init__(self):
        self.NADD_adp = 0
        self.D2RIVOUT_aAVG = np.zeros((NSEQALL, 1))
        self.D2FLDOUT_aAVG = np.zeros((NSEQALL, 1))
        self.D2OUTFLW_aAVG = np.zeros((NSEQALL, 1))
        self.D2RIVVEL_aAVG = np.zeros((NSEQALL, 1))
        self.D2PTHOUT_aAVG = np.zeros((NSEQALL, 1))
        self.D2GDWRTN_aAVG = np.zeros((NSEQALL, 1))
        self.D2RUNOFF_aAVG = np.zeros((NSEQALL, 1))
        self.D2ROFSUB_aAVG = np.zeros((NSEQALL, 1))
        self.D2DAMINF_aAVG = np.zeros((NSEQALL, 1))
        self.D2WEVAPEX_aAVG = np.zeros((NSEQALL, 1))
        self.D1PTHFLW_aAVG = np.zeros((NPTHOUT, NPTHLEV))
        self.D1PTHFLWSUM_aAVG = np.zeros(NPTHOUT)
        self.D2STORGE_aMAX = np.zeros((NSEQALL, 1))
        self.D2OUTFLW_aMAX = np.zeros((NSEQALL, 1))
        self.D2RIVDPH_aMAX = np.zeros((NSEQALL, 1))

    def CMF_DIAG_RESET_ADPSTP(self):
        print("CMF::DIAG_AVERAGE: reset", JYYYYMMDD, JHHMM)
        self.NADD_adp = 0
        self.D2RIVOUT_aAVG.fill(0)
        self.D2FLDOUT_aAVG.fill(0)
        self.D2OUTFLW_aAVG.fill(0)
        self.D2RIVVEL_aAVG.fill(0)
        self.D2PTHOUT_aAVG.fill(0)
        self.D2GDWRTN_aAVG.fill(0)
        self.D2RUNOFF_aAVG.fill(0)
        self.D2ROFSUB_aAVG.fill(0)
        if LDAMOUT:
            self.D2DAMINF_aAVG.fill(0)
        if LWEVAP:
            self.D2WEVAPEX_aAVG.fill(0)
        self.D1PTHFLW_aAVG.fill(0)
        self.D1PTHFLWSUM_aAVG.fill(0)
        self.D2STORGE_aMAX.fill(0)
        self.D2OUTFLW_aMAX.fill(0)
        self.D2RIVDPH_aMAX.fill(0)

    def CMF_DIAG_AVEMAX_ADPSTP(self):
        self.NADD_adp += DT
        for ISEQ in range(NSEQALL):
            self.D2RIVOUT_aAVG[ISEQ, 0] += D2RIVOUT[ISEQ, 0] * DT
            self.D2FLDOUT_aAVG[ISEQ, 0] += D2FLDOUT[ISEQ, 0] * DT
            self.D2RIVVEL_aAVG[ISEQ, 0] += D2RIVVEL[ISEQ, 0] * DT
            self.D2OUTFLW_aAVG[ISEQ, 0] += D2OUTFLW[ISEQ, 0] * DT
            self.D2PTHOUT_aAVG[ISEQ, 0] += D2PTHOUT[ISEQ, 0] * DT - D2PTHINF[ISEQ, 0] * DT
            self.D2GDWRTN_aAVG[ISEQ, 0] += D2GDWRTN[ISEQ, 0] * DT
            self.D2RUNOFF_aAVG[ISEQ, 0] += D2RUNOFF[ISEQ, 0] * DT
            self.D2ROFSUB_aAVG[ISEQ, 0] += D2ROFSUB[ISEQ, 0] * DT
            self.D2OUTFLW_aMAX[ISEQ, 0] = max(self.D2OUTFLW_aMAX[ISEQ, 0], abs(D2OUTFLW[ISEQ, 0]))
            self.D2RIVDPH_aMAX[ISEQ, 0] = max(self.D2RIVDPH_aMAX[ISEQ, 0], D2RIVDPH[ISEQ, 0])
            self.D2STORGE_aMAX[ISEQ, 0] = max(self.D2STORGE_aMAX[ISEQ, 0], D2STORGE[ISEQ, 0])
            if LWEVAP:
                self.D2WEVAPEX_aAVG[ISEQ, 0] += D2WEVAPEX[ISEQ, 0] * DT
        if LDAMOUT:
            for ISEQ in range(NSEQALL):
                self.D2DAMINF_aAVG[ISEQ, 0] += P2DAMINF[ISEQ, 0] * DT
        if LPTHOUT:
            for IPTH in range(NPTHOUT):
                self.D1PTHFLW_aAVG[IPTH, :] += D1PTHFLW[IPTH, :] * DT
        if LSEDOUT:
            if LSEDOUT:
                sadd_riv += DT
                for ISEQ in range(NSEQALL):
                    d2rivout_sed[ISEQ] += D2RIVOUT[ISEQ, 0] * DT
                    d2rivvel_sed[ISEQ] += D2RIVVEL[ISEQ, 0] * DT

    def CMF_DIAG_GETAVE_ADPSTP(self):
        print("CMF::DIAG_AVERAGE: time-average", self.NADD_adp, JYYYYMMDD, JHHMM)
        self.D2RIVOUT_aAVG /= self.NADD_adp
        self.D2FLDOUT_aAVG /= self.NADD_adp
        self.D2OUTFLW_aAVG /= self.NADD_adp
        self.D2RIVVEL_aAVG /= self.NADD_adp
        self.D2PTHOUT_aAVG /= self.NADD_adp
        self.D2GDWRTN_aAVG /= self.NADD_adp
        self.D2RUNOFF_aAVG /= self.NADD_adp
        self.D2ROFSUB_aAVG /= self.NADD_adp
        if LDAMOUT:
            self.D2DAMINF_aAVG /= self.NADD_adp
        if LWEVAP:
            self.D2WEVAPEX_aAVG /= self.NADD_adp
        self.D1PTHFLW_aAVG /= self.NADD_adp
        for ILEV in range(NPTHLEV):
            self.D1PTHFLWSUM_aAVG += self.D1PTHFLW_aAVG[:, ILEV]

    def CMF_DIAG_RESET_OUTPUT(self):
        print("CMF::DIAG_AVERAGE: reset", JYYYYMMDD, JHHMM)
        self.NADD_out = 0
        self.D2RIVOUT_oAVG.fill(0)
        self.D2FLDOUT_oAVG.fill(0)
        self.D2OUTFLW_oAVG.fill(0)
        self.D2RIVVEL_oAVG.fill(0)
        self.D2PTHOUT_oAVG.fill(0)
        self.D2GDWRTN_oAVG.fill(0)
        self.D2RUNOFF_oAVG.fill(0)
        self.D2ROFSUB_oAVG.fill(0)
        if LDAMOUT:
            self.D2DAMINF_oAVG.fill(0)
        if LWEVAP:
            self.D2WEVAPEX_oAVG.fill(0)
        self.D1PTHFLW_oAVG.fill(0)
        self.D2STORGE_oMAX.fill(0)
        self.D2OUTFLW_oMAX.fill(0)
        self.D2RIVDPH_oMAX.fill(0)

    def CMF_DIAG_AVEMAX_OUTPUT(self):
        self.NADD_out += DT
        for ISEQ in range(NSEQALL):
            self.D2RIVOUT_oAVG[ISEQ, 0] += self.D2RIVOUT_aAVG[ISEQ, 0] * DT
            self.D2FLDOUT_oAVG[ISEQ, 0] += self.D2FLDOUT_aAVG[ISEQ, 0] * DT
            self.D2RIVVEL_oAVG[ISEQ, 0] += self.D2RIVVEL_aAVG[ISEQ, 0] * DT
            self.D2OUTFLW_oAVG[ISEQ, 0] += self.D2OUTFLW_aAVG[ISEQ, 0] * DT
            self.D2PTHOUT_oAVG[ISEQ, 0] += self.D2PTHOUT_aAVG[ISEQ, 0] * DT
            self.D2GDWRTN_oAVG[ISEQ, 0] += self.D2GDWRTN_aAVG[ISEQ, 0] * DT
            self.D2RUNOFF_oAVG[ISEQ, 0] += self.D2RUNOFF_aAVG[ISEQ, 0] * DT
            self.D2ROFSUB_oAVG[ISEQ, 0] += self.D2ROFSUB_aAVG[ISEQ, 0] * DT
            self.D2OUTFLW_oMAX[ISEQ, 0] = max(self.D2OUTFLW_oMAX[ISEQ, 0], abs(self.D2OUTFLW_aMAX[ISEQ, 0]))
            self.D2RIVDPH_oMAX[ISEQ, 0] = max(self.D2RIVDPH_oMAX[ISEQ, 0], self.D2RIVDPH_aMAX[ISEQ, 0])
            self.D2STORGE_oMAX[ISEQ, 0] = max(self.D2STORGE_oMAX[ISEQ, 0], self.D2STORGE_aMAX[ISEQ, 0])
            if LWEVAP:
                self.D2WEVAPEX_oAVG[ISEQ, 0] += self.D2WEVAPEX_aAVG[ISEQ, 0] * DT
        if LDAMOUT:
            for ISEQ in range(NSEQALL):
                self.D2DAMINF_oAVG[ISEQ, 0] += self.D2DAMINF_aAVG[ISEQ, 0] * DT
        if LPTHOUT:
            for IPTH in range(NPTHOUT):
                self.D1PTHFLW_oAVG[IPTH, :] += self.D1PTHFLW_aAVG[IPTH, :] * DT

    def CMF_DIAG_GETAVE_OUTPUT(self):
        print("CMF::DIAG_AVERAGE: time-average", self.NADD_out, JYYYYMMDD, JHHMM)
        self.D2RIVOUT_oAVG /= self.NADD_out
        self.D2FLDOUT_oAVG /= self.NADD_out
        self.D2OUTFLW_oAVG /= self.NADD_out
        self.D2RIVVEL_oAVG /= self.NADD_out
        self.D2PTHOUT_oAVG /= self.NADD_out
        self.D2GDWRTN_oAVG /= self.NADD_out
        self.D2RUNOFF_oAVG /= self.NADD_out
        self.D2ROFSUB_oAVG /= self.NADD_out
        if LDAMOUT:
            self.D2DAMINF_oAVG /= self.NADD_out
        if LWEVAP:
            self.D2WEVAPEX_oAVG /= self.NADD_out
        self.D1PTHFLW_oAVG /= self.NADD_out
