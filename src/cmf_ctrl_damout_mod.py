import numpy as np

class CMF_CTRL_DAMOUT_MOD:
    def __init__(self):
        self.CDAMFILE = "./dam_params.csv"
        self.LDAMTXT = True
        self.LDAMH22 = False
        self.LDAMYBY = False
        self.LiVnorm = False
        self.DamID = None
        self.DamName = None
        self.DamIX = None
        self.DamIY = None
        self.DamLon = None
        self.DamLat = None
        self.upreal = None
        self.R_VolUpa = None
        self.Qf = None
        self.Qn = None
        self.DamYear = None
        self.DamStat = None
        self.EmeVol = None
        self.FldVol = None
        self.ConVol = None
        self.NorVol = None
        self.AdjVol = None
        self.Qa = None
        self.DamSeq = None
        self.I1DAM = None

    def CMF_DAMOUT_NMLIST(self):
        print("CMF::DAMOUT_NMLIST: namelist OPEN in unit: ", self.CSETFILE, self.NSETFILE)
        self.CDAMFILE = "./dam_params.csv"
        self.LDAMTXT = True
        self.LDAMH22 = False
        self.LDAMYBY = False
        self.LiVnorm = False
        print("=== NAMELIST, NDAMOUT ===")
        print("CDAMFILE: ", self.CDAMFILE)
        print("LDAMTXT:  ", self.LDAMTXT)
        print("LDAMH22:  ", self.LDAMH22)
        print("LDAMYBY:  ", self.LDAMYBY)
        print("LiVnorm: ", self.LiVnorm)
        print("CMF::DAMOUT_NMLIST: end")

    def CMF_DAMOUT_INIT(self):
        print("CMF::DAMOUT_INIT: initialize dam", self.CDAMFILE)
        self.NDAM = self.read_ndam(self.CDAMFILE)
        print("CMF::DAMOUT_INIT: number of dams", self.NDAM)
        self.DamID = np.zeros(self.NDAM, dtype=int)
        self.DamName = np.zeros(self.NDAM, dtype='U256')
        self.DamIX = np.zeros(self.NDAM, dtype=int)
        self.DamIY = np.zeros(self.NDAM, dtype=int)
        self.DamLon = np.zeros(self.NDAM)
        self.DamLat = np.zeros(self.NDAM)
        self.upreal = np.zeros(self.NDAM)
        self.Qf = np.zeros(self.NDAM)
        self.Qn = np.zeros(self.NDAM)
        self.DamYear = np.zeros(self.NDAM, dtype=int)
        self.DamStat = np.zeros(self.NDAM, dtype=int)
        self.DamSeq = np.zeros(self.NDAM, dtype=int)
        self.FldVol = np.zeros(self.NDAM)
        self.ConVol = np.zeros(self.NDAM)
        self.EmeVol = np.zeros(self.NDAM)
        self.NorVol = np.zeros(self.NDAM)
        self.AdjVol = np.zeros(self.NDAM)
        self.Qa = np.zeros(self.NDAM)
        self.R_VolUpa = np.zeros(self.NDAM)
        self.I1DAM = np.zeros(self.NSEQMAX, dtype=int)
        self.DamSeq.fill(self.IMIS)
        self.DamStat.fill(self.IMIS)
        self.I1DAM.fill(0)
        self.NDAMX = 0
        for IDAM in range(self.NDAM):
            if self.LDAMYBY:
                self.read_dam_params_yby(IDAM)
            else:
                self.read_dam_params(IDAM)
            self.FldVol[IDAM] = self.FldVol_mcm * 1.0e6
            self.ConVol[IDAM] = self.ConVol_mcm * 1.0e6
            self.EmeVol[IDAM] = self.ConVol[IDAM] + self.FldVol[IDAM] * 0.95
            IX = self.DamIX[IDAM]
            IY = self.DamIY[IDAM]
            if IX <= 0 or IX > self.NX or IY <= 0 or IY > self.NY:
                continue
            ISEQ = self.I2VECTOR(IX, IY)
            if self.I1NEXT[ISEQ] == -9999 or ISEQ <= 0:
                continue
            self.NDAMX += 1
            self.DamSeq[IDAM] = ISEQ
            self.DamStat[IDAM] = 2
            self.I1DAM[ISEQ] = 1
            self.I2MASK[ISEQ, 0] = 2
            if self.LDAMH22:
                self.NorVol[IDAM] = self.ConVol[IDAM] * 0.5
                self.R_VolUpa[IDAM] = self.FldVol[IDAM] * 1.0e-6 / self.upreal[IDAM]
            else:
                Vyr = self.Qn[IDAM] * (365.0 * 24.0 * 60.0 * 60.0)
                Qsto = (self.ConVol[IDAM] * 0.7 + Vyr / 4.0) / (180.0 * 24.0 * 60.0 * 60.0)
                self.Qn[IDAM] = min(self.Qn[IDAM], Qsto) * 1.5
                self.AdjVol[IDAM] = self.ConVol[IDAM] + self.FldVol[IDAM] * 0.1
                self.Qa[IDAM] = (self.Qn[IDAM] + self.Qf[IDAM]) * 0.5
            if self.LDAMYBY:
                if self.ISYYYY == self.DamYear[IDAM]:
                    self.DamStat[IDAM] = 1
                elif self.ISYYYY < self.DamYear[IDAM] and self.DamYear[IDAM] > 0:
                    self.DamStat[IDAM] = -1
                    self.I1DAM[ISEQ] = -1
                    self.FldVol[IDAM] = 0.0
                    self.ConVol[IDAM] = 0.0
        print("CMF::DAMOUT_INIT: allocated dams:", self.NDAMX)
        for ISEQ in range(self.NSEQALL):
            if self.I1DAM[ISEQ] <= 0 and self.I1NEXT[ISEQ] > 0:
                JSEQ = self.I1NEXT[ISEQ]
                if self.I1DAM[JSEQ] == 1 or self.I1DAM[JSEQ] == 11:
                    self.I1DAM[ISEQ] = 10
                    self.I2MASK[ISEQ, 0] = 1
            if self.I1DAM[ISEQ] == 1 and self.I1NEXT[ISEQ] > 0:
                JSEQ = self.I1NEXT[ISEQ]
                if self.I1DAM[JSEQ] == 1 or self.I1DAM[JSEQ] == 11:
                    self.I1DAM[ISEQ] = 11
                    self.I2MASK[ISEQ, 0] = 2
        if not self.LRESTART:
            self.P2DAMSTO[:, 0] = 0.0
            for IDAM in range(self.NDAM):
                if self.DamStat[IDAM] == self.IMIS:
                    continue
                ISEQ = self.DamSeq[IDAM]
                if self.DamStat[IDAM] == -1:
                    self.P2DAMSTO[ISEQ, 0] = self.P2RIVSTO[ISEQ, 0] + self.P2FLDSTO[ISEQ, 0]
                else:
                    self.P2DAMSTO[ISEQ, 0] = self.P2RIVSTO[ISEQ, 0] + self.P2FLDSTO[ISEQ, 0]
                    if self.P2DAMSTO[ISEQ, 0] < self.ConVol[IDAM]:
                        self.P2DAMSTO[ISEQ, 0] = self.ConVol[IDAM]
                        self.P2RIVSTO[ISEQ, 0] = self.ConVol[IDAM]
                        self.P2FLDSTO[ISEQ, 0] = 0.0
        else:
            if self.LDAMYBY:
                for IDAM in range(self.NDAM):
                    if self.DamStat[IDAM] == 1:
                        ISEQ = self.DamSeq[IDAM]
                        self.P2DAMSTO[ISEQ, 0] = self.P2RIVSTO[ISEQ, 0] + self.P2FLDSTO[ISEQ, 0]
                        if self.LiVnorm and self.P2DAMSTO[ISEQ, 0] < self.ConVol[IDAM]:
                            self.P2DAMSTO[ISEQ, 0] = self.ConVol[IDAM]
                            self.P2RIVSTO[ISEQ, 0] = self.ConVol[IDAM]
                            self.P2FLDSTO[ISEQ, 0] = 0.0
        for ISEQ in range(self.NSEQALL):
            self.P2DAMINF[ISEQ, 0] = 0.0
        if self.LPTHOUT:
            for IPTH in range(self.NPTHOUT):
                ISEQP = self.PTH_UPST[IPTH]
                JSEQP = self.PTH_DOWN[IPTH]
                if ISEQP <= 0 or JSEQP <= 0:
                    continue
                if self.I1DAM[ISEQP] > 0 or self.I1DAM[JSEQP] > 0:
                    for ILEV in range(self.NPTHLEV):
                        self.PTH_ELV[IPTH, ILEV] = 1.0e20

    def CMF_DAMOUT_CALC(self):
        self.UPDATE_INFLOW()
        for IDAM in range(self.NDAM):
            if self.DamStat[IDAM] <= 0:
                continue
            ISEQD = self.DamSeq[IDAM]
            DamVol = self.P2DAMSTO[ISEQD, 0]
            DamInflow = self.P2DAMINF[ISEQD, 0]
            if self.LDAMH22:
                if DamVol <= self.NorVol[IDAM]:
                    DamOutflw = self.Qn[IDAM] * (DamVol / self.ConVol[IDAM])
                elif self.NorVol[IDAM] < DamVol <= self.ConVol[IDAM]:
                    if self.Qf[IDAM] <= DamInflow:
                        DamOutflw = self.Qn[IDAM] * 0.5 + (DamVol - self.NorVol[IDAM]) / (self.ConVol[IDAM] - self.NorVol[IDAM]) * (self.Qf[IDAM] - self.Qn[IDAM])
                    else:
                        DamOutflw = self.Qn[IDAM] * 0.5 + ((DamVol - self.NorVol[IDAM]) / (self.EmeVol[IDAM] - self.NorVol[IDAM]))**2 * (self.Qf[IDAM] - self.Qn[IDAM])
                elif self.ConVol[IDAM] < DamVol < self.EmeVol[IDAM]:
                    if self.Qf[IDAM] <= DamInflow:
                        DamOutflw = self.Qf[IDAM] + max((1.0 - self.R_VolUpa[IDAM] / 0.2), 0.0) * (DamVol - self.ConVol[IDAM]) / (self.EmeVol[IDAM] - self.ConVol[IDAM]) * (DamInflow - self.Qf[IDAM])
                    else:
                        DamOutflw = self.Qn[IDAM] * 0.5 + ((DamVol - self.NorVol[IDAM]) / (self.EmeVol[IDAM] - self.NorVol[IDAM]))**2 * (self.Qf[IDAM] - self.Qn[IDAM])
                else:
                    DamOutflw = max(DamInflow, self.Qf[IDAM])
            else:
                if DamVol <= self.ConVol[IDAM]:
                    DamOutflw = self.Qn[IDAM] * (DamVol / self.ConVol[IDAM])**0.5
                elif self.ConVol[IDAM] < DamVol <= self.AdjVol[IDAM]:
                    DamOutflw = self.Qn[IDAM] + ((DamVol - self.ConVol[IDAM]) / (self.AdjVol[IDAM] - self.ConVol[IDAM]))**3.0 * (self.Qa[IDAM] - self.Qn[IDAM])
                elif self.AdjVol[IDAM] < DamVol <= self.EmeVol[IDAM]:
                    if DamInflow >= self.Qf[IDAM]:
                        DamOutflw = self.Qn[IDAM] + ((DamVol - self.ConVol[IDAM]) / (self.EmeVol[IDAM] - self.ConVol[IDAM])) * (DamInflow - self.Qn[IDAM])
                        DamOutTmp = self.Qa[IDAM] + ((DamVol - self.AdjVol[IDAM]) / (self.EmeVol[IDAM] - self.AdjVol[IDAM]))**0.1 * (self.Qf[IDAM] - self.Qa[IDAM])
                        DamOutflw = max(DamOutflw, DamOutTmp)
                    else:
                        DamOutflw = self.Qa[IDAM] + ((DamVol - self.AdjVol[IDAM]) / (self.EmeVol[IDAM] - self.AdjVol[IDAM]))**0.1 * (self.Qf[IDAM] - self.Qa[IDAM])
                else:
                    if DamInflow >= self.Qf[IDAM]:
                        DamOutflw = DamInflow
                    else:
                        DamOutflw = self.Qf[IDAM]
            DamOutflw = min(DamOutflw, DamVol / self.DT, float(self.P2RIVSTO[ISEQD, 0] + self.P2FLDSTO[ISEQD, 0]) / self.DT)
            DamOutflw = max(DamOutflw, 0.0)
            self.D2RIVOUT[ISEQD, 0] = DamOutflw
            self.D2FLDOUT[ISEQD, 0] = 0.0

    def CMF_DAMOUT_WATBAL(self):
        GlbDAMSTO = 0.0
        GlbDAMSTONXT = 0.0
        GlbDAMINF = 0.0
        GlbDAMOUT = 0.0
        for IDAM in range(self.NDAM):
            if self.DamStat[IDAM] == self.IMIS:
                continue
            ISEQD = self.DamSeq[IDAM]
            DamInflow = self.D2RIVINF[ISEQD, 0] + self.D2FLDINF[ISEQD, 0] + self.D2RUNOFF[ISEQD, 0]
            DamOutflw = self.D2RIVOUT[ISEQD, 0] + self.D2FLDOUT[ISEQD, 0]
            GlbDAMSTO += self.P2DAMSTO[ISEQD, 0]
            GlbDAMINF += DamInflow * self.DT
            GlbDAMOUT += DamOutflw * self.DT
            self.P2DAMSTO[ISEQD, 0] += DamInflow * self.DT - DamOutflw * self.DT
            GlbDAMSTONXT += self.P2DAMSTO[ISEQD, 0]
        DamMiss = GlbDAMSTO - GlbDAMSTONXT + GlbDAMINF - GlbDAMOUT
        print("CMF::DAM_CALC: DamMiss at all dams:", DamMiss * 1.0e-9)

    def CMF_DAMOUT_WRTE(self):
        if self.LDAMTXT:
            if not self.IsOpen:
                self.IsOpen = True
                self.CYYYY = f"{self.ISYYYY:04d}"
                self.DAMTXT = f"./damtxt-{self.CYYYY}.txt"
                self.LOGDAM = self.INQUIRE_FID()
                self.open_file(self.LOGDAM, self.DAMTXT)
                self.CLEN = f"{self.NDAMX}"
                self.CFMT = f"(i10,{self.CLEN}(a36))"
                self.JDAM = 0
                for IDAM in range(self.NDAM):
                    if self.DamStat[IDAM] == self.IMIS:
                        continue
                    self.JDAM += 1
                    if self.DamStat[IDAM] == -1:
                        self.WriteTxt[self.JDAM] = f"{self.DamID[IDAM]:12d}{-9.0:12.2f}{-9.0:12.2f}"
                        self.WriteTxt2[self.JDAM] = f"{self.upreal[IDAM]:12.2f}{self.Qf[IDAM]:12.2f}{self.Qn[IDAM]:12.2f}"
                    else:
                        self.WriteTxt[self.JDAM] = f"{self.DamID[IDAM]:12d}{(self.FldVol[IDAM] + self.ConVol[IDAM]) * 1.0e-9:12.2f}{self.ConVol[IDAM] * 1.0e-9:12.2f}"
                        self.WriteTxt2[self.JDAM] = f"{self.upreal[IDAM]:12.2f}{self.Qf[IDAM]:12.2f}{self.Qn[IDAM]:12.2f}"
                self.write_to_file(self.LOGDAM, self.CFMT, self.NDAMX, self.WriteTXT)
                self.CFMT = f"(a10,{self.CLEN}(a36))"
                self.write_to_file(self.LOGDAM, self.CFMT, "Date", self.WriteTXT2)
            self.JDAM = 0
            for IDAM in range(self.NDAM):
                if self.DamStat[IDAM] == self.IMIS:
                    continue
                self.JDAM += 1
                ISEQD = self.DamSeq[IDAM]
                DDamInf = self.P2DAMINF[ISEQD, 0]
                DDamOut = self.D2RIVOUT[ISEQD, 0] + self.D2FLDOUT[ISEQD, 0]
                self.WriteTxt[self.JDAM] = f"{self.P2DAMSTO[ISEQD, 0] * 1.0e-9:12.2f}{DDamInf:12.2f}{DDamOut:12.2f}"
            self.CFMT = f"(i10,{self.CLEN}(a36))"
            self.write_to_file(self.LOGDAM, self.CFMT, self.IYYYYMMDD, self.WriteTXT)

    def UPDATE_INFLOW(self):
        for ISEQ in range(self.NSEQALL):
            if self.I1DAM[ISEQ] > 0:
                self.D2RIVOUT[ISEQ, 0] = 0.0
                self.D2FLDOUT[ISEQ, 0] = 0.0
                self.P2DAMINF[ISEQ, 0] = 0.0
        for ISEQ in range(self.NSEQALL):
            if self.I1DAM[ISEQ] == 10 or self.I1DAM[ISEQ] == 11:
                JSEQ = self.I1NEXT[ISEQ]
                self.P2DAMINF[JSEQ, 0] += self.D2RIVOUT_PRE[ISEQ, 0] + self.D2FLDOUT_PRE[ISEQ, 0]
        for ISEQ in range(self.NSEQRIV):
            if self.I1DAM[ISEQ] == 10:
                JSEQ = self.I1NEXT[ISEQ]
                DSLOPE = (self.D2ELEVTN[ISEQ, 0] - self.D2ELEVTN[JSEQ, 0]) * self.D2NXTDST[ISEQ, 0]**(-1.0)
                DSLOPE = max(DSLOPE, self.PMINSLP)
                DVEL = self.D2RIVMAN[ISEQ, 0]**(-1.0) * DSLOPE**0.5 * self.D2RIVDPH[ISEQ, 0]**(2.0 / 3.0)
                DAREA = self.D2RIVWTH[ISEQ, 0] * self.D2RIVDPH[ISEQ, 0]
                self.D2RIVVEL[ISEQ, 0] = DVEL
                self.D2RIVOUT[ISEQ, 0] = DAREA * DVEL
                self.D2RIVOUT[ISEQ, 0] = min(self.D2RIVOUT[ISEQ, 0], float(self.P2RIVSTO[ISEQ, 0]) / self.DT)
                DSLOPE_F = min(0.005, DSLOPE)
                DVEL_F = self.PMANFLD**(-1.0) * DSLOPE_F**0.5 * self.D2FLDDPH[ISEQ, 0]**(2.0 / 3.0)
                DARE_F = self.P2FLDSTO[ISEQ, 0] * self.D2RIVLEN[ISEQ, 0]**(-1.0)
                DARE_F = max(DARE_F - self.D2FLDDPH[ISEQ, 0] * self.D2RIVWTH[ISEQ, 0], 0.0)
                self.D2FLDOUT[ISEQ, 0] = DARE_F * DVEL_F
                self.D2FLDOUT[ISEQ, 0] = min(self.D2FLDOUT[ISEQ, 0] * 1.0, self.P2FLDSTO[ISEQ, 0] / self.DT)

    def read_ndam(self, cdamfile):
        # Implement the read_ndam function here
        pass

    def read_dam_params_yby(self, idam):
        # Implement the read_dam_params_yby function here
        pass

    def read_dam_params(self, idam):
        # Implement the read_dam_params function here
        pass

    def INQUIRE_FID(self):
        # Implement the INQUIRE_FID function here
        pass

    def open_file(self, logdam, damtxt):
        # Implement the open_file function here
        pass

    def write_to_file(self, logdam, cfmt, *args):
        # Implement the write_to_file function here
        pass
