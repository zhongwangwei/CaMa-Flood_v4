import numpy as np

class CMF_CTRL_RESTART_MOD:
    def __init__(self):
        self.CRESTSTO = "restart"
        self.CRESTDIR = "./"
        self.CVNREST = "restart"
        self.LRESTCDF = False
        self.LRESTDBL = True
        self.IFRQ_RST = 0

    def CMF_RESTART_NMLIST(self):
        print("CMF::RESTART_NMLIST: namelist OPEN in unit: ", self.CSETFILE, self.NSETFILE)
        self.CRESTSTO = "restart"
        self.CRESTDIR = "./"
        self.CVNREST = "restart"
        self.LRESTCDF = False
        self.LRESTDBL = True
        self.IFRQ_RST = 0
        print("=== NAMELIST, NRESTART ===")
        print("CRESTSTO:  ", self.CRESTSTO)
        print("CRESTDIR:  ", self.CRESTDIR)
        print("CVNREST:   ", self.CVNREST)
        print("LRESTCDF:  ", self.LRESTCDF)
        print("LRESTDBL:  ", self.LRESTDBL)
        print("IFRQ_RST:  ", self.IFRQ_RST)

    def CMF_RESTART_INIT(self):
        self.P2RIVSTO.fill(0)
        self.P2FLDSTO.fill(0)
        self.D2RIVOUT.fill(0)
        self.D2FLDOUT.fill(0)
        self.D2RIVOUT_PRE.fill(0)
        self.D2FLDOUT_PRE.fill(0)
        self.D2RIVDPH_PRE.fill(0)
        self.D2FLDSTO_PRE.fill(0)
        if self.LPTHOUT:
            self.D1PTHFLW.fill(0)
            self.D1PTHFLW_PRE.fill(0)
        if self.LDAMOUT:
            self.P2DAMSTO.fill(0)
        if self.LLEVEE:
            self.P2LEVSTO.fill(0)
        if self.LGDWDLY:
            self.P2GDWSTO.fill(0)
        if self.LRESTCDF:
            self.READ_REST_CDF()
        else:
            self.READ_REST_BIN()
        if self.LSTOONLY:
            self.D2FLDSTO_PRE = self.P2FLDSTO

    def READ_REST_BIN(self):
        print('READ_REST: read restart binary: ', self.CRESTSTO)
        if self.LRESTDBL:
            with open(self.CRESTSTO, 'rb') as f:
                self.P2RIVSTO = np.fromfile(f, dtype=np.float64).reshape(self.NSEQMAX, 1)
                self.P2FLDSTO = np.fromfile(f, dtype=np.float64).reshape(self.NSEQMAX, 1)
                if not self.LSTOONLY:
                    self.D2RIVOUT_PRE = np.fromfile(f, dtype=np.float64).reshape(self.NSEQMAX, 1)
                    self.D2FLDOUT_PRE = np.fromfile(f, dtype=np.float64).reshape(self.NSEQMAX, 1)
                    self.D2RIVDPH_PRE = np.fromfile(f, dtype=np.float64).reshape(self.NSEQMAX, 1)
                    self.D2FLDSTO_PRE = np.fromfile(f, dtype=np.float64).reshape(self.NSEQMAX, 1)
                if self.LGDWDLY:
                    self.P2GDWSTO = np.fromfile(f, dtype=np.float64).reshape(self.NSEQMAX, 1)
                if self.LDAMOUT:
                    self.P2DAMSTO = np.fromfile(f, dtype=np.float64).reshape(self.NSEQMAX, 1)
                if self.LLEVEE:
                    self.P2LEVSTO = np.fromfile(f, dtype=np.float64).reshape(self.NSEQMAX, 1)
        else:
            with open(self.CRESTSTO, 'rb') as f:
                self.P2RIVSTO = np.fromfile(f, dtype=np.float32).reshape(self.NSEQMAX, 1)
                self.P2FLDSTO = np.fromfile(f, dtype=np.float32).reshape(self.NSEQMAX, 1)
                if not self.LSTOONLY:
                    self.D2RIVOUT_PRE = np.fromfile(f, dtype=np.float32).reshape(self.NSEQMAX, 1)
                    self.D2FLDOUT_PRE = np.fromfile(f, dtype=np.float32).reshape(self.NSEQMAX, 1)
                    self.D2RIVDPH_PRE = np.fromfile(f, dtype=np.float32).reshape(self.NSEQMAX, 1)
                    self.D2FLDSTO_PRE = np.fromfile(f, dtype=np.float32).reshape(self.NSEQMAX, 1)
                if self.LGDWDLY:
                    self.P2GDWSTO = np.fromfile(f, dtype=np.float32).reshape(self.NSEQMAX, 1)
                if self.LDAMOUT:
                    self.P2DAMSTO = np.fromfile(f, dtype=np.float32).reshape(self.NSEQMAX, 1)
                if self.LLEVEE:
                    self.P2LEVSTO = np.fromfile(f, dtype=np.float32).reshape(self.NSEQMAX, 1)
        if self.LPTHOUT and not self.LSTOONLY:
            with open(self.CRESTSTO + '.pth', 'rb') as f:
                if self.LRESTDBL:
                    self.D1PTHFLW_PRE = np.fromfile(f, dtype=np.float64).reshape(self.NPTHOUT, self.NPTHLEV)
                else:
                    self.D1PTHFLW_PRE = np.fromfile(f, dtype=np.float32).reshape(self.NPTHOUT, self.NPTHLEV)

    def READ_REST_CDF(self):
        print('READ_REST: read restart netcdf: ', self.CRESTSTO)
        import netCDF4 as nc
        with nc.Dataset(self.CRESTSTO, 'r') as ds:
            self.P2RIVSTO = ds.variables['rivsto'][:, :, 0].reshape(self.NSEQMAX, 1)
            self.P2FLDSTO = ds.variables['fldsto'][:, :, 0].reshape(self.NSEQMAX, 1)
            if not self.LSTOONLY:
                self.D2RIVOUT_PRE = ds.variables['rivout_pre'][:, :, 0].reshape(self.NSEQMAX, 1)
                self.D2FLDOUT_PRE = ds.variables['fldout_pre'][:, :, 0].reshape(self.NSEQMAX, 1)
                self.D2RIVDPH_PRE = ds.variables['rivdph_pre'][:, :, 0].reshape(self.NSEQMAX, 1)
                self.D2FLDSTO_PRE = ds.variables['fldsto_pre'][:, :, 0].reshape(self.NSEQMAX, 1)
            if self.LGDWDLY:
                self.P2GDWSTO = ds.variables['gdwsto'][:, :, 0].reshape(self.NSEQMAX, 1)
            if self.LDAMOUT:
                self.P2DAMSTO = ds.variables['damsto'][:, :, 0].reshape(self.NSEQMAX, 1)
            if self.LLEVEE:
                self.P2LEVSTO = ds.variables['levsto'][:, :, 0].reshape(self.NSEQMAX, 1)
            if self.LPTHOUT and not self.LSTOONLY:
                self.D1PTHFLW_PRE = ds.variables['pthflw_pre'][:, :, 0].reshape(self.NPTHOUT, self.NPTHLEV)

    def CMF_RESTART_WRITE(self):
        IREST = 0
        if self.IFRQ_RST >= 0 and self.KSTEP == self.NSTEPS:
            IREST = 1
        if self.IFRQ_RST >= 1 and self.IFRQ_RST <= 24:
            if self.JHOUR % self.IFRQ_RST == 0 and self.JMIN == 0:
                IREST = 1
        if self.IFRQ_RST == 30:
            if self.JDD == 1 and self.JHOUR == 0 and self.JMIN == 0:
                IREST = 1
        if IREST == 1:
            print('CMF::RESTART_WRITE: write time: ', self.JYYYYMMDD, self.JHHMM)
            if self.LRESTCDF:
                self.WRTE_REST_CDF()
            else:
                self.WRTE_REST_BIN()

    def WRTE_REST_BIN(self):
        print('WRTE_REST_BIN: restart file:', self.CRESTDIR + self.CVNREST + self.CDATE + self.CSUFBIN)
        if self.LRESTDBL:
            with open(self.CRESTDIR + self.CVNREST + self.CDATE + self.CSUFBIN, 'wb') as f:
                self.P2RIVSTO.tofile(f)
                self.P2FLDSTO.tofile(f)
                if not self.LSTOONLY:
                    self.D2RIVOUT_PRE.tofile(f)
                    self.D2FLDOUT_PRE.tofile(f)
                    self.D2RIVDPH_PRE.tofile(f)
                    self.D2FLDSTO_PRE.tofile(f)
                if self.LGDWDLY:
                    self.P2GDWSTO.tofile(f)
                if self.LDAMOUT:
                    self.P2DAMSTO.tofile(f)
                if self.LLEVEE:
                    self.P2LEVSTO.tofile(f)
        else:
            with open(self.CRESTDIR + self.CVNREST + self.CDATE + self.CSUFBIN, 'wb') as f:
                self.P2RIVSTO.astype(np.float32).tofile(f)
                self.P2FLDSTO.astype(np.float32).tofile(f)
                if not self.LSTOONLY:
                    self.D2RIVOUT_PRE.astype(np.float32).tofile(f)
                    self.D2FLDOUT_PRE.astype(np.float32).tofile(f)
                    self.D2RIVDPH_PRE.astype(np.float32).tofile(f)
                    self.D2FLDSTO_PRE.astype(np.float32).tofile(f)
                if self.LGDWDLY:
                    self.P2GDWSTO.astype(np.float32).tofile(f)
                if self.LDAMOUT:
                    self.P2DAMSTO.astype(np.float32).tofile(f)
                if self.LLEVEE:
                    self.P2LEVSTO.astype(np.float32).tofile(f)
        if self.LPTHOUT:
            with open(self.CRESTDIR + self.CVNREST + self.CDATE + self.CSUFBIN + '.pth', 'wb') as f:
                if self.LRESTDBL:
                    self.D1PTHFLW_PRE.tofile(f)
                else:
                    self.D1PTHFLW_PRE.astype(np.float32).tofile(f)

    def WRTE_REST_CDF(self):
        print('WRTE_REST:create RESTART NETCDF:', self.CRESTDIR + self.CVNREST + self.CDATE + self.CSUFCDF)
        import netCDF4 as nc
        with nc.Dataset(self.CRESTDIR + self.CVNREST + self.CDATE + self.CSUFCDF, 'w', format='NETCDF4') as ds:
            ds.createDimension('time', None)
            ds.createDimension('lat', self.NY)
            ds.createDimension('lon', self.NX)
            if self.LPTHOUT:
                ds.createDimension('NPTHOUT', self.NPTHOUT)
                ds.createDimension('NPTHLEV', self.NPTHLEV)
            lat = ds.createVariable('lat', 'f4', ('lat',))
            lon = ds.createVariable('lon', 'f4', ('lon',))
            time = ds.createVariable('time', 'f8', ('time',))
            lat.long_name = 'latitude'
            lat.units = 'degrees_north'
            lon.long_name = 'longitude'
            lon.units = 'degrees_east'
            time.long_name = 'time'
            time.units = f'seconds since {self.ISYYYY}-{self.ISMM}-{self.ISDD} {self.ISHOUR}:{self.ISMIN}'
            rivsto = ds.createVariable('rivsto', 'f8', ('lon', 'lat', 'time'), zlib=True, fill_value=self.DMIS)
            fldsto = ds.createVariable('fldsto', 'f8', ('lon', 'lat', 'time'), zlib=True, fill_value=self.DMIS)
            if not self.LSTOONLY:
                rivout_pre = ds.createVariable('rivout_pre', 'f8', ('lon', 'lat', 'time'), zlib=True, fill_value=self.DMIS)
                fldout_pre = ds.createVariable('fldout_pre', 'f8', ('lon', 'lat', 'time'), zlib=True, fill_value=self.DMIS)
                rivdph_pre = ds.createVariable('rivdph_pre', 'f8', ('lon', 'lat', 'time'), zlib=True, fill_value=self.DMIS)
                fldsto_pre = ds.createVariable('fldsto_pre', 'f8', ('lon', 'lat', 'time'), zlib=True, fill_value=self.DMIS)
            if self.LGDWDLY:
                gdwsto = ds.createVariable('gdwsto', 'f8', ('lon', 'lat', 'time'), zlib=True, fill_value=self.DMIS)
            if self.LDAMOUT:
                damsto = ds.createVariable('damsto', 'f8', ('lon', 'lat', 'time'), zlib=True, fill_value=self.DMIS)
            if self.LLEVEE:
                levsto = ds.createVariable('levsto', 'f8', ('lon', 'lat', 'time'), zlib=True, fill_value=self.DMIS)
            if self.LPTHOUT and not self.LSTOONLY:
                pthflw_pre = ds.createVariable('pthflw_pre', 'f8', ('NPTHOUT', 'NPTHLEV', 'time'), zlib=True)
            lat[:] = self.D1LAT
            lon[:] = self.D1LON
            time[0] = (self.KMINNEXT - self.KMINSTART) * 60
            rivsto[:, :, 0] = self.P2RIVSTO.reshape(self.NX, self.NY)
            fldsto[:, :, 0] = self.P2FLDSTO.reshape(self.NX, self.NY)
            if not self.LSTOONLY:
                rivout_pre[:, :, 0] = self.D2RIVOUT_PRE.reshape(self.NX, self.NY)
                fldout_pre[:, :, 0] = self.D2FLDOUT_PRE.reshape(self.NX, self.NY)
                rivdph_pre[:, :, 0] = self.D2RIVDPH_PRE.reshape(self.NX, self.NY)
                fldsto_pre[:, :, 0] = self.D2FLDSTO_PRE.reshape(self.NX, self.NY)
            if self.LGDWDLY:
                gdwsto[:, :, 0] = self.P2GDWSTO.reshape(self.NX, self.NY)
            if self.LDAMOUT:
                damsto[:, :, 0] = self.P2DAMSTO.reshape(self.NX, self.NY)
            if self.LLEVEE:
                levsto[:, :, 0] = self.P2LEVSTO.reshape(self.NX, self.NY)
            if self.LPTHOUT and not self.LSTOONLY:
                pthflw_pre[:, :, 0] = self.D1PTHFLW_PRE

