import numpy as np

class CMF_CTRL_NMLIST_MOD:
    def __init__(self):
        self.LADPSTP = True
        self.LFPLAIN = True
        self.LKINE = False
        self.LFLDOUT = True
        self.LPTHOUT = False
        self.LDAMOUT = False
        self.LLEVEE = False
        self.LSEDOUT = False
        self.LTRACE = False
        self.LOUTINS = False
        self.LROSPLIT = False
        self.LWEVAP = False
        self.LWEVAPFIX = False
        self.LWEXTRACTRIV = False
        self.LGDWDLY = False
        self.LSLPMIX = False
        self.LSLOPEMOUTH = False
        self.LMEANSL = False
        self.LSEALEV = False
        self.LRESTART = False
        self.LSTOONLY = False
        self.LOUTPUT = True
        self.LOUTINI = False
        self.LGRIDMAP = True
        self.LLEAPYR = True
        self.LMAPEND = False
        self.LBITSAFE = False
        self.LSTG_ES = False
        self.CDIMINFO = "NONE"
        self.DT = 24 * 60 * 60
        self.IFRQ_INP = 24
        self.PMANRIV = 0.03
        self.PMANFLD = 0.10
        self.PGRV = 9.8
        self.PDSTMTH = 10000.0
        self.PCADP = 0.7
        self.PMINSLP = 1.0e-5
        self.IMIS = -9999
        self.RMIS = 1.0e20
        self.DMIS = 1.0e20
        self.CSUFBIN = ".bin"
        self.CSUFVEC = ".vec"
        self.CSUFPTH = ".pth"
        self.CSUFCDF = ".nc"
        self.NX = 1440
        self.NY = 720
        self.NLFP = 10
        self.NXIN = 360
        self.NYIN = 180
        self.INPN = 1
        self.WEST = -180.0
        self.EAST = 180.0
        self.NORTH = 90.0
        self.SOUTH = -90.0

    def CMF_CONFIG_NMLIST(self):
        print("CMF::CONFIG_NMLIST: namelist opened")
        self.LADPSTP = True
        self.LFPLAIN = True
        self.LKINE = False
        self.LFLDOUT = True
        self.LPTHOUT = False
        self.LDAMOUT = False
        self.LLEVEE = False
        self.LSEDOUT = False
        self.LTRACE = False
        self.LOUTINS = False
        self.LROSPLIT = False
        self.LWEVAP = False
        self.LWEVAPFIX = False
        self.LWEXTRACTRIV = False
        self.LGDWDLY = False
        self.LSLPMIX = False
        self.LSLOPEMOUTH = False
        self.LMEANSL = False
        self.LSEALEV = False
        self.LRESTART = False
        self.LSTOONLY = False
        self.LOUTPUT = True
        self.LOUTINI = False
        self.LGRIDMAP = True
        self.LLEAPYR = True
        self.LMAPEND = False
        self.LBITSAFE = False
        self.LSTG_ES = False
        self.CDIMINFO = "NONE"
        self.DT = 24 * 60 * 60
        self.IFRQ_INP = 24
        self.PMANRIV = 0.03
        self.PMANFLD = 0.10
        self.PGRV = 9.8
        self.PDSTMTH = 10000.0
        self.PCADP = 0.7
        self.PMINSLP = 1.0e-5
        self.IMIS = -9999
        self.RMIS = 1.0e20
        self.DMIS = 1.0e20
        self.CSUFBIN = ".bin"
        self.CSUFVEC = ".vec"
        self.CSUFPTH = ".pth"
        self.CSUFCDF = ".nc"
        self.NX = 1440
        self.NY = 720
        self.NLFP = 10
        self.NXIN = 360
        self.NYIN = 180
        self.INPN = 1
        self.WEST = -180.0
        self.EAST = 180.0
        self.NORTH = 90.0
        self.SOUTH = -90.0
        print("=== NAMELIST, NRUNVER ===")
        print("LADPSTP ", self.LADPSTP)
        print("LFPLAIN ", self.LFPLAIN)
        print("LKINE   ", self.LKINE)
        print("LFLDOUT ", self.LFLDOUT)
        print("LPTHOUT ", self.LPTHOUT)
        print("LDAMOUT ", self.LDAMOUT)
        print("LLEVEE  ", self.LLEVEE)
        print("LSEDOUT ", self.LSEDOUT)
        print("LTRACE  ", self.LTRACE)
        print("LOUTINS ", self.LOUTINS)
        print("LROSPLIT ", self.LROSPLIT)
        print("LWEVAP   ", self.LWEVAP)
        print("LWEVAPFIX", self.LWEVAPFIX)
        print("LWEXTRACTRIV", self.LWEXTRACTRIV)
        print("LGDWDLY  ", self.LGDWDLY)
        print("LSLPMIX  ", self.LSLPMIX)
        print("LSLOPEMOUTH ", self.LSLOPEMOUTH)
        print("LMEANSL: ", self.LSEALEV)
        print("LSEALEV: ", self.LSEALEV)
        print("LRESTART ", self.LRESTART)
        print("LSTOONLY ", self.LSTOONLY)
        print("LOUTPUT  ", self.LOUTPUT)
        print("LOUTINI  ", self.LOUTINI)
        print("LGRIDMAP ", self.LGRIDMAP)
        print("LLEAPYR  ", self.LLEAPYR)
        print("LMAPEND  ", self.LMAPEND)
        print("LBITSAFE ", self.LBITSAFE)
        print("LSTG_ES " , self.LSTG_ES)
        print("=== NAMELIST, NCONF ===")
        print("CDIMINFO  ", self.CDIMINFO)
        print("DT        ", self.DT)
        print("DTIN      ", self.IFRQ_INP * 60 * 60)
        print("IFRQ_INP  ", self.IFRQ_INP)
        print("=== DIMINFO ===")
        print("NX,NY,NLFP     ", self.NX, self.NY, self.NLFP)
        print("NXIN,NYIN,INPN ", self.NXIN, self.NYIN, self.INPN)
        print("WEST,EAST,NORTH,SOUTH ", self.WEST, self.EAST, self.NORTH, self.SOUTH)
        print("=== NAMELIST, NPARAM ===")
        print("PMANRIV  ", self.PMANRIV)
        print("PMANRIV  ", self.PMANFLD)
        print("PGRV     ", self.PGRV)
        print("PDSTMTH  ", self.PDSTMTH)
        print("PCADP    ", self.PCADP)
        print("PMINSLP  ", self.PMINSLP)
        print("IMIS     ", self.IMIS)
        print("RMIS     ", self.RMIS)
        print("DMIS     ", self.DMIS)
        print("CSUFBIN  ", self.CSUFBIN)
        print("CSUFVEC  ", self.CSUFVEC)
        print("CSUFPTH  ", self.CSUFPTH)
        print("CSUFCDF  ", self.CSUFCDF)
        print("CMF::CONFIG_NMLIST: end")

    def CMF_CONFIG_CHECK(self):
        print("CMF::CONFIG_CHECK: check setting conflicts")
        if self.DT < 60 or self.DT % 60 != 0:
            print("DT= ", self.DT)
            print("DT should be multiple of 60. CaMa-Flood controls time by MINUTE")
            print("stop")
            raise Exception("DT should be multiple of 60")
        if self.IFRQ_INP * 60 * 60 % self.DT != 0:
            print("DTIN, DT= ", self.IFRQ_INP * 60 * 60, self.DT)
            print("DTIN should be multiple of DT")
            print("stop")
            raise Exception("DTIN should be multiple of DT")
        if self.LSEALEV and self.IFRQ_INP * 60 * 60 % self.DT != 0:
            print("DTSL, DT= ", self.IFRQ_INP * 60 * 60, self.DT)
            print("DTSL should be multiple of DT")
            print("stop")
            raise Exception("DTSL should be multiple of DT")
        if not self.LFPLAIN and not self.LKINE:
            print("LFPLAIN=.false. & LKINE=.false.")
            print("CAUTION: NO FLOODPLAIN OPTION recommended to be used with kinematic wave (LKINE=.true.)")
        if self.LKINE and self.LADPSTP:
            print("LKINE=.true. & LADPSTP=.true.")
            print("adaptive time step recommended only with local inertial equation (LKINE=.false.)")
            print("Set appropriate fixed time step for Kinematic Wave")
        if self.LKINE and self.LPTHOUT:
            print("LKINE=.true. & LPATHOUT=.true.")
            print("bifurcation channel flow only available with local inertial equation (LKINE=.false.)")
            print("STOP")
            raise Exception("bifurcation channel flow only available with local inertial equation")
        if self.LGDWDLY and not self.LROSPLIT:
            print("LGDWDLY=true and LROSPLIT=false")
            print("Ground water reservoir can only be active when runoff splitting is on")
        if self.LWEVAPFIX and not self.LWEVAP:
            print("LWEVAPFIX=true and LWEVAP=false")
            print("LWEVAPFIX can only be active if LWEVAP is active")
        if self.LWEXTRACTRIV and not self.LWEVAP:
            print("LWEXTRACTRIV=true and LWEVAP=false")
            print("LWEXTRACTRIV can only be active if LWEVAP is active")
        print("CMF::CONFIG_CHECK: end")
