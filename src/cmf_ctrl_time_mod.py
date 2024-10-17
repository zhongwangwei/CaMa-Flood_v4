import numpy as np

class CMF_CTRL_TIME_MOD:
    def __init__(self):
        self.SYEAR = 2000
        self.SMON = 1
        self.SDAY = 1
        self.SHOUR = 0
        self.EYEAR = 2001
        self.EMON = 1
        self.EDAY = 1
        self.EHOUR = 0

    def CMF_TIME_NMLIST(self):
        print("CMF::TIME_NMLIST: namelist OPEN in unit: ", self.CSETFILE, self.NSETFILE)
        self.SYEAR = 2000
        self.SMON = 1
        self.SDAY = 1
        self.SHOUR = 0
        self.EYEAR = 2001
        self.EMON = 1
        self.EDAY = 1
        self.EHOUR = 0
        print("=== NAMELIST, NSIMTIME ===")
        print("SYEAR,SMON,SDAY,SHOUR:", self.SYEAR, self.SMON, self.SDAY, self.SHOUR)
        print("EYEAR,EMON,EDAY,EHOUR:", self.EYEAR, self.EMON, self.EDAY, self.EHOUR)
        self.YYYY0 = self.SYEAR
        self.MM0 = 1
        self.DD0 = 1
        print("TIME_NMLIST: YYYY0 MM0 DD0 set to : ", self.YYYY0, self.MM0, self.DD0)
        print("CMF::TIME_NMLIST: end")

    def CMF_TIME_INIT(self):
        print("CMF::TIME_INIT: initialize time variables")
        self.ISYYYYMMDD = self.SYEAR * 10000 + self.SMON * 100 + self.SDAY
        self.ISHHMM = self.SHOUR * 100
        self.ISYYYY = self.SYEAR
        self.ISMM = self.SMON
        self.ISDD = self.SDAY
        self.ISHOUR = self.SHOUR
        self.ISMIN = 0
        self.IEYYYYMMDD = self.EYEAR * 10000 + self.EMON * 100 + self.EDAY
        self.IEHHMM = self.EHOUR * 100
        self.IEYYYY = self.EYEAR
        self.IEMM = self.EMON
        self.IEDD = self.EDAY
        self.IEHOUR = self.EHOUR
        self.IEMIN = 0
        print('Start Date:', self.ISYYYYMMDD, self.ISHHMM, self.KMINSTART)
        print('End Date:', self.IEYYYYMMDD, self.IEHHMM, self.KMINEND)
        self.KMINSTART = self.DATE2MIN(self.ISYYYYMMDD, self.ISHHMM)
        self.KMINEND = self.DATE2MIN(self.IEYYYYMMDD, self.IEHHMM)
        self.KMIN = self.KMINSTART
        self.KSTEP = 0
        self.NSTEPS = int(((self.KMINEND - self.KMINSTART) * 60) / self.DT)
        print('NSTEPS:', self.NSTEPS)
        self.IYYYYMMDD = self.ISYYYYMMDD
        self.SPLITDATE(self.IYYYYMMDD, self.IYYYY, self.IMM, self.IDD)
        self.IHHMM = self.ISHHMM
        self.SPLITHOUR(self.IHHMM, self.IHOUR, self.IMIN)
        self.KMINNEXT = self.KMIN
        self.JYYYYMMDD = self.IYYYYMMDD
        self.JHHMM = self.IHHMM
        self.SPLITDATE(self.JYYYYMMDD, self.JYYYY, self.JMM, self.JDD)
        self.SPLITHOUR(self.JHHMM, self.JHOUR, self.JMIN)
        print('Initial Time Step Date:Hour :', self.IYYYYMMDD, '_', self.IHOUR, ':', self.IMIN)
        print("CMF::TIME_INIT: end")

    def CMF_TIME_NEXT(self):
        self.KSTEP += 1
        self.KMINNEXT = self.KMIN + int(self.DT / 60)
        print("CMF::TIME_NEXT: ", self.KSTEP, self.KMIN, self.KMINNEXT, self.DT)
        self.MIN2DATE(self.KMINNEXT, self.JYYYYMMDD, self.JHHMM)
        self.SPLITDATE(self.JYYYYMMDD, self.JYYYY, self.JMM, self.JDD)
        self.SPLITHOUR(self.JHHMM, self.JHOUR, self.JMIN)
        print("Strt of Tstep: KMIN, IYYYYMMDD, IHHMM", self.KMIN, self.IYYYYMMDD, self.IHHMM)
        print("End of Tstep: KMINNEXT, JYYYYMMDD, JHHMM", self.KMINNEXT, self.JYYYYMMDD, self.JHHMM)

    def CMF_TIME_UPDATE(self):
        print("CMF_TIME_UPDATE:")
        self.KMIN = self.KMINNEXT
        self.IYYYYMMDD = self.JYYYYMMDD
        self.IYYYY = self.JYYYY
        self.IMM = self.JMM
        self.IDD = self.JDD
        self.IHHMM = self.JHHMM
        self.IHOUR = self.JHOUR
        self.IMIN = self.JMIN
        print("Current time update: KMIN, IYYYYMMDD, IHHMM", self.KMIN, self.IYYYYMMDD, self.IHHMM)

    def DATE2MIN(self, date, hour):
        # Implement the DATE2MIN function here
        pass

    def MIN2DATE(self, kmin, jyyyymmdd, jhhmm):
        # Implement the MIN2DATE function here
        pass

    def SPLITDATE(self, yyyymmdd, yyyy, mm, dd):
        # Implement the SPLITDATE function here
        pass

    def SPLITHOUR(self, hhmm, hour, minute):
        # Implement the SPLITHOUR function here
        pass
