import numpy as np
from mpi4py import MPI

class CMF_CTRL_MPI_MOD:
    def __init__(self):
        self.ierr = None
        self.Nproc = None
        self.Nid = None
        self.iOMP = None
        self.nOMP = None
        self.REGIONALL = 1
        self.REGIONTHIS = 1
        self.MPI_COMM_CAMA = MPI.COMM_WORLD

    def CMF_MPI_INIT(self, ICOMM_CMF=None):
        if ICOMM_CMF is not None:
            self.MPI_COMM_CAMA = ICOMM_CMF
            self.REGIONALL = MPI.COMM_WORLD.Get_size()
            self.REGIONTHIS = MPI.COMM_WORLD.Get_rank() + 1
        else:
            MPI.Init()
            self.MPI_COMM_CAMA = MPI.COMM_WORLD
            self.REGIONALL = MPI.COMM_WORLD.Get_size()
            self.REGIONTHIS = MPI.COMM_WORLD.Get_rank() + 1

    def CMF_MPI_END(self):
        MPI.Finalize()

    def CMF_MPI_AllReduce_R2MAP(self, R2MAP, RMIS):
        R2TMP = np.full_like(R2MAP, RMIS)
        self.MPI_COMM_CAMA.Allreduce(R2MAP, R2TMP, op=MPI.MIN)
        R2MAP[:] = R2TMP[:]

    def CMF_MPI_AllReduce_R1PTH(self, R1PTH, RMIS):
        R1PTMP = np.full_like(R1PTH, RMIS)
        self.MPI_COMM_CAMA.Allreduce(R1PTH, R1PTMP, op=MPI.MIN)
        R1PTH[:] = R1PTMP[:]

    def CMF_MPI_AllReduce_D2MAP(self, D2MAP, DMIS):
        D2TMP = np.full_like(D2MAP, DMIS)
        self.MPI_COMM_CAMA.Allreduce(D2MAP, D2TMP, op=MPI.MIN)
        D2MAP[:] = D2TMP[:]

    def CMF_MPI_AllReduce_P2MAP(self, P2MAP, DMIS):
        P2TMP = np.full_like(P2MAP, DMIS)
        self.MPI_COMM_CAMA.Allreduce(P2MAP, P2TMP, op=MPI.MIN)
        P2MAP[:] = P2TMP[:]

    def CMF_MPI_AllReduce_D1PTH(self, D1PTH, DMIS, NPTHOUT, NPTHLEV, PTH_UPST, PTH_DOWN):
        for IPTH in range(NPTHOUT):
            if PTH_UPST[IPTH] <= 0 or PTH_DOWN[IPTH] <= 0:
                D1PTH[IPTH, :] = DMIS
        D1PTMP = np.full_like(D1PTH, DMIS)
        self.MPI_COMM_CAMA.Allreduce(D1PTH, D1PTMP, op=MPI.MIN)
        D1PTH[:] = D1PTMP[:]

    def CMF_MPI_AllReduce_P1PTH(self, P1PTH, DMIS, NPTHOUT, NPTHLEV, PTH_UPST, PTH_DOWN):
        for IPTH in range(NPTHOUT):
            if PTH_UPST[IPTH] <= 0 or PTH_DOWN[IPTH] <= 0:
                P1PTH[IPTH, :] = DMIS
        P1PTMP = np.full_like(P1PTH, DMIS)
        self.MPI_COMM_CAMA.Allreduce(P1PTH, P1PTMP, op=MPI.MIN)
        P1PTH[:] = P1PTMP[:]

    def CMF_MPI_ADPSTP(self, DT_MIN):
        DT_LOC = DT_MIN
        DT_OUT = np.array(DT_LOC, dtype=np.float64)
        self.MPI_COMM_CAMA.Allreduce(DT_OUT, DT_OUT, op=MPI.MIN)
        DT_MIN = DT_OUT
        print(f"ADPSTP (MPI_AllReduce): DT_LOC->DTMIN {DT_LOC:.2f} {DT_MIN:.2f}")
