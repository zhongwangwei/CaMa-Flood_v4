import numpy as np
from cmf_drv_control_mod import cmf_drv_input, cmf_drv_init, cmf_drv_end
from cmf_drv_advance_mod import cmf_drv_advance
from cmf_ctrl_forcing_mod import cmf_forcing_get, cmf_forcing_put
from cmf_ctrl_tracer_mod import cmf_tracer_forc_get, cmf_tracer_forc_interp
from cmf_ctrl_mpi_mod import cmf_mpi_init, cmf_mpi_end
from cmf_ctrl_sedinp_mod import cmf_sed_forcing
from yos_cmf_input import nxin, nyin, dt, dtin, ltrace, lsedout
from yos_cmf_time import nsteps

def main():
    # Local variables
    istep = 0
    istepadv = 0
    zbuff = None

    # MPI Initialization
    cmf_mpi_init()

    # Namelist handling
    cmf_drv_input()

    # Initialization
    cmf_drv_init()

    # Allocate data buffer for input forcing
    zbuff = np.zeros((nxin, nyin, 2))

    # Main temporal loop / time-step (nsteps calculated by drv_init)
    istepadv = int(dtin / dt)
    for istep in range(1, nsteps + 1, istepadv):
        # Read forcing from file, This is only relevant in Stand-alone mode
        cmf_forcing_get(zbuff)

        # Interpolate runoff & send to CaMa-Flood
        cmf_forcing_put(zbuff)

        if ltrace:
            cmf_tracer_forc_get()
            cmf_tracer_forc_interp()

        # Advance CaMa-Flood model for istepadv
        cmf_drv_advance(istepadv)

        if lsedout:
            cmf_sed_forcing()

    # Finalize CaMa-Flood
    del zbuff
    cmf_drv_end()

    # MPI specific finalization
    cmf_mpi_end()

if __name__ == "__main__":
    main()
