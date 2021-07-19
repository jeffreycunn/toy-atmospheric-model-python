# Import User_Settings

import user_settings as us

# Define Constants

import constants as cs

from griddef import griddef
from base import base
from cmm_init import cmm_init
from integ import integ


def CMM_driver():
    # Define Grid
    # test = "running CMM_driver"
    print("Running griddef...")
    zu, zw, xu, xs = griddef(us.nz, us.nx, us.dz, us.dx, cs.zu, cs.zw, cs.xu, cs.xs)
    print("griddef ran")
    # Initialize the Base State
    tb, qb, qbs, tbv, pisfc, pib, rhou, rhow, rhb, ub, um, u = base(us.profile_method, us.nx,
                                                                    us.nz, us.dz, us.psurf, us.qsurf, us.q4km, us.ztr,
                                                                    us.temptr,
                                                                    us.ttr, us.tsurf, cs.p0, cs.cp, cs.g, cs.rd, cs.xk,
                                                                    cs.c_v, cs.zu, cs.rhow, cs.rhou, cs.tb,
                                                                    cs.tbv, cs.qb, cs.qbs, cs.rhb, cs.pib, cs.ub, cs.um,
                                                                    cs.u, cs.up)

    # Set Initial Conditions
    print("running cmm_init")
    th, thm, pim, pic, pprt = cmm_init(cs.xs, cs.g, us.nx, us.nz, cs.zu, us.dx, us.dz, us.zcnt, us.xcnt, us.radz,
                                       us.radx, cs.trigpi, cs.cp, cs.rhou, cs.tb, cs.qb, us.thermamp, cs.th, cs.thm,
                                       cs.pim, cs.pic, cs.pprt, us.bubble_switch)
    print("cmm_init ran...")
    # Integrate the Model

    print("Running integ...")
    integ(cs.tbv, cs.pib, cs.p0, cs.lv, cs.rd, cs.ub, cs.g, cs.cp, us.c_sound, cs.rhow, cs.rhou, cs.tb, cs.zu, cs.zw,
          cs.xu, cs.xs, us.x_dist_in, us.latdmpcoef, us.raydmpcoef, us.raydmpz, cs.trigpi, cs.qb, cs.um, cs.u, cs.up,
          cs.wm,
          cs.w, cs.wp, cs.thm, cs.th, cs.thp, cs.pim, cs.pic, cs.pip, cs.pprt, cs.qvtot, cs.qvm, cs.qv, cs.qvp, cs.qcm,
          cs.qc, cs.qcp,
          cs.qrainm, cs.qrain, cs.qrainp, us.dx, us.dz, us.nt, us.nz, us.nx, us.dt, us.asscoef, us.cmixh, us.cmixv,
          us.qc0, us.k1,
          us.k2, cs.thvm, cs.thv, cs.thvp, 'blank', 25)

    print("integ ran")
    # Save Model Data to Output File

    # save (modeloutputpathandfilename)

    return zu, zw, xu, xs, u, th
