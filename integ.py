import numpy as np

from diffusion import diffusion
from rayleighTop import rayleighTop
from kessler import kessler
from asselin import asselin
from lateralSponge import lateralSponge
from create_image import create_image

def integ(tbv, pib, p0, lv, rd, ub, g, cp, c_sound, rhow, rhou, tb, zu, zw, xu, xs, x_dist_in,
          latdmpcoef, raydmpcoef, raydmpz, trigpi, qb, um, u, up, wm, w, wp, thm, th, thp, pim, pic, pip, pprt,
          qvtot, qvm, qv, qvp, qcm, qc, qcp, qrainm, qrain, qrainp, dx, dz, nt, nz, nx, dt, asscoef, cmixh, cmixv,
          qc0, k1, k2, thvm, thv, thvp, modeloutputpathandfilename, dtoutput):
    VT = np.zeros((nx, nz))
    ql = np.zeros((nx, nz))
    tv = np.zeros((nx, nz))
    ##
    # Set Variable name for a given timestep
    thisdump = dtoutput
    # Set var_0 equal to var_init
    th_0 = th
    u_0 = u
    w_0 = w
    pi_0 = pic
    pprt_0 = pprt
    qv_0 = qv
    qc_0 = qc
    qrain_0 = qrain

    # th_str = 'th_0'
    # u_str = 'u_0'
    # w_str = 'w_0'
    # pi_str = 'pi_0'
    # pprt_str = 'pprt_0'
    # qv_str = 'qv_0'
    # qc_str = 'qc_0'
    # qrain_str = 'qrain_0'
    # dt_str = 'dt'
    # dx_str = 'dx'
    # dz_str = 'dz'
    # nz_str = 'nz'
    # nx_str = 'nx'
    # nt_str = 'nt'
    # qb_str = 'qb'
    # zu_str = 'zu'
    # zw_str = 'zw'
    # xu_str = 'xu'
    # xs_str = 'xs'
    # tb_str = 'tb'
    #
    #
    # # Append New Timesteps to .mat file.
    # full_str = ['save ' modeloutputpathandfilename ' ' th_str ' ' u_str ' ' w_str ' ' pi_str ' ' pprt_str ' ' qv_str ' ' qc_str ' ' qrain_str ' ' dt_str ' ' dx_str ' ' dz_str ' ' nt_str ' ' nx_str  ' ' nz_str ' ' qb_str ' ' zu_str ' ' zw_str ' ' xu_str ' ' xs_str ' ' tb_str]
    # eval(full_str)

    ##
    # um[:,:] = ub
    # u[:,:] = ub
    wm[:, :] = 0.0
    w[:, :] = 0.0
    # pic[:,:] = 0.0         # Uncomment to make the initial step non-hydrostatic
    # pim = pic
    thm = th
    moist = 2
    # mtncenter = 25000

    u[260:340, 1:3] = 0
    w[260:340, 1:3] = 0
    u[270:330, 4:7] = 0
    w[270:330, 4:7] = 0
    u[280:320, 8:11] = 0
    w[280:320, 8:11] = 0
    u[290:310, 12:14] = 0
    w[290:310, 12:14] = 0

    # u(560:640,1:3) = 0
    # w(560:640,1:3) = 0
    # u(570:630,4:7) = 0
    # w(570:630,4:7) = 0
    # u(580:620,8:11) = 0
    # w(580:620,8:11) = 0
    # u(590:610,12:14) = 0
    # w(590:610,12:14) = 0

    um = u
    wm = w

    # Forecast Loop

    for n in range(0, nt - 1):
        # display(n)

        if n == 1:
            d2t = dt
            dtx = dt / dx
            dtz = dt / dz
        else:
            d2t = 2 * dt
            dtx = 2.0 * dt / dx
            dtz = 2.0 * dt / dz

        ## Leap-Frog Scheme

        #  do loop for u
        for k in range(1, nz - 2):
            for i in range(1, nx - 2):
                up[i, k] = um[i, k] - 0.25 * dtx * ((u[i + 1, k] + u[i, k]) ** 2
                                                    - (u[i - 1, k] + u[i, k]) ** 2)
                - 0.25 * dtz * (rhow[0, k + 1] * (w[i, k + 1] + w[i - 1, k + 1])
                                * (u[i, k + 1] + u[i, k])
                                - rhow[0, k] * (w[i, k] + w[i - 1, k])
                                * (u[i, k - 1] + u[i, k])) / rhou[0, k]
                - dtx * cp * tbv[0, k] * (pic[i, k] - pic[i - 1, k])

        #  do loop for w with buoyancy in Fovell's form in terms of theta' and qv'
        #  first prepare for w loop
        for k in range(0, nz - 1):
            for i in range(0, nx - 1):
                if moist == 2:
                    tv[i, k] = (tb[0, k] + th[i, k]) * (1 + 0.61 * (qb[0, k] + qv[i, k])) - tbv[0, k]
                    ql[i, k] = qc[i, k] + qrain[i, k]
                else:
                    tv[i, k] = th[i, k]
                    ql[i, k] = 0.0

        #  now do w loop
        for k in range(2, nz - 2):
            for i in range(1, nx - 2):
                wp[i, k] = wm[i, k] - 0.25 * dtx * ((u[i + 1, k] + u[i + 1, k - 1])
                                                    * (w[i + 1, k] + w[i, k])
                                                    - (u[i, k] + u[i, k - 1])
                                                    * (w[i, k] + w[i - 1, k]))
                - 0.25 * dtz * (rhou[0, k] * (w[i, k + 1] + w[i, k]) ** 2
                                - rhou[0, k - 1] * (w[i, k] + w[i, k - 1]) ** 2) / rhow[0, k]
                - dtz * cp * 0.5 * (tbv[0, k] + tbv[0, k - 1]) * (pic[i, k] - pic[i, k - 1])
                + d2t * g * (0.5 * (tv[i, k] / tbv[0, k] + tv[i, k - 1] / tbv[0, k - 1]) - 0.5 * (
                            ql[i, k] + ql[i, k - 1]))

        #
        #  do loop for theta using flux form
        #     for k=2:nz-1
        #         for i=2:nx-1
        #             thp[i, k]=thm[i, k] - 0.5*dtx*(u[i+1, k]*(th[i+1, k]+th[i, k])  ...
        #                 -u[i, k]*(th[i, k]+th[i-1, k])) ...
        #                 - 0.5*dtz*( rhow[k + 1]*w[i, k+1]             ...
        #                 *(th[i, k+1]+th[i, k])          ...
        #                 - rhow[k]*w[i, k]             ...
        #                 *(th[i, k]+th[i, k-1]))         ...
        #                 /rhou[k] ...
        #                 - 0.5*dtz*( rhow[k + 1]*w[i, k+1]             ...
        #                 *(tb[k + 1] - tb[k])            ...
        #                 + rhow[k]*w[i, k]             ...
        #                 *(tb[k] - tb[k-1]))           ...
        #                 /rhou[k]
        #
        #         end
        #     end

        #  do loop for theta using flux form except advective form on third term
        for k in range(1, nz - 2):
            for i in range(1, nx - 2):
                thp[i, k] = thm[i, k] - 0.5 * dtx * (u[i + 1, k] * (th[i + 1, k] + th[i, k])
                                                     - u[i, k] * (th[i, k] + th[i - 1, k]))
                - 0.5 * dtz * (rhow[0, k + 1] * w[i, k + 1]
                               * (th[i, k + 1] + th[i, k])
                               - rhow[0, k] * w[i, k]
                               * (th[i, k] + th[i, k - 1])) / rhou[0, k]
                - 0.5 * dtz * (w[i, k + 1] * (tb[0, k + 1] - tb[0, k])
                               + w[i, k] * (tb[0, k] - tb[0, k - 1]))

        #  do loop for pi
        for k in range(1, nz - 2):
            for i in range(1, nx - 2):
                pip[i, k] = pim[i, k] - ((c_sound ** 2) / (rhou[0, k] * cp * (tbv[0, k] ** 2))) * (
                            dtx * rhou[0, k] * tbv[0, k] * (u[i + 1, k] - u[i, k])
                            + dtz * 0.5 *
                            (rhow[0, k + 1] * w[i, k + 1] * (tb[0, k + 1] + tb[0, k])
                             - rhow[0, k] * w[i, k] * (tb[0, k] + tb[0, k - 1])))

        # Dimensional Pressure
        for k in range(1, nz - 2):
            for i in range(1, nx - 2):
                pprt[i, k] = pip[i, k] * cp * rhou[0, k] * tb[0, k] * (1.0 + 0.61 * qb[0, k])

        ## Call Diffusion Scheme

        up, um = diffusion(ub, up, um, cmixh, cmixv, dx, dz, dt, nx, nz)
        wp, wm = diffusion(np.zeros((1, nz)), wp, wm, cmixh, cmixv, dx, dz, dt, nx, nz)
        thp, thm = diffusion(np.zeros((1, nz)), thp, thm, cmixh, cmixv, dx, dz, dt, nx, nz)
        pip, pim = diffusion(np.zeros((1, nz)), pip, pim, cmixh, cmixv, dx, dz, dt, nx, nz)

        ## Apply Rayleigh Damping Top to Dry Variables

        for k in range(1, nz - 2):
            for i in range(1, nx - 2):
                if zu[0, k] > raydmpz:
                    up[i, k] = rayleighTop(zu[0, k], zu[0, nz - 1], up[i, k], ub[0, k], raydmpcoef, raydmpz, trigpi, i,
                                           k)
                    wp[i, k] = rayleighTop(zw[0, k], zw[0, nz - 1], wp[i, k], 0, raydmpcoef, raydmpz, trigpi, i, k)
                    thp[i, k] = rayleighTop(zu[0, k], zu[0, nz - 1], thp[i, k], 0, raydmpcoef, raydmpz, trigpi, i, k)
                    pip[i, k] = rayleighTop(zu[0, k], zu[0, nz - 1], pip[i, k], 0, raydmpcoef, raydmpz, trigpi, i, k)

        ## Leapfrog Advection of Moist Variables

        #  do loop for qv using flux form except advective form on third term
        for k in range(1, nz - 3):
            for i in range(1, nx - 2):
                qvp[i, k] = qvm[i, k] - 0.5 * dtx * (u[i + 1, k] * (qv[i + 1, k] + qv[i, k])
                                                     - u[i, k] * (qv[i, k] + qv[i - 1, k]))
                - 0.5 * dtz * (rhow[0, k + 1] * w[i, k + 1]
                               * (qv[i, k + 1] + qv[i, k])
                               - rhow[0, k] * w[i, k]
                               * (qv[i, k] + qv[i, k - 1])) / rhou[0, k]
                - 0.5 * dtz * (w[i, k + 1] * (qb[0, k + 1] - qb[0, k])
                               + w[i, k] * (qb[0, k] - qb[0, k - 1]))

        #  do loop for qc using flux form except advective form on third term
        for k in range(1, nz - 2):
            for i in range(1, nx - 2):
                qcp[i, k] = qcm[i, k] - 0.5 * dtx * (
                            u[i + 1, k] * (qc[i + 1, k] + qc[i, k]) - u[i, k] * (qc[i, k] + qc[i - 1, k]))
                - 0.5 * dtz * (rhow[0, k + 1] * w[i, k + 1]
                               * (qc[i, k + 1] + qc[i, k])
                               - rhow[0, k] * w[i, k]
                               * (qc[i, k] + qc[i, k - 1])) / rhou[0, k]

        #  do loop for qrain using flux form except advective form on third term
        for k in range(1, nz - 2):
            for i in range(1, nx - 2):
                qrainp[i, k] = qrainm[i, k] - 0.5 * dtx * (u[i + 1, k] * (qrain[i + 1, k] + qrain[i, k])
                                                           - u[i, k] * (qrain[i, k] + qrain[i - 1, k]))
                - 0.5 * dtz * (rhow[0, k + 1] * (w[i, k + 1] - 0.5 * (VT[i, k + 1] + VT[i, k])) * (
                            qrain[i, k + 1] + qrain[i, k])
                               - rhow[0, k] * (w[i, k] - 0.5 * (VT[i, k] + VT[i, k - 1])) * (
                                           qrain[i, k] + qrain[i, k - 1])) / rhou[0, k]

        for k in range(1, nz - 2):
            for i in range(1, nx - 2):
                # Compute qvtotal

                qvtot[i, k] = qb[0, k] + qvp[i, k]

        ## "Fill" any negative water produced by advection of qc, qrain

        qvtot = np.where(qvtot < 0, 0, qvtot)

        #  REMOVE ALL NEGATIVE WATER PRODUCED BY ADVECTION
        for k in range(1, nz - 2):
            for i in range(1, nx - 2):
                qvp[i, k] = max(qb[0, k] + qvp[i, k], 0.0) - qb[0, k]
                qcp[i, k] = max(qcp[i, k], 0.0)
                qrainp[i, k] = max(qrainp[i, k], 0.0)

        ## Call Kessler Mircophysics Scheme

        qvp, qcp, qrainp, thp, VT = kessler(n, qvp, qvtot, qcp, qrainp, dt, qc0, k1, k2, cp, rd, lv, p0, pib, thp, pip,
                                            rhou, rhow, tb, qb, VT, nx, nz)

        ##  Apply computational diffusion to moist variables (qv,qc,qrain)

        qvp, qvm = diffusion(np.zeros((1, nz)), qvp, qvm, cmixh, cmixv, dx, dz, dt, nx, nz)
        qcp, qcm = diffusion(np.zeros((1, nz)), qcp, qcm, cmixh, cmixv, dx, dz, dt, nx, nz)
        qrainp, qrainm = diffusion(np.zeros((1, nz)), qrainp, qrainm, cmixh, cmixv, dx, dz, dt, nx, nz)

        ## Apply Raleigh Top to Moist Variables
        for k in range(1, nz - 2):
            for i in range(1, nx - 2):
                if zu[0, k] > raydmpz:
                    qvp[i, k] = rayleighTop(zu[0, k], zu[0, nz - 1], qvp[i, k], 0, raydmpcoef, raydmpz, trigpi, i, k)
                    qcp[i, k] = rayleighTop(zu[0, k], zu[0, nz - 1], qcp[i, k], 0, raydmpcoef, raydmpz, trigpi, i, k)
                    qrainp[i, k] = rayleighTop(zu[0, k], zu[0, nz - 1], qrainp[i, k], 0, raydmpcoef, raydmpz, trigpi, i, k)

        ## Call Asselin Filter

        um, u, up = asselin(um, u, up, asscoef, nx, nz)
        wm, w, wp = asselin(wm, w, wp, asscoef, nx, nz)
        thm, th, thp = asselin(thm, th, thp, asscoef, nx, nz)
        pim, pic, pip = asselin(pip, pic, pip, asscoef, nx, nz)
        qvm, qv, qvp = asselin(qvm, qv, qvp, asscoef, nx, nz)
        qcm, qc, qcp = asselin(qcm, qc, qcp, asscoef, nx, nz)
        qrainm, qrain, qrainp = asselin(qrainm, qrain, qrainp, asscoef, nx, nz)

        ## Force top and bottom boundary

        # For w-wind

        #     for i = 1 : nx
        wp[:, 0] = 0.0
        wp[:, 1] = 0.0
        wp[:, nz - 1] = 0.0

        #     end

        ## Use zero gradient at boundaries for u, theta, p, qv,qc, and qrain at
        # top and bottom
        for i in range(1, nx - 2):
            up[i, 0] = up[i, 1]
            up[i, nz-1] = up[i, nz - 2]
            thp[i, 0] = thp[i, 1]
            thp[i, nz-1] = thp[i, nz - 2]
            pip[i, 0] = pip[i, 1]
            pip[i, nz-1] = pip[i, nz - 2]
            qvp[i, 0] = qvp[i, 1]
            qvp[i, nz-1] = qvp[i, nz - 2]
            qcp[i, 0] = qcp[i, 1]
            qcp[i, nz-1] = qcp[i, nz - 2]
            qrainp[i, 0] = qrainp[i, 1]
            qrainp[i, nz-1] = qrainp[i, nz - 2]

        #
        #     for k = 1 : nz
        #         for i = 1 : nx
        #             if (zu[k] < (0.3*xs(i)-4500)) && (xs(i) < 25000) && (xs(i) > 15000)
        #                 th[i, k] = th[i, k+1]
        #                 pic[i, k] = pic[i, k+1]
        #             end
        #             if (zu[k] < (-0.3*xs(i)+10500)) && (xs(i) > 25000) && (xs(i) < 35000)
        #                 th[i, k] = th[i, k+1]
        #                 pic[i, k] = pic[i, k+1]
        #             end
        #
        #             if (zu[k] < (0.3*xu(i)-4500)) && (xu(i) < 25000) && (xu(i) > 15000)
        #                 u[i, k] = 0
        #             end
        #             if (zu[k] < (-0.3*xu(i)+10500)) && (xu(i) > 25000) && (xu(i) < 35000)
        #                 u[i, k] = 0
        #             end
        #
        #             if (zw[k] < (0.3*xs(i)-4500)) && (xs(i) < 25000) && (xs(i) > 15000)
        #                 w[i, k] = 0
        #             end
        #             if (zw[k] < (-0.3*xs(i)+10500)) && (xs(i) > 25000) && (xs(i) < 35000)
        #                 w[i, k] = 0
        #             end
        #         end
        #     end

        # up[260:340, 1:3] = 0
        # wp[260:340, 1:3] = 0
        # up[270:330, 4:7] = 0
        # wp[270:330, 4:7] = 0
        # up[280:320, 8:11] = 0
        # wp[280:320, 8:11] = 0
        # up[290:310, 12:14] = 0
        # wp[290:310, 12:14] = 0

        # up(560:640,1:3) = 0
        # wp(560:640,1:3) = 0
        # up(570:630,4:7) = 0
        # wp(570:630,4:7) = 0
        # up(580:620,8:11) = 0
        # wp(580:620,8:11) = 0
        # up(590:610,12:14) = 0
        # wp(590:610,12:14) = 0

        ## Call Lateral Boundary Sponge

        for k in range(0, nz - 1):
            for i in range(0, nx - 1):
                up[i, k] = lateralSponge(up[i, k], ub[0, k], xu[0, i], x_dist_in, dx, nx, trigpi, latdmpcoef, 1)
                wp[i, k] = lateralSponge(wp[i, k], 0, xs[0, i], x_dist_in, dx, nx, trigpi, latdmpcoef, 0)
                thp[i, k] = lateralSponge(thp[i, k], 0, xs[0, i], x_dist_in, dx, nx, trigpi, latdmpcoef, 0)
                pip[i, k] = lateralSponge(pip[i, k], 0, xs[0, i], x_dist_in, dx, nx, trigpi, latdmpcoef, 0)
                qvp[i, k] = lateralSponge(qvp[i, k], 0, xs[0, i], x_dist_in, dx, nx, trigpi, latdmpcoef, 0)
                qcp[i, k] = lateralSponge(qcp[i, k], 0, xs[0, i], x_dist_in, dx, nx, trigpi, latdmpcoef, 0)
                qrainp[i, k] = lateralSponge(qrainp[i, k], 0, xs[0, i], x_dist_in, dx, nx, trigpi, latdmpcoef, 0)

        ## Apply periodic lateral boundaries for all k

        up[0, :] = up[nx - 1, :]
        up[nx-1, k] = up[1, k]
        wp[0, :] = wp[nx - 1, :]
        wp[nx-1, k] = wp[1, k]
        thp[0, :] = thp[nx - 1, :]
        thp[nx-1, k] = thp[1, k]
        pip[0, :] = pip[nx - 1, :]
        pip[nx-1, k] = pip[1, k]
        qvp[0, :] = qvp[nx - 1, :]
        qvp[nx-1, k] = qvp[1, k]
        qcp[0, :] = qcp[nx - 1, :]
        qcp[nx-1, k] = qcp[1, k]
        qrainp[0, :] = qrainp[nx - 1, :]
        qrainp[nx - 1, :] = qrainp[1, :]

        ## Create Output Matrix with all times

        timenow = dt * n

        if timenow >= thisdump:
        #
        #     # Set Variable name for a given timestep
        #     th_str = ['th_' num2str(n*dt)]
        #     u_str = ['u_' num2str(n*dt)]
        #     w_str = ['w_' num2str(n*dt)]
        #     pi_str = ['pi_' num2str(n*dt)]
        #     pprt_str = ['pprt_' num2str(n*dt)]
        #     qv_str = ['qv_' num2str(n*dt)]
        #     qc_str = ['qc_' num2str(n*dt)]
        #     qrain_str = ['qrain_' num2str(n*dt)]
        #
        #     # Set var_n equal to varp
        #
        #     th_str2 = [th_str '=thp']
        #     u_str2 = [u_str '=up']
        #     w_str2 = [w_str '=wp']
        #     pi_str2 = [pi_str '=pip']
        #     pprt_str2 = [pprt_str '=pprt']
        #     qv_str2 = [qv_str '=qvp']
        #     qc_str2 = [qc_str '=qcp']
        #     qrain_str2 = [qrain_str '=qrainp']
        #
        #     eval(th_str2)
        #     eval(u_str2)
        #     eval(w_str2)
        #     eval(pi_str2)
        #     eval(pprt_str2)
        #     eval(qv_str2)
        #     eval(qc_str2)
        #     eval(qrain_str2)
        #
        #     # Append New Timesteps to .mat file.
        #     full_str = ['save ' modeloutputpathandfilename ' ' th_str ' ' u_str ' ' w_str ' ' pi_str ' ' pprt_str ' ' qv_str ' ' qc_str ' ' qrain_str ' ' '-append']
        #     eval(full_str)
        #

            print("Running create_image")
            create_image(u)
            print("create_image ran")

            thisdump = thisdump + dtoutput



        ## Update Variables
        for i in range(0, nx - 1):
            for k in range(0, nz - 1):
                um[i, k] = u[i, k]
                u[i, k] = up[i, k]
                wm[i, k] = w[i, k]
                w[i, k] = wp[i, k]
                pim[i, k] = pic[i, k]
                pic[i, k] = pip[i, k]
                thm[i, k] = th[i, k]
                th[i, k] = thp[i, k]
                #             thvm[i, k] = thv[i, k]
                #             thv[i, k] = thvp[i, k]
                qvm[i, k] = qv[i, k]
                qv[i, k] = qvp[i, k]
                qcm[i, k] = qc[i, k]
                qc[i, k] = qcp[i, k]
                qrainm[i, k] = qrain[i, k]
                qrain[i, k] = qrainp[i, k]

        print("completed time step: " + str(n) + " of " + str(nt))

