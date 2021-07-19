import math

# function [ tb,qb,qbs,tbv,pisfc, pib, rhou, rhow,rhb, ub,um,u ] = ...
def base( profile_method, nx, nz, dz, psurf, qsurf, q4km, ztr, temptr, ttr, tsurf, p0, cp, g,
          rd, xk, c_v, zu, rhow, rhou, tb, tbv, qb, qbs, rhb, pib, ub, um, u, up):


    if profile_method == 1:
        #########################################
        # Define a dry and neutral environment

        # Set base state theta to 300 K and mixing ratio to 0 kg/kg

        for k in range(1, nz - 1):
            tb[0, k] = 300.0
            qb[0, k] = 0.0
            tbv[0, k] = tb[0, k]*(1 + 0.61*qb[0, k])



        # Assign pi
        pisfc = (psurf / p0)**xk
        pib[0, 2] = pisfc - g*0.5*dz/(cp*tbv[0, 2])

        for k in range(2, nz - 1):

            tbvavg[0, k] = 0.5*(tbv[0, k] + tbv[0, k-1])
            pib[0, k] = pib[0, k-1] - g*dz/(cp*tbvavg[0, k])


        for k in range(1, nz - 1):
        #  Assign mean state density at u point

            rhou[k] = (p0*pib[0, k]**(c_v/rd))/(rd*tbv[0, k])

            # Assign mean state density at w point

            rhow[0, k] = 0.5*(rhou[0, k] + rhou[0, k-1])


        # #Ensure there is no gradient at the boundaries
        # tb[0, 0] = tb[0, 1]
        # tb[0, nz] = tb[0, nz-2]
        # pib[0, 0] = pib[0, 1]
        # pib[0, nz] = pib[0, nz-2]
        # rhou[0, 0] = rhou[0, 1]
        # rhou[0, nz] = rhou[0, nz-2]
        # rhow[0, 1] = 1.15
        # rhow[0, 0] = rhow[0, 1]
        # rhow[0, nz] = rhow[0, nz-2]

    elif profile_method == 2:
        ## fovellwk82soundingdef
        # define a sounding on the u/scalar points

        # Assign Values to vertical grid

        for k in range(1, nz - 1):

        # Assign water vapor mixing ratio
            if zu[0, k] <= 4000.0:
                qb[0, k] = qsurf - (qsurf - q4km)*zu[0, k]/4000.0
            elif zu[0, k] <= 8000.0:
                qb[0, k] = q4km - q4km*(zu[0, k]-4000.0)/4000.0
            else:
                qb[0, k] = 0.0


        # Assign potential temperature
            if zu[0, k]<= ztr:
                tb[0, k] = tsurf + (ttr - tsurf)*(zu[0, k]/ztr)**1.25
            else:
                tb[0, k] = ttr*math.exp(g*(zu[0, k] - ztr)/(cp*temptr))


            # Assign virtual potential temperature
            tbv[0, k] = tb[0, k]*(1 + 0.61*qb[0, k])



        # Assign pi
        pisfc = (psurf / p0)**xk
        pib[0, 2] = pisfc - g*0.5*dz/(cp*tbv[0, 2])

        for k in range(2, nz - 1):

            tbvavg[0, k] = 0.5*(tbv[0, k] + tbv[0, k-1])
            pib[0, k] = pib[0, k-1] - g*dz/(cp*tbvavg[0, k])



        #  Assign mean state density at u point
        for k in range(1, nz-1):

            rhou[0, k] = (p0*pib[0, k]**(c_v/rd))/(rd*tbv[0, k])


            # Assign mean water vapor saturation mixing ratio

            t[0, k] = tb[0, k]*pib[0, k]
            p[0, k] = p0*pib[0, k]**(cp/rd)
            qbs[0, k] = (380.0/p[0, k])*math.exp(17.27*(t[0, k]-273.0)/(t[0, k]-36.0))

            # Calculate and Assign RH

            rhb[0, k] = qb[0, k]/qbs[0, k]

            #  define density at the true surface from known information
        pisfc = (psurf/p0)**xk
        rhow[0, 2] = p0*(pisfc**(c_v/rd))/(rd*tsurf)

        #  define density at the other w points by interpolating from u/scalar points
        for k in range(2, nz-1):
            rhow[0, k] = 0.5*(rhou[0, k] + rhou[0, k-1])


    elif profile_method == 3:


        #  define the Weisman Klemp sounding on the u/scalar points

        for k in range(1, nz-1):
            #  assign potential temperatures and relative humidities
            if zu[0, k] <= ztr:
                tb[0, k] = tsurf + (ttr-tsurf)*(zu[0, k]/ztr)**1.25
                rhb[0, k] = 1 - 0.02*(zu[0, k]/ztr)**1.25
            elif zu[0, k] > ztr:
                tb[0, k] = ttr*math.exp(g*(zu[0, k]-ztr)/(cp*temptr))
                rhb[0, k] = 0.98

            #  first guess sans moisture
            tbv[0, k] = tb[0, k]


        #  compute the non-dimensional pressure profile as a first guess
        #    (needed for qvs below)
        pisfc = (psurf / p0)**xk
        pib[0, 2] = pisfc - g*0.5*dz/(cp*tbv[0, 2])

        for k in range(2, nz - 1):
            tbvavg = 0.5*(tbv[0, k] + tbv[0, k-1])
            pib[0, k] = pib[0, k-1] - g*dz/(cp*tbvavg)

        #  convert to mixing ratios
        for k in range(1, nz-1):
            p = p0*pib[0, k]**(cp/rd)
            temp = tb[0, k]*pib[0, k]
            qvs = (380.0/p) * math.exp((17.27*(temp-273.0))/(temp-36.0))
            #             if zu[k] < 20000
            qb[0, k] = min(qsurf, qvs*rhb[0, k])
            #             else
            #             qb[k] = 0
            #             end
            #  define virtual potential temperatures
            tbv[0, k] = tb[0, k]*(1 + 0.61*qb[0, k])

            # if k==98 or k == 99 or k == 100:
            #     display[k]
            #     display(p)
            #     display(temp)
            #     display(qvs)
            #     display('qb')
            #     display(qb[k])
            #     display('tbv')
            #     display(tbv[k])
            #     display('rhb')
            #     display(rhb[k])


        #  go back and recompute the non-dimensional pressure profile now that
        #    we have humidity included
        pisfc = (psurf / p0)**xk
        pib[0, 2] = pisfc - g*0.5*dz/(cp*tbv[0, 2])

        for k in range(2, nz - 1):

            tbvavg = 0.5*(tbv[0, k] + tbv[0, k-1])
            pib[0, k] = pib[0, k-1] - g*dz/(cp*tbvavg)


        #  compute the density profile
        pisfc = (psurf/p0)**xk
        rhow[0, 2] = p0*(pisfc**(c_v/rd))/(rd*tsurf)

        #  define density at the other w points by interpolating from u/scalar points
        #  define density at u/scalar points
        for k in range(1, nz-1):
            rhou[0, k] = p0*(pib[0, k]**(c_v/rd))/(rd*tbv[0, k])


        #  define density at the true surface from known information
        pisfc = (psurf/p0)**xk
        rhow[0, 2] = p0*(pisfc**(c_v/rd))/(rd*tsurf)

        #  define density at the other w points by interpolating from u/scalar points
        for k in range(2, nz-1):
            rhow[0, k] = 0.5*(rhou[0, k] + rhou[0, k-1])

        #  compute the profiles of saturation mixing ratio and relative humidity
        #  define saturation mixing ratio on u/scalar points
        for k in range(1, nz-1):
            p = p0*pib[0, k]**(cp/rd)
            temp = tb[0, k]*pib[0, k]
            #  use Teten's formula (see Pielke's text, 2nd ed., p. 257-258 for more info)
            qbs[0, k] = (380.0/p) * math.exp((17.27*(temp-273.0))/(temp-36.0))
            rhb[0, k] = qb[0, k] / qbs[0, k]



 # End Method If Statement

    #Ensure there is no gradient at the boundaries
    tb[0, 0] = tb[0, 1]
    tb[0, nz-1] = tb[0, nz-2]
    tbv[0, 0] = tbv[0, 1]
    tbv[0, nz-1] = tbv[0, nz-2]
    pib[0, 0] = pib[0, 1]
    pib[0, nz-1] = pib[0, nz-2]
    rhou[0, 0] = rhou[0, 1]
    rhou[0, nz-1] = rhou[0, nz-2]
    rhow[0, 0] = rhow[0, 1]
    rhow[0, nz-1] = rhow[0, nz-2]
    qb[0, 0] = qb[0, 1]
    qb[0, nz-1] = qb[0, nz-2]

    ## Define U base

    # Uniform Flow

    ub[:] = -10
    um[:, :] = -10
    u[:, :] = -10

    ############################################################
    # High Surface Wind decreasing as altitude increases \

                                               #     usfc=-10.0  # u-wind at the surface (m/s)
                                                                                  #     # * for initwinds=1 (only)
    #     deltau=10.0  # u-wind change over the layer (m/s)
    #     # * for initwinds=1 (only)
    #     shrdepth=3000.0  # depth of deltau layer (m)
    #     #  * for initwinds=1 (only)
    #
    #
    #     for k=2:nz-1
                  #         #  assign base state values
                                                 #         if zu[k]<shrdepth
                                                              #             ub[k] = usfc + zu[k]*deltau/shrdepth
    #         else
    #             ub[k] = usfc + deltau
    #         end
              #         for i= 1 :nx
                                  #             #  copy base state values into prognosed u array \
                                                                                           #             u(i,k) = ub[k]
    #             um(i,k) = ub[k]
    #         end
              #     end

              ###################################################################

              # Low-Level Jet at a given z height (given by shrdepth)
    #
    #     usfc=0  # u-wind at the surface (m/s)
                                   #     # * for initwinds=1 (only)
    #     deltau=-10.0  # u-wind change over the layer (m/s)
    #     # * for initwinds=1 (only)
    #     shrdepth=2000.0  # depth of deltau layer (m)
    #     #  * for initwinds=1 (only)
    #
    #     for k=2:nz-1
                  #         #  assign base state values
                                                 #         if zu[k]<=shrdepth
                                                              #             ub[k] = usfc + zu[k]*deltau/shrdepth
    #         elseif zu[k] <= 2*shrdepth \
                     #             ub[k] = 2*deltau  - deltau*(zu[k]/(shrdepth))
    #         else
    #             ub[k] = 0
    #         end
              #     for i= 1 :nx
                              #             #  copy base state values into prognosed u array \
                                                                                       #             u(i,k) = ub[k]
    #             um(i,k) = ub[k]
    #      end
           #      end
           ######################################################################
           #  handle boundaries


    for i in range(1, nx):
        u[i, 0] = ub[0, 1]
        up[i, 0] = ub[0, 1]
        um[i, 0] = ub[0, 1]
        u[i, nz-1] = ub[0, nz-2]
        up[i, nz-1] = ub[0, nz-2]
        um[i, nz-1] = ub[0, nz-2]

    return tb, qb, qbs, tbv, pisfc, pib, rhou, rhow, rhb, ub, um, u
