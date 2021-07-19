import math

def kessler(n, qvp,qvtot,qcp,qrainp,dt,qc0,k1,k2,cp,rd,lv,p0,pib,thp,pip,rhou,rhow,tb,qb,VT,nx,nz ):

    ## Rain Water Conversions

    for k in range(1, nz-2):
        for i in range(1, nx-2):

            autodt=0.0
            accrdt=0.0
            revapdt=0.0

            #  AUTOCONVERSION OF CLOUD TO RAIN

            if(qcp[i, k]>0.0):
                auto = k1*max(qcp[i, k]-qc0, 0.0)
                autodt = auto*2*dt


            # ACCRETION OF CLOUD BY RAIN

            if qcp[i, k]>0.0 and qrainp[i, k] > 0.0:
                accr = rhou[k]*k2*qcp[i, k]*(qrainp[i, k]**0.875)
                accrdt = accr*2*dt


            # limit amount of cloud water we can lose

            if autodt+accrdt>qcp[i, k]:
                qcsink=autodt+accrdt
                autodt=qcp[i, k]*autodt/qcsink
                accrdt=qcp[i, k]*accrdt/qcsink


            # EVAPORATION OF RAIN

            if(qrainp[i, k]> 0.0):
                cvent = 1.6 + 30.39*(rhou[k]*qrainp[i, k])**0.2046
                pprime = pip[i, k]*cp*rhou[k]*tb[k]
                pbar = p0*pib[k]**(cp/rd)
                pres = pbar + pprime
                temp = (thp[i, k]+tb[k])*(pip[i, k]+pib[k])
                # Teten's formula
                qvsat = (380.0/pres) * math.exp((17.27*(temp-273.0))/(temp-36.0))

                if qvsat>(qb[k]+qvp[i, k]):

                    revap = cvent *(1.0-((qvp[i, k]+qb[0, k])/qvsat))*((rhou[k]*qrainp[i, k])**0.525)/\
                            ( (2.03e4+9.584e6/(qvsat*pres)) * rhou[0, k] )

            #  limit amount of rain we can evaporate

                    revapdt = min(revap*2*dt,qrainp[i, k])

                    if (qvp[i, k]+qb[k]+revapdt)>qvsat:
                        revapdt=qvsat -(qvp[i, k]+qb[0, k])


            qvp[i, k] = qvp[i, k] + revapdt
            qcp[i, k] = qcp[i, k] - (autodt+accrdt)
            qrainp[i, k] = qrainp[i, k] + (autodt+accrdt) - revapdt
            thp[i, k] = thp[i, k] - revapdt*lv/(cp*pib[0, k])
            # print(str(i) + ", " + str(k))


            ## Compute Rainwater fall speeds for each grid box
            #(to be used in qr next time step)

    for k in range(0, nz-1):
        for i in range(0, nx-1):

            VT[i, k] = 0.0
            if qrainp[i, k]>1.0*10**(-12):
                qrr = max(qrainp[i, k]*0.001*rhou[0, k], 0.0)
                vtden = math.sqrt(rhou[0, 1]/rhou[0, k])
                VT[i, k] = 36.34*(qrr**0.1364)*vtden


            ## Computation of condensation/evaporation of cloud mass
            #

            # CONDENSATION/EVAPORATION OF CLOUD MASS
    for k in range(1, nz-2):
        for i in range(1, nx-2):
            pprime = pip[i, k]*cp*rhou[0, k]*tb[0, k]
            pbar = p0*pib[0, k]**(cp/rd)
            pres = pbar + pprime
            temp = (thp[i, k]+tb[0, k])*(pip[i, k]+pib[0, k])
            # Teten's formula
            qvsat = (380.0/pres) *math.exp((17.27*(temp-273.0))/(temp-36.0))

            if qcp[i, k] > 0 or (qvp[i, k]+qb[0, k] > qvsat):

                #  compute the adjustment to total evap/cond owing to slight change in
                #    qvsat by the latent heating/chilling
                phi = qvsat*lv*17.27*237.0/(cp * (temp-36.0)**2)
                dqv = (qvsat-qvp[i, k]-qb[0, k])/(1+phi)

                #  if supersaturated, dqv is negative if subsaturated, dqv is positive
                #     for the subsaturated case, only allow evaporation up to qc, and no more
                dqv = min(dqv, qcp[i, k])

                qvp[i, k] = qvp[i, k] + dqv
                qcp[i, k] = qcp[i, k] - dqv
                thp[i, k] = thp[i, k] - dqv*lv/(cp*pib[0, k])


    return qvp, qcp, qrainp, thp, VT
