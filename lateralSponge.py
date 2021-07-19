import math


def lateralSponge(varp, varb, x_loc, x_dist_in, dx, nx, trigpi, latdmpcoef, u_switch):

    #  by default, make the sponging layer 10# of the total domain width
    xwidth=0.1*(nx-2)*dx

    coef=0.0

    #  if we're within the west sponge zone, then apply damping
    if x_loc <= xwidth:
        coef = 0.5*(1.-math.cos(trigpi*(xwidth-x_loc))/xwidth)

    #  if we're within the east sponge zone, then apply damping
    if x_loc >= ((nx-2)*dx-xwidth):
        coef = 0.5*(1.-math.cos(trigpi*(x_loc-(nx-2)*dx-xwidth))/xwidth)

    #  apply sponge to slowly remove perturbations
    if coef > 0.0:
        varp = varp - coef*(varp-varb)

    return varp
