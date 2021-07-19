import math

def rayleighTop(varz,varztop,varp,varb,raydmpcoef,raydmpz,trigpi,i,k):

    #   Set varb to 0 for w,th,and pic

    # Compute Local Damping Coefficient
    coef = raydmpcoef*0.5*(1-math.cos(trigpi*(varz-raydmpz)/(varztop-raydmpz)))

    # Apply Sponge to Slowly Remove Perturbations
    varp = varp - coef*(varp-varb)

    return varp
