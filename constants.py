import numpy as np
import dask.array as da
import user_settings as us

#########################################################################
#  Physical and mathematical constants

g = 9.8   # gravitational acceleration: m/s^2
cp = 1004.0  # specific heat at constant p: J/(kg K)
rd = 287.0  # dry gas constant:  J/(kg K)
c_v = cp-rd  # specific heat at constant V: J/(kg K)
xk = rd/cp  # kappa, the ratio of rd/cp
p0 = 100000.0  # reference pressure:  Pa
lv = 2500000.0  # latent heat of vaporization:  J/kg
trigpi = np.pi # trigonmetric pi
rhol = 1000.0 # Density of Rain drop (1000 kg/m^3)

#########################################################################

## Initialize Variable Arrays

zw = np.zeros((1, us.nz))
zu = np.zeros((1, us.nz))
xu = np.zeros((1, us.nx))
xs = np.zeros((1, us.nx))
rhow = np.zeros((1, us.nz))
rhou = np.zeros((1, us.nz))
tbv = np.zeros((1, us.nz))
ub = np.zeros((1, us.nz))
tb = np.zeros((1, us.nz))
qb = np.zeros((1, us.nz))
pib = np.zeros((1, us.nz))
qbs = np.zeros((1, us.nz))
rhb = np.zeros((1, us.nz))
qvm = np.zeros((us.nx, us.nz))
qv = np.zeros((us.nx, us.nz))
qvp = np.zeros((us.nx, us.nz))
qvtot = np.zeros((us.nx, us.nz))
qcm = np.zeros((us.nx, us.nz))
qc = np.zeros((us.nx, us.nz))
qcp = np.zeros((us.nx, us.nz))
qrainm = np.zeros((us.nx, us.nz))
qrain = np.zeros((us.nx, us.nz))
qrainp = np.zeros((us.nx, us.nz))
thm = np.zeros((us.nx, us.nz))
th = np.zeros((us.nx, us.nz))
thp = np.zeros((us.nx, us.nz))
thvm = np.zeros((us.nx, us.nz))
thv = np.zeros((us.nx, us.nz))
thvp=np.zeros((us.nx, us.nz))
um = np.zeros((us.nx, us.nz))
u = np.zeros((us.nx, us.nz))
up = np.zeros((us.nx, us.nz))
wm = np.zeros((us.nx, us.nz))
w = np.zeros((us.nx, us.nz))
wp = np.zeros((us.nx, us.nz))
pim = np.zeros((us.nx, us.nz))
pic = np.zeros((us.nx, us.nz))
pip = np.zeros((us.nx, us.nz))
pprt = np.zeros((us.nx, us.nz))


