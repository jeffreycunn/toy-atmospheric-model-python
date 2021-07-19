##  Miscellaneous User Settings

nz = 82 	# The number of grid points in the vertical
nx = 100  # The number of grid points in the horizontal
nt = 50  # Number of timesteps
dz = 250.0	# Vertical grid spacing: m
dx = 500.0  # Horizontal Grid Spacing: m
dt = 1  # Time Step secs
c_sound = 50.0  # Speed of sound m/s
asscoef = 0.1 # Asselin coefficient
cmixh = 0.003 # Horizontal diffusion Coefficient
cmixv = 0.003 # Vertical diffusion coefficient
# ub = 0  # Mean Environmental wind
profile_method = 3 # 1 - Dry Neutral Sounding, 2 - Fovell Weisman-Klemp (1982) Sounding, 3 - Moist Stable

## Kessler Microphysics Settings

qc0 = .001 # Autoconversion Threshold
k1 = .001 # Autoconversion Rate (s^-1)
k2 = 2.2 # Accretion Rate (s^-1)

## Initial Condition Bubble Settings
bubble_switch = 0 # 1 - Bubble On, 0 - Bubble Off
radz = 3000.0  # Vertical Radius (m)
radx = 5000.0 # Horizontal Radius (m)
zcnt = 2000.0  # Center of thermal (meters above ground)
xcnt = 82000 # Horizontal Center (m)
thermamp = -6.0 # Thermal Amplitude (K)

## Rayleigh Damping Top Settings
raydmpcoef = 10/100  # Rayleigh Damping Coefficient
raydmpz = 14000 # Damping layer extends down to xx (meters)

## Lateral Sponge Boundary
sponge_dist = .15  # 10 Percent of Domain Width
x_dist_in = sponge_dist*nx*dx
latdmpcoef = 0.9
##  Parameters to describe the model environmental sounding

if profile_method == 2:

    # surface pressure, 96500 N/m^2
    psurf = 96500.0
    # potential temp. at the surface (K), default is 300.0 K
    tsurf = 300.0
    # mixing ratio at the surface (kg/kg), default is 0.0161
    qsurf = 0.0161
    # mixing ratio @ 4 km AGL (kg/kg), default is 0.001
    q4km = 0.001
    # height of the tropopause (m), default is 12000.0 m
    ztr = 12000.0
    # temperature at the tropopause (K), default is 213.0 K
    temptr = 213.0
    # potential temp. at the tropopause (K), default is 343.0 K
    ttr = 343.0

elif profile_method == 3:

    # surface pressure, 96500 N/m^2
    psurf = 96500.0
    # potential temp. at the surface (K), default is 300.0 K
    tsurf = 281.0
    # mixing ratio at the surface (kg/kg), default is 0.0161
    qsurf = 0.006
    # mixing ratio @ 4 km AGL (kg/kg), default is 0.001
    q4km = 0.001
    # height of the tropopause (m), default is 12000.0 m
    ztr = 12000.0
    # temperature at the tropopause (K), default is 213.0 K
    temptr = 218
    # potential temp. at the tropopause (K), default is 343.0 K
    ttr = 345.0


