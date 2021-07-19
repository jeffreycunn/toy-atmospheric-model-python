def diffusion(varb, varp, varm, cmixh, cmixv, dx, dz, dt, nx, nz):
    # Horizontal Diffusion

    kmixh = cmixh * dx * dx / dt

    for k in range(1, nz - 2):
        for i in range(1, nx - 2):
            varp[i, k] = varp[i, k] + 2.0 * dt * kmixh * (varm[i + 1, k] - 2.0 * varm[i, k] + varm[i - 1, k]) / (
                        dx * dx)

    # Vertical Diffusion
    kmixv = cmixv * dz * dz / dt

    for k in range(2, nz - 2):
        for i in range(1, nx - 2):
            varp[i, k] = varp[i, k] + 2 * dt * (kmixv * ((varm[i, k + 1] - varb[0, k + 1])
                                                         - 2.0 * (varm[i, k] - varb[0, k])
                                                         + (varm[i, k - 1] - varb[0, k - 1])) / (dz * dz))

    return varp, varm
