def asselin( varm, var, varp, asscoef, nx, nz):

    # Apply Asselin filtering for 2dt filtering
    for k in range(1, nz-2):
        for i in range(1, nx-2):

            var[i, k] = var[i, k] + asscoef*(varp[i, k]-2*var[i, k]+varm[i, k])

    return varm, var, varp
