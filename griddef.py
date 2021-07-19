from google.cloud import logging
import numpy as np

# Instantiates a client
# logging_client = logging.Client()
#
# # The name of the log to write to
# log_name = "cmm-python-griddef"
#
# # Selects the log to write to
# logger2 = logging_client.logger(log_name)
def griddef(nz, nx, dz, dx, zu, zw, xu, xs):
    # Assign height values to u-scalar heights
    print("running griddef")
    # logger2.log_text("griddef log....")
    # for k in zu:
    for k in range(1, nz):
        zu[0, k] = (k - 1.5) * dz
        # logger2.log_text(zu[0, k])

    # Assign height values to w heights

    for k in range(1, nz):
        zw[0, k] = dz * (k - 2)

    # Assign x-dist values to u on the x grid-pt

    for i in range(1, nx):
        xu[0, i] = dx * (i - 2)
    #
    # # Assign x-dist values to scalars on the x grid-pt
    #
    for i in range(1, nx):
        xs[0, i] = (i - 1.5) * dx

    return zu, zw, xu, xs
