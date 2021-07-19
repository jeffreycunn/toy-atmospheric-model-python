import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def create_image(data_array):

    fig, ax = plt.subplots()
    im = ax.imshow(data_array, interpolation='bilinear', cmap=cm.RdYlGn,
                   origin='lower')

    plt.savefig('./static/test.png')
