
from itertools import chain, dropwhile
from math import ceil, floor

from matplotlib.pyplot import Axes
import numpy as np


def plot_densities(bands: list[np.ndarray], ax: Axes):
    energies = list(chain.from_iterable(bands))

    minE = min(energies)
    maxE = max(energies)
    
    # set up the list of bins, where each bin is "named" after its lower bound
    bin_size = (maxE-minE)/50

    # TODO: pull out bin calculation to its own function?
    # the dividing and multiplying here is to round the bins down to the nearest integer multiple of bin_size
    lowest_bin_value = floor( minE / bin_size ) * bin_size
    highest_bin_value = floor( maxE / bin_size ) * bin_size

    length = round((highest_bin_value - lowest_bin_value) / bin_size)
    bins = np.linspace(lowest_bin_value, highest_bin_value, length + 1) # use length + 1 as both endpoints are included in the linspace
    counts = np.zeros(bins.shape, dtype = np.int32)

    for energy in energies:
        # overriding builtin bin() for getting binary integers; shouldn't cause issues
        # don't multiply by bin_size since we want an index in counts, not the actual bin value
        # subtracting lowest_bin_value makes sure that the lowest bin always has index 0
        bin = floor( (energy - lowest_bin_value) / bin_size )
        counts[bin] += 1

    total = sum(counts)
    density = counts/(total * bin_size)

    ax.plot(density, bins, "k")


# plot_densities([np.array([0, 20, 77, 5, 5.1])], None)
