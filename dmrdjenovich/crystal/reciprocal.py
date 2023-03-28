import numpy as np

def get_reciprocal(lattice):
    """
    Given an input lattice as a matrix with column vectors
    corresponding to lattice generators, returns its
    corresponding reciprocal lattice generators as column
    vectors in another matrix.
    """
    raw = np.linalg.inv(np.transpose(lattice))*2*np.pi
    out = list(raw)
    for i in range(0, len(out)):
        out[i] = list(out[i])
    return out
