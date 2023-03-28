import sys
import numpy as np

class Voronoi(object):
    """
    Class to enable the construction of an obtuse superbasis
    and a wigner seitz cell for any given 3D lattice.
    """
    
    _tol = 1E-7
    
    @staticmethod
    def get_voronoi(lattice):
        """
        Given a lattice defined by 3 generating vectors,
        returns the 7 distinct non-zero voronoi vectors
        in the L/2L equivalence class.
        The 7 vectors returned and their inverses are
        guaranteed to be sufficient to define a wigner
        seitz cell for the lattice.
        
        Generating vectors are column vectors of the given
        argument.
        
        Returns the 7 vectors as column vectors.
        """
        lattice = np.transpose(lattice)
        if lattice is None or isinstance(lattice, list) or len(lattice) != 3:
            raise ValueError("Expecting a 3-vector basis")
        oSBasis = np.transpose(Voronoi.obtuse_superbasis(lattice))
        to_return = [oSBasis[0], oSBasis[1], oSBasis[2], oSBasis[3],
                    [oSBasis[0][i] + oSBasis[1][i] for i in range(0, 3)],
                    [oSBasis[0][i] + oSBasis[2][i] for i in range(0, 3)],
                    [oSBasis[0][i] + oSBasis[3][i] for i in range(0, 3)]]
        return np.transpose(to_return)
        
    @staticmethod
    def obtuse_basis(lattice):
        """
        Given a lattice defined by 3 generating vectors,
        returns an obtuse basis for the lattice given by
        its three shortest vectors.
        
        Generating vectors are column vectors of the given
        argument.
        
        Returns the 3 vectors as column vectors.
        """
        to_return = [None]*3
        pool = np.transpose(Voronoi.obtuse_superbasis(lattice))
        max = -1
        max_val = 0
        for i in range(0, len(pool)):
            val = sum(pool[i][x]**2 for x in range(0, len(pool[i])))
            if val > max_val:
                max = i
                max_val = val
        j = 0
        for i in range(0, len(pool)):
            if i == max:
                continue
            to_return[j] = pool[i]
            j += 1
        return np.transpose(to_return)
        
    @staticmethod
    def obtuse_superbasis(lattice):
        """
        Given a lattice defined by 3 generating vectors,
        returns an obtuse superbasis for the lattice.
        
        Generating vectors are column vectors of the given
        argument.
        
        Returns the 4 vectors as colunmn vectors.
        """
        if len(lattice) != 3:
            raise Exception("Voronoi reduction only implemented for 3D lattices")
        lattice = np.transpose(lattice)
        temp = [None]*(len(lattice)+1)
        sum = [0]*3
        for i in range(0, len(lattice)):
            temp[i+1] = lattice[i]
            for j in range(0, 3):
                sum[j] -= lattice[i][j]
        is_obtuse = False
        while not is_obtuse:
            is_obtuse = True
            for i in range(0, len(temp)):
                if not is_obtuse:
                    break
                for j in range(i + 1, len(temp)):
                    if not is_obtuse:
                        break
                    test = sum((temp[i][k]*temp[j][k] for k in range(0, 3)))
                    if Voronoi.double_equal(test, 0, Voronoi._tol):
                        continue
                    if test > 0:
                        is_obtuse = False
                        cache = temp[0]
                        temp[0] = temp[i]
                        temp[i] = cache
                        cache = temp[2]
                        temp[2] = temp[j]
                        temp[j] = cache
                        temp[1] = [temp[1][k] + temp[0][k] for k in range(0, 3)]
                        temp[3] = [temp[3][k] + temp[0][k] for k in range(0, 3)]
                        temp[0] = [-temp[0][k] for k in range(0, 3)]
        return np.transpose(temp)
                    
    @staticmethod
    def double_equal(d1, d2):
        """
        Custom double equals method based on a tolerance
        specified in the class.
        """
        if d1 == 0 or d2 == 0:
            return abs(d1) < Voronoi._tol and abs(d2) < Voronoi._tol
        if abs(d2) < 1.0 and abs(d1) > abs(d2)*sys.float_info.max:
            return False
        if abs(d2) > 1.0 and abs(d1) < abs(d2)*sys.float_info.min:
            return False
        return ((1 - Voronoi._tol) < (d1 / d2)) and ((d1 / d2) < (1 + Voronoi._tol))
        
