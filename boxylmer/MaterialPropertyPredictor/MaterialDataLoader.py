from mp_api.client import MPRester 
from abc import ABC, abstractmethod
import numpy as np
from sklearn.model_selection import train_test_split

RSEED = 1996

class AbstractDataLoader(ABC):
    test_size = 0.2
    def __init__(self):
        self.structures = None

        # cached data in self.structures. It's possible to lump this all into a dict or just to use it on the fly, 
        #   but for now, explicitly referencing them will assist when we're unsure what the best descriptors are for the target property.  
        self.lattice_vectors = None
        self.lattice_angles = None
        self.volumes = None
        self.extra_data = None

        self.band_gaps = None

        self.formulas = None

        # hmm... https://journals.aps.org/prb/abstract/10.1103/PhysRevB.89.205118 
        #   suggests there are better ways to represent a crystal lattice than coulomb matrices.
        #   This could be a good idea to explore in the future for more complex properties
        #   This paper ALSO tells me that coulomb matrices are a previous industry standard, 
        #   So it stands to reason there should be some libraries to speed this up out there. 
        #   For now, though, this method is slow when doing lattice-consistent distances and needs optimization.
        self.coulomb_matrices = None
        
        # https://par.nsf.gov/servlets/purl/10187380#:~:text=Coulomb%20Matrix%20Eigenvalues%20(CME)%20are,searches%2C%20and%20interpret%20rotational%20spectra.
        #   Very good primer on why the sorted eigenvalues are both useful and computationally efficient.
        #   Still not as good as a graph representation, but this is much faster to implement given I still have papers to grade :(.
        self.coulomb_matrix_eigenvalues = None
        self.input_length = None
        self.input_cache = None

    @abstractmethod
    def load_data(self):
        pass

    @classmethod
    def calculate_distances_fast(cls, structure):
        coords = np.array([site.coords for site in structure.sites])
        n = len(coords)
        distances = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i+1, n):
                dist = np.linalg.norm(coords[i] - coords[j])
                distances[i, j] = dist
                distances[j, i] = dist

        return distances

    @classmethod
    def calculate_distances_accurate(cls, structure):
        n = len(structure)
        distances = np.zeros((n, n))
        for i in range(n):
            for j in range(i+1, n):
                dist = structure.get_distance(i, j)
                distances[i, j] = dist
                distances[j, i] = dist

        return distances

    # coulomb matrix implemented from https://singroup.github.io/dscribe/0.3.x/tutorials/coulomb_matrix.html
    @classmethod
    def calculate_coulomb_matrix(cls, structure, distance_method):
        n = len(structure)
        mat = np.zeros((n, n))
        z = np.array([atom.specie.Z for atom in structure]) # docs: https://pymatgen.org/pymatgen.core.html#pymatgen.core.periodic_table.Element
        if distance_method == 'fast':
            distances = AbstractDataLoader.calculate_distances_fast(structure)
        elif distance_method == 'accurate':
            distances = AbstractDataLoader.calculate_distances_accurate(structure)
        else:
            raise ValueError("Invalid distance calculation method. Use \'fast\' or \'accurate\'.")
        
        z_product = z[:, np.newaxis] * z[np.newaxis, :]
        z_product = z_product.astype('float64')  
        mat = np.divide(z_product, distances, out=np.zeros_like(z_product), where=distances!=0)

        np.fill_diagonal(mat, 0.5 * z ** 2.4) # faster than using for loop directly (cProfile -> 0.181)
        return mat

    @classmethod 
    def calculate_eigenvalues(cls, coulomb_matrix):
        eigenvalues, _ = np.linalg.eigh(coulomb_matrix)
        sorted_evs = np.sort(eigenvalues)[::-1]
        return sorted_evs
    
    def calculate_max_input_size(self):
        # later we may modify these descriptors, so lets tabulate them. We also have extras!
        vol = 1
        abc = 3
        angles = 3
        eigenvals_len = max([len(s) for s in self.structures])
        extras_len = len(self.extra_data)
        return vol + abc + angles + eigenvals_len + extras_len
    
    def _validate(self):
        ndata = len(self.structures)
        for key in self.extra_data:
            assert len(self.extra_data[key] == ndata)
        # given more constraints later on, this may be expanded if the types of input data structures widen

    def _process_parsed_data(self, distance_method):
        # self.coulomb_matrices = [self.calculate_coulomb_matrix(s, distance_method='accurate') for s in self.structures]
        if distance_method == 'accurate':
            print("Distance method is accurate rather than fast, this may take some time.")
        self.coulomb_matrices = []
        total_iterations = len(self.structures) 
        percent_complete = 0
        for i, s in enumerate(self.structures):
            if distance_method == 'accurate':
                new_percent_complete = (i * 100) // total_iterations 
                if new_percent_complete > percent_complete: 
                    percent_complete = new_percent_complete
                    print(f"Coulomb Matrix Progress: {percent_complete}%", end='\r')
            self.coulomb_matrices.append(
                self.calculate_coulomb_matrix(s, distance_method)
            )

        self.coulomb_matrix_eigenvalues = [self.calculate_eigenvalues(cmat) for cmat in self.coulomb_matrices]
        self.input_length = self.calculate_max_input_size()
        self._validate()

    def get_model_inputs(self, padding_value=0):
        """
        Acquire the model input data as a procured ndarray of values in a standard (Sample, Feature) format. 

        Parameters:
        - padding_value (float, optional): Padding. Default is 0.

        Returns:
        - padded_feature_vectors (np.array): Padded inputs of shape (Sample, Feature)
        """
        if self.input_cache is None: 
            num_structures = len(self.structures)
            padded_feature_vectors = np.full((num_structures, self.input_length), padding_value, dtype=np.float32)

            for i in range(num_structures):
                feature_vector = np.hstack(
                    (
                        self.volumes[i], 
                        self.lattice_vectors[i],
                        self.lattice_angles[i], 
                        self.coulomb_matrix_eigenvalues[i]
                        # np.array([atom.specie.Z for atom in self.structures[i]], dtype=np.float32) # what if we use only the atoms? 
                    ), 
                    dtype=np.float32
                )

                padding_size = self.input_length - len(feature_vector)

                if padding_size >= 0:
                    padded_feature_vectors[i, :len(feature_vector)] = feature_vector
                else:
                    raise ValueError("Input feature vector size exceeded the available input size.")
            self.input_cache = padded_feature_vectors

        return self.input_cache
    
    def get_model_outputs(self):
        """
        Acquire the model ouptut data as a procured single dimensional ndarray. 

        Returns:
        - outputs (np.array): Onputs of shape (Float32)
        """
        return self.band_gaps

    def get_train_data(self, train_size=1-test_size, seed=RSEED):
        # these data sets are typically small, calling the function like this is not likely to bottleneck
        inputs = self.get_model_inputs()
        outputs = self.get_model_outputs()
        train_x, _, train_y, _ = train_test_split(inputs, outputs, train_size=train_size, random_state=seed)
        return train_x, train_y

    def get_test_data(self, test_size=test_size, seed=RSEED):
        inputs = self.get_model_inputs()
        outputs = self.get_model_outputs()
        _, test_x, _, test_y = train_test_split(inputs, outputs, test_size=test_size, random_state=seed)
        return test_x, test_y

    def get_input_length(self):
        return self.input_length

    def __len__(self):
        return len(self.structures)

class MPRLoader(AbstractDataLoader):
    def __init__(self):
        super().__init__()

    def load_data(self, api_key, distance_method='fast', **kwargs):
        with MPRester(api_key) as mpr:
            data = mpr.materials.summary.search(
                fields=[
                    "material_id", 
                    "band_gap", 
                    "formula_pretty", 
                    "volume", 
                    "structure"],
                **kwargs
            )

        self.structures = [d.structure for d in data]
        self.lattice_vectors = np.array([(s.lattice.a, s.lattice.b, s.lattice.c) for s in self.structures], dtype=np.float32)
        self.lattice_angles = np.array([s.lattice.angles for s in self.structures], dtype=np.float32)
        self.volumes = np.array([s.volume for s in self.structures], dtype=np.float32)
        self.formulas = [d.formula_pretty for d in data]
        self.band_gaps = np.array([d.band_gap for d in data], dtype=np.float32)


        # when pulling directly from MPR, extra data doesn't make sense as you don't know how much data will be found ahead of time
        # thus for this concrete class, extras are moot (but we could fill them out later, or just lump all our parsed data into a dict)
        self.extra_data = {} 

        self._process_parsed_data(distance_method=distance_method)
