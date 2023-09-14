from mp_api.client import MPRester 
from abc import ABC, abstractmethod
import numpy as np
from sklearn.model_selection import train_test_split
from pymatgen.core import Element 


RSEED = 1996

class AbstractDataLoader(ABC):
    test_size = 0.2
    def __init__(self, n_eigenvals=None):
        self.structures = None

        # cached data in self.structures. It's possible to lump this all into a dict or just to use it on the fly, 
        #   but for now, explicitly referencing them will assist when we're unsure what the best descriptors are for the target property.  
        self.lattice_vectors = None
        self.lattice_angles = None
        self.volumes = None
        self.densities = None
        self.atomic_densities = None
        self.is_metal = None
        self.natoms = None
        self.energies_per_atom = None

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
        self.n_eigenvals = n_eigenvals

    @abstractmethod
    def load_data(self):
        """
        This method fetches materials data / parses it. Implementation arguments will change.
        """

        pass

    @classmethod
    def _calculate_distances_fast(cls, structure):
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
    def _calculate_distances_accurate(cls, structure):
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
        """
        Calculates the Coulomb matrix for a given atomic structure using either a 'fast' 
        or 'accurate' distance calculation method.

        The Coulomb matrix is a representation of an atomic structure, capturing the
        pairwise Coulombic interactions between atoms in a matrix form.
        Implementation is based on: https://singroup.github.io/dscribe/0.3.x/tutorials/coulomb_matrix.html

        Parameters:
            structure: A pymatgen Structure object representing the atomic structure.
            distance_method: The method for distance calculation between atoms. Supply 'fast' or 'accurate'.
                - Fast will use the immediate cartesian distance of the structure, which will lose important 
                information about the periodic nature of the lattice. Use this for debugging.
                - Accurate will use pymatgens internal distance calculation, but will take around 1000 times 
                longer to process when loading. Use this for final model calculations.
        Returns:
            np.ndarray: The Coulomb matrix of shape (n, n), where n is the number of atoms in the structure.

        Raises:
            ValueError: If an invalid distance calculation method is provided. 
        """ # todo should I raise ValueError or some other? 
        n = len(structure)
        mat = np.zeros((n, n))
        z = np.array([atom.specie.Z for atom in structure]) # docs: https://pymatgen.org/pymatgen.core.html#pymatgen.core.periodic_table.Element
        if distance_method == 'fast':
            distances = AbstractDataLoader._calculate_distances_fast(structure)
        elif distance_method == 'accurate':
            distances = AbstractDataLoader._calculate_distances_accurate(structure)
        else:
            raise ValueError("Invalid distance calculation method. Use \'fast\' or \'accurate\'.")
        
        z_product = z[:, np.newaxis] * z[np.newaxis, :]
        z_product = z_product.astype('float64')  
        mat = np.divide(z_product, distances, out=np.zeros_like(z_product), where=distances!=0)

        np.fill_diagonal(mat, 0.5 * z ** 2.4) # faster than using for loop directly (cProfile -> 0.181)
        return mat

    @classmethod 
    def calculate_eigenvalues(cls, coulomb_matrix):
        "Calculate the eigenvalues of a coulomb matrix."
        eigenvalues, _ = np.linalg.eigh(coulomb_matrix)
        sorted_evs = np.sort(eigenvalues)[::-1]
        return sorted_evs

    @classmethod
    def calculate_coordination_stats(cls, structure, radius):
        """
        Calculate the average and standard deviation of 
        the coordination numbers in a given radius for all 
        of the atoms in a pymatgen structure object.
        
        Parameters:
            structure (pymatgen.Structure): structure to analyze.
            radius (float): The radius (in angstroms) to search for coordinating atoms.
            
        Returns:
            tuple: (average, standard deviation) of coordination numbers.

        # todo I have not been able to back this strategy from the literature yet: Search for better ways to encode constant-size structural info. 
        """
        total_coordination = 0
        coordination_sqr = 0
        
        n = len(structure)
        for site in structure:
            neighbors = structure.get_neighbors(site, r=radius)
            
            coordination_number = len(neighbors)
            total_coordination += coordination_number
            coordination_sqr += coordination_number ** 2
        
        avg_coordination = total_coordination / n
        variance = (coordination_sqr / n) - (avg_coordination ** 2)
        std_deviation = np.sqrt(variance)
        
        return avg_coordination, std_deviation

    @classmethod
    def create_element_position_map(cls, structures):
        """
        Given all structures in a dataset

        Parameters:
            structures (list): A list of pymatgen Structures.

        Returns:
            dict: A dictionary where (Element string -> integer position).
                
        Example:
            >>> position_map = create_element_position_map([structure1, structure2])
            >>> print(position_map)
            {'H': 0, 'O': 1, 'Fe': 2, ...}
        
        """
        unique_elements = set()
        for structure in structures:
            for site in structure:
                e = site.specie.symbol
                unique_elements.add(e)

        sorted_elements = sorted(list(unique_elements), key=lambda x: Element(x).Z)
        element_position_map = {element: index for index, element in enumerate(sorted_elements)}

        return element_position_map

    @classmethod
    def create_element_fraction_vector(cls, structure, element_position_map):
        """
        Generate a vector whose indices are mapped by the element_position_map (element -> position).
     
        Parameters:
            structure (pymatgen.Structure): The structure of the material.
            position_map (dict): A dictionary mapping element strings to positions in the feature vector (see `create_element_position_map`).
            
        Returns:
            np.array: ndarray holding the fractional occurrence of each element in the structure, reflecting the mapping.
            
        Example:
            >>> # assume you have 20 of one atom, 80 of another, and nothing else.
            >>> fractional_vector = create_fractional_occurrence_vector(structure, element_position_map)
            >>> print(fractional_vector)
            [0.2, 0.8, 0, ...]
        """
        total_atoms = len(structure)
        element_vector = np.zeros(len(element_position_map))

        for site in structure:
            e = site.specie.symbol
            position = element_position_map[e]
            element_vector[position] += 1

        element_vector /= total_atoms
        return element_vector

    # todo this should actually be parametrized; a loader should have options as to what information is being collected ideally.
    def calculate_max_input_size(self):  
        "Calculate the expected total size of the input feature vector for each sample."
        # later we may modify these descriptors, so lets tabulate them. We also have extras!
        vol = 1
        abc = 3
        angles = 3
        density = 1
        atomic_density = 1
        is_metal = 1
        natoms = 1
        energy_per_atom = 1
        coordination_stats = 2 * 3
        if self.n_eigenvals is None:
            eigenvals_len = max([len(s) for s in self.structures])
        else:
            eigenvals_len = self.n_eigenvals
        elemental_composition_length = len(self.elemental_composition_map)
        
        extras_len = len(self.extra_data)


        n = vol + abc + angles + density + atomic_density + is_metal + natoms + energy_per_atom + coordination_stats + elemental_composition_length + eigenvals_len + extras_len 

        return n
    
    def _validate(self):
        ndata = len(self.structures)
        for key in self.extra_data:
            assert len(self.extra_data[key] == ndata)
        # given more constraints later on, this may be expanded if the types of input data structures widen

    def _process_parsed_data(self, distance_method):
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
        if distance_method == 'accurate': 
            print("Coulomb Matrix Progress: Complete")
        
        self.coulomb_matrix_eigenvalues = [self.calculate_eigenvalues(cmat) for cmat in self.coulomb_matrices]
        self.elemental_composition_map = self.create_element_position_map(self.structures)
        self.elemental_fraction_vectors = [self.create_element_fraction_vector(s, self.elemental_composition_map) for s in self.structures]

        self.input_length = self.calculate_max_input_size()
        self._validate()

    @classmethod
    def resize_eigenvals(cls, arr, new_length):
        "Pad the eigenvalue vector with zeros to a total length."
        current_length = len(arr)
        if current_length > new_length:
            return arr[:new_length]
        elif current_length < new_length:
            return np.pad(arr, (0, new_length - current_length), 'constant', constant_values=0)
        else:
            return arr

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

                if self.n_eigenvals is None:
                    coulomb_matrix_eigenvalues = self.coulomb_matrix_eigenvalues[i]
                else:
                    coulomb_matrix_eigenvalues = self.resize_eigenvals(self.coulomb_matrix_eigenvalues[i], self.n_eigenvals)

                
                feature_vector = np.hstack(
                    (
                        self.volumes[i], 
                        self.lattice_vectors[i],
                        self.lattice_angles[i], 
                        self.densities[i],
                        self.atomic_densities[i],
                        self.is_metal[i],
                        self.natoms[i],
                        self.energies_per_atom[i],
                        self.calculate_coordination_stats(self.structures[i], 3),
                        self.calculate_coordination_stats(self.structures[i], 2),
                        self.calculate_coordination_stats(self.structures[i], 1),
                        self.elemental_fraction_vectors[i],
                        coulomb_matrix_eigenvalues,
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
        - outputs (np.array): Outputs of shape (Float32)
        """
        return self.band_gaps

    def get_train_data(self, train_size=1-test_size, seed=RSEED):
        """
        Get data used in training the model. 

        Parameters:
            train_size (float): Amount of the dataset to include in the training split. Ranges from 0 to 1.
                Defaults to 1 - test_size. (defined in this class)
            seed (int): Random seed used for reproducibility. Defaults to RSEED, a constant defined in this module.
       
        
        Returns:
            tuple: (train_x, train_y) where train_x is the feature matrix (samples, features) for training and
                train_y is the target band gaps for training.
        """
        # these data sets are typically small, calling the function like this is not likely to bottleneck
        inputs = self.get_model_inputs()
        outputs = self.get_model_outputs()
        train_x, _, train_y, _ = train_test_split(inputs, outputs, train_size=train_size, random_state=seed)
        return train_x, train_y

    def get_test_data(self, test_size=test_size, seed=RSEED):
        """
        Get data used in testing the model. 

        Parameters:
            test_size (float): Amount of the dataset to include in the testing split. Ranges from 0 to 1.
                Defaults to test_size. (defined in this class)
            seed (int): Random seed used for reproducibility. Defaults to RSEED, a constant defined in this module.
            
        Returns:
            tuple: (train_x, train_y) where train_x is the feature matrix (samples, features) for training and
                train_y is the target band gaps for training.
        """
        inputs = self.get_model_inputs()
        outputs = self.get_model_outputs()
        _, test_x, _, test_y = train_test_split(inputs, outputs, test_size=test_size, random_state=seed)
        return test_x, test_y

    def get_input_length(self):
        return self.input_length

    def __len__(self):
        return len(self.structures)

class MPRLoader(AbstractDataLoader):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def load_data(self, api_key, distance_method='fast', **kwargs):
        """
        Implements load_data from AbstractDataLoader.
        
        Parameters:
            api_key (str): The API key for accessing MPR.
            distance_method (str, optional): The method used for distance calculations. Defaults to 'fast'. 
                Using 'accurate' instead will cause a progress meter to be printed. Expect very slow parse times for datasets > 500 items in size.
            **kwargs: Search criteria in MPR (Forwarded to MPRester).
            
        Returns: None
        """
        with MPRester(api_key) as mpr:
            data = mpr.materials.summary.search(
                fields=[
                    "material_id", 
                    "band_gap", 
                    "formula_pretty", 
                    "volume",
                    "density",
                    "density_atomic", 
                    "structure",
                    "is_metal",
                    "energy_per_atom"],
                **kwargs
            )

        self.structures = [d.structure for d in data]
        self.lattice_vectors = np.array([(s.lattice.a, s.lattice.b, s.lattice.c) for s in self.structures], dtype=np.float32)
        self.lattice_angles = np.array([s.lattice.angles for s in self.structures], dtype=np.float32)
        self.volumes = np.array([s.volume for s in self.structures], dtype=np.float32)
        self.densities = np.array([d.density for d in data], dtype=np.float32)
        self.atomic_densities = np.array([d.density_atomic for d in data], dtype=np.float32)
        self.is_metal = np.array([d.is_metal for d in data], dtype=np.float32) # implicit bool -> float for ndarray, I think
        self.natoms = np.array([len(s) for s in self.structures], dtype=np.float32)
        self.energies_per_atom = np.array([d.energy_per_atom for d in data], dtype=np.float32)
        

        self.formulas = [d.formula_pretty for d in data]
        
        
        self.band_gaps = np.array([d.band_gap for d in data], dtype=np.float32)


        # when pulling directly from MPR, extra data doesn't make sense as you don't know how much data will be found ahead of time
        # thus for this concrete class, extras are moot (but we could fill them out later, or just lump all our parsed data into a dict)
        self.extra_data = {} 

        self._process_parsed_data(distance_method=distance_method)
        
        return None
