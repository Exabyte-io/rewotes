def write_atom_species(file,atomic_species):
    with open(file, "a") as file_object:
        file_object.write("ATOMIC_SPECIES \n")
        for atom in atomic_species:
            listed = list(atom.values())
            file_object.write(" ".join(listed)+'\n')


def write_atom_positions(file, positions):
    with open(file, "a") as file_object:
        file_object.write("ATOMIC_POSITIONS (crystal) \n")
        for i in positions:
            try:
                file_object.write(" ".join(i)+'\n')
            except:
                try:
                    file_object.write(" ".join(i.astype(str))+'\n')
                except:
                    file_object.write(" ".join(str(i))+'\n')

def write_cell_parameters(file, cell):
    with open(file, "a") as file_object:
        file_object.write("CELL_PARAMETERS (angstrom) \n")
        for i in cell:
            try:
                file_object.write(" ".join(i)+'\n')
            except:
                try:
                    file_object.write(" ".join(i.astype(str))+'\n')
                except:
                    file_object.write(" ".join(str(i))+'\n')
            # print(i.astype(str))


def write_k_points(file, k):
    with open(file, "a") as file_object:
        file_object.write("K_POINTS automatic \n")
        file_object.write(k+'\n')

def write_k_points_bands(file,k):
    with open(file, "a") as file_object:
        file_object.write("K_POINTS crystal_b \n")
        file_object.write(str(len(k))+'\n')
        for point in k:
            listed = list(point.values())
            file_object.write(" ".join(listed)+'\n')

def write_k_points_matdyn(file,k):
    with open(file, "a") as file_object:
        file_object.write(str(len(k))+'\n')
        for point in k:
            listed = list(point.values())
            file_object.write(" ".join(listed)+'\n')


