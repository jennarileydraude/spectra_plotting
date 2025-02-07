import numpy as np
import matplotlib.pyplot as plt
import os

def get_dihedral_angle(atom_1_index, atom_2_index, atom_3_index, atom_4_index, atom_array):
    vector_1 = atom_array[atom_2_index] - atom_array[atom_1_index]
    vector_2 = atom_array[atom_3_index] - atom_array[atom_2_index]
    vector_3 = atom_array[atom_4_index] - atom_array[atom_3_index]
    
    cross_12 = np.cross(vector_1, vector_2)
    cross_23 = np.cross(vector_2, vector_3)
    
    term_1 = np.dot(vector_2, np.cross(cross_12, cross_23))
    term_2 = np.linalg.norm(vector_2) * np.dot(cross_12, cross_23)
    
    dihedral_angle = np.arctan2(term_1, term_2)
    return np.rad2deg(dihedral_angle)

def read_xyz_file(file):
    with open(file, 'r') as xyz_file:
        next(xyz_file)
        next(xyz_file)
        atom_data_list = []
        for line in xyz_file:
            if line.strip():  # Skip empty lines
                atom_data = line.split()
                atom_coordinates = [float(coord) for coord in atom_data[1:]]
                atom_data_list.append(atom_coordinates)
    return np.array(atom_data_list)

def process_xyz_files(directory, atom_indices):
    dihedral_angles = []
    file_list = sorted([os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.xyz')])
    
    for file in file_list:
        atom_array = read_xyz_file(file)
        angle = get_dihedral_angle(*atom_indices, atom_array)
        dihedral_angles.append(angle)
    return dihedral_angles

def plot_dihedral_angles(dihedral_angles1, dihedral_angles2):
    time = np.arange(1, len(dihedral_angles1) + 1)  # 1 Frame = 1 ns ?
    plt.figure(figsize=(10, 6))
    plt.x_grid, y_grid = np.meshgrid(dihedral_angles1, dihedral_angles2)
    dia_flat1 = dihedral_angles1
    dia_flat2 = dihedral_angles2
    plt.scatter(dia_flat1, dia_flat1, c=time, cmap='winter', alpha=0.7)
    plt.colorbar(label="Time (ns)")
    plt.xlabel('Dihedral Angle 1 (Degrees)')
    plt.ylabel('Dihedral Angle 2 (Degrees)')
    plt.minorticks_on()
    plt.title('Dihedral Angles for Length 6 Oligomer')
    plt.savefig('dia2.png')

# Call functions
directory_path = '/home/jennadraude/plot_dihedrals/2dia/buet_6oligomer_500frames'  
atom1_indices = (494, 493, 492, 485)  
atom2_indicies = (480, 479, 478, 401)
dihedral_angles1 = process_xyz_files(directory_path, atom1_indices)
dihedral_angles2 = process_xyz_files(directory_path, atom2_indicies)
plot_dihedral_angles(dihedral_angles1, dihedral_angles2)
