import numpy as np
from Bio.PDB import PDBParser

def extract_features(pdb_path):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    coords = []
    for atom in structure.get_atoms():
        coords.append(atom.coord)

    coords = np.array(coords)

    # Structural descriptors
    centroid = coords.mean(axis=0)
    radius = np.linalg.norm(coords - centroid, axis=1).mean()
    atom_count = len(coords)

    return np.array([atom_count, radius])
