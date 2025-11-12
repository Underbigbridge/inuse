import ase
from ase import *
from ase.build import surface, bulk, fcc111, fcc111_root
from ase.build import molecule, add_adsorbate, sort
from ase.io import read, write
from ase.visualize import view
from ase.constraints import FixAtoms
import numpy as np

pt111=fcc111_root('Pt', 3, size=(2, 2, 3), vacuum=10)#, orthogonal=True)
mask = [9 < atom.position[2] < 13.0 for atom in pt111]  #  
fixed_indices = [i for i, fixed in enumerate(mask) if fixed]  #  
constraint = FixAtoms(indices=fixed_indices)
pt111.set_constraint(constraint)
def add_complete_h_layer(slab, height):
    """
     On bottom site of Pt, add one layer of H
    """
    # find bottom atoms
    bottom_layer_z = np.min(slab.positions[:, 2])
    bottom_atoms = [atom for atom in slab if abs(atom.position[2] - bottom_layer_z) < 0.1]
    
    print(f"on {len(bottom_atoms)} Pt atoms add H")
    
    # addH
    h_positions = []
    for atom in bottom_atoms:
        h_position = atom.position.copy()
        h_position[2] += height   
        h_positions.append(h_position)
    
    #from ase import Atoms
    h_atoms = Atoms('H' * len(h_positions), positions=h_positions)
    slab.extend(h_atoms)
    
    return slab
pt111_with_h = add_complete_h_layer(pt111, height=5.6)

#view(pt111_with_h)
adsorbate1 = molecule('H2O') 
adsorbate1.rotate(38,'x' )
adsorbate1.rotate(0,'z' )
adsorbate2 = molecule('H2O')
adsorbate2.rotate(90,'y' )
adsorbate2.rotate(-90,'z' )

#adsorbate2 = molecule('H2O') 
#adsorbate2.rotate(0,'x')
add_adsorbate(pt111_with_h, adsorbate1, 2.8, position=(9.604,2.783),offset= 0)
add_adsorbate(pt111_with_h, adsorbate2, 2.8, position=(7.204,1.397),offset= 0)
add_adsorbate(pt111_with_h, adsorbate1, 2.8, position=(4.803,2.783,),offset= 0)
add_adsorbate(pt111_with_h, adsorbate2, 2.8, position=(4.823,5.555),offset= 0)
add_adsorbate(pt111_with_h, adsorbate1, 2.8, position=(7.204,6.941),offset= 0)
add_adsorbate(pt111_with_h, adsorbate2, 2.8, position=(9.624,5.555),offset= 0)
add_adsorbate(pt111_with_h, adsorbate1, 2.8, position=(12.02,6.941),offset= 0)
add_adsorbate(pt111_with_h, adsorbate2, 2.8, position=(2.400,1.397),offset= 0)
pt111_with_h = sort(pt111_with_h)
view(pt111_with_h)
print(pt111_with_h)

pt111_with_h.write('PtrootH2O_H.vasp')