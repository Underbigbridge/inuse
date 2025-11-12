import ase
from ase import *
from ase.build import surface, bulk, fcc111, fcc111_root
from ase.build import molecule, add_adsorbate, sort
from ase.io import read, write
from ase.visualize import view
from ase.constraints import FixAtoms
import numpy as np

#pt111=fcc111_root('Pt', 3, size=(2, 2, 3), vacuum=10)#, orthogonal=True)
pt111=fcc111('Pt', size=(3, 3, 4), vacuum=10)#, orthogonal=True)
mask = [9 < atom.position[2] < 13.0 for atom in pt111]  #  
fixed_indices = [i for i, fixed in enumerate(mask) if fixed]  #  
constraint = FixAtoms(indices=fixed_indices)
pt111.set_constraint(constraint)
def add_complete_h_layer(slab, height):

    # find bottom atoms
    fcc_atoms = [atom for atom in slab if abs(atom.position[2] - 12.26) < 0.4]
    
    print(f"on {len(fcc_atoms)} Pt atoms add H")
    
    # addH
    h_positions = []
    for atom in fcc_atoms:
        h_position = atom.position.copy()
        h_position[2] += height   
        h_positions.append(h_position)
    
    #from ase import Atoms
    h_atoms = Atoms('H' * len(h_positions), positions=h_positions)
    slab.extend(h_atoms)
    
    return slab
pt111_with_h = add_complete_h_layer(pt111, height=5.7)

#view(pt111_with_h)
adsorbate1 = molecule('H2O') 
adsorbate1.rotate(38,'x' )
adsorbate1.rotate(150,'z' )
add_adsorbate(pt111_with_h, adsorbate1, 2.8, position=(2.772,0.000),offset= 0)

adsorbate2 = molecule('H2O')
adsorbate2.rotate(-90,'y' )
add_adsorbate(pt111_with_h, adsorbate2, 2.8, position=(1.386,2.400),offset= 0)

adsorbate3 = adsorbate1
adsorbate3.rotate(-120,'z' )
add_adsorbate(pt111_with_h, adsorbate3, 2.8, position=(2.772,4.801),offset= 0)

adsorbate4 = adsorbate2
adsorbate4.rotate(-120,'z' )
add_adsorbate(pt111_with_h, adsorbate4, 2.8, position=(5.544,4.801),offset= 0)

adsorbate5 = adsorbate3
adsorbate5.rotate(-120,'z' )
add_adsorbate(pt111_with_h, adsorbate3, 2.8, position=(6.930,2.400),offset= 0)

adsorbate6 = adsorbate4
adsorbate6.rotate(-120,'z' )
add_adsorbate(pt111_with_h, adsorbate4, 2.8, position=(5.544,0.000),offset= 0)

pt111_with_h = sort(pt111_with_h)
view(pt111_with_h)
print(pt111_with_h)

pt111_with_h.write('PtH2O_H.vasp')