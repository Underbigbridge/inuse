import ase
from ase import *
from ase.build import surface, bulk, fcc111, fcc111_root
from ase.build import molecule, add_adsorbate
from ase.io import read, write
from ase.visualize import view
from ase.constraints import FixAtoms

pt111 = fcc111('Pt',
               size=(3, 2, 3),
               a=3.9905801652880402,
               vacuum=12,
               orthogonal=True)
adsorbate1 = molecule('H2O') 
adsorbate1.rotate(38,'x' )
adsorbate1.rotate(-90,'z' )

adsorbate2 = molecule('H2O')
adsorbate2.rotate(-90,'y' )


add_adsorbate(pt111, adsorbate1, 3.5, position=(0,0),offset= (0.0,1.0))
add_adsorbate(pt111, adsorbate1, 3.5, position=(0,0),offset= (2.0,0.0))
add_adsorbate(pt111, adsorbate2, 3.5, position=(0,0),offset= 0)
add_adsorbate(pt111, adsorbate2, 3.5, position=(0,0),offset= 1)
 
mask = [11 < atom.position[2] < 14.5 for atom in pt111]  
fixed_indices = [i for i, fixed in enumerate(mask) if fixed]  
constraint = FixAtoms(indices=fixed_indices)
pt111.set_constraint(constraint)
view(pt111)
print(pt111)
pt111.write('pt111_3_2.vasp')