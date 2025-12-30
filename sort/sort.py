#!/usr/bin/env python3

import ase
from ase.build import *
from ase.io import read, write
import sys
import numpy as np
from collections import OrderedDict

input_file = sys.argv[1]

outfile =  f"{input_file.split('.')[0]}_sortout.vasp"

inatoms = read(input_file)
atoms = sort(inatoms)

# True: 用分数坐标(推荐周期体系/slab)；False: 用笛卡尔坐标
use_fractional_z = True

# 取 z 坐标
if use_fractional_z:
    z = atoms.get_scaled_positions()[:, 2]
else:
    z = atoms.get_positions()[:, 2]

symbols = atoms.get_chemical_symbols()


species_order = []
seen = set()
for s in symbols:
    if s not in seen:
        species_order.append(s)
        seen.add(s)

new_indices = []
def sort_key(i):
    return z[i]
for sp in species_order:
    idx = []
    for i, sym in enumerate(symbols):
        if sym == sp:
            idx.append(i)
    # 该元素内部按 z 从小到大排序
    idx_sorted = sorted(idx, key=sort_key)  # or reverse=True
    new_indices.extend(idx_sorted)

atoms_sorted = atoms[new_indices]

# 写回 VASP 文件（ASE 会按 atoms_sorted 的顺序输出原子）
write(outfile, atoms_sorted) #format="vasp", direct=use_fractional_z, vasp5=True)

print(f"Done. Wrote: {outfile}")
print("Element order preserved:", species_order)
print("Sorted within each element by z", "(fractional)" if use_fractional_z else "(cartesian)")

