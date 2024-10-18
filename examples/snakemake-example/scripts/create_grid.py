import numpy as np
import sys
import os

targetfile = os.path.join('.', 'results', sys.argv[2], f'rand{sys.argv[1]}', 'atom.in')

GRIDN = int(sys.argv[3])
BUGSIZE = 1e-6
BUGSPACING = float(sys.argv[4])
BUGZ = 1e-6

grid_dist = BUGSPACING*BUGSIZE
margin = grid_dist/2
bug_vert = np.linspace(grid_dist, grid_dist*GRIDN, num=GRIDN)-margin

n_atoms = GRIDN**2
xy_max = bug_vert[-1]+margin
z_max = 2.0e-4


with open(targetfile, 'w') as ofile:
    ofile.write('NUFEB SIMULATION\n')
    ofile.write('\n')
    ofile.write(f'\t\t{n_atoms} atoms\n')
    ofile.write('\t\t25 atom types\n')
    ofile.write('\n')
    ofile.write(f'\t0 {xy_max:.1e} xlo xhi\n')
    ofile.write(f'\t0 {xy_max:.1e} ylo yhi\n')
    ofile.write(f'\t0 {z_max:.1e} zlo zhi\n')
    ofile.write('\n')
    ofile.write('Atoms\n')
    ofile.write('\n')
    i = 0
    for y in bug_vert:
        for x in bug_vert:
            i = i + 1
            bug_str = f'\t{i} {i} {BUGSIZE:.2e} 150 {x:.2e} {y:.2e} {BUGZ:.2e} {BUGSIZE:.2e}\n'
            ofile.write(bug_str)
