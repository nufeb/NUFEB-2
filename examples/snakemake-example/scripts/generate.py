import os
import sys
import shutil
import fileinput
from nufeb_case_utils import sim_volume

def copytree(src, dst, symlinks=False, ignore=None):
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)


targetdir = os.path.join('.', 'results', sys.argv[2], f'rand{sys.argv[1]}')
templatedir = os.path.join('.', 'resources', 'template')
print(f'Random seed is {sys.argv[1]}')
print(f'Working in {os.getcwd()}')
print(f'Target {targetdir}')
print(f'Template {templatedir}')

if not os.path.isdir(templatedir):
    sys.stderr.write(f'Could not find the template directory {templatedir}\n')
    sys.stderr.write(f'\t This script requires a template NUFEB run which\n')
    sys.stderr.write(f'\t 1. Includes an inputscript.nufeb\n')
    sys.stderr.write(f'\t 1. With a het_seed string for the random seed.\n')
    sys.exit(1)

if os.path.isdir(targetdir):
    copytree(templatedir, targetdir)
else:
    shutil.copytree(templatedir, targetdir)

box_vol = sim_volume(float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]),
                     float(sys.argv[6]), float(sys.argv[7]))

fname = os.path.join(targetdir, 'inputscript.nufeb')
with fileinput.FileInput(fname, inplace=True, backup='.bak') as file:
    for line in file:
        het_line = line.replace('het_seed', sys.argv[1])
        print(het_line.replace('template_fullvol', str(box_vol)), end='')
