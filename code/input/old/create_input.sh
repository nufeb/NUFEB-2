#!/bin/bash
# Matlab povray script

matlab -nodisplay -nosplash -nodesktop -r "data = LAMMPSFOAM('nufeb.bubblemd',-1,7);exit;"