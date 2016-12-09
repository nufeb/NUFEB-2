./Allset.sh
decomposePar > log.decomposePar
mpirun -np 2 lammpsFoam -parallel > log.parallel
