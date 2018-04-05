FROM ubuntu:16.04

#add new user
RUN useradd --create-home --shell /bin/bash admin

#install required libs
RUN apt-get update
RUN apt-get -y install build-essential sudo make git openmpi-bin openmpi-common libopenmpi-dev 
#libvtk6-dev

#switch to user admin, download NUFEB
USER admin
WORKDIR /home/admin
RUN git clone https://nufebdocker:nufeb123@github.com/shelllbw/nufeb.git

#compile NUFEB code
WORKDIR /home/admin/nufeb/code/lammps5Nov16/src
RUN make yes-user-nufeb
RUN make yes-granular
RUN make -j4 mpi

#add NUFEB executable to bashrc
RUN echo "export PATH=$PATH:/home/admin/nufeb/code/lammps5Nov16/src/" >> /home/admin/.bashrc

