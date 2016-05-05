function [varargout] = LAMMPSFOAM(varargin)
% Read LAMMPS dump file one timestep at a time
% Input 
%       Dump file name with path 
%       Starting file pointer position in dump file
%       Number of columns in the dump file
% Output is in the form of a structure with following variables
% .timestep     --> Vector containing all time steps
% .Natoms       --> Vector containing number of atoms at each time step
% .x_bound      --> [t,2] array with xlo,xhi at each time step
% .y_bound      --> [t,2] array with ylo,yhi at each time step
% .z_bound      --> [t,2] array with zlo,zhi at each time step
% .atom_data    --> 2 dimensional array with data 
% .position     --> file pointer for reading next time 
%               
% Example
%       data = LAMMPSFOAM('nufeb.bubblemd',-1,7); 
%           Reads the first timestep in the file dump.LAMMPS
%           Specify position = -1 if only the last dump step in the
%               file is needed
%
% See also readdump, scandump
%
%  Author :  Arun K. Subramaniyan
%            sarunkarthi@gmail.com
%            http://web.ics.purdue.edu/~asubrama/pages/Research_Main.htm
%            School of Aeronautics and Astronautics
%            Purdue University, West Lafayette, IN - 47907, USA.

try
    dump = fopen(varargin{1},'r');
catch
    error('Dumpfile not found!');
end
position = varargin{2}; % from beg of file
ncol = varargin{3}; %number of columns

i=1;
  % Put atom type yourself accodring to situations 
Diam = 1e-5;
Dens = 1010;
nm = 4;  % Type of species. 
sub = 0.04;
o2 = 0.005;
no2 = 0.001;
no3 = 0.001;
nh4 = 0.04;

t=0;
done = 0;
last_status = 0;
if position ~= -1
    fseek(dump,position,'bof');
else
    last_status = 1;
end
while done == 0 & last_status == 0
    id = fgetl(dump);
    switch id
        case 'ITEM: TIMESTEP'
            if t == 0
                timestep(i) = str2num(fgetl(dump));
                t=1;
            end
        case 'ITEM: NUMBER OF ATOMS'
            Natoms = str2num(fgetl(dump));
         
        case 'ITEM: BOX BOUNDS ff ff pp'
       % case 'ITEM: BOX BOUNDS ff ff ff'
            pp =12121
            x_bound(1,:) = str2num(fgetl(dump));
            y_bound(1,:) = str2num(fgetl(dump));
            z_bound(1,:) = str2num(fgetl(dump));
        case 'ITEM: ATOMS id type diameter mass x y z '
            atom_data = zeros(Natoms,ncol);%Allocate memory for atom data
            for j = 1 : 1: Natoms
                atom_data(j,:) = str2num(fgetl(dump));
            end
            done = 1;
            p = ftell(dump);
    end
end

% Getting only the last step
if last_status == 1
    % First get the position of the beginning of the last step in the
    % dumpfile
    while ~feof(dump)
        temp = fgetl(dump);
        if length(temp) == 14
            if strcmpi(temp,'ITEM: TIMESTEP')
                p = ftell(dump); % starting position of line next to the header
            end
        end
    end
    fclose(dump);
    dump = fopen(varargin{1},'r');
    fseek(dump,p,'bof');
    % getting Timestep
    timestep = str2num(fgetl(dump));
    
    while ~feof(dump)
        id = fgetl(dump);
        switch id
            case 'ITEM: NUMBER OF ATOMS'
                Natoms = str2num(fgetl(dump));
     
         case 'ITEM: BOX BOUNDS ff ff pp'
       %   case 'ITEM: BOX BOUNDS ff ff ff'   
                pp = 4324243
                x_bound(1,:) = str2num(fgetl(dump));
                y_bound(1,:) = str2num(fgetl(dump));
                z_bound(1,:) = str2num(fgetl(dump));
            case 'ITEM: ATOMS id type diameter mass x y z '
                 C = textscan(dump,'%f %f %f %f %f %f %f',Natoms(i));
     
     % Coordinates and Velocities.
            for r=1:7
                atom_data(:,r) = C{r};
            end  
        end
    end
end

%  Added By Prashant Gupta 
% In order to post process and write a file in the 

% The box data 
  xlo = y_bound(1,1)*Diam;
  xhi = y_bound(1,2)*Diam;
  
  ylo = z_bound(1,1)*Diam;
  yhi = z_bound(1,2)*Diam;
  
  zlo = x_bound(1,1)*Diam;
  zhi = x_bound(1,2)*Diam;
  
  sz =size(atom_data);
  sz(1)
%atom_data_lammps  = atom_data;

atom_data_lammps(:,1) = (1:sz(1))';
atom_data_lammps(:,2) = randi([1 nm],1,sz(1))';
atom_data_lammps(:,3) = atom_data(:,3)*Diam;
%atom_data_lammps(:,4) = (6/pi)*(Dens)*atom_data(:,4);
atom_data_lammps(:,4) = Dens*ones(sz(1),1);
atom_data_lammps(:,ncol-2) = Diam*atom_data(:,ncol-1);
atom_data_lammps(:,ncol-1) = Diam*atom_data(:,ncol); 
atom_data_lammps(:,ncol) = Diam*atom_data(:,ncol-2);
atom_data_lammps(:,ncol+1) = sub*1e-4*randi([0 1e+4],1,sz(1))';
atom_data_lammps(:,ncol+2) = o2*1e-4*randi([0 1e+4],1,sz(1))';
atom_data_lammps(:,ncol+3) = nh4*1e-4*randi([0 1e+4],1,sz(1))';
atom_data_lammps(:,ncol+4) = no2*1e-4*randi([0 1e+4],1,sz(1))';
atom_data_lammps(:,ncol+5) = no3*1e-4*randi([0 1e+4],1,sz(1))';

Numofatomtype = max(atom_data_lammps(:,2))


Num = Natoms;


% atom_data_lammps
% atom_data 
%atom_data_lammps == x,y,z final data in lammps Foam coord axis. 

%----------Outputs-------------

%OUTPUTS IN SAME VARIABLE STRUCTURE
varargout{1}.timestep = timestep;
varargout{1}.Natoms = Natoms;
varargout{1}.x_bound = x_bound;
varargout{1}.y_bound = y_bound;
varargout{1}.z_bound = z_bound;
varargout{1}.atom_data = atom_data;

% save('myfile.txt', 'atom_data_lammps', '-ASCII')
%fprintf(file_1,'Example formatted output \n\n')

fd = fopen('IC350nut2.in', 'wt'); % Open for writing
fprintf(fd,'  Granular Flow Simulation \n\n');
fprintf(fd,'       ');
fprintf(fd,'%d', Num(1,1));
fprintf(fd,'  atoms \n');
fprintf(fd,'          ');
fprintf(fd,'%d', Numofatomtype);
fprintf(fd,'  atom types \n\n');
fprintf(fd,'   ');
fprintf(fd,'%d', xlo);
fprintf(fd,'  ');
fprintf(fd,'%d', xhi);
fprintf(fd,'  xlo ');
fprintf(fd,'xhi \n');

% Y Bound Data 
fprintf(fd,'   ');
fprintf(fd,'%d', ylo);
fprintf(fd,'  ');
fprintf(fd,'%d', yhi);
fprintf(fd,'  ylo ');
fprintf(fd,'yhi \n');
 % _____________________

% Z Bound Data  
fprintf(fd,'   ');
fprintf(fd,'%d', zlo);
fprintf(fd,'  ');
fprintf(fd,'%d', zhi);
fprintf(fd,'  zlo ');
fprintf(fd,'zhi \n\n');
%___________________________

fprintf(fd,' Atoms\n\n');

for i=1:Natoms
 fprintf(fd,'     ');
 fprintf(fd, '%d %d %d  %d %d %d %d %d %d %d %d %d \n', atom_data_lammps(i,:));
end

% for i=10:99
%  fprintf(fd,'    ');
%  fprintf(fd, '%d   %d %d  %d   %d   %d   %d\n', atom_data_lammps(i,:));
% end
% 
% for i=100:999
%  fprintf(fd,'   ');
%  fprintf(fd, '%d   %d %d  %d   %d   %d   %d\n', atom_data_lammps(i,:));
% end
% 
% for i=1000:9999
%  fprintf(fd,'  ');
%  fprintf(fd, '%d   %d %d  %d   %d   %d   %d\n', atom_data_lammps(i,:));
% end
% 
%   for i=10000:Natoms
%    fprintf(fd,' ');
%    fprintf(fd, '%d   %d %d  %d   %d   %d   %d\n', atom_data_lammps(i,:));
%   end

fclose(fd);


%varargout(1).atom_data_lammps = atom_data_lammps;
varargout{1}.position = p; %gives postion of ITEM: TIMESTEP line
%------------------------------

fclose(dump);

end
