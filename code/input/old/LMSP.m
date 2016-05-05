clear
clc



dump = fopen('snapnorm.bubblemd','r');


i=1;

Ncol = 9;
diam = 1e-4;
dens = 1100;

fd = fopen('snapnorm.bubblemd', 'wt'); % Open for writing

while feof(dump) == 0
    id = fgetl(dump);
    switch id
        case 'ITEM: TIMESTEP'
            timestep(i) = str2num(fgetl(dump));
        case 'ITEM: NUMBER OF ATOMS'
            Natoms(i) = str2num(fgetl(dump));
         case 'ITEM: BOX BOUNDS pp pp pp'
            x_bound(1,:) = str2num(fgetl(dump));
            y_bound(1,:) = str2num(fgetl(dump));
            z_bound(1,:) = str2num(fgetl(dump));  
        case 'ITEM: ATOMS id type diameter x y z vx vy vz '

           C = textscan(dump,'%f %f %f %f %f %f %f %f %f',Natoms(i));
     
     % Coordinates and Velocities.
            for r=1:Ncol
                atom_data(:,r) = C{r};
            end    
            
            atom_data(:,3) = atom_data(:,3)/diam ;
            atom_data(:,4) = atom_data(:,4)/diam ;
            atom_data(:,5) = atom_data(:,5)/diam ;
            atom_data(:,6) = atom_data(:,6)/diam ;
            atom_data(:,7) = atom_data(:,7)/diam ;
            atom_data(:,8) = atom_data(:,8)/diam ;
            atom_data(:,9) = atom_data(:,9)/diam ;
            
            x_bound(1,:) = x_bound(1,:)/diam;
            y_bound(1,:) = y_bound(1,:)/diam;
            z_bound(1,:) = z_bound(1,:)/diam;
            
            fprintf(fd,'ITEM: TIMESTEP \n');
            fprintf(fd,'%d \n', timestep(i));
            fprintf(fd,'ITEM: NUMBER OF ATOMS \n');
            fprintf(fd,'%d \n', Natoms(i));
            fprintf(fd,'ITEM: BOX BOUNDS pp pp pp \n');
            fprintf(fd,'%d %d\n', x_bound(1,:));
            fprintf(fd,'%d %d\n', y_bound(1,:));
            fprintf(fd,'%d %d\n', z_bound(1,:));
            fprintf(fd,'ITEM: ATOMS id type diameter x y z vx vy vz \n'); 
            %fprintf(fd,'ITEM: ATOMS id type diameter x y z vx vy vz \n');
           fprintf(fd, '%i %i %f %f %f %f %f %f %f \n', atom_data'); 
            
            
            
           i = i+1
    end
end 

fclose(dump)
fclose(fd);



