%#########Matlab codes for reading LAMMP model
%#cd('E:/Lammp/Long-Run')

clear;
  dump = fopen('snapshot.bubblemd','r');

format long

i=1;
while feof(dump) == 0
    id = fgetl(dump);
    switch id
        case 'ITEM: TIMESTEP'
            timestep(i) = str2num(fgetl(dump));
        case 'ITEM: NUMBER OF ATOMS'
            Natoms(i) = str2num(fgetl(dump))
        case 'ITEM: BOX BOUNDS pp pp pp'
            x_bound(i,:) = str2num(fgetl(dump));
            y_bound(i,:) = str2num(fgetl(dump));
            z_bound(i,:) = str2num(fgetl(dump));
         case 'ITEM: ATOMS id type diameter x y z vx vy vz fx fy fz '
        
           
           C = textscan(dump,'%f %f %f %f %f %f %f %f %f %f %f %f',Natoms(i));
           id = C{1};
           type = C{2};
           Dia = C{3};
           X = C{4};
           Y = C{5};
           Z  = C{6};
           Vx = C{7};
           Vy = C{8};
           Vz = C{9};
           Fx = C{10};
           Fy = C{11};
           Fz = C{12};
           
           data = [id,type,Dia,X,Y,Z,Vx,Vy,Vz,Fx,Fy,Fz];
         
           fname = sprintf('data%d.xlsx', i);
           save(fname)
           csvwrite(fname,data)
           % For each time-step: These arrays can be accessed
           % here. Do all post processing here for a particular
           % time step. 
           
    i=i+1;
end

end 






