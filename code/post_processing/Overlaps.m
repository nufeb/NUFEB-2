clear;
  dump = fopen('snapshot.bubblemd','r');

format long

i=1;
Ncol = 12;
while feof(dump) == 0
    id = fgetl(dump);
    switch id
        case 'ITEM: TIMESTEP'
            timestep(i) = str2num(fgetl(dump));
        case 'ITEM: NUMBER OF ATOMS'
            Natoms(i) = str2num(fgetl(dump));
        case 'ITEM: BOX BOUNDS pp pp pp'
            x_bound(i,:) = str2num(fgetl(dump));
            y_bound(i,:) = str2num(fgetl(dump));
            z_bound(i,:) = str2num(fgetl(dump));
        %case 'ITEM: ATOMS id x y z vx vy vz omegax omegay omegaz c_2_0[4] c_2_1[4] c_2_2[4] c_2_3[4] '
         case 'ITEM: ATOMS id type diameter x y z vx vy vz fx fy fz '
           C = textscan(dump,'%f %f %f %f %f %f %f %f %f %f %f %f',Natoms(i));
           
            
Tp = C{2};
ID = C{1};
D = C{3};
x = C{4};
y = C{5};
z = C{6};
vx = C{7};
vy= C{8};
vz = C{9};
fx= C{10};
fy = C{11};
fz = C{12};



 Avgd(i) = mean(D);
for j=1:length(x)
    sumrad = (D+D(j))*0.5;
 dist = sumrad -((x - x(j)).^2 + (y - y(j)).^2 + (z - z(j)).^2).^(1/2);
 r = dist(dist>0 & dist ~= sumrad);
 r2 = dist(dist<0);
 
 Totad(j) = sum(r)/Avgd(i);
 Totadr2(j) = sum(r2)/Avgd(i);
 Ovlps(j) = length(r);
 %Ar = Ar/Avgd(i);
 %Totad = Totad+ sum(Ar);
 
end 

Nsumovlps(i) = sum(Totad); 
Nsumovlps2(i) = sum(Totadr2); 
Novlps(i) = sum(Ovlps);



i = i+1;
    
    end
    
end 

%  plot(Nsumovlps,'color','red','LineWidth',2,'Marker','o','MarkerSize',8)
%  hold on
%  
%  
%  plot(-Nsumovlps2,'color','blue','LineWidth',2,'Marker','+','MarkerSize',8)
%  hold on
%  
%  
  plot(Novlps./Natoms,'color','blue','LineWidth',2,'Marker','+','MarkerSize',8)
  hold on

Nsumovlps(end)
-Nsumovlps2(end)
Novlps(end)





