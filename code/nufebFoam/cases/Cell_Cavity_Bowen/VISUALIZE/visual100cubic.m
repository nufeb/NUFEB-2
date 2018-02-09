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
            Natoms(i) = str2num(fgetl(dump));
        case 'ITEM: BOX BOUNDS pp pp pp'
            x_bound(i,:) = str2num(fgetl(dump));
            y_bound(i,:) = str2num(fgetl(dump));
            z_bound(i,:) = str2num(fgetl(dump));
        %case 'ITEM: ATOMS id x y z vx vy vz omegax omegay omegaz c_2_0[4] c_2_1[4] c_2_2[4] c_2_3[4] '
         case 'ITEM: ATOMS id type diameter x y z '
        
           
           C = textscan(dump,'%f %f %f %f %f %f',Natoms(i));
           

Tp = C{2};
ID = C{1};
D = C{3}*1e+4;
X = C{4}*1e+4-0.5;
Y = C{5}*1e+4-0.5;
Z = C{6}*1e+4-0.05;

ii=i;

fid=fopen('data_p.pov','w');



%if (mod(ii,1)==0 & ii>180)



    for l=1:length(C{2})
       
        
        if (Tp(l)==1)
     %    RED = 0.057;
      %  GREEN = 0.012;
     %   BLUE = 0.071;
        
       %  X(l)=0.2;
        RED = 0.8;
        GREEN = 0.2;
        BLUE = 0.2;
        
       
        end 
        
         if (Tp(l)==2)
        RED = 0.2;
        GREEN = 0.8;
        BLUE = 0.2;
        
        
        
        
        end 
        
         if (Tp(l)==3)
        RED = 0.0;
        GREEN = 0.0;
        BLUE = 0.0;
        end 
        
         if (Tp(l)==4)
        RED = 0.2;
        GREEN = 0.2;
        BLUE = 0.9;
        
        
	    if (X(l)>0.025) 
        
	   %   X(l)=0.2;
	      
	    end
	    
        
        
        end 
        
         if (Tp(l)==5)
        RED = 0.96;
        GREEN = 0.51;
        BLUE = 0.19;
        end 
        
         if (Tp(l)==6)
        RED = 0.9;
        GREEN = 0.2;
        BLUE = 0.2;
        end 
        
         if (Tp(l)==7)
        RED = 0.27;
        GREEN = 0.94;
        BLUE = 0.94;
        end 
        
         if (Tp(l)==8)
        RED = 0.94;
        GREEN = 0.20;
        BLUE = 0.90;
        end 
        
         if (Tp(l)==9)
        RED = 0.82;
        GREEN = 0.96;
        BLUE = 0.24;
        end 
        
         if (Tp(l)==10)
        RED = 0.98;
        GREEN = 0.75;
        BLUE = 0.75;
        end 
        
        if (Tp(l)==11)
         RED = 0.57;
        GREEN = 0.12;
        BLUE = 0.71;
        end 

      
fprintf(fid,['sphere{<' num2str(X(l)) ',' num2str(Y(l)) ',' num2str(Z(l)) '>,' num2str(1*D(l)/2) ' pigment {rgb <' num2str(RED) ',' num2str(GREEN) ',' num2str(BLUE) '>} finish{reflection 0.2 specular 0.3 ambient 0.42}}\n']);
    
     
    
    end
    fclose(fid);
    j=i+1000000;   
POV_RAY=PvrRotate100cubic(image, i)

generate=['povray +H3000 +W3000 +FJ Display=-D +O0_images/image' num2str(j) '.jpg import.pov'];
system(generate);

%end
    
    i=i+1;
end




end 


Inputdir = './0_images/' ;  % Input dit of images
fname=dir([Inputdir '*.jpg']); % Specify file extension

wobj = VideoWriter([Inputdir 'video1.avi']); % specify video file name
set(wobj,'FrameRate',5); % Specify frame rate
open(wobj)


for i=1:1:length(fname)
    image = imread([Inputdir fname(i).name]);
    A = image;
    B=A(:,:,1:3);
    writeVideo(wobj, B);
    disp(i)
end
close(wobj);




