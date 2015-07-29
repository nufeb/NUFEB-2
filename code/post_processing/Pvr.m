function [ POV_RAY ] = Pvr(image)

%==================================



% WRITE import.pov FILE FOR POVRAY
fid=fopen('import.pov','w');

fprintf(fid,'// Right-handed coordinate system in which the z-axis points upwards\n');
fprintf(fid,'camera {\n');

fprintf(fid,['location <100e-3,100e-3,100e-3>\n']);
fprintf(fid,'sky y\n');
fprintf(fid,'right -0.24*x*image_width/image_height\n');
fprintf(fid,'up 0.24*y\n');
fprintf(fid,'look_at <0.0,0.0,0.0 >\n');
fprintf(fid,'}\n');
fprintf(fid,'// White background\n');
fprintf(fid,'background{rgb 1}\n');
fprintf(fid,'// Two lights with slightly different colors\n');
fprintf(fid,'light_source{<-8,-20,30> color rgb <0.77,0.75,0.75>}\n');
fprintf(fid,'light_source{<25,-12,12> color rgb <0.38,0.40,0.40>}\n');
fprintf(fid,'// Radius of the Voronoi cell network\n');
fprintf(fid,'#declare r=0.05;\n');
fprintf(fid,'// Radius of the particles\n');
fprintf(fid,'#declare s=0.6;\n');
fprintf(fid,'// Particles\n');
fprintf(fid,'union{\n');
fprintf(fid,'#include "data_p.pov"\n');
fprintf(fid,'}\n');



fclose(fid);

POV_RAY=1


end





