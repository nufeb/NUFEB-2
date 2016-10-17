#include "shapes.inc"
// Right-handed coordinate system in which the z-axis points upwards
camera {
location <0, 60e-2, 0>
sky z
right 0.24*x*image_width/image_height
up 0.24*y
look_at <0, 0, 0>
rotate <0, 0, 3*10>
}
// Create simualtion box
#declare b_x=0.05;
#declare b_y=0.05;
#declare b_z=0.05;
object{ 
Wire_Box(<b_x, b_y, b_z>, <-b_x, -b_y, -b_z>, 0.0002, 0) 
 texture{ pigment{ color rgb<1,1,1>}}
}
// White background
background{rgb 1}
// Two lights with slightly different colors
light_source{<-8,-20,30> color rgb <0.77,0.75,0.75>}
light_source{<25,-12,12> color rgb <0.38,0.40,0.40>}
// Radius of the Voronoi cell network
#declare r=0.05;
// Radius of the particles
#declare s=0.6;
// Particles
union{
#include "data_p.pov"
}
