// Right-handed coordinate system in which the z-axis points upwards
camera {
location <2.0e-2,-20.0e-2,2e-2>
sky y
right -0.24*x*image_width/image_height
up 0.24*y
look_at <2.0e-2,2.0e-2,2e-2 >
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
