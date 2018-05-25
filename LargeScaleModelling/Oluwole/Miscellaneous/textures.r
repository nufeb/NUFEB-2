# Define some useful textures.

green_glass=Texture(" pigment { color rgbf<0.7, 1.0, 0.7, 0.7> } finish {specular 1 roughness 0.001 ambient 0 diffuse 0 reflection 0.1}  ")
red_glass=Texture(" pigment { color rgbf<1.0, 0.4, 0.4, 0.7> } finish {specular 1 roughness 0.001 ambient 0 diffuse 0 reflection 0.1}  ")
blue_glass=Texture(" pigment { color rgbf<0.4, 0.4, 1.0, 0.7> } finish {specular 1 roughness 0.001 ambient 0 diffuse 0 reflection 0.1}  ")
clear_glass=Texture(" pigment { color rgbf<1.0, 1.0, 1.0, 0.7> }  ")
black_glass=Texture(" pigment { color rgbf<0.2, 0.2, 0.2, 0.4> } finish {specular 1 roughness 0.001 ambient 0 diffuse 0 reflection 0.1}  ")
glass=Texture("finish {specular 1 roughness 0.001 ambient 0 diffuse 0 reflection 0.1}" )

mirror=Texture("pigment { color rgb <0.95, 0.95, 0.95>} finish {
    ambient 0.15
    brilliance 5
    diffuse 0.6
    metallic
    specular 0.80
    roughness 1/100
    reflection 0.75
}")

gold=Texture("pigment { color rgb <1.00, 0.875, 0.575>} finish {
    brilliance 6
    diffuse 0.33
    ambient rgb <0.12,0.1175,0.1115>
    reflection rgb <0.75,0.68,0.53>
    metallic 1
    specular 0.80
    roughness 1/120
}")

i_glass=Interior("ior 1.5")