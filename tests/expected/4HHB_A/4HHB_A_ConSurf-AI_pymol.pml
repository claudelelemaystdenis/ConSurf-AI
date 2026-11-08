# ConSurf-AI colouring. Run from the folder containing 4HHB_A_ConSurf-AI_structure.cif:
#   pymol <this_script>.pml
load 4HHB_A_ConSurf-AI_structure.cif, consurf_ai
hide everything, consurf_ai
show cartoon, consurf_ai
color grey80, consurf_ai
show sticks, consurf_ai and organic
color green, consurf_ai and organic and elem C
color atomic, consurf_ai and organic and not elem C
show spheres, consurf_ai and inorganic
set sphere_scale, 0.5, consurf_ai and inorganic
set_color consurf_1, [0.063, 0.784, 0.82]
set_color consurf_2, [0.549, 1.0, 1.0]
set_color consurf_3, [0.843, 1.0, 1.0]
set_color consurf_4, [0.918, 1.0, 1.0]
set_color consurf_5, [1.0, 1.0, 1.0]
set_color consurf_6, [0.988, 0.929, 0.957]
set_color consurf_7, [0.98, 0.788, 0.871]
set_color consurf_8, [0.941, 0.49, 0.671]
set_color consurf_9, [0.627, 0.145, 0.376]
set_color consurf_lowconf, [1.0, 1.0, 0.588]
color consurf_1, consurf_ai and polymer and chain A and resi 15+19+78+115
color consurf_3, consurf_ai and polymer and chain A and resi 4+8+9+12+13+18+22+23+24+26+34+50+53+57+64+67+68+71+77+82+111+121
color consurf_4, consurf_ai and polymer and chain A and resi 5+11+16+17+20+21+30+51+56+60+96+100+106+109+112+113+114+116+120+123+125+130+131+134
color consurf_5, consurf_ai and polymer and chain A and resi 28+35+36+38+45+49+54+63+72+73+76+79+81+85+89+90+99+102+103+104+105+108+119+129+132+133+138
color consurf_6, consurf_ai and polymer and chain A and resi 75
color consurf_7, consurf_ai and polymer and chain A and resi 3+10+14+27+29+32+44+46+47+48+55+61+69+70+101+107+110+118+124+126+135+137
color consurf_8, consurf_ai and polymer and chain A and resi 6+7+25+33+41+42+52+62+65+66+74+80+86+92+117+122+128+139
color consurf_9, consurf_ai and polymer and chain A and resi 1+2+31+37+39+40+43+58+59+83+84+87+88+91+93+94+95+97+98+127+136+140+141
bg_color white
orient consurf_ai and polymer and chain A
