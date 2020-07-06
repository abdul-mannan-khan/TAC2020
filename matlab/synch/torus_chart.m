function [x,y,z] = torus_chart(tt1,tt2)
x = cos(tt1)+cos(tt2+tt1)*cos(tt1);
y = sin(tt1)+cos(tt2+tt1)*sin(tt1);
z = sin(tt2+tt1);