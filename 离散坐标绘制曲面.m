clc;clear;close all;
x = [1 2 3];
y = [0 0 0];
z = [2 4 9];
figure(1);
plot3(x,y,z);
X = [x x x];X = sort(X);
Y = [-1 -1 -1 0 0 0 1 1 1];
Z = [z z z];Z = sort(Z);
[X,Y] = meshgrid(X,Y);
Z = repmat(Z,9,1);
figure(2);
mesh(X,Y,Z)