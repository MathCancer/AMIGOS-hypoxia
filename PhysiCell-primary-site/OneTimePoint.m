close all
clear
clc
cd C:\Users\Furkan\Documents\GitHub\AMIGOS-hypoxia-Forked-\PhysiCell-primary-site
cd output\

load('output00000011_microenvironment0.mat');
XPos = multiscale_microenvironment(1,:);
YPos = multiscale_microenvironment(2,:);
Oxygen = multiscale_microenvironment(5,:);

plot3(XPos,YPos,Oxygen,'.')
cd ..

Oxygen(11326)