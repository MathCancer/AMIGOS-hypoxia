close all
clear
clc
cd C:\Users\Furkan\Documents\GitHub\AMIGOS-hypoxia-Forked-\PhysiCell-primary-site
cd output\

%%
load('output00000080_microenvironment0.mat');
XPos = multiscale_microenvironment(1,:);
YPos = multiscale_microenvironment(2,:);
Oxygen = multiscale_microenvironment(5,:);

plot3(XPos,YPos,Oxygen,'.')

Oxygen(11326)
%%

cd C:\Users\Furkan\Documents\GitHub\AMIGOS-hypoxia-Forked-\PhysiCell-primary-site
cd output\
figure()
load('output00000080_vasculature.mat');
XPos = Vascular_Data (1,:);
YPos = Vascular_Data (2,:);
FuncVas = Vascular_Data(5,:);

plot3(XPos,YPos,FuncVas)

cd ..
