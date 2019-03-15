load('output00000008_vasculature.mat')
FunVas = Vascular_Data(5,:);
XPos = Vascular_Data(1,:);
YPos = Vascular_Data(2,:);
close all
plot3(XPos,YPos,FunVas,'.')