load('output00000005_vasculature.mat')
FunVas = Vascular_Data(5,:);
XPos = Vascular_Data(1,:);
YPos = Vascular_Data(2,:);

plot3(XPos,YPos,FunVas)