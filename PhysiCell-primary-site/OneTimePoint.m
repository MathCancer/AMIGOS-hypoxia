close all
clear
clc
cd C:\Users\Furkan\Documents\GitHub\AMIGOS-hypoxia-Forked-\PhysiCell-primary-site
cd output\

%%
% load('output00000080_microenvironment0.mat');
% XPos = multiscale_microenvironment(1,:);
% YPos = multiscale_microenvironment(2,:);
% Oxygen = multiscale_microenvironment(5,:);

% plot3(XPos,YPos,Oxygen,'.')

%%
load('output00000007_vasculature.mat');

figure()
Y_origin = find(Vascular_Data (2,:) == 10);
XPos1 = Vascular_Data(1,:);
FuncVas = Vascular_Data(5,:);
FuncVas = FuncVas(Y_origin);
XPos1=XPos1(Y_origin);
plot(XPos1,FuncVas);
cd ..
%%
% XPos = Vascular_Data (1,:);
% YPos = Vascular_Data (2,:);
% FuncVas = Vascular_Data(5,:);
% figure()
% plot3(XPos,YPos,FuncVas)
% 


%%
% 
% cd C:\Users\Furkan\Desktop\Vascular_Outputs\output_VEGF_ring_no_vascular_ring
% 
% 
% load('output00000056_vasculature.mat');
% hold on
% Y_origin = find(Vascular_Data (2,:) == 10);
% XPos1 = Vascular_Data(1,:);
% FuncVas1 = Vascular_Data(5,:);
% FuncVas1 = FuncVas1(Y_origin);
% XPos1=XPos1(Y_origin);
% plot(XPos1,FuncVas1,'r');
% title('basal extension rate = a')
% %legend({'basal extension rate = a', 'basal extension rate = a*1000*1000'})