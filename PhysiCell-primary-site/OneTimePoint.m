close all
clear
clc
% cd C:\Users\Furkan\Documents\GitHub\AMIGOS-hypoxia-Forked-\PhysiCell-primary-site
% cd output

cd C:\Users\Furkan\Desktop\Vascular_Outputs\Q)output_61_12_no_freeze_pres_1000_small_expa_1e+12_deg

%%
% load('output00000080_microenvironment0.mat');
% XPos = multiscale_microenvironment(1,:);
% YPos = multiscale_microenvironment(2,:);
% Oxygen = multiscale_microenvironment(5,:);

% plot3(XPos,YPos,Oxygen,'.')

%%
s = what;
MatFiles = s.mat;
VasMatFiles = MatFiles(contains(MatFiles,'vasculature'));
load(VasMatFiles{1})
VasMaxValue=max(Vascular_Data(5,:));
VasMinValue=min(Vascular_Data(5,:));
for i = 1:length(VasMatFiles)
    load(VasMatFiles{i})
    if VasMaxValue < max(Vascular_Data(5,:))
        VasMaxValue = max(Vascular_Data(5,:));
    end
    if VasMinValue > min(Vascular_Data(5,:))
        VasMinValue = min(Vascular_Data(5,:));
    end    
end

h = figure('Renderer', 'painters', 'Position', [10 10 900 600]);
filename = 'crossSectionVasculature.gif';
%%
for i = 1:length(VasMatFiles)
    
    load(VasMatFiles{i});
    figure(1)
    Y_origin = find(Vascular_Data (2,:) == 10);
    XPos1 = Vascular_Data(1,:);
    FuncVas = Vascular_Data(5,:);
    FuncVas = FuncVas(Y_origin);
    XPos1=XPos1(Y_origin);
    plot(XPos1,FuncVas);
    
    BigTitle = num2str(i/24*3);
    suptitle(['Time = ',BigTitle,' days']);
    
        frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if i == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end


%%
% XPos = Vascular_Data (1,:);
% YPos = Vascular_Data (2,:);
% FuncVas = Vascular_Data(5,:);
% figure()
% plot3(XPos,YPos,FuncVas)
% 

cd C:\Users\Furkan\Documents\GitHub\AMIGOS-hypoxia-Forked-\PhysiCell-primary-site
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