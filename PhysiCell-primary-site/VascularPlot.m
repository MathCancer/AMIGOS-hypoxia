close all
clear
clc

% 
% cd C:\Users\Furkan\Documents\GitHub\AMIGOS-hypoxia-Forked-\PhysiCell-primary-site
% cd output

cd C:\Users\Furkan\Desktop\Vascular_Outputs\Q)output_61_12_no_freeze_pres_1000_small_expa_1e+12_deg

%%
% Obtaining names of mat files
s = what;
MatFiles = s.mat;
VasMatFiles = MatFiles(contains(MatFiles,'vasculature'));
MEMatFiles = MatFiles(contains(MatFiles,'micro'));
MEMatFiles(1) = [];


%%
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
filename = 'testAnimated.gif';

%%
for i = 1:length(VasMatFiles)
    load(VasMatFiles{i});
    FunVas = Vascular_Data(5,:);
    XPos = Vascular_Data(1,:);
    YPos = Vascular_Data(2,:);

    FunVas = reshape(FunVas, [150 150]);    
    XPos = reshape(XPos, [150 150]);
    YPos = reshape(YPos, [150 150]);
    
    ax1 = subplot(2,2,1);
    contourf(XPos,YPos,FunVas,'linecolor','none');
    caxis(ax1,[VasMinValue VasMaxValue]);
    colorbar('eastoutside');
    title('Vasculature') ;
    axis image;
    
    load(MEMatFiles{i});
    O2 = multiscale_microenvironment(5,:);
    VEGF = multiscale_microenvironment(7,:);
    O2 = reshape(O2, [150 150]);
    VEGF = reshape(VEGF, [150 150]);
    HP = multiscale_microenvironment(8,:);
    HP = reshape(HP, [150 150]);
    
    
    ax4 = subplot(2,2,2);
    contourf(XPos,YPos,HP,'linecolor','none');
    title('Pressure');
    colorbar('eastoutside');
    axis image;
    
    
    
    ax2 = subplot(2,2,3);
    contourf(XPos,YPos,O2,'linecolor','none');
    caxis(ax2,[0 50]);
    title('Oxygen');
    colorbar('eastoutside');
    axis image;
    
    ax3 = subplot(2,2,4);
    contourf(XPos,YPos,VEGF,'linecolor','none');
    caxis(ax3,[-0.001 1]);
    colorbar('eastoutside');
    set(gca,'color',[0.2422,0.1504,0.6603]);
    title('VEGF');
    axis image;
    
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


cd C:\Users\Furkan\Documents\GitHub\AMIGOS-hypoxia-Forked-\PhysiCell-primary-site
