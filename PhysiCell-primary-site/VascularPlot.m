close all
clear
clc
cd output\
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




%%
for i = 1:length(VasMatFiles)
    load(VasMatFiles{i})
    FunVas = Vascular_Data(5,:);
    XPos = Vascular_Data(1,:);
    YPos = Vascular_Data(2,:);

    FunVas = reshape(FunVas, [150 150]);    
    XPos = reshape(XPos, [150 150]);
    YPos = reshape(YPos, [150 150]);
    
    ax1 = subplot(2,2,1.5);
    contourf(XPos,YPos,FunVas,'linecolor','none')
    caxis(ax1,[VasMinValue VasMaxValue])
    colorbar('eastoutside')
    title('Vasculature') 

    
    load(MEMatFiles{i});
    O2 = multiscale_microenvironment(5,:);
    VEGF = multiscale_microenvironment(7,:);
    O2 = reshape(O2, [150 150]);
    VEGF = reshape(VEGF, [150 150]);
    
    ax2 = subplot(2,2,3);
    contourf(XPos,YPos,O2,'linecolor','none')
    caxis(ax2,[0 100])
    title('Oxygen')
    colorbar('eastoutside')

    
    ax3 = subplot(2,2,4);
    contourf(XPos,YPos,VEGF,'linecolor','none')
    caxis(ax3,[-0.001 1])
    colorbar('eastoutside')
    set(gca,'color',[0.2422,0.1504,0.6603])
    title('VEGF')
    
    BigTitle = num2str(i*60/60/24);
    suptitle(['Time = ',BigTitle,' days']);
    pause(0.005)
end



cd ..