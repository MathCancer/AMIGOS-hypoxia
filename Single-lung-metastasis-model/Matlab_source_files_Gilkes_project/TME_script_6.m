function output = TME_script_6(live_cell_density, necrotic_cell_density, apoptotic_cell_density, vascular_density, oxygen, angiogenic_factor, hypoxic_level, hypoxic_o2_threshold, critical_o2_threshold, dr, dr_fine, dt, domain_size, t);  
  
%This file is to run all the tumor microenvironment updates - O2, hypoxia,
%angiogenic factor, and vasculature functioning.  Adding each piece down
%below as a function call to a separate script.  


persistent n;

parameters.max_time = dt*60; % converts hours to minutes for use in diffusion scripts
max_time = 10;

%  define units (for our sanity)
time_units = 'minutes'; 
space_units = 'micron'; 

% mesh information

parameters.coarse_to_fine_ratio = 3;
parameters.mesh_size_scaling_parameter = 60;
parameters.dr_coarse = dr;
parameters.dr = dr_fine; % determines fine grid size from input coarse grid size
parameters.space_size = domain_size; %32 * 50

if mod(domain_size, 60) ~= 0
    sprintf('error')
end

% time information

parameters.dt = 1E-2;


% Counters

parameters.t = 0;
parameters.n = 0;

% Simulation initializations

parameters.tumor_radius = 10 * parameters.mesh_size_scaling_parameter; %32*8;
parameters.tumor_influence_thickness = 4 * parameters.mesh_size_scaling_parameter; %32*3;
parameters.tumor_viable_rim_size = 100;

% Variable passing

tumor_cell_density_full = live_cell_density + necrotic_cell_density + apoptotic_cell_density;
tumor_cell_density_live = live_cell_density;

% set up parameter values 


% Oxygen Field - from PM
% http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0134999

parameters.oxygen_diffusion_coefficient = 1E5; % um^2/min
parameters.oxygen_secretion_rate = 9.9; %1/min
parameters.tumor_cell_O2_uptake = 10.2;%%10; %1/min might need to increase uptake - rim is just too big ... I wonder how this would look on a cartesian grid, sans the invers square law
parameters.oxygen_decay_rate = 0.01 * parameters.tumor_cell_O2_uptake; % 1/min
parameters.blood_oxygen_level = 1;

% Angiogenic Factor Field
% https://doi.org/10.1051/mmnp/20094404  https://www.mmnp-journal.org/articles/mmnp/abs/2009/04/mmnp20094p118/mmnp20094p118.html

parameters.af_diffusion_coefficient = 166.7;
parameters.af_secretion_rate_max = 1;   %  What should this be?
parameters.af_target_level = 1;
parameters.af_decay_rate = 0.01667;  % 100 - 200 um length scale desired
parameters.angiogenic_factor_length_scale = 400;

%Vascular Growth, Movement, Birth, and Death - 

%parameters.max_vascular_extension_rate = 3.5;%166.7;%16.67; %1 micron/hour
parameters.vascular_death_rate = (1/18/60);%0.2/60;%0.05/60;%0.1;  % NOTE: Necrotic cells don't kill vasculature!
parameters.max_vascular_birth_rate = (1/18/60);%1;
parameters.vascular_proliferation_threshold = 0.001;%0.0005; % non-dimensional concentration of AF
parameters.vascular_proliferation_saturation = 0.5; 
parameters.vascular_chemotaxis_threshold = 0.001;%0.0005;  % non-dimensional concentration of AF
parameters.vascular_chemotaxis_saturation = 0.5; % non-dimensional concentration of AF
parameters.vascular_density_max = 1;
%parameters.effective_vascular_cutoff_threshold = 0;%1E-8; %might try 1E-7 % non-dimensional vascular density

% Tumor Cell "Field"

parameters.tumor_cell_motility = 100; %0.1 micron/hour - from Cartesian implentation
parameters.tumor_cell_density_max = 1;
parameters.hypoxic_o2_threshold = hypoxic_o2_threshold;  % DG marker level I think
parameters.critical_o2_threshold = critical_o2_threshold;

% Simulation Initialization

% figure out mesh size

% Fine

r_nodes = ceil( parameters.space_size / parameters.dr )+1; 
r_length = r_nodes*parameters.dr; 
R = 0 : parameters.dr : r_length - 1;

% Coarse

r_nodes_coarse = ceil( parameters.space_size / parameters.dr_coarse )+1;
r_length_coarse = r_nodes_coarse * parameters.dr_coarse; 
R_coarse = 0 : parameters.dr_coarse : r_length_coarse - 1;

%allocate variables 

%tumor_cell_density = zeros (r_nodes_coarse, 1);
% vascular_density = ones (r_nodes_coarse, 1);
% oxygen = ones( r_nodes, 1 );
% angiogenic_factor = zeros(r_nodes, 1);
% hypoxic_level = zeros(r_nodes,1);


% for i=1:r_nodes_coarse
% 
%     r = parameters.dr_coarse * i; 
%  
%     % create tumor
%     if( r <= parameters.tumor_radius )
%         tumor_cell_density(i) = 1; % initializing to max density
%         %vascular_density(i) = 0.5; 
%     end
% 
%     % create vascular hole at tumor
%     modified_r = max( 0, r - parameters.tumor_radius );         
%     gaussian = exp( -modified_r^2 / ( 2* (parameters.tumor_influence_thickness/3)^2 ) ); 
%     vascular_density(i) = 1-gaussian; 
% 
%     % intial oxygen profile
%     
%     if i<r_nodes_coarse
%         
%         oxygen(3*i - 2) = vascular_density(i);
%         oxygen(3*i - 1) = vascular_density(i);
%         oxygen(3*i) = vascular_density(i);
%        
%     elseif i==r_nodes_coarse
%         
%         oxygen(r_nodes) = vascular_density(i);
%         
%     end
% end
% 
% 
% for i=1:r_nodes
%     
%     r = parameters.dr_coarse * i; 
%     
%     modified_r = max( 0, r - (parameters.tumor_radius-parameters.tumor_viable_rim_size) );         
%     gaussian = exp( -modified_r^2 / ( 2* (parameters.angiogenic_factor_length_scale)^2 ) ); 
%     angiogenic_factor(i) = gaussian;
% 
%      %angiogenic_factor(i) = 0.0001;
%     %angiogenic_factor(i) = 0.000101;
%     
%     %oxygen(i) = 0.5 * parameters.blood_oxygen_level;
% end
   
if (t==0)
    
    n=0;
    
end

% Initial Configuration

if (t<0)

save (sprintf('save%08i.mat', parameters.n)); 
    
    n = 0;

    % Oxygen

    figure(1);
    %hold on;
    plot(R,oxygen,'b');
    %hold on;
    %plot(R_coarse,tumor_cell_density,'k')
    %plot(R_coarse,vascular_density,'r')
    title(sprintf('Sim. Time =%0.8f dt=%0.2g', parameters.t, parameters.dt))  
    ylabel(sprintf('Oxygen, normalized to 1.0 with D=%0.2e and lambda=%0.2e', parameters.oxygen_diffusion_coefficient, parameters.oxygen_decay_rate))
    xlabel(sprintf('X-axis (um) with dr=%0.2g', parameters.dr))
    print('-f1', sprintf('oxygen%08i', parameters.n), '-djpeg');
    %hold off;
    
    % HF
    
    figure(2);
    %hold on;
    plot(R,hypoxic_level);
    title(sprintf('Sim. Time =%0.8f dt=%0.2g', parameters.t, parameters.dt))  
    ylabel(sprintf('Hypoxic Fraction'))
    xlabel(sprintf('X-axis (um) with dr=%0.2g', parameters.dr))
    print('-f2', sprintf('hypoxia%08i', parameters.n), '-djpeg');
    
    % AF
    
    figure(3);
    %hold on;
    plot(R,angiogenic_factor);
    title(sprintf('Sim. Time =%0.8f dt=%0.2g', parameters.t, parameters.dt))  
    ylabel(sprintf('Angiogenic Factor'))
    xlabel(sprintf('X-axis (um) with dr=%0.2g', parameters.dr))
    print('-f3', sprintf('angio_factor%08i', parameters.n), '-djpeg');
    
    % Vascular Density
    
    figure(4);
    %hold on;
    plot(R_coarse,vascular_density);
    title(sprintf('Sim. Time =%0.8f dt=%0.2g', parameters.t, parameters.dt))  
    ylabel(sprintf('Vascular Density'))
    xlabel(sprintf('X-axis (um) with dr=%0.2g', parameters.dr))
    print('-f4', sprintf('vascular_density%08i', parameters.n), '-djpeg');
    
    % Tumor Cell Density
    
    figure(5);
    %hold on;
    plot(R_coarse,tumor_cell_density);
    title(sprintf('Sim. Time =%0.8f dt=%0.2g', parameters.t, parameters.dt))  
    ylabel(sprintf('Tumor Cell Density'))
    xlabel(sprintf('X-axis (um) with dr=%0.2g', parameters.dr))
    print('-f5', sprintf('tumor_cell_density%08i', parameters.n), '-djpeg');

end


% run the codes

while( parameters.t < parameters.max_time + parameters.dt/2 )

    % Oxygen

     output = TME_O2_3(tumor_cell_density_full, vascular_density, oxygen, parameters); 
     
     oxygen = output.oxygen;
    %vascular_density_fine = output.vascular_density_fine;
    %tumor_cell_density_fine = output.tumor_cell_density_fine;
    
    % Hypoxia

     output = TME_hypoxic_level_3(hypoxic_level, oxygen, parameters);
     
     hypoxic_level = output.hypoxic_level;

    % Angiogenic Factor

     output = TME_af_3 (angiogenic_factor, hypoxic_level, tumor_cell_density_live, parameters);
     
     angiogenic_factor = output.angiogenic_factor;

    % Vascular Functionality

     output = TME_vasculature_8 (angiogenic_factor, tumor_cell_density_full, vascular_density, parameters);

     vascular_birth_rate = output.vascular_birth_rate;
     vascular_death_rate = output.vascular_death_rate;
     angiogenic_factor_coarse = output.angiogenic_factor_coarse;
     grad_angiogenic_factor = output.grad_angiogenic_factor;
     velocity = output.velocity;
     vascular_chemotaxis_rate_modifier = output.vascular_chemotaxis_rate_modifier;
    
    % Written output
    
    parameters.t = parameters.t + parameters.dt;
    parameters.n = parameters.n + 1;
    n = n+1;
    
%     if (mod(parameters.n,1/(1*parameters.dt))==0)
%         
%         save (sprintf('save%08i.mat', parameters.n)); 
% 
%         figure(1);
%         plot(R,oxygen);
%         %hold on;
%         %plot(R_coarse,tumor_cell_density,'k')
%         %plot(R_coarse,vascular_density,'r')
%         title(sprintf('Sim. Time =%0.8f dt=%0.2g', t, parameters.dt))  
%         ylabel(sprintf('Oxygen, normalized to 1.0 with D=%0.2e and lambda=%0.2e', parameters.oxygen_diffusion_coefficient, parameters.oxygen_decay_rate))
%         xlabel(sprintf('X-axis (um) with dr=%0.2g', parameters.dr))
%         print('-f1', sprintf('oxygen%08i', parameters.n), '-djpeg');
%         %hold off;
%         
%         figure(2);
%         plot(R,hypoxic_level);
%         title(sprintf('Sim. Time =%0.8f dt=%0.2g', t, parameters.dt))  
%         ylabel(sprintf('Hypoxic Fraction'))
%         xlabel(sprintf('X-axis (um) with dr=%0.2g', parameters.dr))
%         print('-f1', sprintf('hypoxia%08i', parameters.n), '-djpeg');
%         
%         figure(3);
%         plot(R,angiogenic_factor);
%         title(sprintf('Sim. Time =%0.8f dt=%0.2g', t, parameters.dt))  
%         ylabel(sprintf('Angiogenic Factor'))
%         xlabel(sprintf('X-axis (um) with dr=%0.2g', parameters.dr))
%         print('-f1', sprintf('angio_factor%08i', parameters.n), '-djpeg');
%         
%         figure(4);
%         plot(R_coarse,vascular_density);
%         
%         title(sprintf('Sim. Time =%0.8f dt=%0.2g', t, parameters.dt))  
%         ylabel(sprintf('Vascular Density'))
%         xlabel(sprintf('X-axis (um) with dr=%0.2g', parameters.dr))
%         print('-f1', sprintf('vascular_density%08i', parameters.n), '-djpeg');
%         
%         figure(5);
%         plot(R_coarse,tumor_cell_density);
%         title(sprintf('Sim. Time =%0.8f dt=%0.2g', t, parameters.dt))  
%         ylabel(sprintf('Tumor Cell Density'))
%         xlabel(sprintf('X-axis (um) with dr=%0.2g', parameters.dr))
%         print('-f1', sprintf('tumor_cell_density%08i', parameters.n), '-djpeg');
%         
%         
%         
%     end
    

    

end

n = n+1;

% Writen output - once per call to the TME script

%save (sprintf('save%08i.mat', n)); 

% Graphic output - once per call to the TME script

% figure(1);
% plot(R,oxygen);
% %hold on;
% %plot(R_coarse,tumor_cell_density,'k')
% %plot(R_coarse,vascular_density,'r')
% title(sprintf('Sim. Time =%0.8f dt=%0.2g', t, parameters.dt))  
% ylabel(sprintf('Oxygen, normalized to 1.0 with D=%0.2e and lambda=%0.2e', parameters.oxygen_diffusion_coefficient, parameters.oxygen_decay_rate))
% xlabel(sprintf('X-axis (um) with dr=%0.2g', parameters.dr))
% print('-f1', sprintf('oxygen%08i', n), '-djpeg');
% %hold off;
% 
% figure(2);
% plot(R,hypoxic_level);
% title(sprintf('Sim. Time =%0.8f dt=%0.2g', t, parameters.dt))  
% ylabel(sprintf('Hypoxic Fraction'))
% xlabel(sprintf('X-axis (um) with dr=%0.2g', parameters.dr))
% print('-f2', sprintf('hypoxia%08i', n), '-djpeg');
% 
% figure(3);
% plot(R,angiogenic_factor);
% title(sprintf('Sim. Time =%0.8f dt=%0.2g', t, parameters.dt))  
% ylabel(sprintf('Angiogenic Factor'))
% xlabel(sprintf('X-axis (um) with dr=%0.2g', parameters.dr))
% print('-f3', sprintf('angio_factor%08i', n), '-djpeg');
% 
% figure(4);
% plot(R_coarse,vascular_density);
% title(sprintf('Sim. Time =%0.8f dt=%0.2g', t, parameters.dt))  
% ylabel(sprintf('Vascular Density'))
% xlabel(sprintf('X-axis (um) with dr=%0.2g', parameters.dr))
% print('-f4', sprintf('vascular_density%08i', n), '-djpeg');
% 
% figure(5);
% plot(R_coarse,tumor_cell_density);
% title(sprintf('Sim. Time =%0.8f dt=%0.2g', t, parameters.dt))  
% ylabel(sprintf('Tumor Cell Density'))
% xlabel(sprintf('X-axis (um) with dr=%0.2g', parameters.dr))
% print('-f5', sprintf('tumor_cell_density%08i', n), '-djpeg');

output.oxygen = oxygen;
output.vascular_density = vascular_density;
output.hypoxic_level = hypoxic_level;
output.angiogenic_factor = angiogenic_factor;
output.angiogenic_factor_coarse = angiogenic_factor_coarse;
output.velocity = velocity;
output.vascular_chemotaxis_rate_modifier = vascular_chemotaxis_rate_modifier;
output.grad_angiogenic_factor = grad_angiogenic_factor;
output.parameters = parameters;

return;
toc;