clear 
% figure(5)
% clf
% figure(7)
% clf 
% figure(1); clf
% figure(2); clf
% figure(3); clf
% figure(4); clf
% figure(5); clf
% close all 
tic;

3.14

%  define units (for our sanity)
time_units = 'hr';
space_units = 'micron'; 

% mesh/volume information

cell_volume = (4/3)*pi*10^3; % (cubic microns) % http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0134999
dr = 60; 
dr_fine = 20;
domain_size = 3000; 
vertical_scale = 1.25;  

% Time information

dt = 0.01; 
t_max = 24*100; 
output_interval = 1; % hours

% Counters

t= 0 ; 

% Parameters

% Previous oxygen and vascular field info - now in TME_script_*

% D_oxygen = 1e5 * 60; % 1e5 micron^2 / min * 60 min / hour 
% tumor_oxygen_uptake = 10 * 60; % 10  1/min * 60 min/hour 
% oxygen_decay = 0.001 * tumor_oxygen_uptake ; %  0.01 * tumor_oxygen_uptake ; % 0.1 1/min * 60 min/hour
% 
% 
% oxygen_source = 999; % goal : oxygen_far_field = 0.9 * oxygen_blood  when no tumor cells, fully vascularized 
% vascular_degradation_rate = 0.05; % 1 percent degraded per hour where there are full tumor cell density 
% 
% oxygen_far_field = 1.0; 
% oxygen_blood = oxygen_far_field * ( oxygen_source + oxygen_decay)/oxygen_source ; 

% Tumor Cell "Field"

mu = 135;%1000; %20000 % 1;%  / (pi*10^2) / 0.1  % a cell moves per a cell's worth of surface area in 6 minutes 
D = 0*0.001 * mu ; 
birth_rate = (1/18); 
apoptotic_death_rate = 0.01 * birth_rate; 
necrotic_death_rate = 0.5;%1/(1*24);%1/(3*24.0); % 1.0/(3 * 24.0 ); % 3 day survival 
apoptotic_clearance_rate = 1/8.6; 
necrotic_clearance_rate = 1/(60*24.0); 
hypoxic_o2_threshold = 5/38; % Matches PC hypoxic threshold for breast cancer, oxygen well vasculartized equals 38 mmHg
critical_o2_threshold = 2.5/38; % Matches PC critical for breast cancer, oxygen well vasculartized equals 38 mmHg
cell_spatial_proliferation_factor = 1.05; % Gives the percentage above maximum cell packing that cells can grow to (on the voxel level)
cell_spatial_mechanics_factor = 0.95; % Gives the percentage of maxium cell packing at which the cells start to spill into neighboring voxels

% Vascular "field"

max_vascular_extension_rate = 3.5; % um/hour
effective_vascular_cutoff_threshold = 1E-8; %might try 1E-7 % non-dimensional vascular density


%%%% setup 

% figure out mesh size for TME

% Fine

r_nodes = ceil( domain_size / dr_fine )+1; 
r_length = r_nodes*dr_fine; 
R = 0 : dr_fine : r_length - 1;

% Coarse

r_nodes_coarse = ceil( domain_size / dr )+1; 
r_length_coarse = r_nodes_coarse * dr; 
R_coarse = 0 : dr : r_length_coarse - 1;

% Figure out volume sizes and elements for tumor growth and mechanics

voxel_edges = 0 : dr: domain_size; % Well do I want a centered FVM or an edge?  How are the fluxes set up? ...
                                     % I have a feeling that I will just
                                     % have to go over all the
                                     % calculatations and see.  
                     
voxel_volumes = dr^3 * ones( size( voxel_edges)); 

voxel_volumes(1) = (4/3)*pi*( dr )^3; 

for i=1:length(voxel_edges)-1
    outer_radius = voxel_edges(i+1); 
    inner_radius = voxel_edges(i) ;
    
    outer_volume = (4/3)*pi*outer_radius^3 ;
    inner_volume = (4/3)*pi*inner_radius^3 ;
    
    voxel_volumes(i) = outer_volume - inner_volume;
    
end

voxel_volumes(length(voxel_edges)) = (4/3)*pi*(dr*length(voxel_edges))^3 - (4/3)*pi*(dr*(length(voxel_edges)-1))^3; % Note, I suppose this could be a problem if the simulation ever got out this far.

max_cells = voxel_volumes / cell_volume; % voxel_volumes / cell_volume ; 
max_densities = max_cells ./ voxel_volumes; 

%Allocate/initialize variables 

tumor_cell_density = zeros (r_nodes_coarse, 1);
vascular_density = ones (r_nodes_coarse, 1);
vascular_density_next = ones (r_nodes_coarse, 1);

oxygen = ones( r_nodes, 1 );
angiogenic_factor = zeros(r_nodes, 1);
hypoxic_level = zeros(r_nodes,1);

initial_tumor_cell_vector_size = size(tumor_cell_density);


%oxygen = ones(size(voxel_centers));
%vascular_density = ones(size(voxel_centers)); 

template_neighbors_list = []; 
neighbors_list = template_neighbors_list; 

template_link.i = 0; 
template_link.j = 0; 
template_link.Area = dr^2; 
template_link.dr = dr; 

for i=1:length(voxel_edges)-1
    links(i) = template_link; 
    links(i).i = i; 
    links(i).j = i+1; 
    
    outer_radius = voxel_edges(i+1); 
    
    links(i).Area = 4*pi*outer_radius^2;   
end

% Subvolumes

% for k=1:length(links);
%     i = links(k).i;
%     j = links(k).j;
%     Area = links(k).Area;
%     dr_ij = links(k).dr;
% 
%     sub_volumes(i) = 4/3*pi*((i*dr_ij)^3-((i-1)*dr_ij)^3); % volume of element i
% 
% end

total_cells = zeros( size(voxel_edges) ); 

live_cells = zeros( size(voxel_edges) ); 
apoptotic_cells = zeros( size(voxel_edges) ); 
necrotic_cells = zeros( size(voxel_edges) ); 

i = floor( length(voxel_edges)/2 );
i = find( abs( voxel_edges ) < dr/2 , 1 );
 live_cells( i ) = 27;
 % live_cells(2) = 702;

% Testing O2, AF, and Vascularization profile

for i=1:6
    
    vascular_density (i) = 0;
   live_cells(i) = cell_spatial_mechanics_factor*max_cells(i);
    
end

% for i=1:13
%     
%     live_cells(i) = 0;
%     necrotic_cells(i) = max_cells(i);
%     
% end
% 
% for i=1:(3*17)
%     
%     oxygen(i) = 0;
%     
% end

total_cells = live_cells + apoptotic_cells + necrotic_cells; 


% Intial output

pic_index = 0; 

figure(7);
clf
plot( voxel_edges, total_cells , 'b' ); 

figure(6);
h = area( voxel_edges ,[ necrotic_cells./max_cells ; apoptotic_cells./max_cells ; live_cells./max_cells ]' ); 
h(1).FaceColor = 'k'; 
h(2).FaceColor = [.25 .5 1]; 
h(3).FaceColor = 'r'; 

hold on 
plot( voxel_edges , vascular_density , 'k-o' , 'marker' , 'o', 'markerfacecolor', 'k' , 'markeredgecolor' , [0,.8,0] , 'markersize', 5 , 'linewidth', 1.5 , 'color', [0,.8,0])
hold off     

% hold on 
% plot( voxel_centers , oxygen_coarse , 'c-o' , 'marker' , 'o', 'markerfacecolor', 'c' , 'markeredgecolor' , [0.8,0,0] , 'markersize', 5 , 'linewidth', 1.5 , 'color', [0.8,0,0])
% hold off 

h1 = legend( 'necrotic', 'apoptotic' , 'live' , 'vasculature');
h1.FontSize = 12; 
h1.Location = 'southeast'; 

axis( [0 , 1.1*domain_size , 0 ,vertical_scale] ); 
set( gca , 'FontSize', 12 )
xlabel( 'position (\mum)' , 'FontSize', 12 )
title( sprintf('t = %3.2f days', t/24) , 'FontSize' , 13 )

% filename = sprintf('tumorfig%08u.fig' , pic_index ); 
% savefig(filename);
filename = sprintf('tumor%08u.jpg' , pic_index ); 
print( '-djpeg' , filename); 
pic_index = pic_index+1; 

total_vascular_mass(pic_index) = dot(vascular_density, voxel_volumes);

while( t < t_max - .1*dt )  % This ends up being used to control written output frequency, not anything else!
    t;
    
    % Call to tumor growth and movement, which subsequently calls the TME,
    % which is run at the correct TME time scale.  
    size(vascular_density);
    tumor_and_vascular_growth_and_movement_4 ;
    
    figure(7)   
    hold on 
    plot( voxel_edges, total_cells ./max_cells  , 'b' ); 
    plot( voxel_edges, oxygen_coarse, 'k' ); 
    plot( voxel_edges, vascular_density , 'r' ); 
    hold off
    
    figure(7) ; plot( voxel_edges, total_cells ./ max_cells , 'b' , voxel_edges , oxygen_coarse, 'k' , voxel_edges, vascular_density , 'r' ) ;
    legend('tumor cells', 'oxygen' , 'vascular density' ) 
    
%     filename = sprintf('TME_and_tumor%08u.jpg' , pic_index ); 
%     print( '-djpeg' , filename); 

%     figure(8); plot( voxel_centers , apoptotic_cells ,'m' , voxel_centers , necrotic_cells , 'k' )
%     legend( 'apoptotic' , 'necrotic' );

    figure(9); plot( voxel_edges , live_cells ./ (total_cells + eps ) , 'r' , voxel_edges, apoptotic_cells ./( total_cells+eps ) , 'b' , voxel_edges, necrotic_cells ./ (total_cells + eps ) , 'k' ); 
    legend( 'live fraction', 'apoptotic fraction' , 'necrotic fraction' ); 

    figure(6);
    h = area( voxel_edges ,[ necrotic_cells./max_cells ; apoptotic_cells./max_cells ; live_cells./max_cells ]' ); 
    h(1).FaceColor = 'k'; 
    h(2).FaceColor = [.25 .5 1]; 
    h(3).FaceColor = 'r'; 

   hold on 
   plot( voxel_edges , vascular_density , 'k-o' , 'marker' , 'o', 'markerfacecolor', 'k' , 'markeredgecolor' , [0,.8,0] , 'markersize', 5 , 'linewidth', 1.5 , 'color', [0,.8,0])
   hold off     

    hold on 
    plot( voxel_edges , oxygen_coarse , 'c-o' , 'marker' , 'o', 'markerfacecolor', 'c' , 'markeredgecolor' , [0.8,0,0] , 'markersize', 5 , 'linewidth', 1.5 , 'color', [0.8,0,0])
    hold off 
    
    hold on 
    plot( voxel_edges , angiogenic_factor_coarse , 'y-o' , 'marker' , 'o', 'markerfacecolor', 'y' , 'markeredgecolor' , [0,0,0.8] , 'markersize', 5 , 'linewidth', 1.5 , 'color', [0,0,0.8])
    hold off 
    
    
    h1 = legend( 'necrotic', 'apoptotic' , 'live' , 'vasculature', 'oxygen', 'angiogenic factor');
    h1.FontSize = 12; 
    h1.Location = 'southeast'; 

    axis( [0 , 1.1*domain_size , 0 ,vertical_scale] ); 
    set( gca , 'FontSize', 12 )
    xlabel( 'position (\mum)' , 'FontSize', 12 )
    title( sprintf('t = %3.2f days', t/24) , 'FontSize' , 13 )

%     filename = sprintf('tumorfig%08u.fig' , pic_index ); 
%     savefig(filename);
    
    filename = sprintf('tumor%08u.jpg' , pic_index ); 
    print( '-djpeg' , filename); 

    filename = sprintf('save%08u.mat' , pic_index ); 
    save(filename); 

    pic_index = pic_index+1; 
    
    total_vascular_mass(pic_index) = dot(vascular_density, voxel_volumes);
    
end


figure(7) ; plot( voxel_edges, total_cells ./ max_cells , 'b' , voxel_edges , oxygen_coarse, 'k' , voxel_edges, vascular_density , 'r' ) ;
legend('tumor cells', 'oxygen' , 'vascular density' ) 

% filename = sprintf('TME_and_tumor%08u.jpg' , pic_index ); 
% print( '-djpeg' , filename); 

figure(8); plot( voxel_edges , apoptotic_cells ,'m' , voxel_edges , necrotic_cells , 'k' )
legend( 'apoptotic' , 'necrotic' );

figure(9); plot( voxel_edges , live_cells ./ (total_cells + eps ) , 'r' , voxel_edges, apoptotic_cells ./( total_cells+eps ) , 'b' , voxel_edges, necrotic_cells ./ (total_cells + eps ) , 'k' ); 
legend( 'live fraction', 'apoptotic fraction' , 'necrotic fraction' ); 


figure(6);
h = area( voxel_edges ,[ necrotic_cells./max_cells ; apoptotic_cells./max_cells ; live_cells./max_cells ]' ); 
h(1).FaceColor = 'k'; 
h(2).FaceColor = [.25 .5 1]; 
h(3).FaceColor = 'r'; 

hold on 
plot( voxel_edges , vascular_density , 'k-o' , 'marker' , 'o', 'markerfacecolor', 'k' , 'markeredgecolor' , [0,.8,0] , 'markersize', 5 , 'linewidth', 1.5 , 'color', [0,.8,0])
hold off     

hold on 
plot( voxel_edges , oxygen_coarse , 'c-o' , 'marker' , 'o', 'markerfacecolor', 'c' , 'markeredgecolor' , [0.8,0,0] , 'markersize', 5 , 'linewidth', 1.5 , 'color', [0.8,0,0])
hold off 

hold on 
plot( voxel_edges , angiogenic_factor_coarse , 'y-o' , 'marker' , 'o', 'markerfacecolor', 'y' , 'markeredgecolor' , [0,0,0.8] , 'markersize', 5 , 'linewidth', 1.5 , 'color', [0,0,0.8])
hold off 

h1 = legend( 'necrotic', 'apoptotic' , 'live' , 'vasculature', 'oxygen', 'angiogenic factor');
h1.FontSize = 12; 
h1.Location = 'southeast'; 

axis( [0 , 1.1*domain_size , 0 ,vertical_scale] ); 
set( gca , 'FontSize', 12 )
xlabel( 'position (\mum)' , 'FontSize', 12 )
title( sprintf('t = %3.2f days', t/24) , 'FontSize' , 13 )

% savefig(filename);
% filename = sprintf('tumor%08u.jpg' , pic_index ); 
print( '-djpeg' , filename); 
pic_index = pic_index+1; 

total_vascular_mass(pic_index) = dot(vascular_density, voxel_volumes);

% total_vascular_mass_gain = total_vascular_mass(pic_index) - total_vascular_mass(1)
% percent_mass_gain = total_vascular_mass(1)/total_vascular_mass(pic_index)
% total_run_time_in_hours = pic_index-1
% average_mass_gained_per_hour = total_vascular_mass_gain/total_run_time_in_hours



toc;



