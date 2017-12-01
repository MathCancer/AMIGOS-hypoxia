function output = TME_O2_3(tumor_cell_density_full, vascular_density, oxygen, parameters)


% Solving diffusion equations using Scheme 1 from Versypt and Braatz. 
% Inner and outer BC are both Neuman, using the special inner BC scheme provided by 
% Versypt and Braatz.  Solution is obtained with an implicit method using either the native matlab solver or the
% Thomas algorithm.  In this version, the tumor is static and only the
% tumor microenvironment (O2, hypoxia, angengenic factor, vascular movement, growth, and decay) is allowed to vary.
% Outputting oxygen field as output.oxygen.

% Mesh Sizes

% Fine

r_nodes = ceil( parameters.space_size / parameters.dr )+1; 
r_length = r_nodes*parameters.dr; 
R = 0 : parameters.dr : r_length - 1;

% Coarse

r_nodes_coarse = ceil( parameters.space_size / parameters.dr_coarse )+1; 
r_length_coarse = r_nodes_coarse * parameters.dr_coarse; 
R_coarse = 0 : parameters.dr_coarse : r_length_coarse - 1;

% Allocations

oxygen_diffusion_coefficients = ones( r_nodes, 1 );
oxygen_next = ones( r_nodes, 1 );
vascular_density_fine =  ones( r_nodes, 1 );
tumor_cell_density_fine =  ones( r_nodes, 1 );

% Update diffusion coefficients

for i=1:r_nodes

    oxygen_diffusion_coefficients(i) = parameters.oxygen_diffusion_coefficient;
    
end

% Update Vascular and tumor densities at the fine resolution

vascular_density_fine = interp1(R_coarse,vascular_density,R);

tumor_cell_density_fine = interp1(R_coarse,tumor_cell_density_full,R);

% vascular_density_fine(1) = vascular_density(1);
 

% for i=4:r_nodes-2
%     
%     for j = 
%     
%     vascular_density_fine(i) = (R(i) - R_coarse(floor(i/parameters.coarse_to_fine_ratio) + 2))...
%         / (R_coarse(floor(i/parameters.coarse_to_fine_ratio)+1) - R_coarse(floor(i/parameters.coarse_to_fine_ratio) + 2)) ...
%         * vascular_density(floor(i/parameters.coarse_to_fine_ratio)+1) ...
%         - (R(i) - R_coarse(floor(i/parameters.coarse_to_fine_ratio)+1)) ...
%         / (R_coarse(floor(i/parameters.coarse_to_fine_ratio) + 2) - R_coarse(floor(i/parameters.coarse_to_fine_ratio)+1))...
%         * vascular_density(floor(i/parameters.coarse_to_fine_ratio) +2);
%    
%     tumor_cell_density_fine(i) = min(1,(R(i) - R_coarse(floor(i/parameters.coarse_to_fine_ratio) + 2))...
%         / (R_coarse(floor(i/parameters.coarse_to_fine_ratio) + 1) - R_coarse(floor(i/parameters.coarse_to_fine_ratio) + 2)) ...
%         * tumor_cell_density(floor(i/parameters.coarse_to_fine_ratio) + 1) ...
%         + (R(i) - R_coarse(floor(i/parameters.coarse_to_fine_ratio) + 1)) ...
%         / (R_coarse(floor(i/parameters.coarse_to_fine_ratio) + 2) - R_coarse(floor(i/parameters.coarse_to_fine_ratio) + 1))...
%         * tumor_cell_density(floor(i/parameters.coarse_to_fine_ratio) + 2));
%     
% end
% 
% for i=r_nodes-1
    
%     vascular_density_fine(i) = (R(i) - R_coarse(floor(i/parameters.coarse_to_fine_ratio) + 1))...
%         / (R_coarse(floor(i/parameters.coarse_to_fine_ratio)) - R_coarse(floor(i/parameters.coarse_to_fine_ratio) + 1)) ...
%         * vascular_density(floor(i/parameters.coarse_to_fine_ratio)) ...
%         + (R(i) - R_coarse(floor(i/parameters.coarse_to_fine_ratio))) ...
%         / (R_coarse(floor(i/parameters.coarse_to_fine_ratio) + 1) - R_coarse(floor(i/parameters.coarse_to_fine_ratio)))...
%         * vascular_density(floor(i/parameters.coarse_to_fine_ratio) + 1);
%     
%     tumor_cell_density_fine(i) = max(0,(R(i) - R_coarse(floor(i/parameters.coarse_to_fine_ratio) + 1))...
%         / (R_coarse(floor(i/parameters.coarse_to_fine_ratio)) - R_coarse(floor(i/parameters.coarse_to_fine_ratio) + 1)) ...
%         * tumor_cell_density(floor(i/parameters.coarse_to_fine_ratio)) ...
%         + (R(i) - R_coarse(floor(i/parameters.coarse_to_fine_ratio))) ...
%         / (R_coarse(floor(i/parameters.coarse_to_fine_ratio) + 1) - R_coarse(floor(i/parameters.coarse_to_fine_ratio)))...
%         * tumor_cell_density(floor(i/parameters.coarse_to_fine_ratio) + 1));
%     
% end
% 
% vascular_density_fine(r_nodes) = vascular_density(r_nodes_coarse);
% 
% tumor_cell_density_fine(r_nodes) = tumor_cell_density(r_nodes_coarse);

% former main loop - the implicit method implemented in radial coordinates
% with spatially variable diffusion coefficients.  
    
    % Special Cases: Boundary Conditions
    
    % Inner - straight from V and B
    
    for i=1
        for j=i
            
            solver_coefficient_matrix(i,j) = 1 + (parameters.dt/parameters.dr^2 * 6 * oxygen_diffusion_coefficients(i))...
                + parameters.dt * (parameters.oxygen_decay_rate + tumor_cell_density_fine(i)/parameters.tumor_cell_density_max...
               * parameters.tumor_cell_O2_uptake + vascular_density_fine(i)/parameters.vascular_density_max * parameters.oxygen_secretion_rate);
            
        end
    end
    
    for i=1
        for j=i+1
            
            solver_coefficient_matrix(i,j) = -parameters.dt * (6 * oxygen_diffusion_coefficients(i)...
                /parameters.dr^2);
            
        end
    end
    
    % Outer
    
    for i=r_nodes
        for j=r_nodes
            
            solver_coefficient_matrix(i,j) = 1;%parameters.dt * oxygen_diffusion_coefficients(i) / (2 * parameters.dr^2);
            
        end
    end
    
    for i=r_nodes
        for j= r_nodes-1
            
            solver_coefficient_matrix(i,j) = -1;%-parameters.dt * oxygen_diffusion_coefficients(i)/ (2 * parameters.dr^2);
            
        end
    end
    
    % General Case
    
    % alphas
    
    for i=2:r_nodes-1
       for j=i-1
           
           ii = i-1;
           
           solver_coefficient_matrix (i,j) = -parameters.dt * [oxygen_diffusion_coefficients(i) ...
               /(2 * i * parameters.dr^2) * (i-2) + oxygen_diffusion_coefficients(i-1) ...
               /(2 * parameters.dr^2)];
          
       end
    end
    
    % betas
    
    for i=2:r_nodes-1
        for j=i
            
            solver_coefficient_matrix(i,j) = 1 + parameters.dt * [oxygen_diffusion_coefficients(i) ...
                /parameters.dr^2 + oxygen_diffusion_coefficients(i+1)/(2 * parameters.dr^2) + oxygen_diffusion_coefficients(i-1)...
                /(2 * parameters.dr^2) + (parameters.oxygen_secretion_rate * vascular_density_fine(i)/parameters.vascular_density_max) ...
                + (parameters.tumor_cell_O2_uptake * tumor_cell_density_fine(i)/ parameters.tumor_cell_density_max) + parameters.oxygen_decay_rate];
            
        end
    end
    
    % gammas
    
    for i=2:r_nodes-1
        for j = i+1
            
            ii=i-1;
            
            solver_coefficient_matrix(i,j) = -parameters.dt * [oxygen_diffusion_coefficients(i) / (2 * i * parameters.dr^2) * (i+2)...
                + oxygen_diffusion_coefficients(i+1) / (2 * parameters.dr^2)];
            
        end
    end
    
    % b vector for use with native Matlab Solver
    
    b(r_nodes) = 0; %oxygen(r_nodes) + parameters.oxygen_secretion_rate * vascular_density_fine(r_nodes)...
            %/parameters.vascular_density_max * parameters.blood_oxygen_level;
    
    for i=1:r_nodes-1
        
        b(i) = oxygen(i) + parameters.dt * parameters.oxygen_secretion_rate * vascular_density_fine(i)...
            /parameters.vascular_density_max * parameters.blood_oxygen_level;
    end
    
%     for i=n
%         condition(i+1) = cond(solver_coefficient_matrix);
%     end

    % Thomas algorithm - from pg 143 of Gilat and Subramanian
    
    % Step one - Pulling out diagnonals:
    
    % Alphas - left side/below diagonal
    
    for i=2:r_nodes
       
        alphas(i) = solver_coefficient_matrix(i,i-1);
        
    end
    
    %Betas - diagonal
    
    for i=1:r_nodes
       
        betas(i) = solver_coefficient_matrix(i,i);
        
    end
    
    %Gammas - right side/above diagonal
    
    for i=1:r_nodes-1
       
        gammas(i) = solver_coefficient_matrix(i,i+1);
        
    end
    
    %Step 2: top row
    
    gammas(1) = gammas(1)/betas(1);
    
    b(1) = b(1)/betas(1);
    
    %     gammas(1) = gamma(1)/betas(1);
%     b(1) = b(1)/betas(1);
    
    %Step 3: Middle rows
    
    for i=2:r_nodes-1
       
        gammas(i) = gammas(i) / (betas(i) - alphas(i) * gammas(i-1));
        b(i) = (b(i) - alphas(i) * b(i-1)) / (betas(i)...
            - alphas(i) * gammas(i-1));
        
    end
%         gammas(i) = gammas(i)/(betas(i) - alphas(i)*gammas(i-1));
%         b(i) = (b(i) - b(i-1) * alphas(i))/(betas(i) - alphas(i)*gammas(i-1));
    %Step 4: Bottom row
    
    for i = r_nodes
    
        b(i) = (b(i) - alphas(i) * b(i-1)) / (betas(i)...
            - alphas(i) * gammas(i-1));
%     b(r_nodes) = (b(r_nodes) - b(r_nodes-1) * alphas(r_nodes))/(betas(r_nodes) - alphas(r_nodes)*gammas(r_nodes-1));        
        A = (betas(i)...
            - alphas(i) * gammas(i-1));
        
    end
    
    %Step 5: Solution
    
    oxygen_next(r_nodes) = b(r_nodes);
    
    for i = r_nodes-1:-1:1
       
      oxygen_next(i) = b(i) - gammas(i) * oxygen_next(i+1);
%     oxygen_next(r_nodes) = b(r_nodes);
%     
%     for i=r_nodes-1:-1:1
%         
%        oxygen_next(i) = b(i) - gammas(i) * oxygen_next(i+1); 
    end

    % Matlab Solver Solution

    %oxygen_next = solver_coefficient_matrix \ b'; 

    % Replace previous solution with current solution
    
    oxygen = oxygen_next;
    
% Workspace output
    
output.oxygen = oxygen;
output.vascular_density_fine = vascular_density_fine;
output.tumor_cell_density_fine = tumor_cell_density_fine;
%output.tumor_cell_density = tumor_cell_density;
%output.vascular_density = vascular_density;
%output.oxygen_s_and_s = oxygen_s_and_s;
%output.R = R;

%output.diffusion_coefficients = oxygen_diffusion_coefficients;

return;
    
    
    
    
    
    
    