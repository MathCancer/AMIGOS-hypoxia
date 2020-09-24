function output = TME_af_3 (angiogenic_factor, hypoxic_level, tumor_cell_density_live, parameters)

% Solving angiogenic factor diffusion equations using Scheme 1 from Versypt and Braatz (as the with the O2 equations).  
% Inner and outer BC are both Neuman, using the special inner BC scheme provided by 
% Versypt and Braatz.  Solution is obtained with an implicit method using either the native matlab solver or the
% Thomas algorithm.  In this version, the tumor is static and only the
% tumor microenvironment (O2, hypoxia, angengenic factor, vascular movement, growth, and decay) is allowed to vary.
% Updated TS using the test TS.  TS now working at interior.

% Mesh Sizes

% Fine

r_nodes = ceil( parameters.space_size / parameters.dr )+1; 
r_length = r_nodes*parameters.dr; 
R = 0 : parameters.dr : r_length - 1;


% Coarse

r_nodes_coarse = ceil( parameters.space_size / parameters.dr_coarse )+1; 
r_length_coarse = r_nodes_coarse * parameters.dr_coarse; 
R_coarse = 0 : parameters.dr_coarse : r_length_coarse - 1;

%allocate variables 

af_diffusion_coefficients = ones( r_nodes, 1);

solver_coefficient_matrix = zeros( r_nodes, r_nodes);
alphas = zeros(r_nodes, 1); % alphas - below diagnoal - note that first element MUST be zero
betas = zeros(r_nodes, 1); % Betas - diagonal
gammas = zeros(r_nodes-1, 1); % Gammas - above diagnonal
b = zeros(r_nodes, 1); % solution vector for implicit method

% Set up diffusion coefficients - currently constant

for i=1:r_nodes
    
    af_diffusion_coefficients(i) = parameters.af_diffusion_coefficient;
    
end

tumor_cell_density_fine = interp1(R_coarse,tumor_cell_density_live,R);
% 
%     % Special Cases: Boundary Conditions
%     
%     % Inner - straight from V and B
%     
%     for i=1
%         for j=i
%             
%             solver_coefficient_matrix(i,j) = 1 + (parameters.dt/parameters.dr^2 * 6 * af_diffusion_coefficients(i))...
%                 + parameters.dt * (parameters.af_decay_rate + tumor_cell_density_fine(i)/parameters.tumor_cell_density_max...
%                * parameters.af_secretion_rate_max);
%             
%         end
%     end
%     
%     for i=1
%         for j=i+1
%             
%             solver_coefficient_matrix(i,j) = -parameters.dt * (6 * af_diffusion_coefficients(i)...
%                 /parameters.dr^2);
%             
%         end
%     end
%     
%     % Outer
%     
%     for i=r_nodes
%         for j=r_nodes
%             
%             solver_coefficient_matrix(i,j) = 1;%parameters.dt * oxygen_diffusion_coefficients(i) / (2 * parameters.dr^2);
%             
%         end
%     end
%     
%     for i=r_nodes
%         for j= r_nodes-1
%             
%             solver_coefficient_matrix(i,j) = -1;%-parameters.dt * oxygen_diffusion_coefficients(i)/ (2 * parameters.dr^2);
%             
%         end
%     end
%     
%     % General Case
%     
%     % alphas
%     
%     for i=2:r_nodes-1
%        for j=i-1
%            
%            ii = i-1;
%            
%            solver_coefficient_matrix (i,j) = -parameters.dt * [af_diffusion_coefficients(i) ...
%                /(2 * i * parameters.dr^2) * (i-2) + af_diffusion_coefficients(i-1) ...
%                /(2 * parameters.dr^2)];
%           
%        end
%     end
%     
%     % betas
%     
%     for i=2:r_nodes-1
%         for j=i
%             
%             solver_coefficient_matrix(i,j) = 1 + parameters.dt * [af_diffusion_coefficients(i) ...
%                 /parameters.dr^2 + af_diffusion_coefficients(i+1)/(2 * parameters.dr^2) + af_diffusion_coefficients(i-1)...
%                 /(2 * parameters.dr^2) ...
%                 + (parameters.af_secretion_rate_max * tumor_cell_density_fine(i)/ parameters.tumor_cell_density_max) + parameters.af_decay_rate];
%             
%         end
%     end
%     
%     % gammas
%     
%     for i=2:r_nodes-1
%         for j = i+1
%             
%             ii=i-1;
%             
%             solver_coefficient_matrix(i,j) = -parameters.dt * [af_diffusion_coefficients(i) / (2 * i * parameters.dr^2) * (i+2)...
%                 + af_diffusion_coefficients(i+1) / (2 * parameters.dr^2)];
%             
%         end
%     end
%     
%     % b vector for use with native Matlab Solver
%     
%     b(r_nodes) = 0; %oxygen(r_nodes) + parameters.oxygen_secretion_rate * vascular_density_fine(r_nodes)...
%             %/parameters.vascular_density_max * parameters.blood_oxygen_level;
%     
%     for i=1:r_nodes-1
%         
%         b(i) = angiogenic_factor(i) + parameters.dt * parameters.af_secretion_rate_max * tumor_cell_density_fine(i)...
%             /parameters.tumor_cell_density_max * parameters.af_target_level;
%     end

% Testing now loop (above).  Old coefficient loop below.

% former main loop - the implicit method implemented in radial coordinates
% with spatially variable diffusion coefficients.  
    
    % Special Cases: Boundary Conditions
    
    % Inner - straight from V and B
    
    for i=1
        for j=i
            
            solver_coefficient_matrix(i,j) = 1 + (parameters.dt/parameters.dr^2 * 6 * af_diffusion_coefficients(i))...
                + parameters.dt * (parameters.af_decay_rate + hypoxic_level(i) * tumor_cell_density_fine(i)/parameters.tumor_cell_density_max * parameters.af_secretion_rate_max);
            
        end
    end
    
    for i=1
        for j=i+1
            
            solver_coefficient_matrix(i,j) = -parameters.dt * (6 * af_diffusion_coefficients(i)...
                /parameters.dr^2);
            
        end
    end
    
    % Outer
    
    for i=r_nodes
        for j=r_nodes
            
            solver_coefficient_matrix(i,j) = 1; %parameters.dt * oxygen_diffusion_coefficients(i) / (2 * parameters.dr^2);
            
        end
    end
    
    for i=r_nodes
        for j= r_nodes-1
            
            solver_coefficient_matrix(i,j) = -1; %parameters.dt * oxygen_diffusion_coefficients(i)/ (2 * parameters.dr^2);
            
        end
    end
    
    % General Case
    
    % alphas
    
    for i=2:r_nodes-1
       for j=i-1
           
           ii = i-1;
           
           solver_coefficient_matrix (i,j) = -parameters.dt * [af_diffusion_coefficients(i) ...
               /(2 * i * parameters.dr^2) * (i-2) + af_diffusion_coefficients(i-1) ...
               /(2 * parameters.dr^2)];
          
       end
    end
    
    % betas
    
    for i=2:r_nodes-1
        for j=i
            
            solver_coefficient_matrix(i,j) = 1 + parameters.dt * [af_diffusion_coefficients(i) ...
                /parameters.dr^2 + af_diffusion_coefficients(i+1)/(2 * parameters.dr^2) + af_diffusion_coefficients(i-1)...
                /(2 * parameters.dr^2) + (parameters.af_secretion_rate_max * tumor_cell_density_fine(i)/parameters.tumor_cell_density_max * hypoxic_level(i)) ...
                 + parameters.af_decay_rate];
            
        end
    end
    
    % gammas
    
    for i=2:r_nodes-1
        for j = i+1
            
            ii=i-1;
            
            solver_coefficient_matrix(i,j) = -parameters.dt * [af_diffusion_coefficients(i) / (2 * i * parameters.dr^2) * (i+2)...
                + af_diffusion_coefficients(i+1) / (2 * parameters.dr^2)];
            
        end
    end
    
    %b vector for use with native Matlab Solver
    
    b(r_nodes) = 0;  % Enforces OBC 0 flux
    
    for i=1:r_nodes-1
        
        b(i) = angiogenic_factor(i) + parameters.dt * parameters.af_secretion_rate_max * tumor_cell_density_fine(i)...
            /parameters.tumor_cell_density_max * hypoxic_level(i) * parameters.af_target_level;
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
    
    %Step 3: Middle rows
    
    for i=2:r_nodes-1
       
        gammas(i) = gammas(i) / (betas(i) - alphas(i) * gammas(i-1));
        b(i) = (b(i) - alphas(i) * b(i-1)) / (betas(i)...
            - alphas(i) * gammas(i-1));
        
    end
    
    %Step 4: Bottom row
    
    for i = r_nodes
    
        b(i) = (b(i) - alphas(i) * b(i-1)) / (betas(i)...
            - alphas(i) * gammas(i-1));
        
        A = (betas(i)...
            - alphas(i) * gammas(i-1));
        
    end
    
    %Step 5: Solution
    
    angiogenic_factor_next(r_nodes) = b(r_nodes);
    
    for i = r_nodes-1:-1:1
       
        
        angiogenic_factor_next(i) = b(i) - gammas(i) * angiogenic_factor_next(i+1);
        
    end


    % Matlab Solver Solution

    %angiogenic_factor_next = solver_coefficient_matrix \ b; 

    % Replace previous solution with current solution
    
    angiogenic_factor = angiogenic_factor_next;

output.angiogenic_factor = angiogenic_factor;
output.tumor_cell_density_fine = tumor_cell_density_fine;

return;