% Contains logsitic growth and mechanics that take into account multple
% cell packing density targets
% Switched to a max for the modified logistic birth

net_change_in_vascular_mass = zeros(r_nodes_coarse,1);

t_temp = 0; 

 while( t_temp < output_interval - .1*dt )  % Sets up output interval for charts and workspace saving done in tumor_script.
    
     % Forming input for TME
     
    live_cell_density=live_cells./max_cells;  % Also includes necrotic tissue to capture necrotic tissue cuasing vascular degradation and oxgen depletion
    
    necrotic_cell_density = necrotic_cells./max_cells;
    
    apoptotic_cell_density = apoptotic_cells./max_cells;
    
    second_tumor_cell_vector_size = size(live_cell_density');

    if (second_tumor_cell_vector_size == initial_tumor_cell_vector_size)

         else sprintf('Error in tumor cell vector size')

    end
    
    % Calling/Updating TME
 
     output = TME_script_6(live_cell_density', necrotic_cell_density', apoptotic_cell_density', vascular_density, oxygen, angiogenic_factor, hypoxic_level, hypoxic_o2_threshold, critical_o2_threshold, dr, dr_fine, dt, domain_size, t);  
  
     % Required TME output
     
     oxygen = output.oxygen;
     vascular_density = output.vascular_density;
     angiogenic_factor = output.angiogenic_factor;
     angiogenic_factor_coarse = output.angiogenic_factor_coarse;
     velocity = output.velocity;
     vascular_chemotaxis_rate_modifier = output.vascular_chemotaxis_rate_modifier;
     grad_angiogenic_factor = output.grad_angiogenic_factor;
     vascular_birth_rate = output.vascular_birth_rate;
     vascular_death_rate = output.vascular_death_rate;
     %parameters = output.parameters;


    % AF: Fine to Coarse - See TME_vasculature for birth and
    % death rate

    % Vascular growth: birth and death - See TME_vasculature for birth and
    % death rates

    
    for i=1:r_nodes_coarse

        vascular_density_next(i) = vascular_density(i) + 60 * dt*vascular_birth_rate(i)*vascular_density(i) - 60 * dt*vascular_death_rate(i)*vascular_density(i); 

    end

    vascular_density = vascular_density_next; 

    % Vascular Chemotaxis/Advection
    
    % Velocity - see TME_vasculature for direction and scaling velocity calculations
    
    velocity = max_vascular_extension_rate*velocity;
    
    % Flux/Mechanics
    
    % Interior/First Node
    
    for k=1;
        i = links(k).i;
        j = links(k).j;
        Area = links(k).Area;
        dr_ij = links(k).dr;
        
        net_change_in_vascular_mass(i) = ...
            - dt * Area * max(0,velocity(j)) * vascular_density(i)...
            - dt * Area * min(0,velocity(j)) * vascular_density(j);

%         ...
%             + dt * max(0, velocity(j))*vascular_density(i) * Area;
        
        % Net change in vascular mass = time * interface velocity on
        % outside interior voxel * vascular_density of neighbor * shared
        % surface area (interior pointing, or negative velocity + the same
        % for a postive, or outflowing velocity from center voxel.  
        
    end
    
    % Middle Nodes
    
    for k=2:length(links);
        i = links(k).i;
        j = links(k).j;
        Outer_Area = links(k).Area;
        Inner_Area = links(k-1).Area;
        dr_ij = links(k).dr;
        
        net_change_in_vascular_mass(i) = ...
              - dt * Outer_Area * max(0,velocity(j)) * vascular_density(i) ...
              - dt * Outer_Area * min(0,velocity(j)) * vascular_density(j) ...
              + dt * Inner_Area * max(0,velocity(i)) * vascular_density(i-1) ...
              + dt * Inner_Area * min(0,velocity(i)) * vascular_density(i);
              
                
%               dt * velocity(j) * vascular_density(j)*Outer_Area...
%             - dt * velocity(i)*vascular_density(i-1) * Inner_Area;        
    
    
    
%               dt * min(0, velocity(j)) * vascular_density(j)*Outer_Area...
%             + dt * max(0, velocity(j))*vascular_density(i) * Outer_Area...
%             + dt * min(0,velocity(i)) * vascular_density(i)*Inner_Area...
%             + dt * max(0, velocity(i))*vascular_density(i-1) * Inner_Area;
        
        % net change in vascular mass at voxel i = outer interface velocity
        % * outer sphere density * outer surfarce area (if velocity is into
        % voxel).  Use density in i if velocity is positive (out of voxel).
        %  Repeat for more interior surface of voxel i.  
       
    end
    
    % Exterios/Final Node
    
    for k=length(links);
        i = links(k).i;
        j = links(k).j;
        Area = links(k).Area;
        dr_ij = links(k).dr;
        
        net_change_in_vascular_mass(j) = ...
              + dt * Area * max(0,velocity(i)) * vascular_density(i) ...
              + dt * Area * min(0,velocity(i)) * vascular_density(j);        
    
   % - dt * max(0, velocity(j))*vascular_density(i-1) * Inner_Area; 
        
%         dt * min(0,velocity(i)) * vascular_density(i)*Inner_Area...
%             + dt * max(0, velocity(i))*vascular_density(i-1) * Inner_Area;
        
    end
    
%     net_change_in_vascular_density = net_change_in_vascular_mass./voxel_volumes;
    
    for i=1:length(links)+1
        
                vascular_density_next(i) = vascular_density(i) + net_change_in_vascular_mass(i)/voxel_volumes(i);
    
    end
    
    vascular_density = vascular_density_next; 


    
    % Vascular Death

    % Apply cut off vascular density under which vasculature is considered
    % broken and unrepairable within particular voxel, such that new EC have to
    % move in

    for i=1:r_nodes_coarse

       if vascular_density(i) < effective_vascular_cutoff_threshold
           vascular_density (i) =0;
       end

    end

    %vascular_density = vascular_density_next; 

    % Oxygen Coaresing
    
    % Fine to Coarse 

    % See page 19 of June notebook

    for i=1

        oxygen_coarse(i) = oxygen(i);

    end

    for i=2:r_nodes

        volume(i) = 4/3 * pi * (R(i)-R(i-1));
        volume_weighted_oxygen(i) = oxygen(i) * volume(i);

    end

    for i=2:r_nodes_coarse-1

        volume_coarse(i) = 4/3 * pi * (R_coarse(i) - R_coarse(i-1));
        volume_weighted_oxygen_coarse(i) = volume_weighted_oxygen(3*i-2) + volume_weighted_oxygen(3*i-3) + volume_weighted_oxygen(3*i-4);
        oxygen_coarse(i) = volume_weighted_oxygen_coarse(i)/volume_coarse(i);

    end

    for i=r_nodes_coarse

        volume_coarse(i) = 4/3 * pi * (R_coarse(i) - R_coarse(i-1));
        volume_weighted_oxygen_coarse(i) = volume_weighted_oxygen(r_nodes) + volume_weighted_oxygen(r_nodes-1) + volume_weighted_oxygen(r_nodes-2);
        oxygen_coarse(i) = volume_weighted_oxygen_coarse(i)/volume_coarse(i);

    end

    % Logistic Tumor Growth - Modified logistic growth targeting
    % proliferation to be greater than maximum voxel carrying capacity
    
    for i=1:length(voxel_edges)
       b = birth_rate * max(0 , (oxygen_coarse(i)-  hypoxic_o2_threshold )/(1-hypoxic_o2_threshold) ) * max(0,(1.0 -  total_cells(i) /  ( cell_spatial_proliferation_factor * max_cells(i))) );  % Testing the cell max ... 
       d_necrotic = necrotic_death_rate * min( 1 , (hypoxic_o2_threshold - oxygen_coarse(i))/(hypoxic_o2_threshold - critical_o2_threshold)) ;
       if( d_necrotic < 0 )
           d_necrotic = 0; 
       end
       
       d = apoptotic_death_rate + d_necrotic; 
        
       LC = live_cells(i) / (1 - (b-d)*dt);  
       apoptotic_cells_next(i) = (apoptotic_cells(i)+dt*apoptotic_death_rate*LC) / ( 1+apoptotic_clearance_rate*dt ); 
       necrotic_cells_next(i) = (necrotic_cells(i)+dt*d_necrotic*LC) / ( 1+necrotic_clearance_rate*dt ); 
       live_cells_next(i) = LC; 
    end

    total_cells_next = live_cells + apoptotic_cells + necrotic_cells; 

    % Post proliferation operator updates
    
    apoptotic_cells = apoptotic_cells_next; 
    live_cells = live_cells_next; 
    necrotic_cells = necrotic_cells_next; 
    total_cells = total_cells_next; 
    cell_density = total_cells_next ./ voxel_volumes; 

    % Mechanics/Flux operator - modified for a mechnaics target
    
    for k=1:length(links);
        i = links(k).i;
        j = links(k).j; 
        Area = links(k).Area; 
        dr_ij = links(k).dr; 

        h_ij = H(total_cells_next(i)-max_cells(i)*cell_spatial_mechanics_factor) * ... 
            H( cell_density(i) - cell_density(j) ); 

        h_ji = H(total_cells_next(j)-max_cells(j)*cell_spatial_mechanics_factor) * ... 
            H( cell_density(j) - cell_density(i) ); 

        f_live_i = live_cells_next(i) / ( total_cells_next(i) + eps ); 
        f_apoptotic_i = apoptotic_cells_next(i) / ( total_cells_next(i) + eps ); 
        f_necrotic_i = necrotic_cells_next(i) / ( total_cells_next(i) + eps ); 

        f_live_j = live_cells_next(j) / ( total_cells_next(j) + eps ); 
        f_apoptotic_j = apoptotic_cells_next(j) / ( total_cells_next(j) + eps ); 
        f_necrotic_j = necrotic_cells_next(j) / ( total_cells_next(j) + eps ); 


        J_ij_live = -( mu*(h_ij*f_live_i + h_ji*f_live_j) + D )*( cell_density(j) - cell_density(i))/dr_ij; 
        J_ij_apoptotic = -( mu*(h_ij*f_apoptotic_i + h_ji*f_apoptotic_j) + D )*( cell_density(j) - cell_density(i))/dr_ij; 
        J_ij_necrotic = -( mu*(h_ij*f_necrotic_i + h_ji*f_necrotic_j) + D )*( cell_density(j) - cell_density(i))/dr_ij; 


        dN_ij = dt*J_ij_live*Area; 
        live_cells(i) = live_cells(i) - dN_ij; 
        live_cells(j) = live_cells(j) + dN_ij; 

        dN_ij = dt*J_ij_apoptotic*Area; 
        apoptotic_cells(i) = apoptotic_cells(i) - dN_ij; 
        apoptotic_cells(j) = apoptotic_cells(j) + dN_ij; 

        dN_ij = dt*J_ij_necrotic*Area; 
        necrotic_cells(i) = necrotic_cells(i) - dN_ij; 
        necrotic_cells(j) = necrotic_cells(j) + dN_ij; 

    %     total_cells(i) = total_cells(i) - dN_ij; 
    %     total_cells(j) = total_cells(j) + dN_ij; 
    end

total_cells = live_cells + apoptotic_cells + necrotic_cells; 
% boundary condition


% main loop

% for i=2:length(voxel_edges)-1
%     mu_left = mu*( H(total_cells(i-1)-max_cells(i-1))*H(total_cells(i-1)-total_cells(i)) + H(total_cells(i)-max_cells(i))*H(total_cells(i)-total_cells(i-1)) ) + D; 
%     mu_right = mu*( H(total_cells(i+1)-max_cells(i+1))*H(total_cells(i+1)-total_cells(i)) + H(total_cells(i)-max_cells(i))*H(total_cells(i)-total_cells(i+1)) ) + D; 
%     
%     total_cells_next(i) = total_cells(i) + (dt/dx^2)*( mu_left*total_cells(i-1) - (mu_left+mu_right)*total_cells(i) + mu_right*total_cells(i+1)  ); 
% end

% total_cells = total_cells_next; 

t_temp = t_temp+dt;
t = t+dt ;
end

return; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Subfunctions

% function out = vascular_birth_rate_function (angiogenic_factor_coarse, vascular_density, parameters)
%     % make a smooth cutoff, then ramp rapidly. Do I still need this?
%     
% %    coefficient = smooth_heaviside( angiogenic_factor - parameters.vascular_proliferation_threshold , 10/parameters.vascular_proliferation_threshold );
% %    out = coefficient * parameters.max_vascular_birth_rate * ( 1 - vascular_density / parameters.vascular_density_max ) ;
% 
% %  out = parameters.max_vascular_birth_rate * angiogenic_factor*( 1 - vascular_density / parameters.vascular_density_max ) ;
% %    return; 
% 
% %    out = 0.0; 
%     if( angiogenic_factor_coarse > parameters.vascular_proliferation_threshold )
%         out = parameters.max_vascular_birth_rate * ( 1 - vascular_density / parameters.vascular_density_max ); 
%     else
%         out = 0;
%     end
%     
% 
%     return; 
% 
% function out = update_vascular_birth_rate (angiogenic_factor_coarse, vascular_density, parameters)
% 
%     [N] = size( vascular_density ); 
%     
%     for i=1:N
%         
%         out(i) = vascular_birth_rate_function( angiogenic_factor_coarse(i), vascular_density(i), parameters ); 
%         
%     end
%     return; 
% 
% function out = vascular_death_rate_function( tumor_cell_density, parameters )    
%     
%      [N] = size( tumor_cell_density );
%      
%      for i=1:N
%         
%         out(i) = parameters.vascular_death_rate * (tumor_cell_density/parameters.tumor_cell_density_max);
%         
%      end
%     
%     return;
% 
% function out = update_vascular_death_rate ( tumor_cell_density, parameters )    
% 
%     [N] = size( tumor_cell_density ); 
%     
%     for i=1:N
%        
%         out(i) = vascular_death_rate_function( tumor_cell_density(i), parameters ); 
%         
%     end
%     return; 
%     
% function out = vascular_chemotaxis_rate_function(angiogenic_factor_coarse, parameters)
% 
% 
%     if angiogenic_factor_coarse<parameters.vascular_chemotaxis_threshold
%         out = 0;
% 
%     else
%         out = (angiogenic_factor_coarse - parameters.vascular_chemotaxis_threshold)/(1-parameters.vascular_chemotaxis_threshold);
%     end
% 
% return;
%     
%     
%     
% function out = update_vascular_chemotaxis_rate (angiogenic_factor_coarse, parameters)
% 
%     [N] = size(angiogenic_factor_coarse);
%     
%     for i=1:N
%         
%         out(i) = vascular_chemotaxis_rate_function (angiogenic_factor_coarse(i), parameters);
%         
%     end
% 
% return;
%     
% % function out = calculate_gradient_1D( input , parameters )
% %     [M,N] = size( input );
% %     
% %     out = zeros(M,N,2); 
% %     
% %     out(2:M-1,:,1) = input(3:M,:) - input(1:M-2,:) ; 
% %     out(:,:,1) = out(:,:,1) / 2; 
% %     out(:,:,1) = out(:,:,1) / parameters.dr_coarse; 
% %     
% %     return
%     
% function out = calculate_gradient_1D( input, parameters )
%         
%         [M] = size( input );
%         
%         N = M(1);
%         
%         out = zeros(N,1);
%         
%         
%         for i=2:N-1
%         
%             out(i) = input(i+1) - input(i-1);
%             out(i) = out(i) / (2 * parameters.dr_coarse);
%         
%         end
%         
%         return




