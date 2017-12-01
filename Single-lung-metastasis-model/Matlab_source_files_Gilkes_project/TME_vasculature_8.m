function output = TME_vasculature_8 (angiogenic_factor, tumor_cell_density_full, vascular_density, parameters)

% Solving diffusion equations using Scheme 1 from Versypt and Braatz. 
% Inner and outer BC are both Neuman, using the special inner BC scheme provided by 
% Versypt and Braatz.  Solution is obtained with an implicit method using either the native matlab solver or the
% Thomas algorithm.  In this version, the tumor is static and only the
% tumor microenvironment (O2, hypoxia, angengenic factor, vascular movement, growth, and decay) is allowed to vary.
% Updated interior boundary condition (dp/dt = 0 strictly enforced).
% 6/28/17 11:30.  Now it is growing in the center.

% Including a vascular cut off such that below a thresholdlevel the
% vasculacutre dies completely.  In v. 3, adding speed of chemotaxis is
% only dependendt on vascular extension rate.  And now adding in a
% concentration dependency on speed ... and saturation to said speed.

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

angiogenic_factor_coarse = ones(r_nodes_coarse, 1);
vascular_density_next = zeros(r_nodes_coarse, 1);
volume = zeros(r_nodes, 1);
volume_weighted_af = zeros(r_nodes, 1);
volume_coarse =  ones(r_nodes_coarse, 1);
volume_weighted_af_coarse = ones(r_nodes_coarse, 1);
velocity = zeros(r_nodes_coarse,1);

% Fine to Coarse 

% See page 19 of June notebook

for i=1
    
    angiogenic_factor_coarse(i) = angiogenic_factor(i);
    
end

for i=2:r_nodes
    
    volume(i) = 4/3 * pi * (R(i)-R(i-1));
    
    volume_weighted_af(i) = angiogenic_factor(i) * volume(i);
    
end

for i=2:r_nodes_coarse-1
    
    volume_coarse(i) = 4/3 * pi * (R_coarse(i) - R_coarse(i-1));
        
    volume_weighted_af_coarse(i) = volume_weighted_af(3*i-2) + volume_weighted_af(3*i-3) + volume_weighted_af(3*i-4);
    
    angiogenic_factor_coarse(i) = volume_weighted_af_coarse(i)/volume_coarse(i);
    
end

for i=r_nodes_coarse

    volume_coarse(i) = 4/3 * pi * (R_coarse(i) - R_coarse(i-1));
    
    volume_weighted_af_coarse(i) = volume_weighted_af(r_nodes) + volume_weighted_af(r_nodes-1) + volume_weighted_af(r_nodes-2);
    
    angiogenic_factor_coarse(i) = volume_weighted_af_coarse(i)/volume_coarse(i);
    
end

% Updates

vascular_birth_rate = update_vascular_birth_rate (angiogenic_factor_coarse, vascular_density, parameters);
vascular_death_rate =  update_vascular_death_rate ( tumor_cell_density_full, parameters );
vascular_chemotaxis_rate_modifier = update_vascular_chemotaxis_rate_modifier (angiogenic_factor_coarse, parameters);

% size(vascular_birth_rate)
% size(vascular_death_rate)
% size(vascular_density_next)
% size(vascular_density)

% Growth and Decay of Vasculature - taken straight from Cartesian model

% Vascular growth: birth and death 

% for i=1:r_nodes_coarse
%     
%     vascular_density_next(i) = vascular_density(i) + parameters.dt*vascular_birth_rate(i)*vascular_density(i) - parameters.dt*vascular_death_rate(i)*vascular_density(i); 
%     
% end
% 
% vascular_density = vascular_density_next; 

% Vascular growth: advective term (chemotaxis-based elongation)
% Upwinding Scheme - see the notes from PM - already have the
% discretization - need to make it implicit (could try the explicit stuff
% just in case ...)  Email PM also - it is time to check in.

grad_angiogenic_factor = calculate_gradient_1D( angiogenic_factor_coarse, parameters ); 

% % Enforcing outer boundary neuman condition
% 
% grad_angiogenic_factor(r_nodes_coarse,1) = grad_angiogenic_factor(r_nodes_coarse-1,1); 


%velocity = zeros(28); % How many directions can velocity have in a
%spherically symetric model?  12.07.17 Why was this here?  Doesn't make
%sense at all ... at least not at the moment ... 

for i=2:r_nodes_coarse-1
    
%grad_angiogenic_factor(i) = 1;   

     velocity(i) = grad_angiogenic_factor(i,1)/abs((grad_angiogenic_factor(i,1)+eps))*vascular_chemotaxis_rate_modifier(i); 
     %velocity(2) = parameters.vascular_extension_rate * grad_angiogenic_factor(i,j,2); 
         
    
end

% Considering that the advection equation may be a little bit harder to
% handle in spherical coordinates than originally anticipated.  Working on
% FVM instead.  Sounds interesting; wish I had done it sooner.  

% Just looking at this seems to indicate that somethign more complicated
% may be needed - https://arxiv.org/pdf/1701.04834.pdf  They don't split up
% the dr(rho * u) term and also use a WENO scheme.  They are tough after
% all, so really not that surprising.  I hope the FVM works better.  

% for i = 1
%     
%    vascular_density_next(i) = vascular_density(i+1);%vascular_density(i) - parameters.dt * 3*(velocity(i+1) - velocity(i))/parameters.dr_coarse*vascular_density(i);
%     
% end
% 
% for i=2:r_nodes_coarse-1
%     
%          vascular_density_next(i) = vascular_density(i) - parameters.dt * (2 * velocity(i)/R_coarse(i) * vascular_density(i) ...
%              + (velocity(i+1) - velocity(i-1))/(2 * parameters.dr_coarse) * vascular_density(i) + max(0, velocity(i)) * (vascular_density(i) - vascular_density(i-1))...
%              /parameters.dr_coarse + min(0,velocity(i)) * (vascular_density(i+1) - vascular_density(i))/parameters.dr_coarse);
%          % NOTE: TESTING OUT A MINUS SIGN ON THE UPWINDING SCHEME.  
% end
% 
% % for i=r_nodes_coarse-1
% %     
% %          vascular_density_next(i) = vascular_density(i) - parameters.dt * (2 * velocity(i)/R_coarse(i) * vascular_density(i) ...
% %              + (0 - velocity(i-1))/(2 * parameters.dr_coarse) * vascular_density(i) + max(0, velocity(i)) * (vascular_density(i) - vascular_density(i-1))...
% %              /parameters.dr_coarse + min(0,velocity(i)) * (vascular_density(i+1) - vascular_density(i))/parameters.dr_coarse);
% %          
% % end
% 
% for i=r_nodes_coarse % gives a rough zero flux apporximation
%     
%     vascular_density_next(i) = vascular_density_next(i-1);
%     
% end

% Vascular Death

% Apply cut off vascular density under which vasculature is considered
% broken and unrepairable within particular voxel, such that new EC have to
% move in

% for i=1:r_nodes_coarse
%     
%    if vascular_density_next(i) < parameters.effective_vascular_cutoff_threshold
%        vascular_density_next (i) =0;
%    end
%     
% end

vascular_density = vascular_density_next; 

% Cartesian

% 
% % this is an alternate if you want to model that vascular networks grow at a constant rate once activated.              
% 
% %             velocity(1) = grad_angiogenic_factor(i,j,1); 
% %             velocity(2) = grad_angiogenic_factor(i,j,2); 
% %             velocity = velocity / (norm(velocity) + eps ); 
% %             velocity = velocity * parameters.vascular_extension_rate; 
% 
% % upwinded scheme -- here's a reference
% %  https://en.wikipedia.org/wiki/Upwind_scheme
% 
%         vascular_density_next(i,j) = vascular_density(i,j) - advection_constant_x*max(0,velocity(1))*( vascular_density(i,j) - vascular_density(i-1,j) ) ...
%             - advection_constant_x*min(0,velocity(1))*( vascular_density(i+1,j) - vascular_density(i,j) ) ... 
%             - advection_constant_y*max(0,velocity(2))*( vascular_density(i,j) - vascular_density(i,j-1) ) ... 
%             - advection_constant_y*min(0,velocity(2))*( vascular_density(i,j+1) - vascular_density(i,j) );
% 
%     end
% end


output.vascular_density = vascular_density;
output.vascular_birth_rate = vascular_birth_rate;
output.vascular_death_rate = vascular_death_rate;
output.vascular_chemotaxis_rate_modifier = vascular_chemotaxis_rate_modifier;
output.angiogenic_factor_coarse = angiogenic_factor_coarse;
output.grad_angiogenic_factor = grad_angiogenic_factor;
output.volume = volume;
output.volume_weighted_af = volume_weighted_af;
output.velocity = velocity;



return;

% Subfunctions

function out = vascular_birth_rate_function (angiogenic_factor_coarse, vascular_density, parameters)

    if( angiogenic_factor_coarse > parameters.vascular_proliferation_threshold )
        out = parameters.max_vascular_birth_rate * ( 1 - vascular_density / parameters.vascular_density_max )...
            *min(1 , max( 0, (angiogenic_factor_coarse - parameters.vascular_proliferation_threshold)/(parameters.vascular_proliferation_saturation ...
            - parameters.vascular_proliferation_threshold))); 
    else
        out = 0;
    end
    

    return; 

function out = update_vascular_birth_rate (angiogenic_factor_coarse, vascular_density, parameters)

    [N] = size( vascular_density ); 
    
    for i=1:N
        
        out(i) = vascular_birth_rate_function( angiogenic_factor_coarse(i), vascular_density(i), parameters ); 
        
    end
    return; 

function out = vascular_death_rate_function( tumor_cell_density_full, parameters )    
    
     [N] = size( tumor_cell_density_full );
     
     for i=1:N
        
        out(i) = parameters.vascular_death_rate * (tumor_cell_density_full/parameters.tumor_cell_density_max);
        
     end
    
    return;

function out = update_vascular_death_rate ( tumor_cell_density_full, parameters )    

    [N] = size( tumor_cell_density_full ); 
    
    for i=1:N
       
        out(i) = vascular_death_rate_function( tumor_cell_density_full(i), parameters ); 
        
    end
    return; 
    
function out = vascular_chemotaxis_rate_modifier_fn(angiogenic_factor_coarse, parameters)


    if angiogenic_factor_coarse<parameters.vascular_chemotaxis_threshold
        out = 0;

    else
        out = min(1 , max( 0, (angiogenic_factor_coarse - parameters.vascular_chemotaxis_threshold)/(parameters.vascular_chemotaxis_saturation - parameters.vascular_chemotaxis_threshold))); 
    end

return;
    
    
    
function out = update_vascular_chemotaxis_rate_modifier (angiogenic_factor_coarse, parameters)

    [N] = size(angiogenic_factor_coarse);
    
    for i=1:N
        
        out(i) = vascular_chemotaxis_rate_modifier_fn (angiogenic_factor_coarse(i), parameters);
        
    end

return;
    
% function out = calculate_gradient_1D( input , parameters )
%     [M,N] = size( input );
%     
%     out = zeros(M,N,2); 
%     
%     out(2:M-1,:,1) = input(3:M,:) - input(1:M-2,:) ; 
%     out(:,:,1) = out(:,:,1) / 2; 
%     out(:,:,1) = out(:,:,1) / parameters.dr_coarse; 
%     
%     return
    
function out = calculate_gradient_1D( input, parameters )
        
        [M] = size( input );
        
        N = M(1);
        
        out = zeros(N,1);
        
        
        for i=2:N-1
        
            out(i) = input(i+1) - input(i-1);
            out(i) = out(i) / (2 * parameters.dr_coarse);
        
        end
        
        return
