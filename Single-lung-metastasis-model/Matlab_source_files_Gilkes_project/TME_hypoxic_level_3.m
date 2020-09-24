function output = TME_hypoxic_level_3 (hypoxic_level, oxygen, parameters)

% Determining the hypoxic status of regions within the computational domain
% and outputting set status as output.hypoxic_level

hypoxic_level = update_hypoxic_level (oxygen, parameters );

output.hypoxic_level = hypoxic_level;

% should I have this at the coarse scale?  Does it do anything weird to
% have hypoxia be not at the cell level?  I mean, it really only exists at
% the cell level?

return;

% Subfunctions

function out = hypoxic_level_function ( oxygen , parameters )
    
    if( oxygen > parameters.hypoxic_o2_threshold )
        out = 0;
    elseif ( oxygen < parameters.critical_o2_threshold)
        out = 1;
    else 
        out = ((parameters.hypoxic_o2_threshold - oxygen)/(parameters.hypoxic_o2_threshold - parameters.critical_o2_threshold));
    end
    
    return;
   
function out = update_hypoxic_level ( oxygen, parameters )

    [R] = size( oxygen ); 
    
    for i=1:R

        out (i) = hypoxic_level_function (oxygen(i), parameters);

    end
    return;