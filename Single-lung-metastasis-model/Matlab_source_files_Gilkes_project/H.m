function out = H(in)
out = 1.0; 
if( in < eps )
    out = 0.0; 
end
return; 