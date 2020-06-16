function shape = geomFunction(geomParams)
% shape = geomFunction(geomParams)
%
% Compute crater geometry 
%
% Inputs:
% geomParams = parameters that specify geometry of crater
%
% Outputs:
% shape = array containing information about crater geometry. The first
% column is the depth at 1m intervals and the second column is the 
% corresponding radius. Note that shape starts from the deepest depths.

    
    depth = ceil(geomParams(1));
    z = 0:1:depth;
    
    N = length(geomParams)-1;
    dz = depth/(N-1);
    
    i = 1:length(geomParams)-1;
    depths = (i-1).*dz;
    
    rs = geomParams(2:end);

    r = interp1(depths, rs, z);
    
    shape = [flipud(z'), flipud(r')];
    
    



