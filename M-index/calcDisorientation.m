function [ theta ] = calcDisorientation(orientation)
%calcDisorientation returns smallest angle for matrix of given symmetry
%   By premultiplying an orientation matrix by each of the 
%   crystallographically related solutions for a given symmetry, the
%   smallest possible angle is found. This is the 'disortientation' angle
%   (REF NEEDED)
%
%   Inputs:  orientation - an orientation matrix
%               symmetry - identifier for crystal symmetry (NEEDS SORTING)
%
%   Outputs:       theta - the disorientation angle
%
%   Usage: [ theta ] = calcDisorientation(orientation,symmetry)

%% Determine symmetry
% * This needs work currently only for OLIVINE *

% symmetry operators for orthrombic crystal symmetry
SymOp{1} = [1, 0, 0; 0, 1, 0; 0, 0, 1];
SymOp{2} = [1, 0, 0; 0, -1, 0; 0, 0, -1];
SymOp{3} = [-1, 0, 0; 0, 1, 0; 0, 0, -1];
SymOp{4} = [-1, 0, 0; 0, -1, 0; 0, 0, 1];


% Multiply the orientation by all symmetry operators and find the minimum
% angle of all of these
for i = 1:length(SymOp)
    
    rot_tmp   = orientation * SymOp{i};
    theta_tmp(i) = acos(((rot_tmp(1,1) + rot_tmp(2,2) + rot_tmp(3,3) - 1)/2));
    
end

theta = min(theta_tmp);

end

