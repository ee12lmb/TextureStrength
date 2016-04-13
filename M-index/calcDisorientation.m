function [ theta ] = calcDisorientation(orientation,crystal)
%CALCDISORIENTATION returns smallest angle for matrix of given symmetry
%
%   CALC_DISORIENTATION takes a rotation matrix as input. It finds all the
%   symmetric equivilents by premultiplying by all the symmetry operators
%   for a given crystal symmetry (e.g. Appendix of Wheeler et al., 2001).
%   Following this, an angle of rotation is extractred from these
%   symmetrically equivilent rotation matrices. The smallest of these
%   angles is taken as the "disorientation" (Grimmer 1979).
%
%   For use in m_indexDisc
%
%   Inputs:  orientation - an orientation matrix
%            symmetry    - identifier for crystal symmetry (NEEDS SORTING)
%
%   Outputs: theta       - the disorientation angle
%
%   Lewis Bailey - University of Leeds, School of Earth and Environment 
%   2015-16 Undergraduate final year project
%
%   References
%
%   Grimmer, H. (1979). "The distribution of disorientation angles if all
%   relative orientations of neighbouring grains are equally probable".
%   Scripta Metallurgica, 13(2), pp. 161 – 164.
%   doi:http://dx.doi.org/10.1016/0036-9748(79)90058-9.
%
%   Wheeler, J., Prior, D., Jiang, Z., Spiess, R. and Trimby, P. (2001).
%   "The petrological significance of misorientations between grains".
%   Contributions to Mineralogy and Petrology, 141(1), pp. 109–124.
%   doi:10.1007/s004100000225.
%
%   Usage: [ theta ] = calcDisorientation(orientation,crystal)
%
%   See also: DISCRETEMDF, M_INDEXDISC

%% Determine symmetry

% Take input crystal and define appropriate symmetry operators
switch lower(crystal)
    
    case {'olivine', 'post-perovskite'}

        % symmetry operators for orthrombic crystal symmetry
        SymOp{1} = [1, 0, 0;...
                    0, 1, 0;...
                    0, 0, 1];
                
        SymOp{2} = [1, 0, 0;...
                    0,-1, 0;...
                    0, 0,-1];
                
        SymOp{3} = [-1, 0, 0;...
                     0, 1, 0;...
                     0, 0,-1];
                 
        SymOp{4} = [-1, 0, 0;...
                     0,-1, 0;...
                     0, 0, 1];
         
%==========================================================================    
    case 'quartz'

          SymOp{1} = [1, 0, 0;...
                      0, 1, 0;...
                      0, 0, 1];
                  
          SymOp{2} = [0.5000   -0.8660    0.0000;...
                     -0.8660   -0.5000    0.0000;...
                      0.0000    0.0000   -1.0000];
                  
          SymOp{3} = [-0.5000   -0.8660        0;...
                       0.8660   -0.5000        0;...
                            0         0   1.0000];
                        
          SymOp{4} = [-1.0000    0.0000    0.0000;...
                       0.0000    1.0000    0.0000;...
                       0.0000    0.0000   -1.0000];
                   
          SymOp{5} = [-0.5000    0.8660         0;...
                      -0.8660   -0.5000         0;...
                            0         0    1.0000];
                        
          SymOp{6} = [0.5000    0.8660    0.0000;...
                      0.8660    0.5000    0.0000;...
                      0.0000    0.0000   -1.0000];
        
    otherwise
        % now use MTEX to extract symmetry each time - this is likely to be
        % much slower than the above methods, but this will work for any
        % symmetry

        % find appropirate symmetry for crystal
        CS = lookupSym(crystal);

        % extract symmetry operators
        extract = CS.matrix;
        [~,~,nops] = size(extract); % find the number of operators
        
        for i = 1:nops
            SymOp{i} = extract(:,:,i); % put operators in cell array
        end
    
end

% Multiply the orientation by all symmetry operators and find the minimum
% angle of all of these
for i = 1:length(SymOp)
    
    rot_tmp = SymOp{i} * orientation;
    theta_tmp(i) = acos(((rot_tmp(1,1) + rot_tmp(2,2) + rot_tmp(3,3) - 1)/2));
    
end

theta = min(theta_tmp);

end

