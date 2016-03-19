function [ theta ] = calcDisorientation(orientation,crystal)
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
%   Usage: [ theta ] = calcDisorientation(orientation,crystal)

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
        
%         % find appropirate symmetry for crystal
%         CS = lookupSym(crystal);
% 
%         % extract symmetry operators
%         extract = CS.matrix;
%         [~,~,nops] = size(extract); % find the number of operators
        
        % symmetry operators for quartz (USING LOOKUPSYMS DEFINTION
%         SymOp{1} = [1, 0, 0;...
%                     0, 1, 0;...
%                     0, 0, 1];
%                 
%         SymOp{2} = [  -0.5, -0.866,   0;...
%                      0.866,   -0.5,   0;...
%                          0,      0,   1];
%                      
%         SymOp{3} = [  -0.5,  0.866,   0;...
%                     -0.866,   -0.5,   0;...
%                          0,      0,   1];
%                      
%         SymOp{4} = [-1, 0, 0;...
%                      0,-1, 0;...
%                      0, 0,-1];
%                  
%         SymOp{5} = [   0.5,  0.866,   0;...
%                     -0.866,    0.5,   0;...
%                          0,      0,  -1];
%                      
%         SymOp{6} = [   0.5, -0.866,   0;...
%                      0.866,    0.5,   0;...
%                          0,      0,  -1];

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

