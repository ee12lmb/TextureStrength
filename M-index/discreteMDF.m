function [ theta, M, theta_max ] = discreteMDF(input_texture,CS)
%DISCRETEMDF calculates misorientation angles for all grains in
%input_texture
%   Takes either VPSC file path or texture cell array/matrix (depending on
%   how many timesteps there are). Misorientation matricies for every
%   combination of grains are calculated and the associated angle.
%
%   Inputs:  input_texture - either VPSC file path or texture array/matrix
%                            of euler angles
%
%   Outputs: theta         - angles of misorientation between all grains
%            M             - misorientation matricies for all grains 
%            theta_max     - max theoretical angle for this symmetry 
%
%   Usage: [ theta, misor ] = discreteMDF(input_texture,CS)


% set up MTEX package
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/mtex-4.1.3/
startup_mtex;

% add path to read files
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/readfiles/


%% Input checks

% check if input is raw VPSC, texture array or single matrix
if (ischar(input_texture) == 1)
    
    % if a file path is given then read in
    [textures,ngrains,strain,blocks] = read_VPSC(input_texture);
    
elseif (iscell(input_texture) == 1)
    
    % if texture is cell array then extract number of timesteps (blocks)
    blocks = length(input_texture);
    
    % extract number of xtls
    for i = 1:blocks
        nxtl(i) = length(input_texture{i});
    end
    
    % check that the number of grains does not change between time steps
    assert((range(nxtl) == 0), ...
        'Number of grains not consitent across time steps!')
    ngrains = nxtl(1); % assign scalar value for ngrains
        
    % assign to textures for consistancy
    textures = input_texture;
    
else % texture now must be single matrix
    
    % assign to textures for consistancy
    textures = input_texture;
    
    % ngrains is now just length of matrix
    ngrains = length(textures);
    
    % only one block of data
    blocks = 1;
    
end

%% Calculate outputs

% find the maximum angle for this crystal symmetry
[~,uniform_angles] = calcAngleDistribution(CS);
theta_max = uniform_angles(length(uniform_angles));

if (blocks == 1) % if texture only has one timestep
    
    eulers_r = textures*degree;
    
    % calculate 3x3 orientation matrix for each grain
    for i = 1:ngrains 
        
        % separate out angles for ease
        phi1 = eulers_r(1,i);
        Phi  = eulers_r(2,i);
        phi2 = eulers_r(3,i);
        
        % calculate elements of orientation matrix for this grain
        % could make this more efficient by calculating common trig
        % multiplications?
        
        g_tmp(1,1) = cos(phi1) * cos(phi2) - sin(phi1) * sin(phi2) * cos(Phi);
        g_tmp(1,2) = sin(phi1) * cos(phi2) + cos(phi1) * sin(phi2) * cos(Phi);        
        g_tmp(1,3) = sin(phi2) * sin(Phi);        
        g_tmp(2,1) = -cos(phi1) * sin(phi2) - sin(phi1) * cos(phi2) * cos(Phi);        
        g_tmp(2,2) = -sin(phi1) * sin(phi2) + cos(phi1) * cos(phi2) * cos(Phi);        
        g_tmp(2,3) = cos(phi2) * sin(Phi);        
        g_tmp(3,1) = sin(phi1) * sin(Phi);        
        g_tmp(3,2) = -cos(phi1) * sin(Phi);        
        g_tmp(3,3) = cos(Phi);
        
        g{i} = g_tmp; % store each orientation matrix in cell a cell in g
         
%          i
%          g{i}
%          inv(g{i})
%          transpose(g{i})
%          det(g{i})
%         assert(all(all((inv(g{i}) - transpose(g{i})) < eps^2)), 'Not a rotation matrix!')
%         assert((det(g{i})-1.0 < eps^2), 'Not a rotation metrix!')
    end
    
    % find misorientation angles for all grains
    Nangles = 1; % initialise counter for angles
    
    for i = 1:ngrains % loop over all grains
        for j = i+1:ngrains % compare ith grain with all others after it
            
            
            % calculate misorientation matrix between two grains
           %  M_tmp = inv(g{i}) * g{j};
            M_tmp = g{j} * transpose(g{i});
            M{Nangles} = M_tmp;
    
            % extract angle MIN angle given crystal symmetry
            theta(Nangles) = calcDisorientation(M_tmp);
                        
            Nangles = Nangles + 1;
            
        end
    end
    
else % now deal with multiple textures (e.g. cell array)
    
    for k = 1:blocks
        
        eulers_r = textures{k}*degree;
    
        % calculate 3x3 orientation matrix for each grain
        for i = 1:ngrains 

            % separate out angles for ease
            phi1 = eulers_r(1,i);
            Phi  = eulers_r(2,i);
            phi2 = eulers_r(3,i);

            % calculate elements of orientation matrix for this grain
            % could make this more efficient by calculating common trig
            % multiplications?

            g_tmp(1,1) = cos(phi1) * cos(phi2) - sin(phi1) * sin(phi2) * cos(Phi);
            g_tmp(1,2) = sin(phi1) * cos(phi2) + cos(phi1) * sin(phi2) * cos(Phi);        
            g_tmp(1,3) = sin(phi2) * sin(Phi);        
            g_tmp(2,1) = -cos(phi1) * sin(phi2) - sin(phi1) * cos(phi2) * cos(Phi);        
            g_tmp(2,2) = -sin(phi1) * sin(phi2) + cos(phi1) * cos(phi2) * cos(Phi);        
            g_tmp(2,3) = cos(phi2) * sin(Phi);        
            g_tmp(3,1) = sin(phi1) * sin(Phi);        
            g_tmp(3,2) = -cos(phi1) * sin(Phi);        
            g_tmp(3,3) = cos(Phi);

            g{i} = g_tmp; % store each orientation matrix in cell a cell in g

        end

        % find misorientation angles for all grains
        Nangles = 1; % initialise counter for angles

        for i = 1:ngrains % loop over all grains
            for j = i+1:ngrains % compare ith grain with all others after it

                % calculate misorientation matrix between two grains
                M_tmp = inv(g{i}) * g{j};
                M_tmp_array{Nangles} = M_tmp;

                % extract MIN angle given crystal symmetry 
                theta_tmp(Nangles) = calcDisorientation(M_tmp);

                Nangles = Nangles + 1;

            end
        end
        
        % now store this timesteps outputs in cell arrays
        M{k}     = M_tmp_array;
        theta{k} = theta_tmp;       
        
    end
        
    
end

end

