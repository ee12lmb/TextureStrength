function [ J ] = j_hist(input_texture,timestep,n,runs,seed)
%J_HIST Plots histogram of J at a certain strain and number of (random) samples. 
%   Detailed explanation goes here


%% Setup & read in data

addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/
setup_env

% if input_texture is file path, read in  
if (ischar(input_texture) == 1) 
    [textures,ngrains,strain,blocks] = read_VPSC(input_texture,n,seed);
    
else % if input is array then reassign for consitency 
    textures = input_texture;
    
    % find number of timesteps in intput texture
    if (iscell(textures) == 1)
        % if input is a cell array then the number of arrays is how many
        % blocks we have
        blocks = length(textures);
    else
        blocks = 1;
    end
    
end 

%% Calculate J-index 
% Loop over calculation 'runs' times, pulling out n samples for a certain
% timestep, but start with a different seed each time.

% check that requested timestep exists
assert((blocks >= timestep), ... 
    'Requested timestep out of range in input texture!')

if (iscell(textures) == 1) % do we have multiple textures (e.g. cell array)
    
    for i = 1:runs

        % call J_ngrains for this run 
        J(i) = J_ngrains(textures{timestep},n,seed);

        % update seed for next run
        seed = seed + runs; % is this random enough?

    end
else
    
    for i = 1:runs

    % call J_ngrains for this run 
    J(i) = J_ngrains(textures,n,seed);

    % update seed for next run
    seed = seed + runs; % is this random enough?

    end
    
end
    
%% Plot histogram

histogram(J)
xlabel('J-index')
ylabel('Frequency')

runOK = 0;
end

