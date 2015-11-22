function [ m,strain ] = m_indexDisc(input_texture,CS,n,seed)
%M_INDEXDISC calculates the M-index using a discrete method
%   Takes either a VPSC file path or a cell array/matrix of a texture that
%   has already been read in (see read_VPSC).
%
%   Inputs:  input_texture - file path/texture array/texture matrix 
%            CS            - crystal symmetry
%            n             - number of samples to pull from texture
%            seed          - allows repeatability, i.e. generates the same
%                            'random' samples if set equal to previous run
%
%   Outputs: m             - the M-index as calculated by the discrete
%                            method as outlined by Skemer (REF)
%            strain        - strain vector if input is file path (otherwise
%                            this is known before function call)
%
%   Usage: [ m, strain ] = m_indexDisc(input_texture,n,seed)

tic;

%% Setup & read data

% % set up MTEX package
% addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/mtex-4.1.3/
% startup_mtex;
% 
% % add path to read files
% addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/readfiles/


addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/
setup_env


% check if input is raw VPSC or texture array
if (ischar(input_texture) == 1)
    
    % if a file path is given then read in
    [textures,ngrains,strain,blocks] = sample_VPSC(input_texture,n,seed);
    
else
    
   % if input is not a file then pass to sample_texture to deal with
   [textures,blocks,~] = sample_texture(input_texture,n,seed);
   
   % strain information cannot be extracted from inputted texture
   %+but should already be known from previous read_VPSC
   strain = 'Input is texture - strain already extracted';
 
end


%% Calculate and bin theoretical random dist

% find theoretic distribution density
[ uniform_density, uniform_angles ] = calcAngleDistribution(CS);

% find max theoretical angle
theta_max = uniform_angles(length(uniform_angles));

% turn to degrees
theta_max = theta_max/degree;
uniform_angles = uniform_angles/degree;

% calculate bin dimensions (here set to one degree - *could take input?*)
bins = linspace(0,theta_max,theta_max+1);

% calcAngleDistribution returns angles and densities as a fraction of 300
% (not sure why, seems arbitrary). To solve, need to sum up density in the
% bins that we have defined

% find number of angles returned from calcAngleDist in our defined bins
Nangles_in_bins = histc(uniform_angles,bins);

% intialise freq sum for all bins
uniform_freq = zeros(1,length(bins));  


j = 1; % initialise loop counter
for i = 1:length(bins) % now loop through and sum relevant densities into each bin

    for null = 1:Nangles_in_bins(i) % sum up correct indicies for this bin
        
        % e.g. if three angles lie in this bin, loop three times and sum
        uniform_freq(i) = uniform_freq(i) + uniform_density(j);

        % increase counter
        j = j +1;

    end 
    
    % have increased height for this bin, now deal with increased bin width
    uniform_freq(i) = uniform_freq(i)/Nangles_in_bins(i);
    
end % end bin loop

% normalise freq
uniform_freq = uniform_freq/sum(uniform_freq);


%% Calculate and bin misorientation angles for each input texture

if (blocks == 1)
    
    % find misorientation angle distribution
    [ disorentation, ~ ] = discreteMDF(textures,CS);

    % turn angle to degrees
    disorentation = disorentation/degree;
    theta_max = theta_max/degree;

    % bin angles
    disor_freq = histc(disorentation,bins);

    % normalise 
    disor_freq = disor_freq/sum(disor_freq);

    %% Calculate M-index
    m = sum((abs(uniform_freq - disor_freq))/2);
    
else

    for i = 1:blocks 

        % find misorientation angle distribution
        [ disorentation, ~ ] = discreteMDF(textures{i},CS);

        % turn angle to degrees
        disorentation = disorentation/degree;
        theta_max = theta_max/degree;

        % bin angles
        disor_freq = histc(disorentation,bins);

        % normalise 
        disor_freq = disor_freq/sum(disor_freq);

        %% Calculate M-index
        m(i) = sum((abs(uniform_freq - disor_freq))/2);

    end

end
toc
end

