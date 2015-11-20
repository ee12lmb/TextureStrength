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



%% Setup & read data

% set up MTEX package
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/mtex-4.1.3/
startup_mtex;

% add path to read files
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/readfiles/

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

figure(1)
bar(uniform_angles,uniform_density,'histc')

% calculate bin dimensions (here set to one degree - *could take input?*)
bins = linspace(0,theta_max,theta_max+1)

% calcAngleDistribution returns angles and densities as a fraction of 300
% (not sure why, seems arbitrary). To solve, need to sum up density in the
% bins that we have defined

uniform_freq = zeros(length(bins),1);  % intialise freq sum for all bins
min_bin = 1;                           % initialise minimum bin counter

for j = 1:length(uniform_angles) % loop over all angles
    
    % if a bin has been passed, we will never find an angle lower than this
    % i.e. the return from calcAngleDistribution in monotonically
    % increasing
    
    for i = min_bin:length(bins)-1 % check over all bins

        % bin in the same way has histc (see 'help histc')
        if ((bins(i) <= uniform_angles(j)) && (uniform_angles(j) < bins(i+1))) 

            disp('************************')
            disp('Checking if statement')
            bins(i)
            bins(i+1)
            uniform_angles(j)

            
            % increase density count for this bin
            disp('input frequency')
            uniform_freq(i)
            disp('desity to add')
            uniform_density(j)
            uniform_freq(i) = uniform_freq(i) + uniform_density(j);
            disp('output frequency')
            uniform_freq(i)
            disp('************************')
            
            %min_bin = min_bin + 1;
            
            break % move on to next angle (don't loop through rest of bins) 
            
        elseif (uniform_angles(j) == bins(i+1)) % deal with final bin
            
            % histc will include anything equal to final value in final bin
            uniform_freq(i+1) = uniform_freq(i+1) + uniform_density(j);
            
        end
    end
end
% 
figure(2)
%  uniform_freq
%  length(uniform_freq)
%  sum(uniform_freq)
sum(uniform_density)
bar(bins,uniform_freq,'histc')
m = uniform_freq;
%% Calculate and bin misorientation angles

% find misorientation angle distribution
[ disorentation, ~ ] = discreteMDF(textures,CS);

% turn angle to degrees
disorentation = disorentation/degree;
theta_max = theta_max/degree;

% bin angles
disor_freq = histc(disorentation,bins);

% normalise 
disor_freq = disor_freq/sum(disor_freq);



end

