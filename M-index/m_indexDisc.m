function [ m,strain ] = m_indexDisc(input_texture,CS,n,seed,varargin)
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

addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/
setup_env


% check if input is raw VPSC or texture array
if (ischar(input_texture) == 1)
    
    % if a file path is given then read in
    [textures,ngrains,strain,blocks] = sample_VPSC(input_texture,n,seed);
    output = 0; % output format (see bottom of function)
    
else
    
   % if input is not a file then pass to sample_texture to deal with
   [textures,blocks,~] = sample_texture(input_texture,n,seed);
   
   % strain information cannot be extracted from inputted texture
   %+but should already be known from previous read_VPSC
   strain = 'Input is texture - strain already extracted';
   output = 1;
 
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


%% Build output to file (if requested *NEEDS WORK*)

time = toc;

% check that we only have two input arguments
assert((length(varargin) < 2),'Too many optional arguments!')

if(isempty(varargin))
    return % no options, do nothing
    
elseif ((ischar(varargin{1}))) % if optional argument is file path
    
    % assume that filepath checked in shell script/matlab can handle this
    fid = fopen(varargin{1},'a'); % open file for writing (append, so can add headers in shell)

    switch output
        case 0     % our input was a file path so we know strain
            
            % build header
            fprintf(fid,'-------------------------------------------------------------\n');
            fprintf(fid,'Output data file from m_indexDisc run...\n');
            fprintf(fid,'Input read from file: %s\n',input_texture);
            fprintf(fid,'Number of grains sampled: %i\tSeed: %i\n',n,seed);
            fprintf(fid,'Elapsed time (s): %f\n\n',time);
            fprintf(fid,'%10s %10s\n','Strain','M-index');
            fprintf(fid,'-------------------------------------------------------------\n');

              for i = 1:length(m)
                  fprintf(fid,'%10.5f %10.5f\n',strain(i),m(i));
              end
                
        case 1     % our input was inputted texture so we don't know strain
            
            % build header
            fprintf(fid,'-------------------------------------------------------------\n');
            fprintf(fid,'Output data file from m_indexDisc run...\n');
            fprintf(fid,'Input read from texture pre-loaded in matlab, strain unknown\n');
            fprintf(fid,'Number of grains sampled: %i\tSeed: %i\n',n,seed);
            fprintf(fid,'Elapsed time (s): %f\n\n',time);
            fprintf(fid,'%10s\n','M-index');
            fprintf(fid,'-------------------------------------------------------------\n');

              for i = 1:length(m)
                  fprintf(fid,'%10.5f\n',m(i));
              end
    end
    fclose(fid);
else
    disp('Could not output data to file...')
    disp('Final argument should be string containing output file path!')
end
    
    




end

