function [ m, strain ] = m_indexCont(input_texture,n,seed,varargin)
%M_INDEXCONT returns m index and strain vector for an input texture 
%   Takes either a VPSC file path or a cell array/matrix of a texture that
%   has already been read in (see read_VPSC).
%
%   Inputs:  input_texture - file path/texture array/texture matrix 
%            n             - number of samples to pull from texture
%            seed          - allows repeatability, i.e. generates the same
%                            'random' samples if set equal to previous run
%
%   Outputs: m             - the M-index as calculated by the continuous
%                            function method outlined in Mainprice (REF)
%            strain        - strain vector if input is file path (otherwise
%                            this is known before function call)
%
%   Usage: [ m, strain ] = m_indexCont(input_texture,n,seed)

tic;
%% Setup & read data

addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/
setup_env

% check if input is raw VPSC or texture array
if (ischar(input_texture) == 1)
    
    % if a file path is given then read in
    [textures,ngrains,strain,blocks] = sample_VPSC(input_texture,n,seed);
    output = 0; % specify output format for file
    
else
    
   % if input is not a file then pass to sample_texture to deal with
   [textures,blocks] = sample_texture(input_texture,n,seed);
   
   % strain information cannot be extracted from inputted texture
   %+but should already be known from previous read_VPSC
   strain = 'Input is texture - strain already extracted'; 
   output = 1;
    
end

% check for optional arguments
iarg = 1;
wantout = 1; % we don't want output unless the 'filename' flag is active

while iarg<(length(varargin))
    switch varargin{iarg}
        case 'filename'
            iarg = iarg + 1; % take next argument as filename 
            outfile = varargin{iarg};
           
            % check that we are not overwriting a file
            check = exist(outfile,'file');
            assert((check == 0),'Output file already exists!')
           
            wantout = 0;  % we do want the output to file
        otherwise
            error('Unknown flag')
    end
end

% Set up symmetry
CS = crystalSymmetry('Pbnm', [4.75, 10.20, 5.98]);
SS = specimenSymmetry('-1');

% get misorientation angle distribution for specified symmetry
[uniform_density,uniform_angles] = calcAngleDistribution(CS);
uniform_density_N = uniform_density/sum(uniform_density); % normalise

%% Calculate M-index 

% check how many textures are given (i.e. do we need to reference cell
% array or just normal matrix)
if (blocks == 1)
  
    % call function to retrive relevant MDF info
    [~, MDFdensity,~] = continuousMDF(textures,CS,SS);
    MDF_density = MDFdensity/sum(MDFdensity); % normalise
    
    % calculate M-index
    m = (sum((abs(uniform_density_N - MDF_density))/2));
    
else % now deal with multiple timesteps (e.g. cell arrays)
    
    for i = 1:blocks
        % call function to retrive relevant MDF info
        [~, MDFdensity,~] = continuousMDF(textures{i},CS,SS);
        MDF_density = MDFdensity/sum(MDFdensity); % normalise
    
        % calculate M-index
        m(i) = (sum((abs(uniform_density_N - MDF_density))/2));

    end
    
end

%% Build output

time = toc;

if (wantout == 0) % if the filepath has been given as an option
    
    % assume that filepath checked in shell script/matlab can handle this
    fid = fopen(outfile,'a'); % open file for writing (append, so can add headers in shell)

    switch output
        case 0     % our input was a file path so we know strain
            
            % build header
            fprintf(fid,'-------------------------------------------------------------\n');
            fprintf(fid,'Output data file from m_indexCont run...\n');
            fprintf(fid,'Input read from file: %s\n',input_texture);
            fprintf(fid,'Number of grains sampled: %i\tSeed: %i\n',n,seed);
            fprintf(fid,'Elapsed time (s): %f\n\n',time);
            fprintf(fid,'Data columns: Strain | M-index\n');
            fprintf(fid,'-------------------------------------------------------------\n\n');
            fprintf(fid,'++DATA++\n');

              for i = 1:length(m)
                  fprintf(fid,'%10.5f %10.5f\n',strain(i),m(i));
              end
                
        case 1     % our input was inputted texture so we don't know strain
            
            % build header
            fprintf(fid,'-------------------------------------------------------------\n');
            fprintf(fid,'Output data file from j_index run...\n');
            fprintf(fid,'Input read from texture pre-loaded in matlab, strain unknown\n');
            fprintf(fid,'Number of grains sampled: %i\tSeed: %i\n',n,seed);
            fprintf(fid,'Elapsed time (s): %f\n\n',time);
            fprintf(fid,'Data columns: Strain | M-index');
            fprintf(fid,'-------------------------------------------------------------\n\n');
            fprintf(fid,'++DATA++\n');

              for i = 1:length(m)
                  fprintf(fid,'%10.5f\n',m(i));
              end
    end
    fclose(fid);
    
    
end
end

