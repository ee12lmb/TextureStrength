function [ index ] = index_repeat(input_texture,step,n,repeat,seed,varargin)
%INDEX_REPEAT Calculates a given texture index a number of times
%
%   A specified texture index (see inputs) is calculated 'repeat' number of
%   times using 'n' number of grains. For each calculation a different
%   seed is used (beginning at the initial 'seed') resulting in a new
%   random 'n' samples - essentially sampling with replacement.
%
%   Inputs:  input_texture - file path of VPSC or EBSD (*.ctf) data
%                            OR matrix/cell array of euler angles
%            step          - specifies which strain step to calculate the
%                            index if multiple steps are present
%            n             - number of grains to pull from texture
%            repeat        - the number of times to calculate the index,
%                            sampling with replacement
%            seed          - allows repeatability, i.e. generates the same
%                            'random' samples if set equal to previous run
%
%   Outputs: index         - a vector of length 'repeat' containing all
%                            index values calculated
%
%   Optional arguments...
%
%   'crystal'  - specifies the crystal symmetry to use, must be followed 
%                by a recognised crystal name e.g. ...'crystal','quartz'.
%                Accepted crystal names: 'olivine','quartz',
%                'post-perovskite'
%
%   'index'    - specifies the index to use in the calculation. Must be
%                followed by a valid index name e.g. ...'index','j'.
%                Possible index names are;
%
%                'j'  - J-index
%                'mc' - Continuous M-index
%                'md' - Discrete M-index
%
%   'bin'      - if 'md' option is selected, this specifies the bin width,
%                in degrees. Must be followed a number > 0 e.g. ...'bin',0.5
%
%   'binning'  - if 'md' option is selected, specifies the binning
%                algorithm, possible options are 'interp' or 'rebin', e.g.
%                ...'binning','interp'
%
%   'outfile'  - will output data to file specified, including meta data in
%                headers e.g. ...'outfile','PATH'
%
%   '-v'       - verbose, will output progress
%
%   Lewis Bailey - University of Leeds, School of Earth and Environment 
%   2015-16 Undergraduate final year project
%
%   Usage: [ index ] = index_repeat(input_texture,step,n,repeat,seed,varargin)
%
%   See also: J_INDEX, M_INDEXDISC, M_INDEXCONT


tic;
t = clock;
%% Setup evironment

addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/
setup_env;

%% check for optional arguments
iarg = 1;
wantout = 1; % we don't want output unless the 'filename' flag is active
method = 0;  % default method to use is J-index

% setup defautl symmetry (olivine)
CS = crystalSymmetry('Pbnm', [4.75, 10.20, 5.98]);
SS = specimenSymmetry('-1');
crystal = 'olivine';
bin = 1;
binning = 'interp';
verbose = 0; % by default dont print out progress

while iarg<=(length(varargin))
    switch lower(varargin{iarg})
        
        % deal with filename flag-----------------------------------------
        case 'outfile'
            iarg = iarg + 1; % take next argument as filename 
            outfile = varargin{iarg};
           
            % check that we are not overwriting a file
            check = exist(outfile,'file');
            assert((check == 0),'Output file already exists!')
           
            wantout = 0;  % we do want the output to file
            
        % deal with index flag--------------------------------------------   
        case 'index'
            iarg = iarg + 1; % take next argument as index type 
            index_check = varargin{iarg};
            
            % set flag to know which index function to call
            switch lower(index_check) % set method for later
                case 'j'
                    method = 0; 
                case 'mc'
                    method = 1;
                case 'md'
                    method = 2;
                otherwise
                    error('Unknown index')
            end
        
         case 'crystal'  % find the appropriate symmetry 
            
            iarg = iarg + 1; % take next argument
            crystal = varargin{iarg};
            
        case 'bin' % for use with discrete M-index (1 degree default)
            
            iarg = iarg + 1;
            bin = varargin{iarg};
            
        case 'binning' % use with discrete M-index
            
            iarg = iarg + 1;
            binning = varargin{iarg};
            
        case '-v' % verbose - for printing out progress
            
            verbose = 1;
           
        %-----------------------------------------------------------------
        otherwise
            error('Unknown flag')
    end
    iarg = iarg + 1; % move to next option
end

%% Read in all texture info
% Want whole texture to allow faster sampling later (rather than reading in
% from file each time).

% check if input is raw VPSC or texture array
if (ischar(input_texture) == 1)
    
    % check file extension to see if EBSD
    [~,~,ext] = fileparts(input_texture);
    
    if (strcmp(ext,'.ctf') == 1) % if input is path to EBSD file
        
        textures = read_EBSD(input_texture,crystal);
        blocks = 1; % only ever one block in .ctf file
        ngrains = length(textures);
        output = 1;  % cannot output strain
        
    else % input must be VPSC

        % if a file path is given then read in
        [textures,ngrains,strain,blocks] = read_VPSC(input_texture);
        output = 0; % output format (see bottom of function)
    
    end
    
else
    
   % use sample_texture to gain relevant info about texture 
   [~,blocks,ngrains] = sample_texture(input_texture,1,1);
   textures = input_texture; % for consistency later
   input_texture = 'Matlab matrix';
   output = 1;
 
end


%% Repeat index calculation for different sampled texture

% assert that inputs are vaild
assert((step <= blocks),'Requested strain step is out of range of input texture')
assert((n <= ngrains(1)),'Number of samples in the file is less than the number requested')

% separate out correct block if texture is cell array (i.e. multiple blocks)
if (iscell(textures) == 1)
    textures = textures{step};
end

index = zeros(repeat,1); % initialise index array 
ini_seed = seed; % capture initial seed for output file 

for i = 1:repeat % loop over index calculation as many times as requested

    if (verbose == 1)
        fprintf('Running iteration %i of %i...',i,repeat)
        tic;
    end
    
    if (method == 0) % J-index
        [ index(i), ~] = j_index(textures,n,seed,'crystal',crystal);
        
    elseif (method == 1) % continuous M-index
        [ index(i), ~] = m_indexCont(textures,n,seed,'crystal',crystal);
        
    elseif (method == 2) % discrete M-index
        [ index(i), ~] = m_indexDisc(textures,n,seed,'crystal',crystal,'bin',bin,'binning',binning);
    end
    
    seed = seed + 1; % now sample with a different seed 

    % case when we are using all grains
    if (n == ngrains(1)) % no need to repeat as output will be the same
        
        if (verbose == 1)
            disp('Max grains requested - all repeat values will be equal!')
        end
        
        for j = 2:repeat          % loop over requested number of times
            index(j) = index(1);  % can assign all repeat values equal
        end
        
        break % break from repeat loop
    end
    
    if (verbose == 1) 
        time = toc;
        fprintf('done (%f seconds)\n',time)
    end
    
end

%% Build output

time = toc;

if (wantout == 0) % if the filepath has been given as an option
    
    % assume that filepath checked in shell script/matlab can handle this
    fid = fopen(outfile,'w'); % open file for writing 

    switch output
        case 0     % our input was a file path so we know strain
            
            % build header
            fprintf(fid,'IR\t%i\n',length(index)); % code for read_texout 
            fprintf(fid,'+Function:\tindex_repeat\n');
            fprintf(fid,'+Time/date:\t%i:%i %i/%i/%i\n',t(4),t(5),t(3),t(2),t(1));
            fprintf(fid,'+Input file:\t%s\n',input_texture);
            
            if (method == 2) % if discrete, useful to print out bin size
                fprintf(fid,'+Index:\t\t%s (bin size: %f degrees, binning: %s)\n',lower(index_check),bin,binning);
            else
                fprintf(fid,'+Index:\t\t%s\n',lower(index_check));
            end
            
            fprintf(fid,'+Strain step:\t%i\n',step);
            fprintf(fid,'+Strain:\t%f\n',strain(step));
            fprintf(fid,'+Crystal:\t%s\n',crystal);
            fprintf(fid,'+Grains:\t%i\n',n);
            fprintf(fid,'+Seed:\t\t%i\n',ini_seed);
            fprintf(fid,'+Repeat:\t%i\n',repeat);
            fprintf(fid,'+Time taken(s):\t%f\n',time);
            fprintf(fid,'+Columns:\tIndex\n');
            fprintf(fid,'Data\n');

            for i = 1:length(index)
                fprintf(fid,'%10.5f\n',index(i));
            end
                
        case 1     % our input was inputted texture so we don't know strain
            
                        % build header
            fprintf(fid,'IR\t%i\n',length(index)); % code for read_texout 
            fprintf(fid,'+Function:\tindex_repeat\n');
            fprintf(fid,'+Time/date:\t%i:%i %i/%i/%i\n',t(4),t(5),t(3),t(2),t(1));
            fprintf(fid,'+Input file:\t%s\n',input_texture);
            
            if (method == 2) % if discrete, useful to print out bin size
                fprintf(fid,'+Index:\t\t%s (bin size: %f degrees, binning: %s)\n',lower(index_check),bin,binning);
            else
                fprintf(fid,'+Index:\t\t%s\n',lower(index_check));
            end
            
            fprintf(fid,'+Strain step:\t%i\n',step);
            fprintf(fid,'+Strain:\tn/a\n');
            fprintf(fid,'+Crystal:\t%s\n',crystal);
            fprintf(fid,'+Grains:\t%i\n',n);
            fprintf(fid,'+Seed:\t\t%i\n',ini_seed);
            fprintf(fid,'+Repeat:\t%i\n',repeat);
            fprintf(fid,'+Time taken(s):\t%f\n',time);
            fprintf(fid,'+Columns:\tIndex\n');
            fprintf(fid,'Data\n');

            for i = 1:length(index)
                fprintf(fid,'%10.5f\n',index(i));
            end

    end
    fclose(fid);
    
    
end

end