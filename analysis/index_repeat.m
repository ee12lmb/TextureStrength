function [ index ] = index_repeat(input_texture,step,n,repeat,seed,varargin)
%INDEX_REPEAT Calculates a given texture index a number of times
%   A specified texture index (see inputs) is calculated 'repeat' number of
%   times using 'n' number of samples. For each calculation a different
%   seed is used (beginning at the initial 'seed') resulting in a new
%   random 'n' samples (with replacement).


tic;
t = clock;
%% Setup evironment

addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/
setup_env

%% check for optional arguments
iarg = 1;
wantout = 1; % we don't want output unless the 'filename' flag is active
method = 0;  % default method to use is J-index

while iarg<(length(varargin))
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
    
    % if a file path is given then read in
    [textures,ngrains,strain,blocks] = read_VPSC(input_texture);
    output = 0; % output format (see bottom of function)
    
else
    
   % use sample_texture to gain relevant info about texture 
   [~,blocks,ngrains] = sample_texture(input_texture,1,1);
   textures = input_texture; % for consistency later
   output = 1;
 
end


%% Repeat index calculation for different sampled texture

% assert that inputs are vaild
assert((step < blocks),'Requested strain step is out of range of input texture')
assert((n <= ngrains(1)),'Number of samples in the file is less than the number requested')

% separate out correct block if texture is cell array (i.e. multiple blocks)
if (iscell(textures) == 1)
    textures = textures{step};
end

index = zeros(repeat,1); % initialise index array 
ini_seed = seed; % capture initial seed for output file 

for i = 1:repeat % loop over index calculation as many times as requested

    if (method == 0) % J-index
        [ index(i), ~] = j_index(textures,n,seed);
        
    elseif (method == 1) % continuous M-index
        [ index(i), ~] = m_indexCont(textures,n,seed);
        
    elseif (method == 2) % discrete M-index
        % need to sort this out
        error('Discrete M-index not ready')
    end
    
    seed = seed + 1; % now sample with a different seed 

end

%% Build output

time = toc;

if (wantout == 0) % if the filepath has been given as an option
    
    % assume that filepath checked in shell script/matlab can handle this
    fid = fopen(outfile,'a'); % open file for writing (append, so can add headers in shell)

    switch output
        case 0     % our input was a file path so we know strain
            
            % build header
            fprintf(fid,'IR\t%i\n',length(index)); % code for read_texout 
            fprintf(fid,'+Function:\tindex_repeat\n');
            fprintf(fid,'+Time/date:\t%i:%i %i/%i/%i\n',t(4),t(5),t(3),t(2),t(1));
            fprintf(fid,'+Input file:\t%s\n',input_texture);
            fprintf(fid,'+Index:\t\t%s\n',lower(index_check));
            fprintf(fid,'+Strain step:\t%i\n',step);
            fprintf(fid,'+Strain:\t%f\n',strain(step));
            fprintf(fid,'+Grains:\t%i\n',n);
            fprintf(fid,'+Seed:\t\t%i\n',ini_seed);
            fprintf(fid,'+Repeat:\t%i\n',repeat);
            fprintf(fid,'+Time taken(s):\t%f\n',time);
            fprintf(fid,'+Columns:\tIndex\n\n');
            fprintf(fid,'Data\n');

            for i = 1:length(index)
                fprintf(fid,'%10.5f\n',index(i));
            end
                
        case 1     % our input was inputted texture so we don't know strain
            
                        % build header
            fprintf(fid,'IR\t%i\n',length(index)); % code for read_texout 
            fprintf(fid,'+Function:\tindex_repeat\n');
            fprintf(fid,'+Time/date:\t%i:%i %i/%i/%i\n',t(4),t(5),t(3),t(2),t(1));
            fprintf(fid,'+Input file:\tn/a\n');
            fprintf(fid,'+Index:\t\t%s\n',lower(index_check));
            fprintf(fid,'+Strain step:\t%i\n',step);
            fprintf(fid,'+Strain:\tn/a\n');
            fprintf(fid,'+Grains:\t%i\n',n);
            fprintf(fid,'+Seed:\t\t%i\n',ini_seed);
            fprintf(fid,'+Repeat:\t%i\n',repeat);
            fprintf(fid,'+Time taken(s):\t%f\n',time);
            fprintf(fid,'+Columns:\tIndex\n\n');
            fprintf(fid,'Data\n');

            for i = 1:length(index)
                fprintf(fid,'%10.5f\n',index(i));
            end

    end
    fclose(fid);
    
    
end

end