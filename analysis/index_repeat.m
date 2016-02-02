function [ index ] = index_repeat(input_texture,step,n,repeat,seed,varargin)
%INDEX_REPEAT Calculates a given texture index a number of times
%   A specified texture index (see inputs) is calculated 'repeat' number of
%   times using 'n' number of samples. For each calculation a different
%   seed is used (beginning at the initial 'seed') resulting in a new
%   random 'n' samples (with replacement).


tic;
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
        case 'filename'
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
assert((step < blocks),'Requested time step is out of range of input texture')
assert((n <= ngrains),'Number of samples in the file is less than the number requested')

% separate out correct block if texture is cell array (i.e. multiple blocks)
if (iscell(textures) == 1)
    textures = textures{step};
end

index = zeros(repeat,1); % initialise index array 

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


end