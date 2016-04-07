function [ uA, lmA, strain ] = calc_anisotropy(input_texture,n,seed,varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%% Setup & read data

addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/
setup_env;


% Setup defaults if no options given
crystal = 'quartz';
CS = lookupSym(crystal);


iarg = 1;
while iarg<(length(varargin))
    switch varargin{iarg}
        case 'outfile'
            
            iarg = iarg + 1; % take next argument as filename 
            outfile = varargin{iarg};
           
            % check that we are not overwriting a file
            check = exist(outfile,'file');
            assert((check == 0),'Output file already exists!')
           
            wantout = 0;  % we do want the output to file
            
         case 'crystal'  % find the appropriate symmetry 
            
            iarg = iarg + 1; % take next argument
            crystal = varargin{iarg};
            
            switch crystal
                case 'post-perovskite'
                    crystal = 'llm_mgsio3ppv';
            end
            
            CS = lookupSym(crystal);
            
        otherwise
            error('Unknown flag')
    end
    iarg = iarg + 1;
end

fprintf('\nReading data...')
% determine input type and extract relevant information
[ textures, strain, blocks, input_texture, output ] = get_inputInfo(input_texture,n,seed,crystal);
fprintf('done\n')

%% Calculate elasticity measures

% look up crystal elasticity and density for single crystal
[C_single, rho_single] = MS_elasticDB(crystal);

if (blocks == 1)  % we only have one strain step

    % copy and roate this elasticity to align with each crystal
    Cs = zeros(6,6,n);
    for i = 1:n
        Cs(:,:,i) = MS_rotEuler(C_single, textures(1,i), textures(2,i), textures(3,i));
    end 

    % the next step is to average the elasticity of each grain:

    rhos = ones(nxtls,1)*rho_single; % All crystals have the same density.
    vfs  = ones(nxtls,1);            % Same volume fraction for each point

    % normalised by MS_VRH.
    [C_poly_av, rh_poly_av] = MS_VRH(vfs, Cs, rhos);

    % The question is how anisotropic is C_poly_av. The two most useful measures are given by:

    [ uA, lmA ] = MS_anisotropy( C_poly_av );
    
else  % we have multiple strain steps
    
    for b = 1:blocks
        
        tic;
        fprintf('Calculating block %i...',b)
        

        % copy and roate this elasticity to align with each crystal
        Cs = zeros(6,6,n);
        for i = 1:n
            Cs(:,:,i) = MS_rotEuler(C_single, textures{b}(1,i), textures{b}(2,i), textures{b}(3,i));
        end 

        % the next step is to average the elasticity of each grain:

        rhos = ones(n,1)*rho_single; % All crystals have the same density.
        vfs  = ones(n,1);            % Same volume fraction for each point

        % normalised by MS_VRH.
        [C_poly_av, rh_poly_av] = MS_VRH(vfs, Cs, rhos);

        % The question is how anisotropic is C_poly_av. The two most useful measures are given by:

        [ uA(b), lmA(b) ] = MS_anisotropy( C_poly_av );
        
        t = toc;
        fprintf('done (%f seconds)\n',t)
        
    end


end

