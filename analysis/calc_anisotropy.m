function [ uA, lmA, strain ] = calc_anisotropy(input_texture,n,seed,varargin)
%CALC_ANISOTROPY calculates two anisotropy indices for an input texture
%   
%   Takes input texture as either a VPSC or EBSD (*.ctf) file path, or a
%   cell array/matrix of a Euler angles that has already been read in,
%   randomly sampling n grains (without replacement).
%
%   CALC_ANISOTROPY calculates the universal anisotropy index (uA,
%   Ranganathan and Ostoja-Starzewski, 2008) and the Ledbetter and Migliori
%   (2006) index (lmA). If the input is VPSC, then a vector containing the
%   strain at each timestep is returned.
%
%   Uses functionality from the MSAT toolkit (available from
%   http://github.com/andreww/MSAT)
%
%   Inputs:  input_texture - file path/texture array/texture matrix 
%            n             - number of grains to use
%            seed          - allows repeatability, i.e. generates the same
%                            'random' samples if set equal to previous run
%
%   Outputs: uA            - Universal anisotropy index
%            lmA           - Ledbetter and Migliori index
%            strain        - Strain at each step (VPSC only)
%
%   Optional arguments...
%
%   'crystal'  - specfies the crystal symmetry to use, must be followed by a
%                supported crystal symmetry, e.g. 'crystal','quartz'
%
%                currently supported symmetries: 'olivine', 'quartz',
%                'post-perovskite'
%
%   'outfile'  - must be followed by a file path indicating where to dump
%                output text file. Text file includes headers giving meta
%                data information e.g. no. grains, input file etc.
%
%   Lewis Bailey - University of Leeds, School of Earth and Environment 
%   2015-16 Undergraduate final year project
%
%   References
%
%   Ledbetter, H. and Migliori, A. (2006). "A general elastic-anisotropy 
%   measure". Journal of applied physics, 100(6)
%
%   Ranganathan, S. I. and Ostoja-Starzewski, M. (2008). "Universal elastic 
%   anisotropy index". Physical Review Letters, 101(5)
%
%   Usage: [ uA, lmA, strain ] = calc_anisotropy(input_texture,n,seed,varargin)
%
%   See also: J_INDEX, M_INDEXCONT, M_INDEXDISC, INDEX_REPEAT



tic;
t = clock;
%% Setup & read data

addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/
setup_env;


% Setup defaults if no options given
crystal = 'quartz';
wantout = 1;

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
    [C_poly_av, ~] = MS_VRH(vfs, Cs, rhos);

    % The question is how anisotropic is C_poly_av. The two most useful measures are given by:
    [ uA, lmA ] = MS_anisotropy( C_poly_av );
    
else  % we have multiple strain steps
    
    for b = 1:blocks
        

        fprintf('Calculating block %i of %i...',b,blocks)
        

        % copy and roate this elasticity to align with each crystal
        Cs = zeros(6,6,n);
        for i = 1:n
            Cs(:,:,i) = MS_rotEuler(C_single, textures{b}(1,i), textures{b}(2,i), textures{b}(3,i));
        end 

        % the next step is to average the elasticity of each grain:
        rhos = ones(n,1)*rho_single; % All crystals have the same density.
        vfs  = ones(n,1);            % Same volume fraction for each point

        % normalised by MS_VRH.
        [C_poly_av, ~] = MS_VRH(vfs, Cs, rhos);

        % The question is how anisotropic is C_poly_av. The two most useful measures are given by:
        [ uA(b), lmA(b) ] = MS_anisotropy( C_poly_av );
        
        fprintf('done\n')
        
    end

%% Build output to file (if requested)

time = toc;

if (wantout == 0) % if the filepath has been given as an option
    
    % assume that filepath checked in shell script/matlab can handle this
    fid = fopen(outfile,'w'); % open file for writing - will overwrite

    switch output
        case 0     % our input was a file path so we know strain
            
            fprintf(fid,['AN3\t%i\n',...  % build header in one fprintf call
                        '+Function:\tcalc_anisotropy\n',...
                        '+Time/date:\t%i:%i %i/%i/%i\n',...
                        '+Input file:\t%s\n',...
                        '+Crystal:\t%s\n'...
                        '+Grains:\t%i\n',...
                        '+Seed:\t\t%i\n',...
                        '+Time taken(s):\t%f\n',...
                        '+Columns:\tstrain,uA,lmA\n',...
                        'Data\n'],...
                        length(uA),t(4),t(5),t(3),t(2),t(1),input_texture,crystal,n,seed,time);
            
            fprintf(fid,'%10.5f %10.5f %10.5f\n',[strain;uA;lmA]); % dump data to file
       
                
        case 1     % our input was inputted texture so we don't know strain
            
            % build header
            fprintf(fid,['AN2\t%i\n',...  % build header in one fprintf call
                        '+Function:\tcalc_anisotropy\n',...
                        '+Time/date:\t%i:%i %i/%i/%i\n',...
                        '+Input file:\tn/a - input in matlab\n',...
                        '+Crystal:\t%s\n'...
                        '+Grains:\t%i\n',...
                        '+Seed:\t\t%i\n',...
                        '+Time taken(s):\t%f\n',...
                        '+Columns:\tuA,lmA\n',...
                        'Data\n'],...
                        length(uA),t(4),t(5),t(3),t(2),t(1),crystal,n,seed,time);
            
            fprintf(fid,'%10.5f %10.5f\n',[uA;lmA]); % dump data to file
 
    end
    fclose(fid);
    
    
end
    
    
end

