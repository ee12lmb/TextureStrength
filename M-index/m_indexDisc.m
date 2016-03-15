function [ m,strain,disor_freq ] = m_indexDisc(input_texture,n,seed,varargin)
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
t = clock;
%% Setup & read data

addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/
setup_env;

% check for optional arguments
iarg = 1;
wantout = 1; % we don't want output unless the 'filename' flag is active

% setup defautlts (olivine, 1 degree bins)
CS = crystalSymmetry('Pbnm', [4.75, 10.20, 5.98]);
SS = specimenSymmetry('-1');
crystal = 'olivine';
binSize = 1;
hist = 0;
binType = 0; % by default use interp binning
binning = 'interp';

while iarg<=(length(varargin))
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
            CS = lookupSym(crystal);
            
        case 'bin' % sets bin size
            
            iarg = iarg + 1; 
            binSize = varargin{iarg};
            
        case 'hist' % will plot histograms at all time steps
            
            hist = 1;
            disp('Plotting histograms...')
            
        case 'binning'
            
            iarg = iarg + 1;
            binning = varargin{iarg};
            
            switch lower(binning) % assign binning flag for later
                case 'interp'
                    binType = 0;
                case 'rebin'
                    binType = 1;
                otherwise
                    error('Unrecognised binning algorithm')
            end
            
             
        otherwise
            error('Unknown flag')
    end
    iarg = iarg + 1;
end


% determine input type and extract relevant information
[ textures, strain, blocks, input_texture, output ] = get_inputInfo(input_texture,n,seed,crystal);
   

%% Calculate and bin theoretical random dist

% find theoretic distribution density
[ uniform_density, uniform_angles ] = calcAngleDistribution(CS);


% find max theoretical angle
theta_max = uniform_angles(length(uniform_angles));

% turn to degrees
theta_max = theta_max/degree;
uniform_angles = uniform_angles/degree;
%uniform_density = uniform_density/sum(uniform_density);


% calculate bin dimensions (here set to one degree - *could take input?*)
bins = linspace(0,theta_max,((1/binSize)*theta_max)+1);
%mid_bin = bins + 0.5*binSize;
%mid_bin(length(mid_bin)) = bins(length(bins));

if (binType == 0) %----------INTERP ALGORITHM--------------------------------
    
    % initialise freq array
    uniform_freq = zeros(1,length(bins));
    
    % deal with first bin(always zero)
    uniform_freq(1) = uniform_density(1);
    
    % loop over each of our new bins
    for i = 2:length(bins)
        
        % find which of the MTEX bins is our upper limit
        for j = 2:length(uniform_angles)
            
            % found 1st MTEX angle that is bigger than current bin
            if (bins(i) < uniform_angles(j)) 
                
                up = uniform_angles(j);
                low = uniform_angles(j-1);
                
                % find the gradient
                grad = (uniform_density(j) - uniform_density(j-1))/...
                        (uniform_angles(j) - uniform_angles(j-1));
                    
                xinc = bins(i) - uniform_angles(j-1);
                yinc = grad*xinc;
                uniform_density(j-1);
                
                % use gradient to find y incriment
                % e.g. dy = m*dx
                uniform_freq(i) = ((uniform_density(j-1) + (grad*(bins(i) - uniform_angles(j-1))))/(1/300*theta_max))*binSize;
                
                break % we've found limits so move to next bin
            end % if
            
        end % loop to find bounds
                      
    end % looping over defined bins
    
    % normalise
    uniform_freq = uniform_freq/(sum(uniform_freq)*binSize);
    uniform_density = uniform_density/sum(uniform_density);
    
    figure(100)
    hold on
    plot(uniform_angles,uniform_density,'r.-')
    plot(bins,uniform_freq,'bo')
    

else % ---------------------REBINNING ALGORITHM-------------------------------

    % calcAngleDistribution returns angles and densities as a fraction of 300
    % (not sure why, seems arbitrary). To solve, need to sum up density in the
    % bins that we have defined

    % find number of angles returned from calcAngleDist in our defined bins
    Nangles_in_bins = histc(uniform_angles,bins);

    % intialise freq sum for all bins
    uniform_freq = zeros(1,length(bins));  


    j = 1; % initialise loop counter
    for i = 1:length(bins) % now loop through and sum relevant densities into each bin
    return
        if (Nangles_in_bins(i) == 0); continue; end

        for null = 1:Nangles_in_bins(i) % sum up correct indicies for this bin

            % e.g. if three angles lie in this bin, loop three times and sum
            uniform_freq(i) = uniform_freq(i) + uniform_density(j);

            % increase counter
            j = j +1;

        end 

        % have increased height for this bin, now deal with increased bin width
        uniform_freq(i) = uniform_freq(i)/(binSize*Nangles_in_bins(i));

    end % end bin loop

    % normalise freq
    %uniform_freq = uniform_freq/300
    uniform_freq = uniform_freq/sum(uniform_freq);
    uniform_density = uniform_density/sum(uniform_density);

end

if (hist == 1)
    figure(1)
    bar(bins,uniform_freq,'histc')
    axis([0 180 0 0.025])
end

%% Calculate and bin misorientation angles for each input texture

if (blocks == 1)
    
    % find misorientation angle distribution
    [ disorentation, ~ ] = discreteMDF(textures,crystal);

    % turn angle to degrees
    disorentation = disorentation/degree;

    % bin angles
    disor_freq = histc(disorentation,bins);

    % normalise 
    disor_freq = disor_freq/(sum(disor_freq)*binSize);
    
    if (hist == 1)
        figure(i+1)
        bar(bins,disor_freq,'histc')
        axis([0 180 0 0.025]) 
    end


    %% Calculate M-index
    %m = (theta_max/(length(bins))*sum((abs(uniform_freq - disor_freq))/2);
    %m = (theta_max/(2*length(bins)))*sum(abs(uniform_freq - disor_freq));
    m = (binSize/2)*sum(abs(uniform_freq - disor_freq));

else

    for i = 1:blocks 

        % find misorientation angle distribution
        [ disorentation, ~ ] = discreteMDF(textures{i},crystal);

        % turn angle to degrees
        disorentation = disorentation/degree;

        % bin angles
        disor_freq = histc(disorentation,bins);

        % normalise 
        disor_freq = disor_freq/(sum(disor_freq)*binSize);
        
        if (hist == 1)
            figure(i+1)
            bar(bins,disor_freq,'histc')
            axis([0 180 0 0.025]) 
        end

        %% Calculate M-index

        %m(i) = (theta_max/(length(bins))*sum((abs(uniform_freq - disor_freq))/2);
        %m(i) = (theta_max/(2*length(bins)))*sum(abs(uniform_freq - disor_freq));
        m(i) = (binSize/2)*sum(abs(uniform_freq - disor_freq));

    end

end


%% Build output to file (if requested *NEEDS WORK*)

time = toc;

if (wantout == 0) % if the filepath has been given as an option
    
    % assume that filepath checked in shell script/matlab can handle this
    fid = fopen(outfile,'w'); % open file for writing - will overwrite

    switch output
        case 0     % our input was a file path so we know strain
            
            fprintf(fid,['MD2\t%i\n',...  % build header in one fprintf call
                        '+Function:\tm_indexDisc\n',...
                        '+Time/date:\t%i:%i %i/%i/%i\n',...
                        '+Input file:\t%s\n',...
                        '+Crystal:\t%s\n'...
                        '+Grains:\t%i (bin size: %f degrees, binning: %s)\n',...
                        '+Seed:\t\t%i\n',...
                        '+Time taken(s):\t%f\n',...
                        '+Columns:\tStrain,M-index\n',...
                        'Data\n'],...
                        length(m),t(4),t(5),t(3),t(2),t(1),input_texture,crystal,n,binSize,binning,seed,time);
            
            fprintf(fid,'%10.5f %10.5f\n',[strain;m]); % dump data to file
       
                
        case 1     % our input was inputted texture so we don't know strain
            
            % build header
            fprintf(fid,['MD1\t%i\n',...
                         '+Function:\tm_indexDisc\n',...
                         '+Time/date:\t%i:%i %i/%i/%i\n',...
                         '+Input file:\t%s\n',...
                         '+Crystal:\t%s\n'...
                         '+Grains:\t%i (bin size: %f degrees, binning: %s)\n',...
                         '+Seed:\t\t%i\n','+Time taken(s):\t%f\n',...
                         '+Columns:\tM-index\n','Data\n'],...
                         length(m),t(4),t(5),t(3),t(2),t(1),input_texture,crystal,n,binSize,binning,seed,time);  

            fprintf(fid,'%10.5f\n',m);
 
    end
    fclose(fid);
    
    
end
   
end

