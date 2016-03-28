function [CS, phase_names] = get_symmetry(filename)
% GET_SYMMETRY Function to correllate crystal symmetry to phase order,
% adding 'notIndexed' to unused phases.
% 
%   The function reads the EBSD file to find how many and what phases are
%   present in the sample. It also amends the names of phases to the ones 
%   accepted in the melt model, accessory phases are also discounted at 
%   this stage. 
%     
%   The outputs of the get_symmetry function are a list of crystal 
%   symmetries (CS) for each nominated phase present in the sample and a 
%   list of those phase names. If the phase is not present in the database,
%   an error message will be displayed and the code will fail.
%    
%
%   INPUT
%   filename - ctf file
%   
%   OUTPUT
%   CS - Crystal symmetries for the file
%   phase_names - List of phases within the sample
% 
%
%% ***********************************************************************
% Opens ctf file
fileID = fopen(filename, 'r');
formatSpec = '%s %s %s %s %s %s %s %s %s %s %s %s %s';
header = textscan(fileID,formatSpec,30,'CommentStyle','##',...
    'Delimiter','\t');


% Finding the number of phases present from the ctf header
number_of_phases = str2double(cell2mat(header{1,2}(14)));

if (isnan(number_of_phases)) % input file is different from standard .ctf
   
    number_of_phases = str2double(cell2mat(header{1,2}(13)));

    % phases_present represents all phases indexed via ebsd
    phase_end = 13 + number_of_phases;
    phases_present = header{1,3}(14:phase_end);
    
else
    % phases_present represents all phases indexed via ebsd
    phase_end = 14 + number_of_phases;
    phases_present = header{1,3}(15:phase_end);
    
end

fclose(fileID);

% Changing ebsd phase names to more appropriate ones
all_phases = cell(1,number_of_phases);
for i = 1:number_of_phases
    phase = (char(phases_present(i)));
    switch(phase)
        case 'Quartz-new'
            all_phases{i} = 'Quartz';
        case 'An38 Andesine'
            all_phases{i} = 'An38 And';
        case 'K_Feldspar'
            all_phases{i} = 'KFeldspar';
        case 'Ferropargasite amph'
            all_phases{i} = 'Hornblende';
        case 'Diopside CaMgSi2O6'
            all_phases{i} = 'Diopside';
        case 'Biotite - C 2/c'
            all_phases{i} = 'Biotite';
        case 'hyperst'
            all_phases{i} = 'Hypersthene';
        otherwise
            all_phases{i} = phase;
    end 
end

% Correlating crystal symmetries to phases in the same order as ctf file
CS = cell(1,(number_of_phases));
CS{1} = 'notIndexed';
for z = 1:number_of_phases
    if (strcmp(all_phases(z), 'Almandine'))
        CS{z+1} =  crystalSymmetry('m-3m', [11.531 11.531 11.531],...
            'mineral', 'Almandine', 'color', 'magenta');
    elseif (strcmp(all_phases(z), 'An0-10 Alb'))
        CS{z+1} =  crystalSymmetry('-1', [8.161 12.875 7.11],...
            [93.53,116.46,90.24]*degree, 'X||a*', 'Z||c', 'mineral',...
            'An0-10 Alb', 'color', 'blue');
    elseif (strcmp(all_phases(z), 'An25 Olig'))
        CS{z+1} =  crystalSymmetry('-1', [8.169 12.851 7.124],...
            [93.63,116.4,89.46]*degree, 'X||a*', 'Z||c', 'mineral',...
            'An25 Olig', 'color', 'blue');
    elseif (strcmp(all_phases(z), 'An38 And'))
        CS{z+1} = crystalSymmetry('-1', [8.151 12.829 7.103],...
            [93.62,116.21,89.7]*degree, 'X||a*', 'Z||c', 'mineral',...
            'An38 And', 'color', 'blue');
    elseif (strcmp(all_phases(z), 'An67 Lab'))
        CS{z+1} = crystalSymmetry('-1', [8.17 12.86 7.11],...
            [93.6,116.3,89.8]*degree, 'X||a*', 'Z||c', 'mineral',...
            'An67 Lab', 'color', 'blue');
    elseif (strcmp(all_phases(z), 'An85 Byt'))
        CS{z+1} = crystalSymmetry('-1', [8.188 12.882 14.196],...
            [93.37,116.04,90.87]*degree, 'X||a*', 'Z||c', 'mineral',...
            'An85 Byt', 'color', 'blue');
    elseif (strcmp(all_phases(z), 'An90-100 Ano'))
        CS{z+1} = crystalSymmetry('-1', [8.173 12.869 14.165],...
            [93.11,115.91,91.261]*degree, 'X||a*', 'Z||c', 'mineral',...
            'An90-100 Ano', 'color', 'blue');
    elseif (strcmp(all_phases(z), 'Andalusite'))
        CS{z+1} =  crystalSymmetry('2/m', [7.798 7.9031 5.5566],...
            [90,90,90]*degree, 'X||a*', 'Y||b', 'Z||c', 'mineral',...
            'Andalusite', 'color', 'light blue');
    elseif (strcmp(all_phases(z), 'Augite'))
        CS{z+1} = crystalSymmetry('12/m1', [9.707 8.858 5.274],...
            [90,106.52,90]*degree, 'X||a*', 'Y||b', 'Z||c', 'mineral',...
            'Augite', 'color', 'green');
    elseif (strcmp(all_phases(z), 'Biotite'))
        CS{z+1} =  crystalSymmetry('12/m1', [5.3570  9.2450 20.2340],...
            [90,94.98,90]*degree, 'X||a*', 'Y||b', 'Z||c', 'mineral',...
            'Biotite', 'color', 'orange');
    elseif (strcmp(all_phases(z), 'Calcite'))
        CS{z+1} =  crystalSymmetry('-3m1', [ 4.99 4.99 17.064],...
            [90,90,120]*degree, 'X||a*', 'Y||b', 'Z||c', 'mineral',...
            'Calcite', 'color', 'light blue');
    elseif (strcmp(all_phases(z), 'Diopside'))
        CS{z+1} =  crystalSymmetry('12/m1', [9.746 8.99 5.251],...
            [90,90,90]*degree, 'X||a*', 'Y||b', 'Z||c', 'mineral',...
            'Diopside', 'color', 'dark green');
    elseif (strcmp(all_phases(z), 'Enstatite'))
        CS{z+1} =  crystalSymmetry('2/m', [18.227 8.819 5.179],...
            [90,90,90]*degree, 'X||a*', 'Y||b', 'Z||c', 'mineral',...
            'Enstatite', 'color', 'dark green');
    elseif (strcmp(all_phases(z), 'Forsterite'))
        CS{z+1} =  crystalSymmetry('2/m', [4.7540 10.1971 5.9806],...
            [90,105.63,90]*degree, 'X||a*', 'Y||b', 'Z||c', 'mineral',...
            'Forsterite', 'color', 'dark green');
    elseif (strcmp(all_phases(z), 'Garnet'))
        CS{z+1} =  crystalSymmetry('12/m1', [11.459 11.459 11.459],...
            [90,90,90]*degree, 'X||a*', 'Y||b', 'Z||c', 'mineral',...
            'Garnet', 'color', 'red');
    elseif (strcmp(all_phases(z), 'Hornblende'))
        CS{z+1} = crystalSymmetry('12/m1', [9.895 18.119 5.332],...
            [90,105.17,90]*degree, 'X||a*', 'Y||b', 'Z||c', 'mineral',...
            'Hornblende', 'color', 'dark blue');
    elseif (strcmp(all_phases(z), 'Hypersthene'))
        CS{z+1} =  crystalSymmetry('mmm', [18.2 8.86 5.22],...
            [90,90,90]*degree, 'X||a*', 'Y||b', 'Z||c', 'mineral',...
            'Hypersthene', 'color', 'dark green');
    elseif (strcmp(all_phases(z), 'KFeldspar'))
        CS{z+1} = crystalSymmetry('12/m1', [8.572 12.961 7.222],...
            [90,115.93,90]*degree, 'X||a*', 'Y||b', 'Z||c', 'mineral',...
            'KFeldspar', 'color', 'cyan');
    elseif (strcmp(all_phases(z), 'Muscovite'))
        CS{z+1} =  crystalSymmetry('12/m1', [ 5.189 8.995 20.097],...
            [90,95.18,90]*degree, 'X||a*', 'Y||b', 'Z||c', 'mineral',...
            'Muscovite', 'color', 'yellow');
    elseif (strcmp(all_phases(z), 'Quartz'))
        CS{z+1} = crystalSymmetry('pointId', 18, [4.913 4.913 5.504],...
            [90,90,120]*degree,'X||a*', 'Y||b', 'Z||c', 'mineral',...
            'Quartz','color', 'red');           
    elseif (strcmp(all_phases(z), 'Scapolite'))
        CS{z+1} = crystalSymmetry('4/m', [12.069 12.069 7.581],...
            'mineral', 'Scapolite', 'color', 'yellow');
    elseif (strcmp(all_phases(z), 'Sillimanite'))
        CS{z+1} =  crystalSymmetry('2/m', [7.798 7.9031 5.5566],...
            [90,90,90]*degree, 'X||a*', 'Y||b', 'Z||c', 'mineral',...
            'Sillimanite', 'color', 'light blue');
    elseif (strcmp(all_phases(z), 'Staurolite'))
        CS{z+1} =  crystalSymmetry('2/m', [7.871 16.620 5.656],...
            [90,90,90]*degree, 'X||a*', 'Y||b', 'Z||c', 'mineral',...
            'Staurolite', 'color', 'light blue');
    elseif (strcmp(all_phases(z), 'Tourmaline'))
        CS{z+1} = crystalSymmetry('-3', [15.992 15.992 7.19],...
            [90,90,120]*degree, 'X||a*', 'Y||b', 'Z||c', 'mineral',...
            'Tourmaline','color', 'blue');   
    else
        CS{z+1} = 'notIndexed';
    end
end

% Selecting only the phases present in 'get_phase_data' to plot pole
% figures and perform seismic anlaysis on
phase_present = cell(1,9);
for x = 1:number_of_phases
    if (strmatch(all_phases(x), {'Almandine','An0-10 Alb','An25 Olig','An38 And',...
            'An67 Lab','An85 Byt','An90-100 Ano','Andalusite','Augite',...
            'Biotite','Calcite','Diopside','Enstatite','Forsterite','Garnet',...
            'Hornblende','Hypersthene','KFeldspar','Muscovite','Quartz','Scapolite'...
            'Sillimanite','Staurolite','Tourmaline'}))
        phase_present{x} = all_phases{x};
    else
        phase_present{x} = [];
    end         
end

phase_names = phase_present(~cellfun('isempty',phase_present));

end