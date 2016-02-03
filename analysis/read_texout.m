function [ data,type ] = read_texout(file)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%% setup environment
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/
setup_env

%% Read in file 

% open file for reading
fid = fopen(file,'r');

% determine type of file from first line
info   = textscan(fid,'%s %d',1);

% unpack info (bit clunky)
type   = info{1};
type   = type{1};
length = info{2};


switch type
    case 'IR' % reading in index repeat file
        
        for i = 1:14 % skip headers
            fgetl(fid)
        end

        data = fscanf(fid, '%f', [length 1]);
        
    otherwise 
        error('Input file type not recognised')
end


end

