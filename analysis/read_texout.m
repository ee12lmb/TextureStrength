function [ data,type ] = read_texout(file)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


%% setup environment
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/dev/
setup_env

%% Read in file 

% open file for reading
fid = fopen(file,'r');

% determine type of file and length of data from first line
info = textscan(fid,'%s %d',1);

% unpack info (bit clunky)
type   = info{1};
type   = type{1};
length = info{2};


switch type
    case 'IR' % reading in index repeat file
        
        for i = 1:14 % skip headers
            fgetl(fid);
        end

        % read in data
        data = fscanf(fid, '%g', [length 1]);
        
    case {'MC2','MD2','J2'} % reading in two column index output file
        
        for i = 1:10 % skip headers
            fgetl(fid);
        end
        
        % read in data
        data = fscanf(fid, '%g %g', [2 length]);
        data = data';
        
    case {'MC1','MD1','J1'} % reading in single column index output file
        
        for i = 1:10
            fgetl(fid);
        end
        
        % read in data
        data = fscanf(fid, '%g', [length 1]);
        
    otherwise 
        error('Input file type not recognised')
end


end

