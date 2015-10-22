function [ runOK ] = jstrain(infile,outfile)
%JSTRAIN calculates and plots J-index as a function of strain
%   Detailed explanation goes here

% set up MTEX package
addpath /nfs/see-fs-01_teaching/ee12lmb/project/source/mtex-4.1.3/
startup_mtex

% Read whole data file
[textures,ngrains,blocks] = MVT_read_VPSC_file(infile);

% Set up symmetry
CS = crystalSymmetry('Pbnm', [4.75, 10.20, 5.98]);
SS = specimenSymmetry('-1');

% loop to calculate J for each 'block'/timestep
for i = 1:blocks-1
    
    % extract euler angles and turn to radians
    eulers_r = textures{i}*degree;
    
    % creatd object of orientations for all angles in this block
    g = orientation('Euler', eulers_r(1,:), eulers_r(2,:), eulers_r(3,:), CS, SS, 'Bunge');
    
    % calculate ODF for this block
    odf = calcODF(g,'HALFWIDTH', 10*degree, 'silent');
    
    % save J-index value for plotting against strain
    J(i) = textureindex(odf) ; 
    
    i
end

%% Plots
size = 19 ;
% create strain values vector
strain = linspace(0.02,0.5,25) ;

% plot
plot(strain,J,'-vk','MarkerSize',10,'LineWidth',2)
xlabel('Strain','FontWeight','bold','FontSize',size+2)
ylabel('J-index','FontWeight','bold','FontSize',size+2)
axis([0 0.6 1 3])

% format figure
set(gcf,'PaperPositionMode', 'manual', ...
        'PaperOrientation','landscape', ...
        'PaperUnits','normalized', ...
        'PaperPosition',[0.03 0.03 0.9 0.9])% ...
        %'Position',[0.25 0.25 1 1])

           
% format axis
set(gca,'fontsize',size)

print(outfile,'-dpdf','-r800')

runOK = 0;

end

