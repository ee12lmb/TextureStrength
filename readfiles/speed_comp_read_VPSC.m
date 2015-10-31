%% Script to analise read in speed of two versions of read_VPSC
clear;clc;close all;

% read_VPSC:        Uses loop of fgetl and sscanf to read in Eulers
% read_VPSC_old:    Uses fgetl loop to return pointer to correct position


% time taken to read an individual block for the old and new read_VPSC
new = [0.809148,0.771807,0.796950,0.768232,0.790912,0.751279,0.756557,0.766800,0.733605, ...
    0.768388,0.741968,0.740386,0.774137,0.731706,0.760248,0.748203,0.737556,0.773499, ...
    0.733017,0.751363,0.750588,0.731464,0.766767,0.752732,0.746390,0.755381];

old = [0.433511,0.780844,1.150376,1.552762,1.946433,2.579460,2.553688,2.928495,3.273509, ...
    3.642660,4.005739,4.365827,4.744245,5.096555,5.447208,5.837268,6.328803,6.743401, ...
    6.937615,7.222874,7.684128,7.963349,8.235077,8.588553,9.046859,9.566501];

% find the cumulative times
cumulOld = cumsum(old);
cumulNew = cumsum(new);

% create x axis array for number of blocks
nblocks = linspace(1,26,26);

%% Plots

subplot(1,2,1)
plot(nblocks,old,'bo')
hold on
plot(nblocks,new,'ro')
title('Time taken for each timestep (10,000 grains)')
xlabel('Timestep')
ylabel('Time (s)')
legend('fscanf','fgetl and sscanf','Location','northwest')
legend boxoff

subplot(1,2,2)
plot(nblocks,cumulOld,'bo')
hold on 
plot(nblocks,cumulNew,'ro')
title('Cululative time taken for each timestep (10,000 grains)')
xlabel('Timestep')
ylabel('Cumulative time (s)')
legend('fscanf','fgetl and sscanf','Location','northwest')
legend boxoff



