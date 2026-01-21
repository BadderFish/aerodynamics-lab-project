%% Convert files from XFLR5 to Matlab
clear
close all
clc

extractFirstTwoColumns('CL-alpha_complete.txt', 'CL-alpha_complete.txt');
convertAirfoilToDat('CL-alpha_complete.txt','CL-alpha_complete.dat');

extractFirstTwoColumns('CDi-alpha_complete.txt', 'CDi-alpha_complete.txt');
convertAirfoilToDat('CDi-alpha_complete.txt','CDi-alpha_complete.dat');

% extractFirstTwoColumns('induced_angles.txt', 'induced_angles.txt');
% convertAirfoilToDat('induced_angles.txt','induced_angles.dat');