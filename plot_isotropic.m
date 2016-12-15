function [] = plot_isotropic( filename )
%PLOT_PATTERN plots the isotropic calculated intensities 

set(0,'DefaultFigureWindowStyle','docked')

filename = char(filename); 
file = sprintf('OUTPUTS/%s_isotropic.csv',filename); % gets the 'filename.csv' in the folder SAMPLES
IsotropicMatrix = dlmread(file); % reads the file 'filename.csv' and saves it as an array PatternMatrix

thetaVector = IsotropicMatrix(1:end,1); 
bkg = IsotropicMatrix(1:end,2); 
IsotropicIntensityMatrix = IsotropicMatrix(1:end,3);

figure
plot(thetaVector,bkg,thetaVector,IsotropicIntensityMatrix+bkg)
end
