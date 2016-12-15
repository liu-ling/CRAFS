function [] = plot_ODF( filename )
%PLOT_PATTERN plots the isotropic calculated intensities 

set(0,'DefaultFigureWindowStyle','docked')

filename = char(filename); 
file = sprintf('OUTPUTS/%s_ODF.csv',filename); % gets the 'filename.csv' in the folder SAMPLES
ODFMatrix = dlmread(file); % reads the file 'filename.csv' and saves it as an array PatternMatrix

mu = ODFMatrix(1:end,1); 
ODF = ODFMatrix(1:end,2); 
ODF_sin_mu = ODFMatrix(1:end,3);

figure
plot(mu,ODF_sin_mu)
axis tight
%title('ODF(\mu)sin(\mu)');
xlabel('\mu')
ylabel('ODF(\mu)sin(\mu)')
axis([0 90 0 max(ODF_sin_mu)])

figure
plot(mu,ODF)
axis tight
%title('ODF(\mu)');
xlabel('\mu')
ylabel('ODF(\mu)')
axis([0 90 0 max(ODF)])

end

