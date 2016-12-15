function [] = plot_pattern( filename )
%PLOT_PATTERN plots the calculated and experimental intensities 

set(0,'DefaultFigureWindowStyle','docked')

filename = char(filename); 
file = sprintf('OUTPUTS/%s_pattern.csv',filename); % gets the 'filename.csv' in the folder SAMPLES
PatternMatrix = dlmread(file); % reads the file 'filename.csv' and saves it as an array PatternMatrix

thetaVector = PatternMatrix(1,2:end)'; 
etaVector = PatternMatrix(2:end,1)'; 

n=strfind(filename,'(')-1;
if n~=0
    filenameEx=textscan(filename,sprintf('%%%ds',n));     
    filenameEx=filenameEx{1,1};
    fileEx=filenameEx{1,1};
    file3 = sprintf('SAMPLES/%s.csv',fileEx); % gets the 'filename.csv' in the folder SAMPLES
    ExperimentalMatrix = dlmread(file3); % reads the file 'filename.csv' and saves it as an array ExperimentalMatrix
else

    file3 = sprintf('SAMPLES/%s.csv',filename); % gets the 'filename.csv' in the folder SAMPLES
    ExperimentalMatrix = dlmread(file3); % reads the file 'filename.csv' and saves it as an array ExperimentalMatrix
end  

ExperimentalIntensityMatrix = ExperimentalMatrix(2:end,2:end)'; 

CalculatedIntensityMatrix = PatternMatrix(2:end,2:end)';

figure



subplot(3,1,1);
surf(etaVector,thetaVector,ExperimentalIntensityMatrix);
title('\bf\itExperimental');
view(90,90);
axis tight
shading interp
colorbar
set(gca,'XTick', -90:90:90);
[a b]=caxis;
xlabel('\eta')

colormap(jet(800))
subplot(3,1,2);
surf(etaVector,thetaVector,CalculatedIntensityMatrix);
title('\bf\itFit');
view(90,90); 
axis tight
shading interp
colorbar
set(gca,'XTick', -90:90:90);
caxis([a b])
xlabel('\eta')

subplot(3,1,3)
surf(etaVector,thetaVector,ExperimentalIntensityMatrix-CalculatedIntensityMatrix);
title('\bf\itResidual');
view(90,90);
shading interp
axis tight
set(gca,'XTick', -90:90:90);
colorbar
xlabel('\eta')
ylabel('2\theta')

end

