function [] = plot_hkl(filename,eta)
%PLOT_PATTERN plots the calculated and experimental intensities 
set(0,'DefaultFigureWindowStyle','docked')

filename = char(filename); 
file = sprintf('OUTPUTS/%s_output.txt',filename);
fid = fopen(file);
Parameters=textscan(fid,'%*s%*s%s');
fclose(fid);
Parameters=Parameters{1,1};

Xinput=size(1,33);

for i=1:1:33
    Xinput(i,1)=str2double(Parameters{i,1});
end

Xinput=Xinput';
x0=Xinput;

file2 = sprintf('OUTPUTS/%s_pattern.csv',filename); % gets the 'filename.csv' in the folder SAMPLES
PatternMatrix = dlmread(file2); % reads the file 'filename.csv' and saves it as an array PatternMatrix

thetaVector = PatternMatrix(1,2:end)'; 
etaVector = PatternMatrix(2:end,1)'; 

lambda = PatternMatrix(1,1); 
sVector(:,1) = 2*sind((thetaVector(:,1))./2)./lambda; % Calculates the vector s (sVector) from the experimental 2theta


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

load('StandardData.mat') % load the file 'CelluloseIbData.mat' that contain the arrays: CarbonCoefficients,HydrogenCoefficients,OxygenCoefficients,CarbonCoordinates,HydrogenCoordinates,OxygenCoordinates,GaussianCoefficients,PeakMatrix,StartingGuess,LowerBound and UpperBound 

StructureFactors = structure_factors(sVector,PeakMatrix,HydrogenCoefficients,CarbonCoefficients,OxygenCoefficients,HydrogenCoordinates,CarbonCoordinates,OxygenCoordinates); % Calculate the Structure Factors using the function structure_factors
LorentzFactor = lorentz_factor(sVector,lambda); % Calculates the Lorentz Factor using the function lorentz_factor

[cagl0,cagl1,cagl2,a,b,c,gamma,L200,LDiag,LDelta,L004,p200,pDiag,pDelta,p004,K,A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,C02,C04,C06,C08,muf,Gammaf,Af] = set_parameters(x0,Xinput); % uses 'x0' and 'NewXinput' to set the value of each parameter
[hklRandomIntensitiesMatrix,PeakMatrix] = random_intensities(sVector,PeakMatrix,StructureFactors,LorentzFactor,lambda,cagl0,cagl1,cagl2,a,b,c,gamma,L200,LDiag,LDelta,L004,p200,pDiag,pDelta,p004);
[PreferentialOrientationMatrix] = preferential_orientation(PeakMatrix,etaVector,GaussianCoefficients,C02,C04,C06,C08,muf,Gammaf,Af);

PeakMatrix=real(PeakMatrix);

[hklIntensitiesMatrix] = hkl_intensities(etaVector,sVector,hklRandomIntensitiesMatrix,PreferentialOrientationMatrix);

BKG = background(sVector,A0,A1,A2,A3,A4,A5,A6,A7,A8,A9); 
%--------------------------------------------------------------------------
n=etaVector(1,:)==eta;
hklIntensitiesPlot(:,:)=hklIntensitiesMatrix(:,n,:);
hklIntensitiesPlot(:,:)=K.*hklIntensitiesPlot(:,:);

hklTheta=(180/pi).*PeakMatrix(:,7);
MaxhklIntensities=max(hklIntensitiesPlot);

figure

plot(thetaVector,ExperimentalIntensityMatrix(:,n),'b',thetaVector,CalculatedIntensityMatrix(:,n),'r',thetaVector,BKG,'k',thetaVector,hklIntensitiesPlot(:,1:1:end),':',hklTheta,MaxhklIntensities,'ok')
title(sprintf('%s = %d%s','\eta',eta,'\circ'));
xlabel('2\theta')
ylabel('Intensity')
axis([10 45 0 max(ExperimentalIntensityMatrix(:,n))+100])
legend ('Experimental','Calculated','Background')
text(12,max(ExperimentalIntensityMatrix(:,n)),filename)
%--------------------------------------------------------------------------

end



