function [time,Rwp] = crafs(filename,Xinput)
% CRAFS  Cellulose Rietveld Analysis for Fine Structure
%   [time,Rwp] = crafs(filename,Xinput)
% Input
%   filename : name of the .csv file that contains the matrix with the experimental pattern and wavelenght
%   Xinput : the input list of 33 model parameters separated by commas. When a parameter x in Xinput is x = 99, the algorithm understands that x is a 
%            free fitting coefficient that is going to be refined. When x ~= 99, the algorithm understands that x is fixed and is going to be kept equal 
%            to the input number. The parameters must be given in the following order;
%            [cagl0,cagl1,cagl2,a,b,c,gamma,L200,LDiag,LDelta,L004,p200,pDiag,pDelta,p004,K,A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,C02,C04,C06,C08,muf,Gammaf,Af]
% Output
%   time : runtime in minutes
%   Rwp : estimate of fitting quality
%   filename_pattern.csv : has the same format of “filename.csv”, but contains calculated intensities
%   filename_isotropic.csv : contains three columns, thetaVector,bkg and IsotropicCalculatedIntensityMatrix  
%   filename_ODF.csv : contains three columns, muVector, ODF and ODFsin(muVector)
%   filename_output.txt : contains two columns: model parameters and their Xoutput values. If applicable, model parameters are preceded by “*” to 
%                         indicate parameters kept fixed during refinement. In addition to Xoutput, this file contains the residue Rwp and 
%                         integrated intensities Qcr and Qspl. 
%----------------------------------------------------------------------------------------------------------------------------------------------------

tic; % start clock

%----------------------------------------------------------------------------------------------------------------------------------------------------
%   Reading the input file
%----------------------------------------------------------------------------------------------------------------------------------------------------
filename = char(filename); 
file = sprintf('SAMPLES/%s.csv',filename); % gets the 'filename.csv' in the folder SAMPLES
ExperimentalMatrix = dlmread(file); % reads the file 'filename.csv' and saves it as an array ExperimentalMatrix

lambda = ExperimentalMatrix(1,1); 
thetaVector = ExperimentalMatrix(1,2:end)'; 
etaVector = ExperimentalMatrix(2:end,1)'; 

sVector(:,1) = 2*sind((thetaVector(:,1))./2)./lambda; % Calculates the vector s (sVector) from the experimental 2theta

ExperimentalIntensityMatrix = ExperimentalMatrix(2:end,2:end)'; 
WeightedExperimentalIntensityMatrix = ExperimentalIntensityMatrix./(ExperimentalIntensityMatrix.^0.5); 

load('StandardData.mat') % load the file 'StandardData.mat' that contain the arrays: CarbonCoefficients,HydrogenCoefficients,OxygenCoefficients,CarbonCoordinates,HydrogenCoordinates,OxygenCoordinates,GaussianCoefficients,PeakMatrix,StartingGuess,LowerBound and UpperBound 

StructureFactors = structure_factors(sVector,PeakMatrix,HydrogenCoefficients,CarbonCoefficients,OxygenCoefficients,HydrogenCoordinates,CarbonCoordinates,OxygenCoordinates); % Calculate the Structure Factors using the function structure_factors
LorentzFactor = lorentz_factor(sVector,lambda); % Calculates the Lorentz Factor using the function lorentz_factor

%----------------------------------------------------------------------------------------------------------------------------------------------------
%    Performing Refinement
%----------------------------------------------------------------------------------------------------------------------------------------------------
RefinedParameters=zeros(1,34);

for step=1:1:3 % refinement routine using the lsqcurvefit is done in three steps
    
    [x0,lb,ub,NewXinput] = refinement_step(StartingGuess,LowerBound,UpperBound,RefinedParameters,Xinput,RefinementStepMatrix,step); % uses the function refinement_step to generate the NewXinput according to the refinement step
    f = @(x0,sVector)calculated_weighted_intensities(x0,NewXinput,sVector,PeakMatrix,StructureFactors,LorentzFactor,etaVector,lambda,ExperimentalIntensityMatrix,GaussianCoefficients); % f is the calculated_weighted_intensities function that returns the calculated intensities    
    options = optimset('TolX',1e-8);
    [Xoutput] = lsqcurvefit(f,x0,sVector,WeightedExperimentalIntensityMatrix,lb,ub,options);% the lsqcurvefit returns the vector Xoutput with the refined parameters 
    
    x0=Xoutput; % x0 is formed by the values found in previous refinement
    [cagl0,cagl1,cagl2,a,b,c,gamma,L200,LDiag,LDelta,L004,p200,pDiag,pDelta,p004,K,A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,C02,C04,C06,C08,muf,Gammaf,Af] = set_parameters(x0,NewXinput); % uses 'x0' and 'NewXinput' to set the value of each parameter
    RefinedParameters = [cagl0,cagl1,cagl2,a,b,c,gamma,L200,LDiag,LDelta,L004,p200,pDiag,pDelta,p004,K,A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,C02,C04,C06,C08,muf,Gammaf,Af]; 
    
end   

%----------------------------------------------------------------------------------------------------------------------------------------------------
%    Calculating Outputs
%----------------------------------------------------------------------------------------------------------------------------------------------------
CalculatedIntensitiesMatrix = (calculated_weighted_intensities(x0,Xinput,sVector,PeakMatrix,StructureFactors,LorentzFactor,etaVector,lambda,ExperimentalIntensityMatrix,GaussianCoefficients)).*(ExperimentalIntensityMatrix.^0.5); 
BKG = background(sVector,A0,A1,A2,A3,A4,A5,A6,A7,A8,A9);
hklRandomIntensitiesMatrix = random_intensities(sVector,PeakMatrix,StructureFactors,LorentzFactor,lambda,cagl0,cagl1,cagl2,a,b,c,gamma,L200,LDiag,LDelta,L004,p200,pDiag,pDelta,p004);
IsotropicCalculatedIntensitiesMatrix = K.*sum(hklRandomIntensitiesMatrix,3);

muVector = (0:1:90)';
ODF = odf(C02,C04,C06,C08,muVector,muf,Gammaf,Af,GaussianCoefficients);

Rwp = (sum(sum((1./ExperimentalIntensityMatrix).*(ExperimentalIntensityMatrix-CalculatedIntensitiesMatrix).^2))/sum(sum((1./ExperimentalIntensityMatrix).*ExperimentalIntensityMatrix.^2)))^0.5;

    smin1=find((sVector>0.105)&(sVector<0.115),1,'first');
    smin2=find((sVector>0.105)&(sVector<0.115),1,'last');
    smin3=smin1+round((smin2-smin1)/2);
        
    smax1=find((sVector>0.495)&(sVector<0.505),1,'first');
    smax2=find((sVector>0.495)&(sVector<0.505),1,'last');
    smax3=smax1+round((smax2-smax1)/2);
    
    IntegralLimits=sVector(smin3:smax3,1); % s range defined by integral limits

if isempty(IntegralLimits) 
    Qcr=0;
    Qexp=0;   
else   
    Qcr = trapz(IntegralLimits,(IsotropicCalculatedIntensitiesMatrix(smin3:smax3,1)).*(IntegralLimits.^2)); % Trapezoidal numerical integration
    Qexp = trapz(IntegralLimits,(IsotropicCalculatedIntensitiesMatrix(smin3:smax3,1)+BKG(smin3:smax3,1)).*(IntegralLimits.^2));
end

%----------------------------------------------------------------------------------------------------------------------------------------------------
%    Saving Outputs
%----------------------------------------------------------------------------------------------------------------------------------------------------
if exist('OUTPUTS','dir') == 0   % Creates the folder 'OUTPUTS'
    mkdir('OUTPUTS')
end

if exist(sprintf('OUTPUTS/%s_pattern.csv',filename),'file') == 2  % Checks if a file already exists with the same name in the folder 'OUTPUTS'
   
    FileList=ls('OUTPUTS'); % list the files in the folder 'OUTPUTS'
    cont=1;
    N=0;
    for j=1:1:size(FileList,1)
        number=sscanf(FileList(j,:),sprintf('%s(%%d)',filename));
        if number ~= 0  
            N(cont,1)=number;    
            cont=cont+1;
        end    
    end
    n=max(N)+1; 
    filename = char(sprintf('%s(%d)',filename,n)); % the file name becomes 'filename(n)'   
end


% saves the calculated intensity in a .csv file
OutputCalculated = ExperimentalMatrix; 
OutputCalculated(2:end,2:end) = CalculatedIntensitiesMatrix';
dlmwrite(sprintf('OUTPUTS/%s_pattern.csv',filename),OutputCalculated,';');

% saves the isotropic calculated intensity in a .csv file
OutputIsotropic = [thetaVector BKG IsotropicCalculatedIntensitiesMatrix];
dlmwrite(sprintf('OUTPUTS/%s_isotropic.csv',filename),OutputIsotropic,';')   

% saves the odf in a .csv file
OutputODF = [muVector ODF (ODF.*sind(muVector))];
dlmwrite(sprintf('OUTPUTS/%s_ODF.csv',filename),OutputODF,';')   

% saves the refined parameters in a .csv file
OutputX = char('cagl0','cagl1','cagl2','a','b','c','gamma','L200','LDiag','LDelta','L004','p200','pDiag','pDelta','p004','K','A0','A1','A2','A3','A4','A5','A6','A7','A8','A9','C02','C04','C06','C08','muf','Gammaf','Af','Rwp','Qcr','Qexp') ;
OutputRietveldParameters = [RefinedParameters Rwp Qcr Qexp];
ControlVector = ones(size(OutputRietveldParameters)); % the ControlVector indicates which parameters were kept fixed during the refinement
    for j=1:1:size(Xinput,2)

        if Xinput(1,j)==99        
            ControlVector(1,j)=' ';       
        else
            ControlVector(1,j)='*';
        end    
    end
    ControlVector(1,j+1:1:j+3)=' ';
    fileID = fopen(sprintf('OUTPUTS/%s_output.txt',filename),'w');
    for i=1:1:size(OutputRietveldParameters,2)

        fprintf(fileID,'%c%c%c%c%c%c%c = %d\n' ,ControlVector(1,i),OutputX(i,1:1:6),OutputRietveldParameters(1,i));

    end
    fclose(fileID);

%--------------------------------------------------------------------------------------------------------------------------------------------------------------------
time=toc/60; %stop clock

end