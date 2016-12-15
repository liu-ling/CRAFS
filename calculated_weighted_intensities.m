function [CalculatedWeightedIntensitiesMatrix] = calculated_weighted_intensities(x0,Xinput,sVector,PeakMatrix,StructureFactors,LorentzFactor,etaVector,lambda,ExperimentalIntensityMatrix,GaussianCoefficients)
% CALCULATED_INTENSITIES  Returns the set of calculated intensities
%   [CalculatedWeightedIntensitiesMatrix] = calculated_weighted_intensities(x0,Xinput,sVector,PeakMatrix,StructureFactors,LorentzFactor,etaVector,lambda,ExperimentalIntensityMatrix,GaussianCoefficients)
% Input 
%   x0 : row vector containing the initial guesses for the parameters that will be refined 
%   Xinput : the input list of 33 model parameters separated by commas. When a parameter x in Xinput is x = 99, the algorithm understands that x is a 
%            free fitting coefficient that is going to be refined. When x ~= 99, the algorithm understands that x is fixed and is going to be kept equal 
%            to the input number. The parameters must be given in the following order;
%            [cagl0,cagl1,cagl2,a,b,c,gamma,L200,LDiag,LDelta,L004,p200,pDiag,pDelta,p004,K,A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,C02,C04,C06,C08,muf,Gammaf,Af]
%   sVector : column vector formed by the list of values for s 
%   PeakMatrix : Matrix formed by the Miller indices and multiplicities [h k l m]
%   StructureFactors :  Matrix formed by the structure factors. Each StructureFactor column corresponds to an hkl reflection, and each line to an s value. 
%   LorentzFactor: column vector with the lorentz factor values for each s
%   etaVector: row vector formed by the list of values for eta
%   lambda : X-Ray wavelength 
%   ExperimentalIntensityMatrix : matrix formed by the experimental intensities. Each ExperimentalIntensityMatrix column corresponds to an eta value, and each line to an s value.
%   GaussianCoefficients : matrix formed by the set of coefficients that approximates the Gaussian function 
% Output
%   CalculatedWeightedIntensitiesMatrix : matrix formed by the calculated intensities. Each CalculatedWeightedIntensitiesMatrix column corresponds to an eta value, and each line to an s value.
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
[cagl0,cagl1,cagl2,a,b,c,gamma,L200,LDiag,LDelta,L004,p200,pDiag,pDelta,p004,K,A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,C02,C04,C06,C08,muf,Gammaf,Af] = set_parameters(x0,Xinput); % uses 'x0' and 'Xinput' to set the value of each parameter

[hklRandomIntensitiesMatrix,PeakMatrix] = random_intensities(sVector,PeakMatrix,StructureFactors,LorentzFactor,lambda,cagl0,cagl1,cagl2,a,b,c,gamma,L200,LDiag,LDelta,L004,p200,pDiag,pDelta,p004);
[PreferentialOrientationMatrix] = preferential_orientation(PeakMatrix,etaVector,GaussianCoefficients,C02,C04,C06,C08,muf,Gammaf,Af);

[hklIntensitiesMatrix] = hkl_intensities(etaVector,sVector,hklRandomIntensitiesMatrix,PreferentialOrientationMatrix);

[BKG] = background(sVector,A0,A1,A2,A3,A4,A5,A6,A7,A8,A9);
CalculatedWeightedIntensitiesMatrix = (K.*sum(hklIntensitiesMatrix,3)+repmat(BKG,1,size(etaVector,2)))./(ExperimentalIntensityMatrix.^0.5); 

end

