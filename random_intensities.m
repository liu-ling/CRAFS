function [hklRandomIntensitiesMatrix,PeakMatrix] = random_intensities(sVector,PeakMatrix,StructureFactors,LorentzFactor,lambda,cagl0,cagl1,cagl2,a,b,c,gamma,L200,LDiag,LDelta,L004,p200,pDiag,pDelta,p004)
% RANDOM_INTENSITIES  Returns the set of calculated hkl random intensities
%   [hklRandomIntensitiesMatrix,PeakMatrix] = random_intensities(sVector,PeakMatrix,StructureFactors,LorentzFactor,lambda,cagl0,cagl1,cagl2,a,b,c,gamma,L200,LDiag,LDelta,L004,p200,pDiag,pDelta,p004)
% Input 
%   sVector : column vector formed by the list of values for s 
%   PeakMatrix : Matrix formed by the Miller indices and multiplicities [h k l m]
%   StructureFactors :  Matrix formed by the structure factors. Each StructureFactor column corresponds to an hkl reflection, and each line to an s value. 
%   LorentzFactor: column vector with the lorentz factor values for each s
%   lambda : X-Ray wavelength
%   [cagl0,cagl1,cagl2,a,b,c,gamma,L200,LDiag,LDelta,L004,p200,pDiag,pDelta,p004]: model parameters
% Output
%   hklRandomIntensitiesMatrix : matrix formed by the calculated hkl random intensities. Each hklRandomIntensitiesMatrix column corresponds to an eta value, and each line to an s value.
%   PeakMatrix : [h k l m d sin(theta) 2theta s nu xi L p IntegralBreadthIns IntegralBreadthSpec LorentzIntegralBreadthSpec GaussianIntegralBreadthSpec GaussianIntegralBreadth VoigtIntegralBreadth ObsIntegralBreadth Gamma]
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

%PeakMatrix(:,1)=h
%PeakMatrix(:,2)=k
%PeakMatrix(:,3)=l
%PeakMatrix(:,4)=multiplicity

PeakMatrix(:,5)=((1./sind(gamma)).*sqrt((PeakMatrix(:,1)./a).^2+(PeakMatrix(:,2)./b).^2-...  
        (2.*PeakMatrix(:,1).*PeakMatrix(:,2).*cosd(gamma))/(a.*b)+((PeakMatrix(:,3)./c).^2).*...
            (sind(gamma)).^2)).^-1; % stores the values of d-spacing in the column 5 of the PeakMatrix  
PeakMatrix(:,6)=(lambda)./(2.*PeakMatrix(:,5)); % stores the values of Sintheta in the column 6 of the PeakMatrix
PeakMatrix(:,7)=2.*asin(PeakMatrix(:,6)); % stores the values of 2theta in the column 7 of the PeakMatrix
PeakMatrix(:,8)=2.*PeakMatrix(:,6)./lambda; % stores the values of s(hkl) in the column 8 of the PeakMatrix
PeakMatrix(:,9)=acos((PeakMatrix(:,3).*PeakMatrix(:,5))./c);% stores the values of angle nu in column 9 of the PeakMatrix

%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   Crystallite Shape
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------

for i=1:1:size(PeakMatrix,1)               

if PeakMatrix(i,1)==0 && PeakMatrix(i,2)==0  % PeakMatrix(:,1)=h and PeakMatrix(:,2)=k

    PeakMatrix(i,10) = 1; % stores the values of xi(00l) in the column 10 of the PeakMatrix  

    PeakMatrix(i,11) = L004; % stores the values of L(00l) in the column 11 of the PeakMatrix 
    PeakMatrix(i,12) = p004; % stores the values of p(00l) in the column 12 of the PeakMatrix 

else

    PeakMatrix(i,10) = (cos(2*acos(PeakMatrix(i,1)./((PeakMatrix(i,1).^2+PeakMatrix(i,2).^2).^0.5)))).^2; % stores the values of xi(hkl) in the column 10 of the PeakMatrix

    if PeakMatrix(i,2)==0 % PeakMatrix(:,2)=k

        PeakMatrix(i,11) = (1/sin(PeakMatrix(i,9)))*(PeakMatrix(i,10)*L200); % stores the values of L(h0l) in the column 11 of the PeakMatrix
        PeakMatrix(i,12) = PeakMatrix(i,10)*p200; % stores the values of p(h0l) in the column 12 of the PeakMatrix
    else
        PeakMatrix(i,11) = (1/sin(PeakMatrix(i,9)))*(PeakMatrix(i,10)*L200+(1-PeakMatrix(i,10))*(LDiag+0.5*(PeakMatrix(i,2)/abs(PeakMatrix(i,2)))*LDelta)); % stores the values of L(hkl) in the column 11 of the PeakMatrix
        PeakMatrix(i,12) = PeakMatrix(i,10)*p200+(1-PeakMatrix(i,10))*(pDiag+0.5*(PeakMatrix(i,2)/abs(PeakMatrix(i,2)))*pDelta); % stores the values of p(hkl) in the column 12 of the PeakMatrix
        
        if PeakMatrix(i,12)>1 
            PeakMatrix(i,12) = 1;
        elseif PeakMatrix(i,12)<0 
            PeakMatrix(i,12) = 0;
        end    
    end
    
end
end

%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   Peak Broadening
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
PeakMatrix(:,13) = (((cos(PeakMatrix(:,7)./2)).^2).*(cagl0+cagl1.*tan(PeakMatrix(:,7)./2)+cagl2.*(tan(PeakMatrix(:,7)./2).^2))).^0.5; % stores the values of the IntegralBreadthIns in the column 13 of the PeakMatrix
PeakMatrix(:,14) = (1+((pi*0.07)^2).*(PeakMatrix(:,1).^2+PeakMatrix(:,2).^2+(PeakMatrix(:,3).^2)./2))./PeakMatrix(:,11);% stores the values of the IntegralBreadthSpec in the column 14 of the PeakMatrix
PeakMatrix(:,15) = PeakMatrix(:,14).*((PeakMatrix(:,12)+(1-PeakMatrix(:,12)).*(pi*log(2))^0.5)); % stores the values of the LorentzIntegralBreadthSpec in the column 15 of the PeakMatrix
PeakMatrix(:,16) = (PeakMatrix(:,14)./(pi*log(2))^0.5).*((PeakMatrix(:,12)+(1-PeakMatrix(:,12)).*(pi*log(2))^0.5));% stores the values of the GaussianIntegralBreadthSpec in the column 16 of the PeakMatrix
PeakMatrix(:,17) = (PeakMatrix(:,13).^2+PeakMatrix(:,16).^2).^0.5; % stores the values of the GaussianIntegralBreadth in the column 17 of the PeakMatrix 
    k(:,1) = PeakMatrix(:,15)./(PeakMatrix(:,13).*(pi^0.5));
    for j=1:size(PeakMatrix,1)
        if k(j,1)>27
            k(j,1)=27.1;
        end
    end
PeakMatrix(:,18) = PeakMatrix(:,13).*(exp(-k(:,1).^2))./erfc(k(:,1)); % stores the values of the VoigtIntegralBreadth in the column 18 of the PeakMatrix
PeakMatrix(:,19) = ((PeakMatrix(:,18).*PeakMatrix(:,17))./((PeakMatrix(:,12).*PeakMatrix(:,17))+(1-PeakMatrix(:,12)).*PeakMatrix(:,18))); % stores the values of the ObsIntegralBreadth in the column 19 of the PeakMatrix
PeakMatrix(:,20) = (2.*PeakMatrix(:,19)./pi).*((PeakMatrix(:,12)+(1-PeakMatrix(:,12)).*(pi*log(2))^0.5)); % stores the values of Gamma(hkl) in the column 20 of the PeakMatrix
    
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
%   Peak-Shape
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
PseudoVoigt = zeros(size(sVector(:,1),1),size(PeakMatrix(:,1),1));
hklRandomIntensitiesMatrix = zeros(size(sVector(:,1),1),1,size(PeakMatrix(:,1),1));

for i=1:size(PeakMatrix,1)

    Lorentzian(:,1) = ((2./(pi.*PeakMatrix(i,20)))./(1+(4.*(sVector(:,1)-PeakMatrix(i,8)).^2)./(PeakMatrix(i,20).^2)));
    Gaussian(:,1) = ((((4*log(2)).^0.5)./((pi^0.5).*PeakMatrix(i,20))).*(exp((-4.*log(2).*(sVector(:,1)-PeakMatrix(i,8)).^2)./(PeakMatrix(i,20).^2))));
    PseudoVoigt(:,i) = ((PeakMatrix(i,12).*Lorentzian(:,1)+(1-PeakMatrix(i,12)).*Gaussian(:,1)));
    % Calculates the hkl random intensities
    hklRandomIntensitiesMatrix(:,1,i) = PeakMatrix(i,4).*LorentzFactor(:,1).*StructureFactors(:,i).*PseudoVoigt(:,i);

end

end

