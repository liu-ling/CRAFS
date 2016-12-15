function [StructureFactors] = structure_factors(sVector,PeakMatrix,HydrogenCoefficients,CarbonCoefficients,OxygenCoefficients,HydrogenCoordinates,CarbonCoordinates,OxygenCoordinates)
% STRUCTURE_FACTORS   Calculates the Structure Factors
%   [StructureFactors] = structure_factors(sVector,PeakMatrix,HydrogenCoefficients,CarbonCoefficients,OxygenCoefficients,HydrogenCoordinates,CarbonCoordinates,OxygenCoordinates)
% Input
%   sVector : column vector formed by the list of values for s
%   PeakMatrix : Matrix formed by the Miller indices and multiplicities [h k l m]
%   HydrogenCoefficients : row vector with the coefficients for analytical approximation to the hydrogen scattering factors [a1,a2,a3,14,b1,b2,b3,b4,c1] 
%   CarbonCoefficients : row vector with the coefficients for analytical approximation to the carbon scattering factors [a1,a2,a3,14,b1,b2,b3,b4,c1]
%   OxygenCoefficients : row vector with the coefficients for analytical approximation to the oxygen scattering factors [a1,a2,a3,14,b1,b2,b3,b4,c1]
%   HydrogenCoordinates : matrix formed by the fractional atomic coordinates for hydrogen [x y z] 
%   CarbonCoordinates : matrix formed by the fractional atomic coordinates for carbon [x y z]
%   OxygenCoordinates : matrix formed by the fractional atomic coordinates for oxygen [x y z]
% Output
%   StructureFactors : Matrix formed by the structure factors, each StructureFactor column corresponds to an hkl reflection, and each line to an s value.
%------------------------------------------------------------------------------------------------------------

StructureFactorsHydrogen=atomic_scaterring_factors(HydrogenCoefficients,sVector); % calculates the atomic scattering factors for Hydrogen
StructureFactorsCarbon=atomic_scaterring_factors(CarbonCoefficients,sVector); % calculates the atomic scattering factors for Carbon
StructureFactorsOxygen=atomic_scaterring_factors(OxygenCoefficients,sVector); % calculates the atomic scattering factors for Oxygen

StructureFactors=zeros(size(sVector,1),size(PeakMatrix,1));

    for n=1:1:size(PeakMatrix,1) % routine to calculate the structure factors

         Sum=0;

         h=PeakMatrix(n,1);
         k=PeakMatrix(n,2);
         l=PeakMatrix(n,3);
        
        
            for j=1:1:size(HydrogenCoordinates,1)
                      
                S=StructureFactorsHydrogen.*exp(2.*pi.*1i.*(h.*HydrogenCoordinates(j,1)+k.*HydrogenCoordinates(j,2)+l.*HydrogenCoordinates(j,3)));
                Sum = S + Sum;
            end
        
           
            for j=1:1:size(CarbonCoordinates,1)
                       
                S=StructureFactorsCarbon.*exp(2.*pi.*1i.*(h.*CarbonCoordinates(j,1)+k.*CarbonCoordinates(j,2)+l.*CarbonCoordinates(j,3)));
                Sum = S + Sum;

            end
        
            for j=1:1:size(OxygenCoordinates,1)
            
                S=StructureFactorsOxygen.*exp(2.*pi.*1i.*(h.*OxygenCoordinates(j,1)+k.*OxygenCoordinates(j,2)+l.*OxygenCoordinates(j,3)));
                Sum = S + Sum;
            end
        
  StructureFactors(:,n)=Sum;
  StructureFactors(:,n)=StructureFactors(:,n).*conj(StructureFactors(:,n));
  
    end

end

                   
                   
                   
                   
                   
                   

