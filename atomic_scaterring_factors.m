function [AtomicScatteringFactor] = atomic_scaterring_factors(AtomCoefficients,sVector)
% ATOMIC_SCATERRING_FACTORS   Analytical aproximation to the atomic scattering factors
%   [AtomicScatteringFactor] = atomic_scaterring_factors(AtomCoefficients,sVector)
% Input 
%   AtomCoefficients : row vector with the coefficients for analytical approximation to the atomic scattering factors [a1,a2,a3,14,b1,b2,b3,b4,c]
%   sVector : column vector formed by the list of values for s
% Output
%   AtomicScatteringFactor : column vector with the list of atomic scatterring factors for a given sVector  
%---------------------------------------------------------------------------------------
ss = (sVector(:,1)./2).^2;

a1 = AtomCoefficients(1,1);
a2 = AtomCoefficients(1,2);
a3 = AtomCoefficients(1,3);
a4 = AtomCoefficients(1,4);

b1 = AtomCoefficients(1,5);
b2 = AtomCoefficients(1,6);
b3 = AtomCoefficients(1,7);
b4 = AtomCoefficients(1,8);

c = AtomCoefficients(1,9);

      AtomicScatteringFactor = a1.*exp(-b1.*ss(:,1))+...
          a2.*exp(-b2.*ss(:,1))+...
           a3.*exp(-b3.*ss(:,1))+...
            a4.*exp(-b4.*ss(:,1))+c;
end        
        
        