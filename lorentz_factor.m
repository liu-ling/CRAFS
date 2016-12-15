function [LorentzFactor] = lorentz_factor(sVector,lambda)
% LORENTZ_FACTOR  Calculates the Lorentz Factor
%   [LorentzFactor] = lorentz_factor(sVector,lambda)
% Input
%   sVector : column vector formed by the list of values for s 
%   lambda : X-Ray wavelength
% Output
%   LorentzFactor : column vector with the lorentz factor values for each s
%--------------------------------------------------------------------------

theta = asin((sVector(:,1).*lambda)./2);
LorentzFactor(:,1) = 1./((sin(theta(:,1))).^2);

end

