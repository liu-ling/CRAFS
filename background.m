function [BKG] = background(sVector,A0,A1,A2,A3,A4,A5,A6,A7,A8,A9)
% BACKGROUND   Calculates the Background using Chebyshev Polynomials
%   [BKG] = background(sVector,BackgroundCoefficients)
% Input
%   sVector : column vector formed by the list of values for s 
%   [A0,A1,A2,A3,A4,A5,A6,A7,A8,A9] : BackgroundCoefficients 
% Output
%   BKG : column vector with the background intensities for a given sVector and background coefficients
%------------------------------------------------------------------------------------------------------------------
   
x(:,1) = (2.*(sVector(:,1)-min(sVector(:,1)))./(max(sVector(:,1))-min(sVector(:,1))))-1;

t0 = 1;
t1 = x(:,1);
t2 = 2.*x(:,1).^2-1;
t3 = 4.*x(:,1).^3-3.*x(:,1);
t4 = 8.*x(:,1).^4-8.*x(:,1).^2+1;
t5 = 16.*x(:,1).^5+-20.*x(:,1).^3+5.*x(:,1);
t6 = 32.*x(:,1).^6-48.*x(:,1).^4+18.*x(:,1).^2-1;
t7 = 64.*x(:,1).^7-112.*x(:,1).^5+56.*x(:,1).^3-7.*x(:,1);
t8 = 128.*x(:,1).^8-256.*x(:,1).^6+160.*x(:,1).^4-32.*x(:,1).^2+1;
t9 = 256.*x(:,1).^9-576.*x(:,1).^7+432.*x(:,1).^5-120.*x(:,1).^3+9.*x(:,1);

BKG=A0.*t0+A1.*t1+A2.*t2+A3.*t3+A4.*t4+A5.*t5+A6.*t6+A7.*t7+A8.*t8+A9.*t9;

end

