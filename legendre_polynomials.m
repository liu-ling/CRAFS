function [LegendrePolynomials] = legendre_polynomials(n,x)
% LEGENDRE_POLYNOMIALS   Generates the Legendre Polynomials (not normalized, only even terms)
%   [LegendrePolynomials] = legendre_polynomials(n,x)
% Input
%   n : maximum degree of the legendre polynomials
%   x : argument of the Legendre Polynomials
% Output
%   LegendrePolynomials : column vector with the set of even terms of the not normalized Legendre polynomials of x
%--------------------------------------------------------------------------
P=zeros(n+1,1);

P(1,1)=1;
P(2,1)=x;

for j=2:n    
    i=j+1;
    P(i,1) = (((2*j-1)*x*P(i-1,1)-(j-1)*P(i-2,1))/j); % recurrence relation   
end
LegendrePolynomials=P(1:2:n+1,1);
end

