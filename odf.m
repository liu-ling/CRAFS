function [ODF] = odf(C02,C04,C06,C08,mu,muf,Gammaf,Af,GaussianCoefficients)
% ODF   Calculates the Orientation Distribution Function
%   [ODF] = odf(C02,C04,C06,C08,mu,muF,GammaF,Af,GaussianCoefficients)
% Input
%   C02,C04,C06,C08 : expansion coefficients which describe preferential orientation
%   mu : column vector formed by the list of values for mu
%   muF : mean of the Gaussian component
%   GammaF : full-width at half maximum of the Gaussian component
%   Af : multiplying factor of the Gaussian component
%   GaussianCoefficients : matrix formed by the set of coefficients that approximates the Gaussian function
% Output
%   ODF : column vector formed by the list of values for ODF
%------------------------------------------------------------------------------------------------------------
Cf=zeros(41,1);

if Af~=0  
    for i=1:1:41;
        Cf(i,1)=interp2(0:1:45,5:1:30,GaussianCoefficients(:,:,i),muf,Gammaf); % calculates the set of coefficients that approximates the Gaussian function for selected muf and Gammaf  
    end
end

C0=[1;C02;C04;C06;C08;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
Cf(1,1)=0;

ODF = zeros(size(mu));

for i=1:size(mu,1) 
    
    Sum=0;
    P1 = legendre_polynomials(80,cosd(mu(i,1)));

    for l=1:size(C0,1)

        A=(C0(l,1)+Af*Cf(l,1))*P1(l,1);
        Sum=A+Sum;   

    end

    ODF(i,1)=Sum;
            
end
end
