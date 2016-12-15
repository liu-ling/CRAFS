function [PreferentialOrientationMatrix] = preferential_orientation(PeakMatrix,etaVector,GaussianCoefficients,C02,C04,C06,C08,muf,Gammaf,Af)
% PREFERENTIAL_ORIENTATION  Returns the set of preferential orientation coefficients
%   [PreferentialOrientationMatrix] = preferential_orientation(PeakMatrix,etaVector,GaussianCoefficients,C02,C04,C06,C08,muf,Gammaf,Af)
% Input 
%   PeakMatrix : Matrix formed by the Miller indices and multiplicities [h k l m]
%   etaVector: row vector formed by the list of values for eta
%   GaussianCoefficients : matrix formed by the set of coefficients that approximates the Gaussian function 
%   [C02,C04,C06,C08,muf,Gammaf,Af] : model parameters
% Output
%   PreferentialOrientationMatrix : matrix formed by the preferential orientation coefficients. Each PreferentialOrientationMatrix column corresponds to an eta value, and each line to a hkl plane.
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
PreferentialOrientationMatrix=zeros(1,size(etaVector(1,:),2),size(PeakMatrix(:,1),1));
Cf=zeros(41,1);

if Af~=0  
    for i=1:1:41;
        Cf(i,1) = interp2(0:1:45,5:1:30,GaussianCoefficients(:,:,i),muf,Gammaf); % calculates the set of coefficients that approximates the Gaussian function for selected muf and Gammaf  
    end
end

C0=[1;C02;C04;C06;C08;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
Cf(1,1)=0;

for j=1:size(etaVector,2) 

    for i=1:size(PeakMatrix,1)

        Sum=0;
        P1 = legendre_polynomials(80,cos(PeakMatrix(i,9)));
        P2 = legendre_polynomials(80,(cos(PeakMatrix(i,7)/2))*(cosd(etaVector(1,j))));

        for l=1:size(C0,1)

            R=(C0(l,1)+Af*Cf(l,1))*P1(l,1)*P2(l,1);
            Sum=R+Sum;   

        end

        PreferentialOrientationMatrix(1,j,i)=Sum;
        
    end
end



