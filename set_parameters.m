function [cagl0,cagl1,cagl2,a,b,c,gamma,L200,LDiag,LDelta,L004,p200,pDiag,pDelta,p004,K,A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,C02,C04,C06,C08,muf,Gammaf,Af] = set_parameters(x0,Xinput)
% SET_PARAMETERS  Set parameters according to Xinput
%   [cagl0,cagl1,cagl2,a,b,c,gamma,L200,LDiag,LDelta,L004,p200,pDiag,pDelta,p004,K,A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,C02,C04,C06,C08,muf,Gammaf,Af] = set_parameters(x0,Xinput)
% Input
%   x0 : row vector containing the initial guesses for the parameters that will be refined
%   Xinput : the input list of 33 model parameters separated by commas. When a parameter x in Xinput is x = 99, the algorithm understands that x is a 
%            free fitting coefficient that is going to be refined. When x ~= 99, the algorithm understands that x is fixed and is going to be kept equal 
%            to the input number. The parameters must be given in the following order;
%            [cagl0,cagl1,cagl2,a,b,c,gamma,L200,LDiag,LDelta,L004,p200,pDiag,pDelta,p004,K,A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,C02,C04,C06,C08,muf,Gammaf,Af]
% Output
%   [cagl0,cagl1,cagl2,a,b,c,gamma,L200,LDiag,LDelta,L004,p200,pDiag,pDelta,p004,K,A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,C02,C04,C06,C08,muf,Gammaf,Af] : 33 model parameters separated by commas
%--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

RefinedParameters=zeros(size(Xinput));

j=1;    
for i=1:size(Xinput,2)
    
    if Xinput(1,i)==99
        
        RefinedParameters(1,i)=x0(1,j);
        j=j+1;
        
    else
        
        RefinedParameters(1,i)=Xinput(1,i);
        
    end
end
        
cagl0=RefinedParameters(1,1);
cagl1=RefinedParameters(1,2);
cagl2=RefinedParameters(1,3);

a=RefinedParameters(1,4);
b=RefinedParameters(1,5);
c=RefinedParameters(1,6);
gamma=RefinedParameters(1,7);

L200=RefinedParameters(1,8);
LDiag=RefinedParameters(1,9);
LDelta=RefinedParameters(1,10);
L004=RefinedParameters(1,11);

p200=RefinedParameters(1,12);
pDiag=RefinedParameters(1,13);
pDelta=RefinedParameters(1,14);
p004=RefinedParameters(1,15);

K=RefinedParameters(1,16);

A0=RefinedParameters(1,17);
A1=RefinedParameters(1,18);
A2=RefinedParameters(1,19);
A3=RefinedParameters(1,20);
A4=RefinedParameters(1,21);
A5=RefinedParameters(1,22);
A6=RefinedParameters(1,23);
A7=RefinedParameters(1,24);
A8=RefinedParameters(1,25);
A9=RefinedParameters(1,26);

C02=RefinedParameters(1,27);
C04=RefinedParameters(1,28);
C06=RefinedParameters(1,29);
C08=RefinedParameters(1,30);

muf=RefinedParameters(1,31);
Gammaf=RefinedParameters(1,32);
Af=RefinedParameters(1,33);

end

