function [x0,lb,ub,NewXinput] = refinement_step(StartingGuess,LowerBound,UpperBound,RefinedParameters,Xinput,RefinementStepMatrix,step)
% REFINEMENT_STEP   Generates the Xinput according to the refinemnet step
%   [x0,lb,ub,NewXinput] = refinement_step(StartingGuess,LowerBound,UpperBound,RefinedParameters,Xinput,step)  
% Input
%   StartingGuess : complete set of reasonable starting guesses.
%   LowerBound : complete set of reasonable lower bounds. 
%   UpperBound : complete set of reasonable upper bounds.
%   RefinedParameters : list with the values of the 33 model parameters separated by commas. The parameters must be given in the following order;
%            [cagl0,cagl1,cagl2,a,b,c,gamma,L200,LDiag,LDelta,L004,p200,pDiag,pDelta,p004,K,A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,C02,C04,C06,C08,muf,Gammaf,Af]
%   Xinput : the input list of 33 model parameters separated by commas. When a parameter x in Xinput is x = 99, the algorithm understands that x is a 
%            free fitting coefficient that is going to be refined. When x ~= 99, the algorithm understands that x is fixed and is going to be kept equal 
%            to the input number. The parameters must be given in the following order;
%            [cagl0,cagl1,cagl2,a,b,c,gamma,L200,LDiag,LDelta,L004,p200,pDiag,pDelta,p004,K,A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,C02,C04,C06,C08,muf,Gammaf,Af]
%   RefinementStepMatrix : Each line of RefinementStepMatrix indicates a refinement step and each column indicates a model parameter. 
%            The number '1' indicates a free fitting coefficient that is going to be refined, the number 0 indicates a fixed coefficient.
%   step : refinement step must be an integer between 1 and 3 
% Output
%   x0 : row vector containing the initial guesses for the parameters that will be refined in the current refinement step
%   lb : row vector containing the lower bounds for the parameters that will be refined in the current refinement step
%   ub : row vector containing the upper bounds for the parameters that will be refined in the current refinement step
%   NewXinput : the new Xinput for the current refinement step
%------------------------------------------------------------------------------------------------------------------
if step==1    
    Guess = [StartingGuess;LowerBound;UpperBound]; 
    else
    Guess = [RefinedParameters;LowerBound;UpperBound];                
end
NewXinput = Xinput;
j=1;
for i=1:size(Xinput,2)
  
    if Xinput(1,i)==99 
        
        if RefinementStepMatrix(step,i)==1
        
            x0(1,j) = Guess(1,i);
            lb(1,j) = Guess(2,i);
            ub(1,j) = Guess(3,i);
            j=j+1;
            else
            NewXinput(1,i) = Guess(1,i);
        end
    end
end
end

                

       
    
    

