


function SignalOut = performGaussianDiffusion(SignalIn, G, NumSteps)


SignalOut = zeros(length(SignalIn), NumSteps);
SignalOut(:,1) = SignalIn;

for i = 1 : NumSteps - 1
   
    SignalOut(:,i+1) = G*SignalOut(:,i);
    
end
    
    
    
    
    
    
    
    
end


