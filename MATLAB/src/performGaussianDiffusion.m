


function [SignalOut, IterCount] = performGaussianDiffusion(SignalIn, G, NumSteps)


% Scales of interest
S = 2.^[0:NumSteps-1]';

NumVertex = length(SignalIn);

SignalOut = zeros(NumVertex, NumSteps+1);
SignalOut(:,1) = SignalIn;

IterCount = 0;  

for n = 1 : NumSteps
    
    if n == 1

        SignalOut(:,n+1) = G * SignalOut(:,n);
        IterCount = IterCount + 1;   
        
    else      
        
        M = S(n) - S(n-1);
        u = zeros(NumVertex, M+1);
        u(:,1) = SignalOut(:,n);

        for m = 1 : M
            u(:,m+1) = G * u(:,m);
            IterCount = IterCount + 1;            
        end

        SignalOut(:,n+1) = u(:, M+1);
        
    end
    
    
end
    
 
    
    
    
    
    
    
end


