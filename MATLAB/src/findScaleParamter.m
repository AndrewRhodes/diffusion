% Andrew Rhodes
% WVU
% Jan. 2018
%
% input[Step]:
% input[NumLevels]:
% input[Method]: Either 1 (Natural) or 2 (Cutoff Frequency)
% input[Dimension]: Either 2 or 3
%
% output[ScaleParameter]:




function ScaleParameter = findScaleParamter(Step, NumLevels, Method, Dimension)



ScaleParameter = zeros(NumLevels,1);

if Method == 1 % Natural Scale Growth
    
    if Dimension == 2
        
    elseif Dimension == 3
        
        for i = 1 : NumLevels - 1
            
            ScaleParameter(i+1,1) = sqrt(2 * i * Step);
            
        end
    else
        error('Dimension must be 1 or 2')
    end
    
elseif Method == 2 % Cutoff Frequency Method
    
    maxsample = 5;
    ws = 0 : 0.01 : maxsample;
    NumSample = length(ws);
    cut = sqrt(log(2));
    db3 = 1/sqrt(2);
    ws2 = (ws.^2)';
    H = ones(NumSample,1) ;
    h = 1 ./ ( ones(NumSample,1) + Step * ws2 );
    
    
    for i = 1 : NumLevels -1
        
        % Transfer function
        H = H .* h;
        % Find the frequency at the cutoff values
        [uH, indH] = unique(H);
        CutoffFrequencyImplicit = interp1(uH, ws(indH), db3);
        % Change cutoff frequency to scale parameter
        ScaleParameterFrequencyImplicit = CutoffFrequencyImplicit / cut;
        ScaleParameter(i+1,1) =  1 / ScaleParameterFrequencyImplicit;
        
    end
    
else
    error('Method must be 1 or 2')
end


end















