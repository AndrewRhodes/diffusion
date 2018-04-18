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

% May need to change dimension:2D,3D to type:Gaussian,Laplacian


function ScaleParameter = findScaleParamter(StepSize, Alpha, NumLevels, Method, Dimension)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error Testing

if ~isnumeric(StepSize) && (StepSize <= 0) && (length(StepSize) ~= 1)
    error('Expected StepSize > 0, but StepSize = %.3f', StepSize)
end

if ~isnumeric(Alpha) && (Alpha <= 0) && (length(Alpha) ~= 1)
    error('Expected Alpha > 0, but Alpha = %.3f', Alpha)
end

if ~isnumeric(NumLevels) && (NumLevels <= 0) && (length(NumLevels) ~= 1)
    error('Expected NumLevels > 0, but NumLevels = %.3f', NumLevels)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Scale Parameter

ScaleParameter = zeros(NumLevels,1);

if strcmpi( Method, 'Natural' )
    
    if strcmpi( Dimension, '2D' )
        
        ScaleParameter = sqrt((0:NumLevels)' * StepSize^2);
        
    elseif strcmpi( Dimension, '3D' )
        
        ScaleParameter = sqrt((0:NumLevels)' * 2 * Alpha * StepSize);
        
    else
        
        error('Unrecognized Method ''%s''. Must be ''Natural'' or ''Cutoff''.', Method)
        
    end
    
elseif strcmpi( Method, 'Cutoff' )
    
    ws = 0 : 0.01 : 5;
    NumSample = length(ws);
    cut = sqrt(log(2));
    db3 = 1/sqrt(2);
    ws2 = (ws.^2)';
    H = ones(NumSample,1) ;
    h = 1 ./ ( ones(NumSample,1) + Alpha * StepSize * ws2 );
    
    
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
    
    error('Unrecognized Dimension ''%s''. Must be ''3D'' or ''2D''.', Dimension)
    
end


end















