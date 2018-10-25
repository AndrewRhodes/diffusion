




function [DoG, AbsDoG] = buildDoGNew(SignalIn, ScaleParameterIn)

[NumVertices, NumLevels] = size(SignalIn);

NumScales = length(ScaleParameterIn);

if NumLevels ~= NumScales
    error('Expected 2nd dimension of SignalIn to equal length of ScaleParameterIn.')
end


DoG = zeros(NumVertices, NumLevels-1);
AbsDoG = zeros(NumVertices, NumLevels-1);

for i = 1 : NumLevels - 1
%         DoG(:,i) = (SignalIn(:,i+1) - SignalIn(:,i));
%         DoG(:,i) = (SignalIn(:,i+1) - SignalIn(:,i)) / (ScaleParameterIn(i+1)/ScaleParameterIn(i) - 1) ;
        DoG(:,i) = (SignalIn(:,i+1) - SignalIn(:,i)) * (ScaleParameterIn(i,1)^2 / ( ScaleParameterIn(i+1,1)^2 - ScaleParameterIn(i,1)^2 ) ) ;
        
%         AbsDoG(:,i) = abs(SignalIn(:,i+1) - SignalIn(:,i)) ;
        AbsDoG(:,i) = abs(SignalIn(:,i+1) - SignalIn(:,i)) * (ScaleParameterIn(i,1)^2 / ( ScaleParameterIn(i+1,1)^2 - ScaleParameterIn(i,1)^2 ) ) ;
        
end







end







    