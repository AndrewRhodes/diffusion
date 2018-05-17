










function SignalOut = performEplotProjection(CPSignalIn, Eplot)

[NumCP, NumSteps] = size(CPSignalIn);

[NumVertices, NumCP2] = size(Eplot);

if NumCP ~= NumCP2
    error('The 1st dimension of ''CPSignal'' must match the 2nd dimension of ''Eplot''.')
end

SignalOut = zeros( NumVertices, NumSteps );

for i = 1 : NumSteps
    
    SignalOut(:,i) = Eplot * CPSignalIn(:,i);
    
end




end