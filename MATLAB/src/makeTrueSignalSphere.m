



function TrueSignal = makeTrueSignalSphere(TrueSignalModel, NumSteps, ScaleParameter, varargin)

NumModelInputs = nargin(TrueSignalModel);

NumInputs = length(varargin);

if NumInputs ~= ( NumModelInputs - 1 )
   error('Must be enough inputs for truth model.')
end
    
Input = varargin;



for i = 1 : NumSteps 
    
    
    TrueSignal(:,i) = TrueSignalModel(ScaleParameter(i+1), Input{1,1};
    
end












end