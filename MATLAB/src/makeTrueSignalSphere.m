



function TrueSignal = makeTrueSignalSphere(TrueSignalModel, NumSteps, ScaleParameter, Input)

% NumModelInputs = nargin(TrueSignalModel);

% NumInputs = length(varargin);

% if NumInputs ~= ( NumModelInputs - 1 )
%    error('Must be enough inputs for truth model.')
% end
    
% Input = varargin;


TrueSignal = zeros(length(Input), NumSteps);
% TrueSignal(:,1) = TrueSignalModel(0, Input);

for i = 1 : NumSteps 
    
    
    TrueSignal(:,i) = TrueSignalModel(ScaleParameter(i), Input);
    
end












end