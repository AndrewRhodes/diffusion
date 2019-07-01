% Andrew Rhodes
% WVU
% Jan. 2019
%
% input[]:
% input[]:
% input[]:
% input[]:
% input[]:
%
% output[]:

function [tn, tau_n, sigma_n] = findScaleStep(k, t_0, alpha, N)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error Testing

if ~isnumeric(t_0) && (t_0 <= 0) && (length(t_0) ~= 1)
    error('Expected StepSize > 0, but StepSize = %.3f', StepSize)
end

if ~isnumeric(alpha) && (alpha <= 0) && (length(alpha) ~= 1)
    error('Expected Alpha > 0, but Alpha = %.3f', alpha)
end

if ~isnumeric(N) && (N <= 0) && (length(N) ~= 1)
    error('Expected NumLevels > 0, but NumLevels = %.3f', NumSteps)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Scale Parameter


tn = k .^ (2*(0:N)') * t_0;
tau_n = diff(tn);

sigma_n = sqrt(2 * alpha * tn);



end