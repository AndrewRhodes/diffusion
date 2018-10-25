




function [ItL, LBM] = makeCotangentLaplaceBeltrami( FileName, OrderBDF, StepSize, Alpha)


if Alpha < 0
    error('Expected Alpha >= 0, but input is ''%0.3f''', Alpha)
end

if ~strcmp( FileName(end-2:end), 'off')
    error('Expected Filename to end in ''.off'', but was ''%s''', Filename(end-3:end))
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % Cotangent Laplacian
[II, JJ, SS, AA] = cotlpmatrix(FileName);
W=sparse(II, JJ, SS);
A=AA;
Alength = length(A);
A1 = sparse(1:Alength, 1:Alength, 1./A);
LBM = A1 * W;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ItL = cell(OrderBDF,1);

ItL{1,1} = speye(Alength, Alength) - Alpha * StepSize * LBM;


if OrderBDF == 2 || OrderBDF == 3 || OrderBDF == 4
%     I23tL
    ItL{2,1} = speye(Alength, Alength) - (2/3) * Alpha * StepSize * LBM;
end

if OrderBDF == 3 || OrderBDF == 4
%     I611tL
    ItL{3,1} = speye(Alength, Alength) - (6/11) * Alpha * StepSize * LBM;
end

if OrderBDF == 4
%     I1225tL
    ItL{4,1} = speye(Alength, Alength) - (12/25) * Alpha * StepSize * LBM;
end




end











