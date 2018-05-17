




function Error = findDiffusionError(TrueSignal, DiffusedSignal, NumSteps, ModelInput, ShowPlot)

[m1Row, m1Col] = size(TrueSignal);

[m2Row, m2Col] = size(DiffusedSignal);

if m2Row~=m1Row || m1Col~=m2Col
    error('True and Diffused Signal Matrice must be same size.')
end

Error = zeros(NumSteps,1);

if ShowPlot
    figure(1)
end
for i = 1 : NumSteps - 1
    
    Error(i+1,1) = norm(TrueSignal(:,i+1) - DiffusedSignal(:,i+1), inf);
    
    
    if ShowPlot
        clf
        plot(ModelInput, TrueSignal(:,1), 'ko')
        hold on
        plot(ModelInput, TrueSignal(:,i), 'gd')
        plot(ModelInput, DiffusedSignal(:,i),'r.')
        legend('Original', 'Truth', 'Numeric')
        drawnow
    end           
    
    
end
if ShowPlot
    close(figure(1))
end


end