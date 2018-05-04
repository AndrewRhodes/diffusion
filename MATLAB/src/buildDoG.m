




function DoG = buildDoG(SignalIn, ScaleParameterIn, Normalize)

[NumVertices, NumLevels] = size(SignalIn);

NumScales = length(ScaleParameterIn);

if NumLevels ~= NumScales
    error('Expected 2nd dimension of SignalIn to equal length of ScaleParameterIn.')
end


DoG = zeros(NumVertices, NumLevels-1);    


for i = 1 : NumLevels - 1
    
    if Normalize
        DoG(:,i) = 2 * ( SignalIn(:,i+1) - SignalIn(:,i) ) * (ScaleParameterIn(i,1)^2 / ( ScaleParameterIn(i+1,1)^2 - ScaleParameterIn(i,1)^2 ) );
    else
        DoG(:,i) = ( SignalIn(:,i+1) - SignalIn(:,i) ) / ( ScaleParameterIn(i+1,1)^2 - ScaleParameterIn(i,1)^2 );
    end
    

    % Scale Invariant Laplace
%     Fbar = (sum(SurfaceCurvatures(:,1,i)) / PointCloud.LocationCount) * ones(PointCloud.LocationCount, 1) ;
%     
%     SigmaL = sqrt(sum((LaplacianScale(:,1,i)-Fbar).^2)) / sqrt(PointCloud.LocationCount);
%     %  
%     LaplaceScaleInvariant(:,1,i) = ( LaplacianScale(:,1,i) - Fbar ) / SigmaL;
    
end







end