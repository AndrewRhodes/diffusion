    




function DoG = buildDoG(SignalIn, ScaleParameterIn, DoGNormalize)

[NumVertices, NumLevels] = size(SignalIn);

NumScales = length(ScaleParameterIn);

if NumLevels ~= NumScales
    error('Expected 2nd dimension of SignalIn to equal length of ScaleParameterIn.')
end


DoG = zeros(NumVertices, NumLevels-1);


for i = 1 : NumLevels - 1
    
    if strcmpi(DoGNormalize, 'DoG')
        
        DoG(:,i) = SignalIn(:,i+1) - SignalIn(:,i) ;
        
    elseif strcmpi(DoGNormalize, 'AbsDoG')
        
        DoG(:,i) = abs(SignalIn(:,i+1) - SignalIn(:,i)) ;
        
    elseif strcmpi(DoGNormalize, 'NLoG')
%                        
%         DoG(:,i) = ( SignalIn(:,i+1) - SignalIn(:,i) ) ...
%                 * (ScaleParameterIn(i+1,1)^2 ...
%                 / ( ScaleParameterIn(i+1,1)^2 - (i/(i+1))*ScaleParameterIn(i,1)^2 ) );
            

        DoG(:,i) = ( SignalIn(:,i+1) - SignalIn(:,i) ) ...
            * (ScaleParameterIn(i+1,1)^2 ...
            / ( ScaleParameterIn(i+1,1)^2 - ScaleParameterIn(i,1)^2 ) );


%         DoG(:,i) = (i+1)*( SignalIn(:,i+1) - SignalIn(:,i) ) / ( ScaleParameterIn(i+1,1)^2 - ScaleParameterIn(i,1)^2 - 1 );

%         DoG(:,i) = ( SignalIn(:,i+1) - SignalIn(:,i) );  

        
    elseif strcmpi(DoGNormalize, 'AbsNLoG')
        
        DoG(:,i) = abs( SignalIn(:,i+1) - SignalIn(:,i) ) ...
                * (ScaleParameterIn(i+1,1)^2 ...
                / ( ScaleParameterIn(i+1,1)^2 - ((i+1)/i)*ScaleParameterIn(i,1)^2 ) );
            
    else
        
        error('buildDoG::Normalization Structure not recognized.')
    end
    
    
    % Scale Invariant Laplace
    %     Fbar = (sum(SurfaceCurvatures(:,1,i)) / PointCloud.LocationCount) * ones(PointCloud.LocationCount, 1) ;
    %
    %     SigmaL = sqrt(sum((LaplacianScale(:,1,i)-Fbar).^2)) / sqrt(PointCloud.LocationCount);
    %     %
    %     LaplaceScaleInvariant(:,1,i) = ( LaplacianScale(:,1,i) - Fbar ) / SigmaL;
    
end







end






    