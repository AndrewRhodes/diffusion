% Andrew Rhodes
% ASEL
% March 2018


close all
clear
clc

global ProjectRoot; % Additional Paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MCDivs = [3,4,5,6,7,8];
MCSigFrac = [1,2,3,4,5,6,7,8,9];

MCError = zeros(length(MCSigFrac), 2, length(MCDivs));
MCErrorAll = cell(length(MCSigFrac), length(MCDivs));


for MCd = 1 : length(MCDivs)
    
    for MCsf = 1 : length(MCSigFrac)
        
        clearvars -except MCDivs MCd MCSigFrac MCsf MCError MCErrorAll ProjectRoot
        
        
        % NumberDivisions = 5;
        NumberDivisions = MCDivs(MCd);
        
        NumSigma = 5;
        Model = 'Icosphere';
        
        alpha = 1;
        ShowPlot = 0;
        maxTauNumer = 1;
        sigmaFraction = 10;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Model File Location
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        FileLocationModel = strcat(ProjectRoot,'/models/sphere/');
        FileNameModelPly = strcat(Model,num2str(NumberDivisions),'.ply');
        FileNameModelOff = strcat(Model,num2str(NumberDivisions),'.off');
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make the Sphere and Signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        PointCloud = getIcosphere( fullfile(FileLocationModel, FileNameModelOff), NumberDivisions);
        
        
        PointCloud.LocationCount = size(PointCloud.Location,1);
        PointCloud.FaceCount = size(PointCloud.Face, 1);
        PointCloud.FaceArea = findFaceArea(PointCloud.Location,PointCloud.Face);
        PointCloud = findMeshResolution(PointCloud, 'Model');
        PointCloud.VertexArea = accumarray( reshape(PointCloud.Face, 3*PointCloud.FaceCount ,1), repmat(PointCloud.FaceArea,3,1), [PointCloud.LocationCount, 1] )/3;
        
        % % % % % % % % % %
        Sigma = sigmaFraction * PointCloud.Resolution;
        MaxTau = maxTauNumer / PointCloud.Resolution;
        NumSteps = round(MaxTau);
        % % % % % % % % % %
        
        
        [PointCloud.Theta, PointCloud.Phi, PointCloud.Radius] = cart2sph(PointCloud.Location(:,1) ,PointCloud.Location(:,2), PointCloud.Location(:,3));
        
        TrueSignalModel = @(sigma, Phi) exp(-sigma^2)*sin(Phi);
        
        
        SignalOriginal = TrueSignalModel(0, PointCloud.Phi);
        
        PointCloud.Signal = SignalOriginal;
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup the Explicit Gaussian
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        G = makeExplicitGaussian(PointCloud, Sigma, NumSigma);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Scale Parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        ScaleParameter = findScaleParamter(Sigma, alpha, NumSteps, 'Gaussian', 'Natural');
        
        ScaleParameterAbsolute = bsxfun(@plus, ScaleParameter, PointCloud.Resolution);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Diffusion of Signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Signal = performGaussianDiffusion(PointCloud.Signal, G, NumSteps);
        
        
        TrueSignal = makeTrueSignalSphere(TrueSignalModel, NumSteps, ScaleParameter, PointCloud.Phi);
        
        
        Error = findDiffusionError(TrueSignal, Signal, NumSteps, PointCloud.Phi, ShowPlot);
        
        
        MCError(MCsf, 1:2, MCd) = [NumSteps, Error(NumSteps - 1)];
        MCErrorAll{MCsf, MCd} = [(1:NumSteps)', Error];
        
        
    end
end




