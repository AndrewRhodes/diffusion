% Andrew Rhodes
% ASEL
% March 2018


close all
clear
clc

addpath('../src/')
ProjectRoot = setupprojectpaths % Additional Paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

alpha = 1;
NumberDivisions = 3; % For building Icosphere

% Spherical Harmonics are for defined for l = 1 ... inf, but the 
% coefficients on higher terms quickly go to zero. 
MaxDegreeL = 50;

options.rho = 3;
options.dtype = 'geodesic';
% options.dtype = 'euclidean';

FileLocation = '../models/Sphere/';
FileName = strcat('Icosphere',num2str(NumberDivisions),'.off');

ShowPlot = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the Sphere and Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fullfile(FileLocation, FileName), 'file')
    
    [Location, Faces] = icosphere(NumberDivisions);
    [VerticesOut, FacesOut] = clearMeshDuplicates(Location, Faces );
    
    save_off(VerticesOut, FacesOut, fullfile(FileLocation, FileName))
    
    Sphere.Face = FacesOut;
    Sphere.Location = VerticesOut;
    
else
    
    [Sphere.Location, Sphere.Face] = read_off( fullfile(FileLocation, FileName) );

    [m, n] = size(Sphere.Location);
    if m < n
        Sphere.Location = Sphere.Location';
    end
    
    [m, n] = size(Sphere.Face);
    if m < n
        Sphere.Face = Sphere.Face';
    end
    
end


Sphere.LocationCount = size(Sphere.Location,1);
Sphere.FaceCount = size(Sphere.Face, 1);
Sphere.FaceArea = findFaceArea(Sphere.Location,Sphere.Face);
Sphere = findMeshResolution(Sphere, 'Model');


% % % % % % % % % % 
tauExplicit = Sphere.Resolution / 8;
MaxTauExplicit = 10 / Sphere.Resolution;
NumStepsExplicit = round(MaxTauExplicit);
% % % % % % % % % % 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Sphere.Theta, Sphere.Phi, Sphere.Radius] = cart2sph(Sphere.Location(:,1) ,Sphere.Location(:,2), Sphere.Location(:,3));

SphericalHarmonic = makeRealSphericalHarmonic( MaxDegreeL, Sphere.Theta, Sphere.Phi );

% Define the signal

ExactSignal = @(sigma, SignalOriginal, MaxDegreeL) sum(cell2mat(cellfun(@times, num2cell( (exp(-(sigma^2/2)*(1:MaxDegreeL).*((1:MaxDegreeL)+1))')), SignalOriginal, 'UniformOutput', 0)),1);

% % ExactSignal = @(sigma, SignalOriginal, MaxDegreeL) sum(cell2mat(cellfun(@times, cellfun(@times, num2cell( (exp(-(1:MaxDegreeL).*((1:MaxDegreeL)+1)*sigma^2/2)')), num2cell( (exp(-(1:MaxDegreeL).^2/9))' ), 'UniformOutput', 0), SignalOriginal, 'UniformOutput', 0)),1);

SignalOriginal = ExactSignal(0, SphericalHarmonic, MaxDegreeL);
Sphere.Signal = SignalOriginal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the Laplace Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FileLocationMeshLP = strcat(ProjectRoot,'/models/Sphere/meshLP/');
FileNameLapMat = strcat('LapMatMeshWeights','_Div',num2str(NumberDivisions),'_NumSigma',num2str(options.rho),'_',options.dtype,'.mat');
FileNameArea = strcat('Area','_Div',num2str(NumberDivisions),'_NumSigma',num2str(options.rho),'_',options.dtype,'.mat');
FileNamehEdge = strcat('hEdge2','_Div',num2str(NumberDivisions),'_NumSigma',num2str(options.rho),'_',options.dtype,'.mat');

if ~exist(fullfile(FileLocationMeshLP, FileNameLapMat), 'file')
    
    [LapMatMeshWeights, Area, hEdge2] = symmshlp_matrix(fullfile(FileLocation, FileName), options);
    
    save(fullfile(FileLocationMeshLP, FileNameLapMat), 'LapMatMeshWeights')
    save(fullfile(FileLocationMeshLP, FileNameArea), 'Area')
    save(fullfile(FileLocationMeshLP, FileNamehEdge), 'hEdge2')
else
    
    load(fullfile(FileLocationMeshLP, FileNameLapMat))
    load( fullfile(FileLocationMeshLP, FileNameArea) )
    load( fullfile(FileLocationMeshLP, FileNamehEdge) )
    
end

hEdge = (hEdge2/2);

Alength = length(Area);

A1 = sparse(1:Alength, 1:Alength, 1./Area);

LBM = A1 * LapMatMeshWeights;

ItL = speye(Alength, Alength) - alpha*tauExplicit * LBM;

I23tL = speye(Alength, Alength) - (2/3)*alpha*tauExplicit * LBM;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale Parameter Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ScaleParameter = findScaleParamter(tauExplicit, alpha, NumStepsExplicit, 1, 3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SignalExplicit = zeros(Sphere.LocationCount, NumStepsExplicit);
SignalExplicit(:,1) = Sphere.Signal;
Truth = zeros(Sphere.LocationCount, NumStepsExplicit);
Truth(:,1) = Sphere.Signal;

AbsErr = zeros(NumStepsExplicit,1);


WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumStepsExplicit-1));
if ShowPlot
    figure
end

for i = 1 : NumStepsExplicit - 1
    
% %     if i ==1
% %     SignalExplicit(:,i+1) = ItL \ SignalExplicit(:,i);
        [SignalExplicit(:,i+1), flag] = bicg(ItL, SignalExplicit(:,i), 1e-10, 60);
% %     [SignalExplicit(:,i+1), flag] = gmres(ItL, SignalExplicit(:,i), [], 1e-10, 100);
% %     else
% %         [SignalExplicit(:,i+1), flag] = gmres(I23tL, (4/3)*SignalExplicit(:,i) - (1/3)*SignalExplicit(:,i-1), [], 1e-10, 100);
% %     end
    
    if flag
        flag
    end
    
    Truth(:,i+1) = ExactSignal(ScaleParameter(i+1), SphericalHarmonic, MaxDegreeL);
	AbsErr(i+1,1) = norm(Truth(:,i+1) - SignalExplicit(:,i+1), inf);
    
    if ShowPlot
        clf
        plot(Sphere.Theta, SignalOriginal,'ko')
        hold on
        plot(Sphere.Theta, Truth(:,i), 'gd')
        plot(Sphere.Theta, SignalExplicit(:,i),'r.')
    end
    
    waitbar(i/NumStepsExplicit, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumStepsExplicit-1));
    
end

waitbar(i/NumStepsExplicit, WaitBar, sprintf('Diffusion Complete'));
close(WaitBar)
if ShowPlot
    close(figure)
end


% MCError(MCs, 1:2, MCp) = [NumStepsImplicit, AbsErr(NumStepsImplicit - 1)];
% MCErrorAll{MCs, MCp} = [(1:NumStepsImplicit)', AbsErr];









