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
NumberDivisions = 4; % For building Icosphere

options.rho = 6;
options.dtype = 'geodesic';
% options.dtype = 'euclidean';

FileLocation = '../models/sphere/';
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

% NumStepsExplicit = 1000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Sphere.Theta, Sphere.Phi, Sphere.Radius] = cart2sph(Sphere.Location(:,1) ,Sphere.Location(:,2), Sphere.Location(:,3));


% SphericalHarmonic = makeRealSphericalHarmonic( MaxDegreeL, Sphere.Theta, Sphere.Phi );
% SphericalHarmonic = sum(cell2mat(SphericalHarmonic),1);
% Define the signal

% ExactSignal = @(sigma, SignalOriginal, MaxDegreeL) sum(cell2mat( cellfun(@times, num2cell( (exp(-(sigma^2/2)*(1:MaxDegreeL).*((1:MaxDegreeL)+1))')), SignalOriginal, 'UniformOutput', 0)),1);

% ExactSignal = @(sigma, SignalOriginal, MaxDegreeL) sum(cell2mat(cellfun(@times, num2cell( sqrt( (2*(1:MaxDegreeL)'+1) / (4*pi) ) ),  cellfun(@times, num2cell( exp(-(sigma^2/2)*(1:MaxDegreeL).*((1:MaxDegreeL)+1))' ), SignalOriginal, 'UniformOutput', 0), 'UniformOutput', 0)),1);
% ExactSignal = @(sigma, Phi) pi*exp(-sigma^2)*sin(bsxfun(@minus,Phi,pi/2));
ExactSignal = @(sigma, Phi) exp(-sigma^2 )*sin(Phi);

% ExactSignal = @(sigma, OriginalSignal) exp(-sigma^2/2)*OriginalSignal;

% SignalOriginal = sum(cell2mat(cellfun(@times, num2cell(sqrt( bsxfun(@plus, 2*(1:MaxDegreeL), 1)' / (4*pi) ) ), SphericalHarmonic, 'UniformOutput', 0) ),1)';


% SignalOriginal = findExactSignal(0, Sphere.Theta, Sphere.Phi, MaxDegreeL);
% SignalOriginal = ExactSignal(0, SphericalHarmonic, MaxDegreeL);
SignalOriginal = ExactSignal(0, Sphere.Phi);
% SignalOriginal = cos(Sphere.Theta);
% SignalOriginal = findExactSignal(0, Sphere.Theta, Sphere.Phi, 1);
% SignalOriginal = zeros(Sphere.LocationCount,1);
% SignalOriginal(end) = 1; 
% SignalOriginal = findExactSignal(0, Sphere.Theta, Sphere.Phi, MaxDegreeL);
% SignalOriginal = sum(cell2mat(SphericalHarmonic),1)';
Sphere.Signal = SignalOriginal;


% 
% DistanceFromImpulse = acos( Sphere.Location(1,:) * Sphere.Location')';
% EuclidDistanceFromImpulse = sqrt(sum(bsxfun(@minus, Sphere.Location, Sphere.Location(1,:)).^2,2));
% DistanceFromImpulse = asin( (EuclidDistanceFromImpulse/2) .* sqrt(4-EuclidDistanceFromImpulse.^2));
% 
% 
% ExactSignal = @(sigma, OriginalSignal)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the Laplace Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FileLocationMeshLP = strcat(ProjectRoot,'/models/sphere/meshLP/');
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

I611tL = speye(Alength, Alength) - (6/11)*alpha*tauExplicit * LBM;

I1225tL = speye(Alength, Alength) - (12/25)*alpha*tauExplicit * LBM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale Parameter Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ScaleParameter = findScaleParamter(tauExplicit, alpha, NumStepsExplicit, 'Laplacian', 'natural');
% ScaleParameter = findScaleParamter(tauExplicit, alpha, NumStepsExplicit, 'Laplacian', 'cutoff');

% ScaleParameter = sqrt((0:NumStepsExplicit) * 2 * 2 * alpha * tauExplicit)';


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
    figure(1)
end

for i = 1 : NumStepsExplicit - 1
    
    if i == 1
        [SignalExplicit(:,i+1), flag] = bicg(ItL, SignalExplicit(:,i), 1e-10, 100);
    elseif i == 2
        [SignalExplicit(:,i+1), flag] = bicg(I23tL, (4/3)*SignalExplicit(:,i) - (1/3)*SignalExplicit(:,i-1), 1e-10, 100);
    elseif i == 3
        [SignalExplicit(:,i+1), flag] = bicg(I611tL, (18/11)*SignalExplicit(:,i) - (9/11)*SignalExplicit(:,i-1) + (2/11)*SignalExplicit(:,i-2), 1e-10, 100);
    else
        [SignalExplicit(:,i+1), flag] = bicg(I1225tL, (48/25)*SignalExplicit(:,i)-(36/25)*SignalExplicit(:,i-1) + (16/25)*SignalExplicit(:,i-2) - (3/25)*SignalExplicit(:,i-3), 1e-10, 100);
    end
    
    
    
    if flag
        flag
    end
    
%     TruthGauss = exp(-DistanceFromImpulse.^2 / (2 * ScaleParameter(i+1)^2) );
%     TruthGauss = TruthGauss * max(SignalExplicit(:,i+1));
%     Truth(:,i+1) = TruthGauss;
    
%         Truth(:,i+1) = findExactSignal(ScaleParameter(i+1), Sphere.Theta, Sphere.Phi, MaxDegreeL);
%     Truth(:,i+1) = ExactSignal(ScaleParameter(i+1), SphericalHarmonic);
%     Truth(:,i+1) = ExactSignal(ScaleParameter(i+1), SphericalHarmonic, MaxDegreeL);
    Truth(:,i+1) = ExactSignal(ScaleParameter(i+1), Sphere.Phi);
	AbsErr(i+1,1) = norm(Truth(:,i+1) - SignalExplicit(:,i+1), inf);
    
    if ShowPlot
        clf
        plot(Sphere.Phi, SignalOriginal,'ko')
        hold on
        plot(Sphere.Phi, Truth(:,i), 'gd')
        plot(Sphere.Phi, SignalExplicit(:,i),'r.')
    end
    
    waitbar(i/NumStepsExplicit, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumStepsExplicit-1));
end

waitbar(i/NumStepsExplicit, WaitBar, sprintf('Diffusion Complete'));
close(WaitBar)
if ShowPlot
    close(figure(1))
end

figure
loglog(1:NumStepsExplicit, AbsErr)

% MCError(MCs, 1:2, MCp) = [NumStepsImplicit, AbsErr(NumStepsImplicit - 1)];
% MCErrorAll{MCs, MCp} = [(1:NumStepsImplicit)', AbsErr];


% Mean3 = mean( SignalExplicit3 );
% Mean4 = mean( SignalExplicit4 );
% Mean5 = mean( SignalExplicit5 );
% 
% figure
% plot(1:length(Mean3), Mean3);
% hold on
% plot(1:length(Mean4), Mean4);
% plot(1:length(Mean5), Mean5);
% legend('3','4','5')




% for i = 1 : NumStepsExplicit
%     
%     V(i,1) = var( bsxfun(@minus, SignalExplicit(:,i), sum(SignalExplicit(:,i))/Sphere.LocationCount) );
%     M(i,1) = abs(mean(bsxfun(@minus, SignalExplicit(:,i),  sum(SignalExplicit(:,i))/Sphere.LocationCount )));
%     
% end
% 
% 
% loglog(1:NumStepsExplicit, V)
% loglog(1:NumStepsExplicit, M)
