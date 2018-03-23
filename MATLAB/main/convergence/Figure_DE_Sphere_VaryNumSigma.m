% Andrew Rhodes
% ASEL
% March 2018


close all
clear
clc

ProjectRoot = addprojectpaths % Additional Paths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MCdiv = [2, 3, 4, 5, 6, 7];
MCrho = [1, 2, 3, 4, 5, 6, 7, 8, 9];

MCError = zeros(length(MCdiv), 2, length(MCrho));
MCErrorAll = cell(length(MCdiv), length(MCrho));

for MCr = 1 : length(MCrho)
   
for MCd = 1 : length(MCdiv)
    
    clearvars -except MCr MCd MCdiv MCrho MCErrorAll MCError

    NumberDivisions = MCdiv(MCd)
    options.rho = MCrho(MCr)



%     options.rho = 5;
% options.dtype = 'geodesic';
options.dtype = 'euclidean';

alpha = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup File Name Directions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FileLocation = strcat(ProjectRoot,'/models/Sphere/');
FileName = strcat('Icosphere',num2str(NumberDivisions),'.off');

FileLocationMeshLP = strcat(ProjectRoot,'/models/Sphere/meshLP/');
FileNameLapMat = strcat('LapMatMeshWeights','_Div',num2str(NumberDivisions),'_NumSigma',num2str(options.rho),'_',options.dtype,'.mat');
FileNameArea = strcat('Area','_Div',num2str(NumberDivisions),'_NumSigma',num2str(options.rho),'_',options.dtype,'.mat');
FileNamehEdge = strcat('hEdge2','_Div',num2str(NumberDivisions),'_NumSigma',num2str(options.rho),'_',options.dtype,'.mat');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the Sphere and Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~exist(fullfile(FileLocation, FileName), 'file')
    
    [Sphere.Location, Sphere.Face] = icosphere(NumberDivisions);
    
    save_off(Sphere.Location, Sphere.Face, fullfile(FileLocation, FileName))
    
else
    
    [Sphere.Location, Sphere.Face] = read_off(fullfile(FileLocation, FileName));

    [m, n] = size(Sphere.Location);
    if m < n
        Sphere.Location = Sphere.Location';
    end
    
    [m, n] = size(Sphere.Face);
    if m < n
        Sphere.Face = Sphere.Face';
    end
    
end



Sphere.FaceCount = size(Sphere.Face, 1);
Sphere.LocationCount = size(Sphere.Location,1);
Sphere.FaceArea = findFaceAreas(Sphere.Location,Sphere.Face);

Sphere.FaceArea = findFaceArea(Sphere.Location, Sphere.Face);
Sphere.Signal = zeros(Sphere.LocationCount,1);
Sphere = findMeshResolution(Sphere, 'Model');


% % % % % % % % % % 
tauExplicit = Sphere.Resolution / 8;
MaxTauExplicit = 1 / Sphere.Resolution;
NumStepsExplicit = round(MaxTauExplicit);
% % % % % % % % % % 




[Theta, Phi, Radius] = cart2sph(Sphere.Location(:,1) ,Sphere.Location(:,2), Sphere.Location(:,3));


% [Ymn,THETA,PHI,Xm,Ym,Zm]=spharm(10,10,[length(Theta), length(Phi)],0);


L = 10;
M = 10;

Lmn=legendre(L,cos(Phi));

if L~=0
  Lmn=squeeze(Lmn(M+1,:,:));
end

a1=((2*L+1)/(4*pi));
a2=factorial(L-M)/factorial(L+M);
C=sqrt(a1*a2);

Ymn=C*Lmn'.*exp(i*M*Theta);





% Define the signal
% SignalOriginal = cos(1.5*Phi + pi/2) + sin(2.5*Theta - pi/2);
% SignalOriginal = cos(Theta);
SignalOriginal = reshape(real(Ymn),[],1);

Sphere.Signal = SignalOriginal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the Laplace Beltrami
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if ~exist(fullfile(FileLocationMeshLP, FileNameLapMat), 'file')
    
    [LapMatMeshWeights, Area, hEdge2] = symmshlp_matrix(fullfile(FileLocation, FileName), options);
    
    save(fullfile(FileLocationMeshLP, FileNameLapMat), 'LapMatMeshWeights', '-v7.3')
    save(fullfile(FileLocationMeshLP, FileNameArea), 'Area', '-v7.3')
    save(fullfile(FileLocationMeshLP, FileNamehEdge), 'hEdge2', '-v7.3')
else
    
    load( fullfile(FileLocationMeshLP, FileNameLapMat))
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

% ScaleParameter1 = findScaleParamter(tauExplicit, alpha, NumStepsExplicit, 1, 3);

ScaleParameter = zeros(NumStepsExplicit,1);
for i = 1 : NumStepsExplicit
   
    ScaleParameter(i+1,1) = sqrt(2*i^2*alpha*tauExplicit);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Diffusion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SignalExplicit = zeros(Sphere.LocationCount, NumStepsExplicit);
SignalExplicit(:,1) = Sphere.Signal;

AbsErr = zeros(NumStepsExplicit,1);


WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumStepsExplicit-1));

figure(1)
for i = 1 : NumStepsExplicit - 1
    
%     if i == 1
    %     [SignalExplicit(:,i+1), flag] = bicg(ItL, SignalExplicit(:,i), [], 1e-10, 60);
        [SignalExplicit(:,i+1), flag] = gmres(ItL, SignalExplicit(:,i), [], 1e-10, 100);
%     else
%         [SignalExplicit(:,i+1), flag] = gmres(I23tL, (4/3)*SignalExplicit(:,i) - (1/3)*SignalExplicit(:,i-1), [], 1e-10, 100);
%     end
    
    if flag
        disp(flag)
    end
    
    clf
    plot(Theta, SignalOriginal,'b.')
    hold on
    plot(Theta, SignalExplicit(:,i),'ko')
    
    Truth = exp(-ScaleParameter(i)^2/2) * SignalOriginal;
% % %     Truth = exp(-ScaleParameter(i)^2/2) .* cos(Theta) + exp(-9*ScaleParameter(i)^2/2) .* sin(3*Theta);
    plot(Theta, Truth,'rs')
	AbsErr(i,1) = norm(Truth - SignalExplicit(:,i), inf);
   
    
    waitbar(i/NumStepsExplicit, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', i, NumStepsExplicit-1));
    if i ==1
        pause
    else
        pause%(0.2)
    end
end

waitbar(i/NumStepsExplicit, WaitBar, sprintf('Diffusion Complete'));
close(WaitBar)
close(figure(1))


MCError(MCd, 1:2, MCr) = [NumStepsExplicit, AbsErr(NumStepsExplicit - 1)];
MCErrorAll{MCd, MCr} = [(1:NumStepsExplicit)', AbsErr];


end

end




colors = ['c','m','k','b','r','c','m','g','k'];
shapes = {'','','.','^','.','*','+','o','s'};
lin  = {'-','-',':','-','-','-','-.',':','--'};
msizes = [8, 8, 8, 8, 8, 8, 10, 12, 15];

Points = zeros(size(MCErrorAll,1),2, size(MCErrorAll,2));
for j = 1 : size(MCErrorAll,2)
    for i = 1 : size(MCErrorAll,1)
        Points(i,1:2,j) = [MCErrorAll{i,j}(end,1), MCErrorAll{i,j}(end,2)];
    end
end


figure('units','normalized','outerposition',[0 0 1 1])
for i = 1 : size(Points,3)
    loglog( Points(:,1,i), Points(:,2,i), strcat(colors(i), shapes{i}, lin{i}), 'linewidth', 3, 'markersize', msizes(i) )
    hold on
end





















