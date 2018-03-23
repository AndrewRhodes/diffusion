% Andrew Rhodes
% October 2017
% Diffusion on a mesh using Colin B. Macdonald's closest point algorithm,
% Bertalmio's level set approach, and my implicit Euler method
% On a Sphere


close all
clear
clc



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Defined variables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spacing of implicit surface points
spacing = 0.1; 

% Max Time for Iteration Implicit
tauImplicit = spacing / 4;

MaxLevelImplicit = 5;

NumLevelsImplicit = MaxLevelImplicit / tauImplicit;


% degree of interpolation
porder = 3; % I think this can be 2 for Bertalmio method

% dimension of the surface
dim = 3;

% Order of numerical differncing. 1 for 2nd order, 2 for 4th order
Lorder = 2;

% define bandwidth around surface
bandwidth = 1.00001*spacing*sqrt( (dim-1)*((porder+1)/2)^2 + ((Lorder/2+(porder+1)/2)^2) );

IcosphereDiffisions = 3;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional Paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('~/AFOSR/MATLAB/')
addpath('~/Desktop/Ashish/CS3 Code/')
addpath('~/GitProjects/pose/MATLAB_PointCloudDescriptors/OURCVFH/models/')
addpath('~/GitProjects/matlab-utilities/')
addpath('~/GitProjects/mesh_resampling_toolbox/MATLAB_Modules/')
addpath(genpath('~/Documents/Software/cp_matrices/'))
addpath('~/Documents/Software/')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make Explicit Surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[vertices, faces] = icosphere(IcosphereDiffisions);

PointCloud.Location = vertices;
PointCloud.Face = faces;
PointCloud.LocationCount = size(PointCloud.Location,1);
PointCloud.FaceCount = size(PointCloud.Face, 1);
PointCloud.FaceNormal = findFaceNormals(PointCloud.Location, PointCloud.Face);
PointCloud.FaceArea = findFaceAreas(PointCloud.Location, PointCloud.Face);
PointCloud = findMeshNormals(PointCloud);

clear faces vertices

save_off(PointCloud.Location, PointCloud.Face, strcat('Icosphere',num2str(IcosphereDiffisions),'.off'));
[LapMatMeshWeights, Area, hEdge] = symmshlp_matrix(strcat('Icosphere',num2str(IcosphereDiffisions),'.off'));

tauExplicit = (hEdge/2)^2;

MaxTauExplicit  = 5;
NumStepsExplcit = round(MaxTauExplicit / tauExplicit);

A1 = sparse(1:length(Area),1:length(Area), 1./Area);

LBM = A1 * LapMatMeshWeights;

ItL = speye(length(Area),length(Area)) - tauExplicit * LBM;

% Find the scale parameters
maxsample = 3; 
ws = 0 : 0.0001 : maxsample;
NumSample = length(ws);
ScaleParameterSpatial = zeros(NumStepsExplcit,1);
cut = sqrt(log(2));
db3 = 1/sqrt(2);
ws2 = (ws.^2)';
H = 1 ./ ( ones(NumSample,1) + 1 * ws2 );
h = 1 ./ ( ones(NumSample,1) + 1 * ws2 );

for i = 1 : NumStepsExplcit
    % Transfer function
    % 1st order
    H = H .* h;
    % Find the frequency at the cutoff values
    CutoffFrequency = interp1(H, ws, db3);
    % Change cutoff frequency to scale parameter
    ScaleParameterFrequency = CutoffFrequency / cut;   
    ScaleParameterSpatial(i,1) =  sqrt(tauExplicit) / ScaleParameterFrequency;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make starting verticies for implicit grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MinPoint = round(min(PointCloud.Location) - bandwidth - spacing, 1);
MaxPoint = round(max(PointCloud.Location) + bandwidth + spacing, 1);

x1d = (MinPoint(1):spacing:MaxPoint(1))';
y1d = (MinPoint(2):spacing:MaxPoint(2))';
z1d = (MinPoint(3):spacing:MaxPoint(3))';

% Find Closet Point

% XYZ = MinPoint - spacing + IJK * spacing

% CP is closest point on triangular surface, not necessarily a vertex,
% could exist on a face.

[IJK,DIST,CP,XYZ,CPFACE] = tri2cp(PointCloud.Face, PointCloud.Location, spacing, MinPoint, porder, 1);


NumGridPoints = length(IJK);

BandSearchSize = [length(x1d), length(y1d), length(z1d)];

Band = sub2ind(BandSearchSize, IJK(:,1), IJK(:,2), IJK(:,3));



% Interpolate vertex values to closest points on surface
FaceInterpolateWeights = interpBarycenterTriangle(PointCloud, CP, CPFACE);


% Create L, E, M
L = laplacian_3d_matrix(x1d,y1d,z1d, Lorder, Band, Band);

Esphere = interpLagrange3D(BandSearchSize, MinPoint, PointCloud.Location, porder, Band, spacing);
Ecp = interpLagrange3D(BandSearchSize, MinPoint, CP, porder, Band, spacing);

M = lapsharp(L, Ecp);

ItM = speye(size(M)) - tauImplicit * M;


% Extrapolate Signal
[th, phi, r] = cart2sph(PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3));

SignalExplicit = zeros(PointCloud.LocationCount, NumStepsExplcit);
SignalExplicit(:,1) = cos(phi + pi/2);

SignalImplicit = zeros(PointCloud.LocationCount, NumLevelsImplicit);
SignalImplicit(:,1) = cos(phi + pi/2);

Signal = FaceInterpolateWeights * SignalImplicit(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform Diffusion 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Implcit
for i = 1 : NumLevelsImplicit - 1
    
    
    [SignalNew, flag] = bicg(ItM, Signal, 1e-12, 30);
    if flag
        flag
    end
    
    Signal = SignalNew;
    
    SignalImplicit(:,i+1) = Esphere * Signal;
    
    t = i*tauImplicit;
    
    sphplot = SignalImplicit(:,i+1);
    
    err = norm(exp(-2*t)*cos(phi + pi/2)-sphplot,inf) / norm(exp(-2*t)*cos(phi + pi/2),inf);
    [t tauImplicit spacing err]
    
    sphplot = reshape(sphplot, PointCloud.LocationCount, 1);
    trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), sphplot);
    title( ['soln at time ' num2str(t) ', i= ' num2str(i)] );
    xlabel('x'); ylabel('y'); zlabel('z');
    %caxis([-1.05 1.05]);   % lock color axis
    axis equal; shading interp;
    colorbar;
    drawnow(); pause(0);
    
    
end




for i = 1 : NumStepsExplcit - 1
    
    
    [SignalExplicit(:,i+1), flag] = bicg(ItL, SignalExplicit(:,i), 1e-12, 30);
    if flag
        flag
    end
    
    t = ScaleParameterSpatial(i,1);
    
    sphplot = SignalExplicit(:,i+1);
    
    err = norm(exp(-2*t)*cos(phi + pi/2)-sphplot,inf) / norm(exp(-2*t)*cos(phi + pi/2),inf);
    [t tauExplicit hEdge err]
    
     sphplot = reshape(sphplot, PointCloud.LocationCount, 1);
    trisurf(PointCloud.Face, PointCloud.Location(:,1), PointCloud.Location(:,2), PointCloud.Location(:,3), sphplot);
    title( ['soln at time ' num2str(t) ', i= ' num2str(i)] );
    xlabel('x'); ylabel('y'); zlabel('z');
    %caxis([-1.05 1.05]);   % lock color axis
    axis equal; shading interp;
    colorbar;
    drawnow(); pause(0);
    
    
    pause(0.2)
end























