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

% Max Time for Iteration
MaxLevel = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program Defined variables 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% degree of interpolation
porder = 3; % I think this can be 2 for Bertalmio method

% dimension of the surface
dim = 3;

% Order of numerical differncing. 1 for 2nd order, 2 for 4th order
Lorder = 1;

% define bandwidth around surface
bandwidth = 1.00001*spacing*sqrt( (dim-1)*((porder+1)/2)^2 + ((Lorder+(porder+1)/2)^2) );


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
% Load Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MinPoint = [-2,-2,-2];
MaxPoint = [2,2,2];

[xMesh, yMesh, zMesh] = meshgrid(MinPoint(1):spacing:MaxPoint(1),MinPoint(2):spacing:MaxPoint(2),MinPoint(3):spacing:MaxPoint(3));


[cpx, cpy, cpz, dist] = cpSphere(xMesh, yMesh, zMesh, 1, [0,0,0]);

PointsInBandLogic = (abs(dist) <= bandwidth );

band = find(PointsInBandLogic);

CP = [cpx(:), cpy(:), cpz(:)];
CP = CP(band,:);

% XYZ = MinPoint - spacing + IJK * spacing

% Ind = round((bsxfun(@minus, CP, MinPoint) + spacing) ./ spacing);

xMeshInd = (xMesh - MinPoint(1) + spacing)./spacing ;
yMeshInd = (yMesh - MinPoint(2) + spacing)./spacing ;
zMeshInd = (zMesh - MinPoint(3) + spacing)./spacing ;

IJK = [xMeshInd(:), yMeshInd(:), zMeshInd(:)];
IJK = IJK(band,:);

XYZ = [xMesh(:), yMesh(:), zMesh(:)];
XYZ = XYZ(band,:);

NumGridPoints = length(IJK);


Ecp = interp3_matrix(MinPoint(1):spacing:MaxPoint(1), MinPoint(2):spacing:MaxPoint(2), MinPoint(3):spacing:MaxPoint(3), CP(:,1), CP(:,2), CP(:,3), porder, band);


[xs,ys,zs] = sphere(100);
[th_plot, phi_plot, r] = cart2sph(xs(:), ys(:), zs(:));

Eplot = interp3_matrix(MinPoint(1):spacing:MaxPoint(1), MinPoint(2):spacing:MaxPoint(2), MinPoint(3):spacing:MaxPoint(3), xs(:), ys(:), zs(:), porder, band);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate surface normal vectors for each grid point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% direction between closest point and grid point, normalized to unit length

CPNormal = bsxfun(@rdivide, XYZ, sqrt(sum(XYZ.^2,2)));
CPNormal(dist(band)<0,:) = - CPNormal(dist(band)<0,:);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Nearest Neighbor in axis directions to construct connectivity matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BandSearchSize = [length(cpx), length(cpy), length(cpz)];

GridPoint = (1:NumGridPoints)';
GridPointOnes = ones(NumGridPoints,1);

ConnectivityLocationsWeights = cell(6,1);


for PosNegSign = 1 : 2
    for ax = 1 : 3
        
        NeighIJK = IJK;
        
        NeighIJK(:,ax) = NeighIJK(:,ax) + (-1)^PosNegSign * GridPointOnes;
        
        NeighInd = sub2ind(BandSearchSize, NeighIJK(:,1), NeighIJK(:,2), NeighIJK(:,3));
        
        [BandNeighInd, IdxBand, IdxNeighInd] = intersect(band, NeighInd);
        
        ConnectivityLocationsWeights{PosNegSign^2 + ax - 1, 1} = [IdxBand, IdxNeighInd];
    end
end



% Connectivity Axis (x,y,z) Direction (Negative, Positive)
Cxn = sparse([ConnectivityLocationsWeights{1,1}(:,1);GridPoint], [ConnectivityLocationsWeights{1,1}(:,2);GridPoint], [-ones(length(ConnectivityLocationsWeights{1,1}(:,1)),1);GridPointOnes], NumGridPoints, NumGridPoints);
Cyn = sparse([ConnectivityLocationsWeights{2,1}(:,1);GridPoint], [ConnectivityLocationsWeights{2,1}(:,2);GridPoint], [-ones(length(ConnectivityLocationsWeights{2,1}(:,1)),1);GridPointOnes], NumGridPoints, NumGridPoints);
Czn = sparse([ConnectivityLocationsWeights{3,1}(:,1);GridPoint], [ConnectivityLocationsWeights{3,1}(:,2);GridPoint], [-ones(length(ConnectivityLocationsWeights{3,1}(:,1)),1);GridPointOnes], NumGridPoints, NumGridPoints);

Cxp = sparse([ConnectivityLocationsWeights{4,1}(:,1);GridPoint], [ConnectivityLocationsWeights{4,1}(:,2);GridPoint], [ones(length(ConnectivityLocationsWeights{4,1}(:,1)),1);-GridPointOnes], NumGridPoints, NumGridPoints);
Cyp = sparse([ConnectivityLocationsWeights{5,1}(:,1);GridPoint], [ConnectivityLocationsWeights{5,1}(:,2);GridPoint], [ones(length(ConnectivityLocationsWeights{5,1}(:,1)),1);-GridPointOnes], NumGridPoints, NumGridPoints);
Czp = sparse([ConnectivityLocationsWeights{6,1}(:,1);GridPoint], [ConnectivityLocationsWeights{6,1}(:,2);GridPoint], [ones(length(ConnectivityLocationsWeights{6,1}(:,1)),1);-GridPointOnes], NumGridPoints, NumGridPoints);



% Normals Axis (x,y,z) Diagonal
Nxd = sparse(GridPoint, GridPoint, reshape(CPNormal(:,1),[],1), NumGridPoints, NumGridPoints);
Nyd = sparse(GridPoint, GridPoint, reshape(CPNormal(:,2),[],1), NumGridPoints, NumGridPoints);
Nzd = sparse(GridPoint, GridPoint, reshape(CPNormal(:,3),[],1), NumGridPoints, NumGridPoints);



SparseIdentity = speye(NumGridPoints, NumGridPoints);

dt = spacing/4;

SumAlphaPos = Nxd*Cxp + Nyd*Cyp + Nzd*Czp;
SumBetaPos = Cxn * ( Cxp - Nxd*SumAlphaPos ) + Cyn * ( Cyp - Nyd*SumAlphaPos ) + Czn * ( Czp - Nzd*SumAlphaPos );


SumAlphaNeg = Nxd*Cxn + Nyd*Cyn + Nzd*Czn;
SumBetaNeg = Cxp * ( Cxn - Nxd*SumAlphaNeg ) + Cyp * ( Cyn - Nyd*SumAlphaNeg ) + Czp * ( Czn - Nzd*SumAlphaNeg );



ISBPie = (SparseIdentity - dt*SumBetaPos);
ISBNie = (SparseIdentity - dt*SumBetaNeg);

ISBPee = (SparseIdentity + dt*SumBetaPos);
ISBNee = (SparseIdentity + dt*SumBetaNeg);

% IFullPie = SparseIdentity - dt * (Cxn*Cxp + Cyn*Cyp + Czn*Czp)*E;
% IFullNie = SparseIdentity - dt * (Cxp*Cxn + Cyp*Cyn + Czp*Czn)*E;



% Need to remove the positive eigen-values similar to MacDonald. 


NumLevels = round(MaxLevel/dt);

[th, phi, r] = cart2sph(xMesh, yMesh, zMesh);
u = cos(phi + pi/2);

u = u(band);



Signal = zeros(length(u), NumLevels);
Signal(:,1) = u;



UseImplicit = 0;
UseExplicit = 1;

if UseExplicit && UseImplicit
    error('Can only use one Euler method.')
end


figure(1)
WaitBar = waitbar(0, sprintf('Implicit Euler Diffusion %i of %i', 0, NumLevels-1));

for level = 2 : NumLevels
    
% Implicit Euler
    if UseImplicit
        
        if mod(level-1,2)
            
            [Signal(:,level), flag] = bicg(ISBPie, Signal(:,level-1), 1e-10, 50);
            
            
        elseif mod(level,2)
            
            [Signal(:,level), flag] = bicg(ISBNie, Signal(:,level-1), 1e-10, 50);
            
        else
            sprintf('Mistake')
        end
        
        
        if flag
            flag
        end
        
        
    end
    
    % Explicit Euler
    if UseExplicit 
        
        if mod(level-1,2)
            
            Signal(:,level) = ISBPee * Signal(:,level-1);
            
            
        elseif mod(level,2)
            
            Signal(:,level) = ISBNee * Signal(:,level-1);

        else
            sprintf('Mistake')
        end
        
        Signal(:,level) = Ecp * Signal(:,level);
    end
    

    
    t = level*dt;
    
    sphplot = Eplot*Signal(:,level);
    err = norm(exp(-2*t)*cos(phi_plot + pi/2)-sphplot,inf) / norm(exp(-2*t)*cos(phi_plot + pi/2),inf);
    [t dt spacing err]
    surf(xs, ys, zs, reshape(sphplot, size(xs)));
    xlabel('x'); 
    ylabel('y');
    zlabel('z');
    axis equal; 
    shading interp;
    colorbar;
    drawnow(); 
    pause(0);
      
      
    waitbar(level/NumLevels, WaitBar, sprintf('Implicit Euler Diffusion %i of %i', level, NumLevels-1));
    
    
end


waitbar(level/NumLevels, WaitBar, sprintf('Diffusion Complete'));
close(WaitBar)














