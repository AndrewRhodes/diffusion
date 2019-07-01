








% MEX meshGauss.cpp
mex -I~/Documents/Eigen/ -I~/GitAndrew/diffusion/MATLAB/src/MeshLP/ ...
    'GCC=/home/andrew/gcc-6.3.0/bin/gcc' ...
    ~/GitAndrew/diffusion/MATLAB/src/MeshGaus/meshGauss.cpp ...
    ~/GitAndrew/diffusion/MATLAB/src/MeshGaus/meshGausMatrix.cpp ...        
    ~/GitAndrew/diffusion/MATLAB/src/MeshLP/tmesh.cpp ...
    ~/GitAndrew/diffusion/MATLAB/src/MeshLP/point.cpp ...
    ~/GitAndrew/diffusion/MATLAB/src/MeshLP/matrix.cpp ...
    ~/GitAndrew/diffusion/MATLAB/src/MeshLP/offobj.cpp ...    
    -outdir ~/GitAndrew/diffusion/MATLAB/src/MeshGaus/






opt.hs = 2;
opt.rho = 10;
opt.htype = 0;
opt.dtype = 0;
opt.atype = 0;

global ProjectRoot;
FileLocationModel = strcat(ProjectRoot,'/models/object/');
ModelFolder = 'armadillo/';
Model = 'Armadillo_e1_100000';
FileNameModelOff = strcat(ModelFolder,Model,'.off');
filename = fullfile( FileLocationModel, FileNameModelOff );
% mesh_Gauss(filename, opt)
[I, J, S, h] = meshGauss(filename, opt);

W = sparse(I,J,S);


Wsum = sparse(1:length(W), 1:length(W), 1./sum(W,2));    
W = Wsum * W;
    
    
    


opt.hs = 2;
opt.rho = 3;
opt.htype = 'ddr';
opt.dtype = 'euclidean';
    
[II, JJ, SS, AA, h] = symmshlpmatrix(filename, opt);
WW = sparse(II,JJ,SS);

