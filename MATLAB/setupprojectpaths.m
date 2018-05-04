% Andrew Rhodes
% Add all subfolders of project root directory


function setupprojectpaths

global ProjectRoot;

% find root of project folder
CurDir = pwd;
ProjectRootStart = strfind(CurDir, 'MATLAB');
ProjectRoot = CurDir(1:ProjectRootStart+5);

% Add all directories in project root


addpath(genpath(ProjectRoot));



%% Setup Directoryies for Plane
if ~exist( strcat(ProjectRoot, '/models/plane/'), 'dir')
    mkdir( ProjectRoot, 'models/plane' )
end

if ~exist( strcat(ProjectRoot, '/models/plane/CPLaplace'), 'dir')
    mkdir( ProjectRoot, 'models/plane/CPLaplace' )
end

if ~exist( strcat(ProjectRoot, '/models/plane/meshLP'), 'dir')
    mkdir( ProjectRoot, 'models/plane/meshLP' )
end

%% Setup Directories for circle
if ~exist( strcat(ProjectRoot, '/models/circle/'), 'dir')
    mkdir( ProjectRoot, 'models/circle/' )
end

if ~exist( strcat(ProjectRoot, '/models/circle/CPGauss/'), 'dir')
    mkdir( ProjectRoot, 'models/circle/CPGauss/' )
end

if ~exist( strcat(ProjectRoot, '/models/circle/CPLaplace/'), 'dir')
    mkdir( ProjectRoot, 'models/circle/CPLaplace/' )
end


%% Setup Directories for sphere
if ~exist( strcat(ProjectRoot, '/models/sphere'), 'dir')
    mkdir( ProjectRoot, 'models/sphere' )
end

if ~exist( strcat(ProjectRoot, '/models/sphere/CPLaplace'), 'dir')
    mkdir( ProjectRoot, 'models/sphere/CPLaplace' )
end

if ~exist( strcat(ProjectRoot, '/models/sphere/meshLP'), 'dir')
    mkdir( ProjectRoot, 'models/sphere/meshLP' )
end

if ~exist( strcat(ProjectRoot, '/data/'), 'dir')
    mkdir( ProjectRoot, 'data' )
end
    




end