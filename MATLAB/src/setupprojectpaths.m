% Andrew Rhodes
% Add all subfolders of project root directory


function ProjectRoot = setupprojectpaths

% find root of project folder
CurDir = pwd;
ProjectRootStart = strfind(CurDir, 'MATLAB');
ProjectRoot = CurDir(1:ProjectRootStart+5);

% Add all directories in project root


addpath(genpath(ProjectRoot))



if ~exist(strcat(ProjectRoot, '/models/Plane/'), 'dir')
    mkdir(ProjectRoot, 'models/Plane')
end

if ~exist(strcat(ProjectRoot, '/models/Plane/CPLaplace'), 'dir')
    mkdir(ProjectRoot, 'models/Plane/CPLaplace')
end

if ~exist(strcat(ProjectRoot, '/models/Plane/meshLP'), 'dir')
    mkdir(ProjectRoot, 'models/Plane/meshLP')
end



if ~exist(strcat(ProjectRoot, '/models/Sphere'), 'dir')
    mkdir(ProjectRoot, 'models/Sphere')
end

if ~exist(strcat(ProjectRoot, '/models/Sphere/CPLaplace'), 'dir')
    mkdir(ProjectRoot, 'models/Sphere/CPLaplace')
end

if ~exist(strcat(ProjectRoot, '/models/Sphere/meshLP'), 'dir')
    mkdir(ProjectRoot, 'models/Sphere/meshLP')
end

if ~exist(strcat(ProjectRoot, '/data/'), 'dir')
    mkdir(ProjectRoot, 'data')
end
    




end