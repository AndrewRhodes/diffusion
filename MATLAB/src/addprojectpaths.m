% Andrew Rhodes
% Add all subfolders of project root directory


function addprojectpaths

% find root of project folder
CurDir = pwd;
ProjectRootStart = strfind(CurDir, 'MATLAB');
ProjectRoot = CurDir(1:ProjectRootStart-1);

% Add all directories in project root


addpath(genpath(ProjectRoot))



end