% % ipath =['-I' fullfile(matlabroot,'extern', 'include') ' -I/usr/include/CGAL'];
% % 
% % Lpath =['-L' fullfile(matlabroot,'bin','glnx64') ...
% %        ' -L' fullfile(matlabroot,'sys','os','glnx64')...
% %        ' -L/usr/lib'];
% %    %-lXmu -lXi -lXt -lXext -lICE
% % lpath =['-lmwlapack -lpthread -lCGAL -lgmp -lboost_thread'];
% % 
% % flag = '-largeArrayDims';
% % %mex('-v', ipath, 'graphlpmatrix.C', 'comp_llpmatrix.C', 'lpmatrix.C', 'point_cloud.C', Lpath, lpath);
% % 
% % mex('-v', ipath, 'pcdlpmatrix.cpp', 'comp_llpmatrix.cpp', 'lpmatrix.cpp', 'point_cloud.cpp', Lpath, lpath);
% % 
% % % mex('-v', ipath, 'pcdlp.cpp', 'comp_llpmatrix.cpp', 'point_cloud.cpp', Lpath, lpath);



% MEX pcdlpmatrix.cpp
mex 'GCC=/home/andrew/gcc-6.3.0/bin/gcc' pcdlpmatrix.cpp lpmatrix.cpp comp_llpmatrix.cpp point_cloud.cpp ...
    -L/usr/local/MATLAB/R2018b/bin/glnx64 ...
    -L/usr/local/MATLAB/R2018b/sys/os/glnx64 -L/usr/lib -lmwlapack -lpthread ...
    -lCGAL -lgmp -lboost_thread -lboost_system




