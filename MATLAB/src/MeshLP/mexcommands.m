










% MEX symmshlpmatrix.C
mex -v 'GCC=/usr/bin/gcc-4.7' symmshlpmatrix.cpp tmesh.cpp point.cpp ...
    matrix.cpp offobj.cpp comp_meshlpmatrix.cpp meshlpmatrix.cpp



% MEX cotlpmatrix.C
mex 'GCC=/usr/bin/gcc-4.7' cotlpmatrix.cpp tmesh.cpp point.cpp ...
    matrix.cpp offobj.cpp comp_meshlpmatrix.cpp meshlpmatrix.cpp