


#ifndef _MESHGAUSMATRIX_H_
#define _MESHGAUSMATRIX_H_

#include <vector>
#include "tmesh.h"
#include "geodesics/geodesic_algorithm_exact.h"


using namespace std;


void generate_mesh_gaus_matrix(char* filename, unsigned int htype, double hs, double rho, double &h, vector<unsigned int>& Iv, vector<unsigned int>& Jv, vector<double>& Sv, unsigned int &nv);

void generate_mesh_gaus_matrix(char* filename, unsigned int htype, double hs, double rho, double &h, vector<unsigned int>& Iv, vector<unsigned int>& Jv, vector<double>& Sv, vector<double>& Av, unsigned int &nv);




void generate_mesh_gaus_geod_matrix(char* filename, unsigned int htype, double hs, double rho, double &h, vector<unsigned int>& Iv, vector<unsigned int>& Jv, vector<double>& Sv, unsigned int &nv);

void generate_mesh_gaus_geod_matrix(char* filename, unsigned int htype, double hs, double rho, double &h, vector<unsigned int>& Iv, vector<unsigned int>& Jv, vector<double>& Sv, vector<double>& Av, unsigned int &nv);




void compute_one2part_Euclidean_vdist(unsigned int vid_start, TMesh& mesh, vector<pair<unsigned int, double> >& vgdists, double maxdist);

void compute_one2part_Geodesic_vdist(unsigned int vid_start, geodesic::Mesh& geod_mesh, geodesic::GeodesicAlgorithmExact& algorithm, vector<pair<unsigned int, double> >& vgdists, double maxdist);

#endif // _MESHGAUSMATRIX_H_
