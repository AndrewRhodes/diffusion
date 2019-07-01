#include <mex.h>
#include "meshGausMatrix.h"
#include <string>
#include <Eigen/Eigen>
#include <engine.h>
#include "geodesics/geodesic_algorithm_exact.h"

void generate_mesh_gaus_matrix(char* filename, unsigned int htype, double hs, double rho, double &h, vector<unsigned int>& Iv, vector<unsigned int>& Jv, vector<double>& Sv, unsigned int &nv){

  TMesh mesh;
  if( !(mesh.ReadOffFile(filename)) ){
    mexPrintf("Failed to open file %s\n", filename);
  }

  double maxs, mins, avers;
  mesh.MeshSize(maxs, mins, avers);

  mexPrintf("avers: %f\n", avers);

  if(htype == 0){
    h = avers * hs;
  }
  else{
    h = hs;
  }


  nv = mesh.v_count();
  unsigned int nf = mesh.f_count();
  unsigned int vid;
  vector<pair<unsigned int, double> >vgdists;
  // vector<double> totalweight(nv, 0);

  double hh = 2.0 * 2.0 * h * h;
  double weight;
  // mexPrintf("hh: %0.3f\n",hh);
  for(unsigned int i = 0; i < nv; i++){

    vgdists.clear();
    compute_one2part_Euclidean_vdist(i, mesh, vgdists, h*rho);

    for(unsigned int j = 0; j < vgdists.size(); j++){
      vid = vgdists[j].first;

      weight = (4.0 / (M_PI * hh * hh)) * exp(-vgdists[j].second * vgdists[j].second / hh);

      Iv.push_back(i + 1);
      Jv.push_back(vid + 1);
      Sv.push_back(weight);

    }
  }
}




void generate_mesh_gaus_matrix(char* filename, unsigned int htype, double hs, double rho, double &h, vector<unsigned int>& Iv, vector<unsigned int>& Jv, vector<double>& Sv, vector<double>& Av, unsigned int &nv){

  TMesh mesh;
  if( !(mesh.ReadOffFile(filename)) ){
    mexPrintf("Failed to open file %s\n", filename);
  }

  double maxs, mins, avers;
  mesh.MeshSize(maxs, mins, avers);

  mexPrintf("avers: %f\n", avers);

  if(htype == 0){
    h = avers * hs;
  }
  else{
    h = hs;
  }


  nv = mesh.v_count();
  unsigned int nf = mesh.f_count();
  unsigned int vid;

  Av.clear();
  Av.resize(nv,0);

  unsigned int vid0, vid1, vid2;
  double farea;
  for(unsigned int f = 0; f < nf; f++){
    vid0 = mesh.facet(f).vert(0);
		vid1 = mesh.facet(f).vert(1);
		vid2 = mesh.facet(f).vert(2);
		VECTOR3 vv = cross( mesh.vertex(vid1).coord() - mesh.vertex(vid0).coord(),
								  mesh.vertex(vid2).coord() - mesh.vertex(vid0).coord() );
		farea = norm(vv) / 2.0;

		Av[vid0] += farea / 3;
		Av[vid1] += farea / 3;
		Av[vid2] += farea / 3;
  }




  vector<pair<unsigned int, double> >vgdists;
  vector<double> totalweight(nv, 0);
  vector<double> areaweight(nv, 0);

  double hh = h * h;
  double weight;

  for(unsigned int i = 0; i < nv; i++){

    vgdists.clear();
    compute_one2part_Euclidean_vdist(i, mesh, vgdists, h*rho);

    for(unsigned int j = 0; j < vgdists.size(); j++){
      vid = vgdists[j].first;

      weight = (4.0 / (M_PI * hh * hh)) * exp(-vgdists[j].second * vgdists[j].second / hh);

      areaweight[i] *= Av[vid] * Av[i];
      areaweight[vid] *= Av[vid] * Av[i];

      Iv.push_back(i + 1);
      Jv.push_back(vid + 1);
      Sv.push_back(weight);

    }
  }
}









void generate_mesh_gaus_geod_matrix(char* filename, unsigned int htype, double hs, double rho, double &h, vector<unsigned int>& Iv, vector<unsigned int>& Jv, vector<double>& Sv, unsigned int &nv){

  TMesh mesh;
  if( !(mesh.ReadOffFile(filename)) ){
    mexPrintf("Failed to open file %s\n", filename);
  }

  double maxs, mins, avers;
  mesh.MeshSize(maxs, mins, avers);

  mexPrintf("avers: %f\n", avers);

  if(htype == 0){
    h = avers * hs;
  }
  else{
    h = hs;
  }


  nv = mesh.v_count();
  unsigned int nf = mesh.f_count();


  //for geodesic computation
  vector<double> vertices;
  vector<unsigned int> faces;
  for(unsigned int i = 0; i < nv; i ++){
    vertices.push_back((mesh.vertex(i).coord())[0]);
    vertices.push_back((mesh.vertex(i).coord())[1]);
    vertices.push_back((mesh.vertex(i).coord())[2]);
  }
  for(unsigned int i = 0; i < nf; i ++){
    faces.push_back(mesh.facet(i).vert(0));
    faces.push_back(mesh.facet(i).vert(1));
    faces.push_back(mesh.facet(i).vert(2));
  }

  geodesic::Mesh geod_mesh;
	geod_mesh.initialize_mesh_data(vertices, faces);    //create internal mesh data structure including edges
	geodesic::GeodesicAlgorithmExact algorithm(&geod_mesh); //create exact algorithm for the mesh


  vector<pair<unsigned int, double> >vgdists;
  vector<double> totalweight(nv, 0);
  unsigned int vid;
  double hh = h * h;

  for(unsigned int i = 0; i < nv; i++){

    vgdists.clear();
    compute_one2part_Geodesic_vdist(i, geod_mesh, algorithm, vgdists, h * rho);

    for(unsigned int j = 0; j < vgdists.size(); j++){
      vid = vgdists[j].first;

      double weight = (4.0 / (M_PI * hh * hh)) * exp(-vgdists[j].second * vgdists[j].second / hh);

      Iv.push_back(i+1);
      Jv.push_back(vid + 1);
      Sv.push_back(weight);

    }
  }
}





void generate_mesh_gaus_geod_matrix(char* filename, unsigned int htype, double hs, double rho, double &h, vector<unsigned int>& Iv, vector<unsigned int>& Jv, vector<double>& Sv, vector<double>& Av, unsigned int &nv){

  TMesh mesh;
  if( !(mesh.ReadOffFile(filename)) ){
    mexPrintf("Failed to open file %s\n", filename);
  }

  double maxs, mins, avers;
  mesh.MeshSize(maxs, mins, avers);

  mexPrintf("avers: %f\n", avers);

  if(htype == 0){
    h = avers * hs;
  }
  else{
    h = hs;
  }


  nv = mesh.v_count();
  unsigned int nf = mesh.f_count();

  Av.clear();
  Av.resize(nv,0);

  unsigned int vid0, vid1, vid2;
  double farea;
  for(unsigned int f = 0; f < nf; f++){
    vid0 = mesh.facet(f).vert(0);
		vid1 = mesh.facet(f).vert(1);
		vid2 = mesh.facet(f).vert(2);
		VECTOR3 vv = cross( mesh.vertex(vid1).coord() - mesh.vertex(vid0).coord(),
								  mesh.vertex(vid2).coord() - mesh.vertex(vid0).coord() );
		farea = norm(vv) / 2.0;

		Av[vid0] += farea / 3;
		Av[vid1] += farea / 3;
		Av[vid2] += farea / 3;
  }



  //for geodesic computation
  vector<double> vertices;
  vector<unsigned int> faces;
  for(unsigned int i = 0; i < nv; i ++){
    vertices.push_back((mesh.vertex(i).coord())[0]);
    vertices.push_back((mesh.vertex(i).coord())[1]);
    vertices.push_back((mesh.vertex(i).coord())[2]);
  }
  for(unsigned int i = 0; i < nf; i ++){
    faces.push_back(mesh.facet(i).vert(0));
    faces.push_back(mesh.facet(i).vert(1));
    faces.push_back(mesh.facet(i).vert(2));
  }

  geodesic::Mesh geod_mesh;
	geod_mesh.initialize_mesh_data(vertices, faces);    //create internal mesh data structure including edges
	geodesic::GeodesicAlgorithmExact algorithm(&geod_mesh); //create exact algorithm for the mesh


  vector<pair<unsigned int, double> >vgdists;
  vector<double> totalweight(nv, 0);
  vector<double> areaweight(nv, 0);
  unsigned int vid;
  double hh = h * h;

  for(unsigned int i = 0; i < nv; i++){

    vgdists.clear();
    compute_one2part_Geodesic_vdist(i, geod_mesh, algorithm, vgdists, h * rho);

    for(unsigned int j = 0; j < vgdists.size(); j++){
      vid = vgdists[j].first;


      double weight = (4.0 / (M_PI * hh * hh)) * exp(-vgdists[j].second * vgdists[j].second / hh);

      areaweight[i] *= Av[vid] * Av[i];
      areaweight[vid] *= Av[vid] * Av[i];

      Iv.push_back(i+1);
      Jv.push_back(vid + 1);
      Sv.push_back(weight);

    }
  }
}



















void compute_one2part_Euclidean_vdist(unsigned int vid_start, TMesh& mesh, vector<pair<unsigned int, double> >& vgdists, double maxdist){
	for(unsigned int j = 0; j < mesh.v_count(); j ++){
		double d = sqrt( dot(mesh.vertex(vid_start).coord() - mesh.vertex(j).coord(),
									mesh.vertex(vid_start).coord() - mesh.vertex(j).coord()) );
		if(d <= maxdist){
		 	vgdists.push_back( make_pair(j, d) );
		}
	}
}




void compute_one2part_Geodesic_vdist(unsigned int vid_start, geodesic::Mesh& geod_mesh, geodesic::GeodesicAlgorithmExact& algorithm, vector<pair<unsigned int, double> >& vgdists, double maxdist)
{
	geodesic::SurfacePoint source(&geod_mesh.vertices()[vid_start]);      //create source
	vector<geodesic::SurfacePoint> all_sources(1,source);  //in general, there could be multiple sources, but now we have only one
	algorithm.propagate(all_sources, maxdist);   //cover the whole mesh

	double distance;
	for(unsigned int i = 0; i < geod_mesh.vertices().size(); ++ i){
   	geodesic::SurfacePoint p(&geod_mesh.vertices()[i]);
		unsigned int best_source = algorithm.best_source(p,distance);      //for a given surface point, find closets source and distance to this source
		if(distance <=  maxdist){
			vgdists.push_back( make_pair(i, distance) );
			//cout << distance <<endl;    //print geodesic distance for every vertex
		}
	}
}
