#ifndef _MESH_CUH_
#define _MESH_CUH_

//namespace SPH{

#include "defs.h"

#include <stdio.h>

#ifdef  CUDA_BUILD
#include <cuda.h>
#endif

#include "utils.h"

//namespace MetFEM{

struct elem_data{
    double3             *centroid,*normal;
};
//Element is not anymore here, is everything flattened in mesh_d class
class TriMesh_d{
	
	public:

	//Element 						elem_data;
	double3 						*node,*node_v; //Positions and veloc, 
	int									*elnode;			//3 per element
	double 							*pplane; //Only for testing
  //Element data, TODO: PASS TO ELEMDATA
  double3             *centroid,*normal;
  int                 *nfar;
  double3             m_v,m_w;
  int nodecount, elemcount;
  int                 id;
  
	
	//double							v;						//Constant Uniform v
	TriMesh_d(){}
	void AxisPlaneMesh(const int &axis, bool positaxisorent, const double3 p1, const double3 p2,  const int &dens);
  void SetVel(const double3 &v) {m_v = v;} //Like in WeldForm CPU version
	inline void ApplyConstVel(const double3 &v);
	inline void CalcCentroidVelFromNodes();
	inline dev_t void UpdatePlaneCoeff();
	inline void UpdatePos(const double &dt);
	inline dev_t void Move(double dt);
  inline dev_t void CalcNormals();
	inline dev_t void CalcSpheres();
	inline dev_t void CalcCentroids();
  inline dev_t void CheckNormals();
  
  
};

#ifdef  CUDA_BUILD

__global__ inline void MeshUpdateKernel(TriMesh_d *mesh_d, double dt);
__global__ inline void CalcSpheresKernel(TriMesh_d *mesh_d);

__global__ inline void CheckNormalsKernel(TriMesh_d *mesh_d);
#endif

//};//SPH

#include "Mesh.C"

#endif
