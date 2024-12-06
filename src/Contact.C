//TODO:!!! Pending contact
//contactforce = 0
//In global initialize 
// domain max_contact_force
//#include "Mesh.cuh"

namespace MetFEM{
#define MAX_NB_COUNT    20

#include "Mesh.h"

#define HTOL 1.0e-6
#define DFAC 0.0

/*
//TODO: CHANGE TO SEVERAL CONTACT SURFACES 
inline void __device__ Domain_d::UpdateContactParticles(int mesh_id){
  
  int e = threadIdx.x + blockDim.x*blockIdx.x;	
  //printf("elemn %d\n", e);
  
  if (e < trimesh[mesh_id]->elemcount) {
    //printf("UPDATING e %d\n",e);
    //int e = element[i];
    double3 vv = make_double3(0.);
    for (int en = 0; en<3; en++){
      //printf("particle %d \n",i);
      //printf ("node %d \n",trimesh->elnode[3*e+en]);
      // if (trimesh->elnode[3*e+en] < trimesh->nodecount )
      vv += trimesh[mesh_id] -> node_v[trimesh[mesh_id]->elnode[3*e+en]];
      // else 
        // printf("error \n");
    }
    

    v [first_fem_particle_idx[mesh_id] + e] = vv/3.;
    a [first_fem_particle_idx[mesh_id] + e] = make_double3(0.);
    //printf(" particle %d , v %f %f %f \n", e, vv.x, vv.y, vv.z);
    // if (length(normal[e])<1.0e-3)
      // printf("UPDATING ERROR ZERO mesh normal, %f %f %f\n", trimesh -> normal[e].x , trimesh -> normal[e].y, trimesh -> normal[e].z);
    normal[first_fem_particle_idx[mesh_id] + e] = trimesh[mesh_id] -> normal[e];
    //printf("mesh normal, %f %f %f\n", trimesh -> normal[e].x , trimesh -> normal[e].y, trimesh -> normal[e].z);
  }
}

void __global__ CalcContactForcesKernel(Domain_d *dom_d,	const uint *particlenbcount,
																	const uint *neighborWriteOffsets,
																	const uint *neighbors
																	){
	dom_d->CalcContactForcesWang(
	particlenbcount,
	neighborWriteOffsets,
	neighbors);

}
*/

//WANG ET AL, FRICTION COEFFICIENT APPROACH    
//Inputs
//max_contact_force
//Neighbours
//vectors v
//////////////////////////////// 
//// From Wang: Simulating frictional contact in smoothed particle hydrodynamics
//// https://link.springer.com/article/10.1007/s11431-013-5262-x
//// Wang, Wu, GU, HUA, Science China 2013
////////////////////////////////
inline void dev_t Domain_d::CalcContactForcesWang(){
	//int i = threadIdx.x + blockDim.x*blockIdx.x;	
  
	double min_force_ts_=1000.;
// https://stackoverflow.com/questions/10850155/whats-the-difference-between-static-and-dynamic-schedule-in-openmp
			
  //max_contact_force = 0.;
	double min_contact_force = 1000.;
	int inside_pairs = 0;
  //printf("test\n");
  #ifdef BUILD_GPU
  par_loop(i,ext_nodes_count)
  #else
  for  (int i=0; i < ext_nodes_count;i++ )   //i particle is from SOLID domain, j are always rigid 
  #endif
  {
    
    /*
    contforce[i] = make_double3(0.); //RESET
    // CONTACT OFFSET IS FIX BY NOW
    int neibcount = contneib_count[i];
    
    int test = 0; //Should be once per nb

    for (int k=0;k < neibcount;k++) { //Or size
      int j = contneib_part[i*MAX_NB_COUNT+k];

      int mid = mesh_id[j];      
      
      double3 xij;
      double K;
      //int e = element[j]; //Index of mesh Element associated with node
      
      int e = j - first_fem_particle_idx[mid];

      double3 vr = v[i] - v[j];		//Fraser 3-137

      double3 x_pred = x[i] + v[i] * deltat + a[i] * deltat * deltat/2.;

      normal[j] = trimesh[mid]->normal[e];

      double dist = dot (normal[j],x_pred)  - trimesh[mid]->pplane[e];

      if (dist < h[i] ) {

          double3 Qj = x[i] - dist * normal[j];

          bool inside = true;
          int l=0,n;		   
          //printf("Entering while \n");
          while (l<3 && inside){
            n = l+1;	if (n>2) n = 0;
            // double crit = dot (cross ( *trimesh->node[e -> node[j]] - *trimesh->node[e -> node[i]],
                                                                // Qj  - *trimesh->node[e -> node[i]]),
                              // normal[j]);
            double crit = dot (cross ( trimesh[mid]->node[trimesh[mid]->elnode[3*e+n]] - trimesh[mid]->node[trimesh[mid]->elnode[3*e+l]],
                                                                Qj  - trimesh[mid]->node[trimesh[mid]->elnode[3*e+l]]),
                              normal[j]);
            if (crit < 0.0) inside = false;
            l++;
          }
          //printf("Outside while\n");
          
          if (inside ) { //Contact point inside element, contact proceeds
            // //Calculate penetration depth (Fraser 3-49)
            double delta = h[i] - dist;
            double delta_ = - dot( normal[j] , vr);	//Penetration rate, Fraser 3-138
            //printf("delta: %f\n", delta);
            // // DAMPING
            // //Calculate SPH and FEM elements stiffness (series)
            // //Since FEM is assumed as rigid, stiffness is simply the SPH one 
            double kij = 2.0 * m[i] / (deltat * deltat);


						double omega = sqrt (kij/m[i]);
						double psi_cont = 2. * m[i] * omega * DFAC; // Fraser Eqn 3-158

            //Normal Force
            //contforce[i] = (0.05*kij * delta - psi_cont * delta_) * normal[j]; // NORMAL DIRECTION, Fraser 3-159
            contforce[i] = 0.4 * kij * delta  * normal[j]; // NORMAL DIRECTION, Fraser 3-159
            
            //contforce[i].x = contforce[i].y = 0.0; ///// TO TEST BAD CONTACT
            a[i] += (contforce[i] / m[i]);
            //a[i].x = a[i].y = 0.0;
            // //NORMALS NOT RIGHT. IF REPLACING a[i] BY THIS IS OK
             // a[i].x = a[i].y = 0.0; ///// TO TEST BAD CONTACT
             // a[i].z = -1000;
           
            test++;
  
            if (friction_sta > 0.){
              double3 du = x_pred - x[i] - v[j] * deltat ;  
              //printf ("vj %f %f %f\n",v[j].x,v[j].y,v[j].z);
              double3 delta_tg = du - dot(du, normal[j])* normal[j];
              double3 tg_force = kij * delta_tg;
              
              //double dS = pow(m[i]/rho[i],0.33333); //Fraser 3-119
              if (length(tg_force) < friction_sta * length(contforce[i]) ){ //STATIC; NO SLIP
                a[i] -= tg_force / m[i];   
              } else {
                double3 tgforce_dyn = friction_dyn * length(contforce[i]) * tg_force/length(tg_force);
                contforce[i] -= tgforce_dyn;
                a[i] -= tgforce_dyn / m[i];
              }
            }

        }// if inside

      } //dist <h
    }//neibcount	for (int k=0;k < neibcount;k++) { //Or size
    
  */
  } //i<first fem index
	//Correct time step!
//	std::min(deltat,dt_fext)
} //Contact Forces

}; //Namespace
//};//SPH
