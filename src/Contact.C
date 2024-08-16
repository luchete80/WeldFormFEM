//TODO:!!! Pending contact
//contactforce = 0
//In global initialize 
// domain max_contact_force
//#include "Mesh.cuh"

#define MAX_NB_COUNT    20

#include "Mesh.cuh"

#define HTOL 1.0e-6
#define DFAC 0.0

//namespace SPH{
__global__ inline void CalcContactNbKernel(Domain_d *dom_d,
const uint *particlenbcount,
const uint *neighborWriteOffsets,
const uint *neighbors) {
 	
  dom_d->CalcContactNb(
	particlenbcount,
	neighborWriteOffsets,
	neighbors); 
  
}

//Up until now max nb count is 20, this is done in order to fix contneib_part
//and to not reallocate dynamically each nb search
__device__ inline void Domain_d::CalcContactNb(const uint *particlenbcount,
const uint *neighborWriteOffsets,
const uint *neighbors){
  //printf("searching nbs less than index %d\n",first_fem_particle_idx[0] );
  int i = threadIdx.x + blockDim.x*blockIdx.x;	
  if (i < first_fem_particle_idx[0] ) {
    contneib_count[i] = 0;
    int neibcount = particlenbcount[i];
    const uint writeOffset = neighborWriteOffsets[i];
    
    for (int k=0;k < neibcount;k++) { //Or size
      int j = neighbors[writeOffset + k];
      double3 xij = x[i]-x[j];
      //double h = h[i] + h[j];  //not necessary to take average
      for (int mc=0;mc<trimesh_count;mc++){
        if ( (ID[i] == id_free_surf && ID[j] == contact_surf_id[mc]) /*||
             (ID[j] == id_free_surf && ID[i] == contact_surf_id) */) {
          /////if ( norm (Particles[P1]->x - Particles[P2]->x) < ( Particles[P1]->h + Particles[P2]->h ) ){ //2* cutoff being cutoff (h1+h2)/2
          if ( length(xij) < 2.*h[i] ){   
            contneib_part[MAX_NB_COUNT*i + contneib_count[i]] = j;
            contneib_count[i]++;
            if (contneib_count[i]>MAX_NB_COUNT) printf("ERROR, MAX CONTACT NB COUNT REACHED, %d counted\n",contneib_count);
          }
        } //IDs OK
      } //mesh idx
    }//for k neighbours
  }// i < fem index
}

inline void __global__ UpdateContactParticlesKernel(Domain_d *dom, int mesh_id){
  dom->UpdateContactParticles(mesh_id);
}

//inline void __device__ Domain_d::UpdateContactParticles(){
  // for (int m=0; m<trimesh.size();m++){
    // for ( int e = 0; e < trimesh[mid]->element.Size(); e++ ){
      // Vec3_t v = 0.;
      // for (int en = 0;en<3;en++)
        // v += *trimesh[mid] -> node_v[trimesh[mid]->element[e] ->node[en]];
      // Particles[first_fem_particle_idx[m] + e] -> v = 
      // Particles[first_fem_particle_idx[m] + e] -> va = 
      // Particles[first_fem_particle_idx[m] + e] -> vb = v/3.;
      // Particles[first_fem_particle_idx[m] + e] -> a = 0.; 
      // Particles[first_fem_particle_idx[m] + e] -> normal  = trimesh[mid]->element[e] -> normal;
      // //cout << "v "<< v/3.<<", n "<<Particles[first_fem_particle_idx + e] -> normal<<endl;
    // } 
  // }
//}

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
																	/*const int &id, const double &totmass*/){
	dom_d->CalcContactForcesWang(
	particlenbcount,
	neighborWriteOffsets,
	neighbors);

}
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
void __device__ inline Domain_d::CalcContactForcesWang(const uint *particlenbcount,
                                                  const uint *neighborWriteOffsets,
                                                  const uint *neighbors){
	int i = threadIdx.x + blockDim.x*blockIdx.x;	

	double min_force_ts_=1000.;
// https://stackoverflow.com/questions/10850155/whats-the-difference-between-static-and-dynamic-schedule-in-openmp
			
	max_contact_force = 0.;
	double min_contact_force = 1000.;
	int inside_pairs = 0;
  //printf("test\n");
  if (i < first_fem_particle_idx[0] ) {  //i particle is from SOLID domain, j are always rigid 
    
    contforce[i] = make_double3(0.); //RESET
    // CONTACT OFFSET IS FIX BY NOW
    int neibcount = contneib_count[i];
    
    int test = 0; //Should be once per nb
    //printf("i, first fem part, neibcount %d\n",neibcount);
    // printf("Nb indexed,i:%d\n",i);
    // In this Weldform GPU version, is clear than i is SPH particle and 2 is RIGID PARTICLE
    //if (i==11301) {printf("Partcile 11301, nb count %d\n",neibcount);}
    for (int k=0;k < neibcount;k++) { //Or size
      int j = contneib_part[i*MAX_NB_COUNT+k];
      //printf("j neib %d: \n",j);
      //int j = neighbors[writeOffset + k];
      // if (ID[j] != contact_surf_id)
        // printf("CONTACT, j particle, %d on CONTACT SURFACE\n", i);
      
      int mid = mesh_id[j];      
      
      double3 xij;
      double K;
      //int e = element[j]; //Index of mesh Element associated with node
      
      //TODO: MAKE ONE PER MESH
      int e = j - first_fem_particle_idx[mid];
      //printf ("Element e %d\n", e);
      // Summing the smoothed pressure, velocity and stress for fixed particles from neighbour particles
      //IT IS CONVENIENT TO FIX SINCE FSMPairs are significantly smaller
      //cout << "Contact pair size: "<<ContPairs[k].Size()<<endl;

      double3 vr = v[i] - v[j];		//Fraser 3-137
      //printf("particle 2 vr %f %f %f\n",v[j].x, v[j].y,v[j].z);
      // if (k==0)
        // printf("particle 2 (%d) x %f %f %f\n",j, x[j].x, x[j].y, x[j].z);
      //cout << "Particle P1v: "<<v[i]<<endl;
      //cout << "Particle P2v: "<<Particles[P2]->v<<endl;
      
      //ok, FEM Particles normals can be calculated by two ways, the one used to
      //calculate SPH ones, and could be given by mesh input
      //delta_ Is the projection of relative velocity 
      // double delta_ = - dot( normal[j] , vr);	//Penetration rate, Fraser 3-138
      
      double3 x_pred = x[i] + v[i] * deltat + a[i] * deltat * deltat/2.;
      //printf("xpred %f %f %f\n", x_pred.x, x_pred.y, x_pred.z);
      //printf ("pplane: %f", trimesh->pplane[e]);
      
      normal[j] = trimesh[mid]->normal[e];
      // normal[j] = make_double3(0.0,0.0,-1.0);
      
      // if (length(normal[j])<1.0e-3)
        // printf("NORMAL CALC ERROR in particle %d. ZERO, Normal j %.6e %.6e %.6e\nNormal e %.6e %.6e %.6e\n",j,
      // normal[j].x,normal[j].y,normal[j].z,
      // trimesh->normal[e].x,trimesh->normal[e].y,trimesh->normal[e].z);
      //normal[j] = make_double3(0.,0.,-1.);
      
      double dist = dot (normal[j],x_pred)  - trimesh[mid]->pplane[e];
      //double dist = dot (normal[j],x[i])  - trimesh->pplane[e];
      //printf("normal j %d %f %f %f\n", j, normal[j].x, normal[j].y, normal[j].z);
      //printf("OUTSIDE part %d pplane %f,dist %.5e, h %f\n", i, trimesh->pplane[e], dist, h[i]);
      if (dist < h[i] ) {
        //printf ("INSIDE part %d pplane %f, dist: %f, h %f \n", i,trimesh->pplane[e],dist, h[i]);
        //printf("INSIDE part %d dist %\n", i,dist);
        
        // double deltat_cont = ( h[i] + pplane - dot (normal[j],	x[i]) ) / (-delta_);								//Eq 3-142 
        // double3 Ri = x[i] + deltat_cont * vr;	//Eq 3-139 Ray from SPH particle in the rel velocity direction
        // //printf("delta %f \n",delta_);
        
        //cont_forces[i] = make_double3(0.,0.,0.); //RESET

          //cout << "Inside dt contact" <<endl;
          /////Find point of contact Qj
          /////This is in Fraser algorithm
          double3 Qj = x[i] - dist * normal[j];
          // //Check if it is inside triangular element
          // //Find a vector 
          // //Fraser 3-147
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
            //printf("deltat, kij %f %f\n", deltat, kij); 
            //double omega = sqrt (kij/m[i]);
            //double psi_cont = kij * delta; // Fraser Eqn 3-158
            
            //printf("Normal j %f %f %f \n", normal[j].x,normal[j].y,normal[j].z);
                      
            // // TANGENTIAL COMPONENNT
            // // Fraser Eqn 3-167
            // // TODO - recalculate vr here too!
            // double3 tgvr = vr + delta_ * normal[j];  // -dot(vr,normal) * normal
            // double3 tgdir = tgvr / length(tgvr);
            

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
            // printf ("CONTACT!! part %d, pos %.3e %.3e %.6e\n",i,x[i].x,x[i].y,x[i].z);
            // printf ("INSIDE part %d pplane %f, dist: %f, h %f \n", i,trimesh->pplane[e],dist, h[i]);
            
          // if (abs(contforce[i].x) > 1.0e-10 || abs(contforce[i].y) > 1.0e-10)
            // printf("NORMAL CALC ERROR in particle %d. ZERO, Normal j %.6e %.6e %.6e\nNormal e %.6e %.6e %.6e\n",j,
          // normal[j].x,normal[j].y,normal[j].z,
          // trimesh->normal[e].x,trimesh->normal[e].y,trimesh->normal[e].z);
            //printf("step: %d Particle %i, x_pred %f %f %f, dist %f h %f delta %.4e pplane %f kij %f contforce %f %f %f \n", step, i, x_pred.x, x_pred.y,x_pred.z,dist, h[i],delta, trimesh->pplane[e], kij, contforce[i].x, contforce[i].y, contforce[i].z);
            // if (abs(contforce[i].x)>1.0e-3 || abs(contforce[i].y)>1.0e-3 ){
              // printf("CONTACT FORCE x != 0!!!\n");
              // printf("Normal j %f %f %f \n", normal[j].x,normal[j].y,normal[j].z);
            // }
            //printf("dist %f, delta: %.4e, h %f kij %f abscf %.4e \n", dist,delta, h[i], kij, length(contforce[i]));
            //if (length(contforce[i])>0.0 && length(contforce[i]) < 30.){
            //if (i==11301) {
            //printf("step: %d Particle %i, x_pred %f %f %f, dist %f h %f delta %.4e pplane %f kij %f contforce %f %f %f \n", step, i, x_pred.x, x_pred.y,x_pred.z,dist, h[i],delta, trimesh->pplane[e], kij, contforce[i].x, contforce[i].y, contforce[i].z);
            //}
            // TANGENTIAL FORCE //////    
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
    
    //if(test == 0 && contneib_count[i] > 0) printf ("NO CONTACT FORCE APPLIED part %d, pos %.3e %.3e %.6e\n",i,x[i].x,x[i].y,x[i].z);
    //if(test > 1)  printf ("ERROR. MORE THAT ONCE PER PARTICLE");  
    //else if (test == 0)printf ("ERROR. NOT PARTICLES");  
  } //i<first fem index
	//Correct time step!
//	std::min(deltat,dt_fext)
} //Contact Forces

//};//SPH
