//TODO:!!! Pending contact
//contactforce = 0
//In global initialize 
// domain max_contact_force
//#include "Mesh.h"
#include "Domain_d.h"
#include "Mesh.h"

namespace MetFEM{
#define MAX_NB_COUNT    20



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
void dev_t Domain_d::CalcContactForcesWang(){
	//int i = threadIdx.x + blockDim.x*blockIdx.x;	
  
	double min_force_ts_=1000.;
// https://stackoverflow.com/questions/10850155/whats-the-difference-between-static-and-dynamic-schedule-in-openmp
			
  //max_contact_force = 0.;
	double min_contact_force = 1000.;
	int inside_pairs = 0;
  //printf("test\n");
  /*
  #ifdef BUILD_GPU
  par_loop(i,ext_nodes_count)
  #else
  for  (int i=0; i < ext_nodes_count;i++ )   //i particle is from SOLID domain, j are always rigid 
  #endif
  {
  */
  
  for  (int i=0; i < m_dim*m_node_count;i++ ) 
    contforce[i]=0.0;
  
  #ifdef BUILD_GPU
  par_loop(i,m_node_count)
  #else
  for  (int i=0; i < m_node_count;i++ )   //i particle is from SOLID domain, j are always rigid 
  #endif
  {
    //printf("Node %d\n",i);
    //printf("Element count %d\n", trimesh->elemcount);
    double min_dist = 1e10;
    double delta;

    if (ext_nodes[i]){ 
      //Search NEGATIVE DISTANCE node ON MESH
      int j=0;
      bool end = false;
      while (!end){//----------------------------------------------------------
      //CHCK isNodeinElement WITH getShapeFunctionAtPoint
      
      //AND isNodeinSideFace (redYnEla sideface.C)
      
      //for (int j=0;j<trimesh->elemcount;j++){
        //FROM WeldForm Contact Wang
        //dist =  dot (Particles[P2]->normal, x_pred ) - trimesh[m]-> element[Particles[P2]->element] -> pplane;
        //qj = Particles[P1]->x - dist * Particles[P2]->normal;
        
        //Original (wrong)
        //double3 dist = getNodePos3(i) - trimesh->centroid[j];
        //double d = norm2(dist);
        
        //IN wang is predicted
        double d = dot(trimesh->normal[j],getNodePos3(i))  - trimesh->pplane[j];
        //double d = norm2(dist);
        //Original (wrong
        //delta = dot(dist,trimesh->normal[j]);
        delta = dot(getNodePos3(i) - trimesh->centroid[j],trimesh->normal[j]);
        
        //ACCORDING TO WANG
			//double3 vr = Particles[P1]->v - Particles[P2]->v;		//Fraser 3-137
			//delta = - dot(trimesh->normal[j], vr);	//Penetration rate, Fraser 3-138        
        
        if (delta <0 /*&& dist < CERTAIN ELEMENT DISTANCE*/){
          //printf ("ELEMENT %d DELTA <0------\n", e);
          
          //ORIGINAL
          double3 Qj = getPosVec3(i) - d * trimesh->normal[j];
          if (i==254)
            printf("node 254 qj %f %f %f \n", Qj.x,Qj.y,Qj.z);

          bool inside = true;
          int l=0,n;		   
          //printf("Entering while \n");
          while (l<3 && inside){
            n = l+1;	if (n>2) n = 0;
            //double crit;
            double crit = dot (cross ( trimesh->node[trimesh->elnode[3*j+n]] - trimesh->node[trimesh->elnode[3*j+l]],
                                                                Qj  - trimesh->node[trimesh->elnode[3*j+l]]),
                               trimesh->normal[j]);
            if (crit < 0.0) inside = false;
            l++;
          }
          // if (i==254 && !inside)
            // printf("Node 254 not inside, Normal distance %f\n",d);
          // if (i==220)
            // printf("Node 220 inside, Normal distance %f\n",d);
          
            if (inside ){
              //printf("delta: %.3e\n",delta);
              //printf("dist %f %f %f\n",dist.x,dist.y,dist.z);
              //printf("Node: %d, Mesh Element %d INSIDE!--------------------\n",i, e);
              
                //fn = 2. * node->mass * delta / SQ(timeStep);
    //REDYNELAfiniteELement/contact.C
  //computetangentialForce(fn, Ft);
              double3 cf =  - 2.0 * m_mdiag[i] * delta * trimesh->normal[j]/(dt*dt);
              printf("Node %d CF %f %f %f\n",i, cf.x,cf.y,cf.z);
              contforce[m_dim*i] = cf.x;contforce[m_dim*i+1] = cf.y;contforce[m_dim*i+2] = cf.z;
              
              end = true;//JUST ONE MASTER ELEMENT PER SLAVE NODE
            }
          
          
          
          
        }
        j++;
        if (j==trimesh->elemcount)
          end = true;
      } // WHILE j mesh elements
     //printf("MinDist Node%d %.3e, on element %d\n", i, min_dist, e );
        
         
    
    //contforce[i] = make_double3(0.,0.,0.); //RESET
    // CONTACT OFFSET IS FIX BY NOW

      

/*
      double3 vr = v[i] - v[j];		//Fraser 3-137

      double3 x_pred = x[i] + v[i] * deltat + a[i] * deltat * deltat/2.;

      normal[j] = trimesh->normal[e];

      double dist = dot (normal[j],x_pred)  - trimesh->pplane[e];

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
            double crit = dot (cross ( trimesh->node[trimesh->elnode[3*e+n]] - trimesh->node[trimesh->elnode[3*e+l]],
                                                                Qj  - trimesh->node[trimesh->elnode[3*e+l]]),
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
  
    }//external nodes
  } //i<first fem index
	//Correct time step!
//	std::min(deltat,dt_fext)
} //Contact Forces


void dev_t Domain_d::CalcContactForces(){
	//int i = threadIdx.x + blockDim.x*blockIdx.x;	
  
	double min_force_ts_=1000.;
// https://stackoverflow.com/questions/10850155/whats-the-difference-between-static-and-dynamic-schedule-in-openmp
			
  //max_contact_force = 0.;
	double min_contact_force = 1000.;
	int inside_pairs = 0;
  //printf("test\n");
  /*
  #ifdef BUILD_GPU
  par_loop(i,ext_nodes_count)
  #else
  for  (int i=0; i < ext_nodes_count;i++ )   //i particle is from SOLID domain, j are always rigid 
  #endif
  {
  */
  
  for  (int i=0; i < m_dim*m_node_count;i++ ) 
    contforce[i]=0.0;
  int max_cf = 0.0;
  #ifdef BUILD_GPU
  par_loop(i,m_node_count)
  #else
  for  (int i=0; i < m_node_count;i++ )   //i particle is from SOLID domain, j are always rigid 
  #endif
  {
    //printf("Node %d\n",i);
    //printf("Element count %d\n", trimesh->elemcount);
    double min_dist = 1e10;
    double delta;
    
    if (ext_nodes[i]){ 
      //Search NEGATIVE DISTANCE node ON MESH
      int j=0;
      bool end = false;
      while (!end){//----------------------------------------------------------
      //CHCK isNodeinElement WITH getShapeFunctionAtPoint
        //printf("Node i %d Element  j %d\n", i,j);      

        double d = dot(trimesh->normal[j],getNodePos3(i))  - trimesh->pplane[j];
        //double d = norm2(dist);
        //Original (wrong
        //delta = dot(dist,trimesh->normal[j]);
        //delta = dot(getNodePos3(i) - trimesh->centroid[j],trimesh->normal[j]);
        
        //ACCORDING TO WANG
			//double3 vr = Particles[P1]->v - Particles[P2]->v;		//Fraser 3-137
			//delta = - dot(trimesh->normal[j], vr);	//Penetration rate, Fraser 3-138        
        //printf("Node %d d %f\n",i, d);
        if (d < 0 /*&& dist < CERTAIN ELEMENT DISTANCE*/){
          //printf ("ELEMENT %d DELTA <0------\n", e);
          
          //ORIGINAL
          double3 Qj = getPosVec3(i) - d * trimesh->normal[j];
          // if (i==254)
            // printf("node 254 qj %f %f %f \n", Qj.x,Qj.y,Qj.z);

          bool inside = true;
          int l=0,n;		   
          //printf("Entering while \n");
          while (l<3 && inside){
            n = l+1;	if (n>2) n = 0;
            //double crit;
            double crit = dot (cross ( trimesh->node[trimesh->elnode[3*j+n]] - trimesh->node[trimesh->elnode[3*j+l]],
                                                                Qj  - trimesh->node[trimesh->elnode[3*j+l]]),
                               trimesh->normal[j]);
            if (crit < 0.0) inside = false;
            l++;
          }
          // if (i==254 && !inside)
            // printf("Node 254 not inside, Normal distance %f\n",d);
          // if (i==220)
            // printf("Node 220 inside, Normal distance %f\n",d);
          
            if (inside ){
              //printf("delta: %.3e\n",delta);
              //printf("dist %f %f %f\n",dist.x,dist.y,dist.z);
              //printf("Node: %d, Mesh Element %d INSIDE!--------------------\n",i, e);
              
                //fn = 2. * node->mass * delta / SQ(timeStep);
    //REDYNELAfiniteELement/contact.C
  //computetangentialForce(fn, Ft);
              //fn = 2. * node->mass * delta / SQ(timeStep);
              //double3 cf =  - 2.0 * m_mdiag[i] * d * trimesh->normal[j]/(dt*dt);
              double3 cf =  - 0.2 * m_mdiag[i] * d * trimesh->normal[j]/(dt*dt);
              //printf("Node %d CF %f %f %f, dist %f mass %f\n",i, cf.x,cf.y,cf.z, d,m_mdiag[i]);
              contforce[m_dim*i] = cf.x;contforce[m_dim*i+1] = cf.y;contforce[m_dim*i+2] = cf.z;
              
              ////FRICTION 
              //1. Calculare slave nodal tg vel 
              //v tan=vs −dot(vs,normal)⋅normal
              double3 vtan = getVelVec(i) - dot(getVelVec(i),trimesh->normal[j]) * trimesh->normal[j];
              //2. Compute Friction Force Magnitude
              //ft = -mu |fn| vtan/(|vtan| + eps)
              double3 ft = - 0.4 * norm(cf) * vtan/(norm(vtan)+1.0e-5);
              //FROM REDYNELA
              //  // Force tangentielle
              // Ft = -(node->mass / Global_Structure->domains.current()->/*times.timeStep*/ currentSolver->getTimeStep()) * Vt;
              contforce[m_dim*i] += ft.x;contforce[m_dim*i+1] += ft.y;contforce[m_dim*i+2] += ft.z;
              
              end = true;//JUST ONE MASTER ELEMENT PER SLAVE NODE
            }
          
          
          
          
        }

        j++;
        if (j==trimesh->elemcount)
          end = true;
      } // WHILE j mesh elements
     //printf("MinDist Node%d %.3e, on element %d\n", i, min_dist, e );
        
         
    
    //contforce[i] = make_double3(0.,0.,0.); //RESET
    // CONTACT OFFSET IS FIX BY NOW

      

/*
      double3 vr = v[i] - v[j];		//Fraser 3-137

      double3 x_pred = x[i] + v[i] * deltat + a[i] * deltat * deltat/2.;

      normal[j] = trimesh->normal[e];

      double dist = dot (normal[j],x_pred)  - trimesh->pplane[e];

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
            double crit = dot (cross ( trimesh->node[trimesh->elnode[3*e+n]] - trimesh->node[trimesh->elnode[3*e+l]],
                                                                Qj  - trimesh->node[trimesh->elnode[3*e+l]]),
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
  
    }//external nodes
  } //i<first fem index
	//Correct time step!
//	std::min(deltat,dt_fext)
} //Contact Forces

}; //Namespace
//};//SPH
