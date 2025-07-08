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


void dev_t Domain_d::CalcContactForces(){
  
  #ifndef CUDA_BUILD
	//int i = threadIdx.x + blockDim.x*blockIdx.x;	
  double pxa[m_node_count];
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
  //~ par_loop(i,m_dim*m_node_count)
    //~ contforce[i]=0.0;
  par_loop(i,m_node_count){
    m_mesh_in_contact[i]=-1;
    pxa[i]=0.0;
  }
  par_loop(i,m_node_count){
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

        double3 x_pred = getPosVec3(i) + getVelVec(i) * dt + getAccVec(i) * dt * dt/2.0;
        
        double dalt =  dot (trimesh->normal[j], x_pred ) - trimesh-> pplane[j];


        //double d = norm2(dist);
        //Original (wrong
        //delta = dot(dist,trimesh->normal[j]);
        //delta = dot(getNodePos3(i) - trimesh->centroid[j],trimesh->normal[j]);

        //~ x_pred = Particles[P1]->x + Particles[P1]->v * deltat + Particles[P1]->a * deltat * deltat/2.0;
        //~ vr_pred = Particles[P1]->v + Particles[P1]->a * deltat - Particles[P2]->v;
        
        //~ dist =  dot (Particles[P2]->normal, x_pred ) - trimesh[m]-> element[Particles[P2]->element] -> pplane;
        
        //ACCORDING TO WANG
			//double3 vr = Particles[P1]->v - Particles[P2]->v;		//Fraser 3-137
			//delta = - dot(trimesh->normal[j], vr);	//Penetration rate, Fraser 3-138        
        //printf("Node %d d %f\n",i, d);
        if (d < 0 /*&& dist < CERTAIN ELEMENT DISTANCE*/){
          //printf ("ELEMENT %d DELTA <0------\n", e);
         
          //ORIGINAL
          double3 Qj = getPosVec3(i) - d * trimesh->normal[j];
          //~ if (i==40)
            //~ printf("Element %d node 40 qj %f %f %f \n", j,Qj.x,Qj.y,Qj.z);

          bool inside = true;
          int l=0,n;		   
          //printf("Entering while \n");
          while (l<3 && inside){
            n = l+1;	if (n>2) n = 0;
            //double crit;
            double crit = dot (cross ( trimesh->node[trimesh->elnode[3*j+n]] - trimesh->node[trimesh->elnode[3*j+l]],
                                                                Qj  - trimesh->node[trimesh->elnode[3*j+l]]),
                               trimesh->normal[j]);
            double3 test = cross ( trimesh->node[trimesh->elnode[3*j+n]] - trimesh->node[trimesh->elnode[3*j+l]],
                                                                Qj  - trimesh->node[trimesh->elnode[3*j+l]]);
                                      
            if (crit < 0.0) inside = false;
            l++;
          }


  if (inside ){
 
              double nodlen = 0.0;
              double nodevol = 0.0;
              // TODO: MAKE FUNCTION FOR THIS
              for (int e=0; e<m_nodel_count[i];e++) { 
                int eglob   = m_nodel     [m_nodel_offset[i]+e]; //Element
                nodlen +=m_elem_length[eglob];
                nodevol+=vol[eglob];
              }
              nodlen /= m_nodel_count[i]; 
              nodevol /=m_nodel_count[i];

              
              
              double beta_d = 0.1;
              //printf("NODLEN %.4e\n", nodevol);
              //double kcont = 0.8 * mat[0]->Elastic().E() * node_area[i]/nodlen;       
              //double kcont = 0.2 * mat[0]->Elastic().E() * pow(nodevol,0.3333);           
              //double kcont = 0.2 * m_mdiag[i] /(dt*dt);
                            
              // THIS IS USED ALSO FOR CONTACT DAMPING
              double3 v_rel = getVelVec(i); // master assumed static or you can use relative motion     
              double v_reln = dot(v_rel, trimesh->normal[j]);         
              //double F_damp = 0.0;

              //double kcont_geo = 4.0*mat[0]->Elastic().E() * node_area[i] / nodlen;
              //double kcont_geo = mat[0]->Elastic().BulkMod()*nodevol;
              double kcont_mass = 0.2 * m_mdiag[i] / (dt * dt);
              //double kcont = std::min(kcont_geo, kcont_mass);
              double kcont = kcont_mass;
              double damping_ratio = 0.01;
              double omega = sqrt(kcont / m_mdiag[i]);
              double ccrit = 2.0 * m_mdiag[i] * omega;
              double F_damp = damping_ratio * ccrit * v_reln;

              double F_normal = m_contPF * kcont * d;

              //double F_damp = beta_d * sqrt(kcont * m_mdiag[i]) * v_reln;
              //printf("KCON geo %f KCON mass %f: \n", kcont, kcont2;
              double3 cf =  - (F_normal  + F_damp) * trimesh->normal[j];
              //~ if (getPosVec3(i).z<0.0001)
                //~ printf("Node %d CF %f %f %f, dist %f mass %f\n",i, cf.x,cf.y,cf.z, d,m_mdiag[i]);
              contforce[m_dim*i] = cf.x;contforce[m_dim*i+1] = cf.y;contforce[m_dim*i+2] = cf.z;
              pxa[i] = p_node[i] * node_area[i];
              //printf("Nodal pressure %.4e\n",p_node[i]);
              //printf("Cont Force %f %f %f \n",cf.x,cf.y,cf.z);
              
              //printf("MESHIN CONTACT %d\n",trimesh->ele_mesh_id[j]);
              m_mesh_in_contact[i]=trimesh->ele_mesh_id[j];
              //m_mesh_in_contact[i]=0;
              //~ if (m_mesh_in_contact[i]>trimesh->mesh_count)
                //~ printf("ERROR!!! \n");
              ////

              double damping = 1.0;  // Tune this (0.01–0.5)
              ////////WANG 2013     
              double3 du = x_pred - getPosVec3(i) - getVelVec(i) * dt ;
              double3 delta_tg = du - dot(du, trimesh->normal[j])*trimesh->normal[j];
              double3 tgforce = kcont * delta_tg;

              
              
               // 1. Relative velocity at contact in tangential direction
              //double3 v_rel = getVelVec(i); // master assumed static or you can use relative motion
              double3 v_tan = v_rel - dot(v_rel, trimesh->normal[j]) * trimesh->normal[j];

              // 2. Incremental tangential displacement (Δut = vt * dt)
              double3 du_tangent = v_tan * dt;

              // 3. Accumulate tangential displacement
              double3 ut_acc = make_vector_t(ut_prev[3*i],ut_prev[3*i+1],ut_prev[3*i+2]);
              ut_acc = ut_acc + du_tangent;
              ut_prev[m_dim*i] += du_tangent.x; ut_prev[m_dim*i+1] += du_tangent.y;ut_prev[m_dim*i+2] += du_tangent.z;
              // 4. Trial tangential force from accumulated displacement
              double3 Ft_trial = -kcont * ut_acc;

              // 5. Coulomb cap

              double Ft_mag = length(Ft_trial);
              double3 Fn = dot(cf,trimesh->normal[j])*trimesh->normal[j];
              double Ft_max = trimesh->mu_sta[0] * norm(Fn); // norm(Fn) is magnitude of normal force
              

              double3 Ft;

              double Ft_max_static = trimesh->mu_sta[0] * norm(Fn);
              if (Ft_mag <= Ft_max_static) {
                  // 6. Cap force, rescale
                  // Still sticking: Use trial force
                   //~ if (i==0){
                      //~ printf("STATIC FORCE %f %f %f\n",tgforce.x,tgforce.y,tgforce.z);
                     //~ printf("DISP NODE 0 %f %f %f\n", u[3*i],u[3*i+1],u[3*i+2]);
                     //~ print
                     //~ }
                  Ft = Ft_trial;
              }
              else {
                  double Ft_max_dynamic =  trimesh->mu_dyn[0] * norm(Fn);
                  Ft = -Ft_max_dynamic * v_tan/norm(v_tan);  // Use velocity direction now
                  for (int d=0;d<3;d++)ut_prev[m_dim*i+d] = 0;  // Optional: reset accumulation when sliding
              }              


              // 7. Apply friction force
              contforce[m_dim * i + 0] += Ft.x;
              contforce[m_dim * i + 1] += Ft.y;
              contforce[m_dim * i + 2] += Ft.z;                     


              q_cont_conv[i] = trimesh->heat_cond * node_area[i]*(20.0-T[i]);

              
              end = true;//JUST ONE MASTER ELEMENT PER SLAVE NODE
            
            } //INSIDE
                    
          
          
          
          
          
        }

        j++;
        if (j==trimesh->elemcount)
          end = true;
      }
  
    }//external nodes
  } //i<first fem index

   
    
  #endif 
	//Correct time step!
//	std::min(deltat,dt_fext)
} //Contact Forces

void Domain_d::calcContactForceFromPressure(){
  
   bool is_elem_sum[m_elem_count];
   bool is_node_sum[m_node_count];
   //double pxa_el[m_elem_count];
   double cfsum=0.0;
   double area = 0.0;
    double cfnsum = 0.0;
    
   for (int e=0;e<m_elem_count;e++)is_elem_sum[e]=false;
   
    for (int i=0;i<m_node_count;i++){
      if(m_mesh_in_contact[i]>-1 && m_mesh_in_contact[i]<1){
        cfnsum += p_node[i]*node_area[i];
        for (int ne=0; ne<m_nodel_count[i];ne++) {
          int e   = m_nodel     [m_nodel_offset[i]+ne]; //Element
          if (!is_elem_sum[e]){
            //pxa_el[e]+=p[e]*m_elem_area[e];
            is_elem_sum[e]=true;
            cfsum += p[e]*m_elem_area[e];
            area+=m_elem_area[e];
          }
            //~ if (!is_node_sum[i]){
              //~ area+=node_area[i];
              //~ }
        }//nodel
      }//mesh in contact
    }
    
  
  printf("Area %.3e\n", area);
  trimesh->react_p_force[0] = cfsum;
  trimesh->react_force[0].z = cfnsum;
   trimesh->cont_area = area;
}

}; //Namespace
//};//SPH
