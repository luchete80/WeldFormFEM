#include <stdio.h>
#include <math.h>
#include <iostream>

#include <string.h> // memcpy
#include "defs.h"

#include "Matrix_temp.h"
#include  "Domain_d.h"

#include "tensor.cuh"

using namespace std;

//// IF YOU WANYT 2D
//#define BIDIM 


#ifdef BIDIM
#define m_dim 2
#define m_nodxelem 4
#define m_gp_count 1

#else

#define m_dim 3
#define m_nodxelem 8
#define m_gp_count 1

#endif

double gauss_points[1][3] = {{0,0}};
double gauss_weights[m_gp_count] = {8.};

//double tf = 1.0e-3;

double str_rate[m_gp_count][3][3];
double rot_rate[m_gp_count][3][3];
double strain[m_gp_count][3][3];
Matrix tau_;
double tau[m_gp_count][3][3];

double stress[m_gp_count][3][3];


using namespace MetFEM;
class domain_test:
public Domain_d {
public:

double mat_G;
double K_mod;
int red_int = 0;
double element_length = 1.0; // Length of the element
int axi_symm = 0;            //FALSE: PLAIN STRAIN
  

tensor3 Sigma_tst;

void impose_bc() {
#ifdef BIDIM

    // vel[2][1] = vel[3][1] = -1.0;
    // vel[0][0] = vel[0][1] = vel[1][1] = 0.0;

    // accel[2][1] = accel[3][1] = 0.0;
    // accel[0][0] = accel[0][1] = accel[1][1] = 0.0;

#else
    //LOCAL CASE 
    //ELEMENT NODES; NOT GLOBAL
    
    v  [m_dim*0+0] = v  [m_dim*0+1] = v  [m_dim*0+2] = 0.0;
    a[m_dim*0+0] = a[m_dim*0+1] = a[m_dim*0+2] = 0.0;
    
    v[m_dim*1+1] = 0.0;     a[m_dim*1+1] = 0.0;
    v[m_dim*1+2] = 0.0;     a[m_dim*1+2] = 0.0;
    
    v[m_dim*2+2] = a[m_dim*2+2] = 0.0;
    
    v[m_dim*3+0] = a[m_dim*3+0] =0.0;
    v[m_dim*3+2] = 0.0; a[m_dim*3+2] = 0.0;
    
    
    for (int i=4;i<8;i++) {
      v[m_dim*i+2] = -1.0;
      a[m_dim*i+2] = 0.0;}

#endif
}



double calc_vol() {
    double vol = 0.0;
    for (int gp = 0; gp < m_gp_count; gp++) {
        vol += m_detJ[gp] * gauss_weights[gp];
    }
    return vol;
}

void calc_strain(double str_rate[m_gp_count][3][3], double dt, double strain[m_gp_count][3][3]) {

    for (int gp = 0; gp < m_gp_count; gp++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                strain[gp][i][j] = dt * str_rate[gp][i][j];
                printf("strain %e",strain[gp][i][j]);
            }
          printf("\n");
        }
        
    }
}

void calc_pressure(double K_, double dstr[m_gp_count][3][3], double stress[m_gp_count][3][3]) {
    double pi_ = 0.0;
    for (int gp = 0; gp < m_gp_count; gp++) {
        for (int i = 0; i < 3; i++) {
            pi_ += dstr[gp][i][i];
        }
    }
    printf ("K %e",K_);
    pi_ = -pi_ / m_gp_count;
    for (int gp = 0; gp < m_gp_count; gp++) {
        p[gp] = -1.0 / 3.0 * (stress[gp][0][0] + stress[gp][1][1] + stress[gp][2][2]) + K_ * pi_;
        //printf("pres %e ",pres[gp]);
    }
    
}


void calc_hg_forces(double rho, double vol, double cs,double *fhg){

#ifdef BIDIM
  double Sig [4][4] = {{1.,-1.,1.,-1.}};
  double hmod[m_dim][4]={{0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0}};
  int jmax = 1;
#else
  double Sig [4][8] = {{ 1., 1.,-1.,-1.,-1.,-1., 1., 1.}, 
                       { 1.,-1.,-1., 1.,-1., 1., 1.,-1.},
                       { 1.,-1., 1.,-1., 1.,-1., 1.,-1.}, 
                       {-1., 1.,-1., 1., 1.,-1., 1.,-1.}};
  double hmod[m_dim][4]={{0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0}};
  int jmax = 8;
#endif

  for (int i = 0; i < m_nodxelem; i++) 
    for (int d = 0; d < m_dim; d++) 
      fhg[i*m_dim+d] = 0.0;
      
  for (int j = 0; j < jmax; j++) 
    for (int i = 0; i < m_nodxelem; i++) 
      for (int d = 0; d < m_dim; d++) {
        hmod[d][j] +=v[i*m_dim+d]*Sig[j][i];
        //printf("hmod: %6e", hmod[d][j]);
      }

  double ch = 0.1 * pow (vol,0.6666666) * rho * 0.25 * cs;
  
      for (int j = 0; j < jmax; j++) 
        for (int i = 0; i < m_nodxelem; i++) 
          for (int d = 0; d < m_dim; d++) 
            fhg[i*m_dim+d] -=ch*hmod[d][j]*Sig[j][i];
  

  
  // for gp in range(len(gauss_points)):
    // for j in range(jmax):
      // for n in range(m_nodxelem):
        // hmod[:,j] +=v[n,:]*Sig[j,n]
  // for j in range(jmax):
    // for n in range(m_nodxelem):
      // f_[n,:] -=hmod[:,j]*Sig[j,n] 
  // ch = 0.06 * pow (vol,0.66666) * rho * 0.25 * cs

  // f_ *= ch
}


  void Solve() {
    double t = 0.0;
      dt = 0.8e-5;
    //double tf = 0.8e-5;
    double tf = 1.0e-3;


    
    printf("Imposing bc..\n");
    impose_bc();
    printf("Done");
  AssignMatAddress();
  
    calcElemJAndDerivatives();
  CalcElemInitialVol(); //ALSO CALC VOL

  CalcElemVol();
  calcElemDensity();  
    printf("deriv %e", getDerivative(0,0,0,0));

    //calc_jacobian(x_, J);
    printf("m_m_detJ %6e\n",m_detJ[0]);
    
    double vol_0 = calc_vol();
    cout << "vol 0 "<<vol_0<<endl;
    double nod_mass = vol_0 * rho[0] / m_nodxelem;
    cout << "nodal mass "<< nod_mass <<endl;
    
        for (int i = 0; i < m_nodxelem; i++) 
            for (int j = 0; j < m_dim; j++)   
               prev_a[m_dim*i+j] = 0.0;

    double rho_b = 0.8182;
    m_alpha = (2.0 * rho_b - 1.0) / (1.0 + rho_b);
    m_beta = (5.0 - 3.0 * rho_b) / ((1.0 + rho_b) * (1.0 + rho_b) * (2.0 - rho_b));
    m_gamma = 1.5 - m_alpha;
    cout << "t, tf "<<t<<", " <<tf<<endl;
    while (t < tf) {
      printf ("Time-------------------------------------------------: %.6e\n", t);

        // PREDICTION PHASE
        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                //u_[i][j] = dt * (v_[i][j] + (0.5 - m_beta) * dt * prev_a_[i][j]);
                //NEW; global
                int ig = i*m_dim + j;
                u_dt[ig] = dt * (v[ig] + (0.5 - m_beta) * dt * prev_a[ig]);
            }
        }

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                //v_[i][j] += (1.0 - m_gamma) * dt * prev_a_[i][j];
                
                v[m_dim*i+j] += (1.0 - m_gamma) * dt * prev_a[m_dim*i+j];                
                printf("v %e",v[m_dim*i+j] );
            }
        }

        impose_bc();

        calcElemJAndDerivatives();


        // tensor3 Sigma = FromFlatSym(m_sigma,          0 );    
        // cout << "Sigma Tensor\n"<<endl;
        // print(Sigma); 
    
        //calc_jacobian(x_, J);
        
        //calc_str_rate(str_rate,rot_rate);
        calcElemStrainRates();
        
        double str_inc[m_nodxelem][3][3];
        calc_strain(str_rate, dt, str_inc);
        K_mod = mat[0]->Elastic().BulkMod();
        calc_pressure(K_mod, str_inc, stress);
        calcElemPressure(); //CRASHES IN 2D
                

        //printf("pressure %e\n",p[0]);
        CalcStressStrain(dt);
        
        //NOT WORKING
        calcElemForces();
        calcElemHourglassForces();
 
        
         //calc_forces(stress, a);
         //calc_forces2(stress, a);
         calc_hg_forces(rho[0], vol_0, mat[0]->cs0, m_f_elem_hg);
        
        // cout << "rho "<<rho<<"cs "<<mat_cs<<endl;
        // calc_hg_forces(rho, vol_0, mat_cs, f_hg);

        // for (int i = 0; i < m_nodxelem; i++) 
          // for (int j = 0; j < m_dim; j++) 
            // m_f_elem[i*m_dim + j] = -a_[i][j] / nod_mass + f_hg[i][j]; //LOCAL
          
        //ASSUMING LOCAL;
        //assemblyForces(); 

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                int ig = i*m_dim + j;
                //a[ig] = (-a[ig] + m_f_elem_hg[ig])/ nod_mass  - m_alpha * prev_a[ig];
                // a_[i][j] /= (1.0 - m_alpha);
                // v_[i][j] += m_gamma * dt * a_[i][j];
  

                printf ("force ELEMENT %6e ",m_f_elem[ig]);
                a[ig] = (-m_f_elem[ig] + m_f_elem_hg[ig])/ nod_mass  -m_alpha * prev_a[ig]; //GLOBAL
                
                printf ("hg f: %e ", m_f_elem_hg[ig]);
                a[ig] /= (1.0-m_alpha);
                v[ig] += m_gamma * dt * a[ig];
                
            }
        }

        impose_bc();

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                // u_[i][j] += m_beta * dt * dt * a_[i][j];
                // x_[i][j] += u_[i][j];
                
                int ig = i*m_dim + j;
                u_dt[ig] += m_beta * dt * dt * a[ig];
                x[ig] += u_dt[ig];
            }
        }

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                //prev_a_[i][j] = a_[i][j];
                
                prev_a[m_dim*i+j] = a[m_dim*i+j];
            }
        }

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                //u_tot_[i][j] += u_[i][j];
                u[m_dim*i+j] += u_dt[m_dim*i+j];
            }
        }

        t += dt;
        cout << "Time "<<t<< ", dt" <<dt<<endl;
    }

     printf("\n -- DISP --\n");
    for (int i = 0; i < m_nodxelem; i++) {
        for (int j = 0; j < m_dim; j++) {
            printf("%.6e ", u[i*m_dim+j]);
        }
      printf("\n");
    }

    printf("\n -- GLOB DISP --\n");
    for (int i = 0; i < m_nodxelem; i++) {
        for (int j = 0; j < m_dim; j++) {
            printf("%.6e ", u[m_dim*i+j]);
        }
      printf("\n");
    }


}

}; //CLASS

  
int main() {
    printf("Begin..\n");


    double t = 0.0;
        
    

	domain_test *dom_d;

  #ifdef CUDA_BUILD
	report_gpu_mem();
	gpuErrchk(cudaMallocManaged(&dom_d, sizeof(MetFEM::Domain_d)) );
	report_gpu_mem();
  #else
    dom_d = new domain_test;
  #endif
  
//  dom_d->m_red_int = false;
  
  double3 V = make_double3(0.0,0.0,0.0);
  double dx = 0.1;
	double3 L = make_double3(dx,dx,0.0);
  
  #ifdef BIDIM
  L.z = 0.0;
  #else
  L.z = dx;    
  #endif
	double r = 0.05;
	
	dom_d->AddBoxLength(V,L,r,true);
  // for (int i = 0; i < m_nodxelem; i++) {
      // for (int j = 0; j < m_dim; j++) {  
        // dom_d->x[ = 

  double E = 206e9;  // Young's modulus in Pa
  double nu = 0.3;   // Poisson's ratio
  double rho = 7850.0;

  dom_d->setDensity(rho);

  
  cout << "Creating Material..:"<<endl;
  Material_ *mat_h = (Material_ *)malloc(dom_d->getElemCount() * sizeof(Material_ *)); 
  Elastic_ el(E,nu);
  // cout << "Mat type  "<<mattype<<endl;

  Material_ *material_h;
  double Ep, c[6];
  // MATERIAL
  //TODO: MATERIALS SHOULD BE A VECTOR

  
  // !ONLY FOR ELASTIC
  // dom%mat_E = young
  // dom%mat_nu = poisson
  
  double mat_modK= E / ( 3.0*(1.0 -2.0*nu) );
  double mat_G= E / (2.0* (1.0 + nu));
  
  // dom%mat_K = mat_modK !!!TODO CREATE MATERIAL
  
  double mat_cs = sqrt(mat_modK/rho);
    
    material_h  = new Material_(el);
    material_h->cs0 = sqrt(material_h->Elastic().BulkMod()/rho); //TODO: INSIDE MATERIAL 
    cout << "CS_0: "<<material_h->cs0<<endl;
    material_h->Ep = Ep;
    material_h->Material_model = BILINEAR;
    
  dom_d->AssignMaterial(material_h);
  
  dom_d->Solve();
}