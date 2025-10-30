/*************************************************************************/
/*  1el_2d_sig_dev_pe.cpp                                        */
/*  WeldformFEM - High-Performance Explicit & Implicit FEM Solvers     */
/*  (CPU/GPU, C++/CUDA)                                                  */
/*                                                                       */
/*  weldform.sph@gmail.com                                                */
/*  ('https://www.opensourcemech.com',)                                    */
/*                                                                       */
/*  Copyright (c) 2023-2025 Luciano Buglioni          */
/*                                                                       */
/*  This file is part of the WeldformFEM project.                     */
/*  Licensed under the GNU General Public License v3.0 or later. */ 
/*  See the LICENSE file in the project    */
/*  root for full license information.                                   */
/*************************************************************************/




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
#define BIDIM 


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
double gauss_weights[m_gp_count] = {1.};

//double tf = 1.0e-3;


using namespace MetFEM;
class domain_test:
public Domain_d {
public:

double E = 206e9;  // Young's modulus in Pa
double nu = 0.3;   // Poisson's ratio
double rho = 7850.0;
double mat_G;
double K_mod;
int red_int = 0;
double element_length = 1.0; // Length of the element
int axi_symm = 0;            //FALSE: PLAIN STRAIN
  


#ifdef BIDIM
double x_[m_nodxelem][m_dim] = {{0.0,0.0},{0.1,0.0},{0.1,0.1},{0.0,0.1}};
#else
double x_[m_nodxelem][m_dim] = {{0.0,0.0,0.0},{0.1,0.0,0.0},{0.1,0.1,0.0},{0.0,0.1,0.0},
                                {0.0,0.0,0.1},{0.1,0.0,0.1},{0.1,0.1,0.1},{0.0,0.1,0.1}};
#endif

// VARS ENDING WITH _ are element local (for test w/o assembly)
double v_[m_nodxelem][m_dim];
double a_[m_nodxelem][m_dim];
double u_[m_nodxelem][m_dim];
double u_tot_[m_nodxelem][m_dim] ={{0.,0.},{0.,0.},{0.,0.},{0.,0.}};
double prev_a_[m_nodxelem][m_dim];
double f_hg[m_nodxelem][m_dim];

// double gauss_points[m_nodxelem][2]={{-0.577350269, -0.577350269},
                                    // {0.577350269, -0.577350269},
                                    // { 0.577350269,  0.577350269},
                                    // {-0.577350269,  0.577350269}};
                                    
// double gauss_weights[m_gp_count] = {1.,1.,1.,1.};



double invJ[m_gp_count][2][2];
double dNdX[m_gp_count][m_dim][m_nodxelem];
double str_rate[m_gp_count][3][3];
double rot_rate[m_gp_count][3][3];
double strain[m_gp_count][3][3];
Matrix tau_;
double tau[m_gp_count][3][3];
double pres[m_gp_count];
double stress[m_gp_count][3][3];
double radius[m_gp_count];
    double J[m_gp_count][2][2];


void impose_bc(double vel[m_nodxelem][m_dim], double accel[m_nodxelem][m_dim]) {
#ifdef BIDIM
//if dim == 2
    vel[2][1] = vel[3][1] = -1.0;
    vel[0][0] = vel[0][1] = vel[1][1] = 0.0;

    accel[2][1] = accel[3][1] = 0.0;
    accel[0][0] = accel[0][1] = accel[1][1] = 0.0;

#else
    //LOCAL CASE 
    //ELEMENT NODES; NOT GLOBAL
    vel  [0][0] = vel  [0][1] = vel  [0][2] = 0.0;
    accel[0][0] = accel[0][1] = accel[0][2] = 0.0;
    
    vel[1][1] = 0.0;     accel[1][2] = 0.0;
    vel[1][2] = 0.0;     accel[1][2] = 0.0;
    
    vel[2][2] = accel[2][2] = 0.0;
    
    vel[3][0] = accel[3][0] =0.0;
    vel[3][2] = 0.0; accel[3][2] = 0.0;
    
    
    //// GLOBAL CASE (nodes 1 and 2 inverted)
    // vel[2][1] = 0.0;     accel[2][2] = 0.0;
    // vel[2][2] = 0.0;     accel[2][2] = 0.0;
    
    // vel[1][2] = accel[1][2] = 0.0;
    
    

    for (int i=0;i<4;i++) {
      vel[4+i][2] = -1.0;
      accel[4+i][2] = 0.0;}

#endif
}


double calc_vol() {
    double vol = 0.0;
    for (int gp = 0; gp < m_gp_count; gp++) {
        vol += m_detJ[gp] * gauss_weights[gp];
    }
    return vol;
}

void velocity_gradient_tensor(double dNdX[m_gp_count][m_dim][m_nodxelem], double vel[m_nodxelem][m_dim], double grad_v[m_nodxelem][m_dim][m_dim]) {
    for (int gp = 0; gp < m_gp_count; gp++) {
        for (int I = 0; I < m_dim; I++) {
            for (int J = 0; J < m_dim; J++){ 
                grad_v[gp][I][J] = 0.0;
                for (int k = 0; k < m_nodxelem; k++) {
                    //grad_v[gp][I][J] += dNdX[gp][J][k] * vel[k][I];
                    grad_v[gp][I][J] += getDerivative(0,gp,J,k) * vel[k][I]/m_detJ[gp];
                    //printf ("deriv %e " , getDerivative(0,gp,J,k)/m_detJ[gp]);
                }

            }
        }
    }
    
}

void calc_str_rate(double dNdX[m_gp_count][m_dim][m_nodxelem], double v[m_nodxelem][m_dim], double str_rate[m_gp_count][3][3],double rot_rate[m_gp_count][3][3]) {
    double grad_v[m_nodxelem][m_dim][m_dim];
    velocity_gradient_tensor(dNdX, v, grad_v);
    for (int gp = 0; gp < m_gp_count; gp++) {
        for (int i = 0; i < m_dim; i++) {
            for (int j = 0; j < m_dim; j++) {
                str_rate[gp][i][j] = 0.5 * (grad_v[gp][i][j] + grad_v[gp][j][i]);
                rot_rate[gp][i][j] = 0.5 * (grad_v[gp][i][j] - grad_v[gp][j][i]);
                //printf("str rate %e", str_rate[gp][i][j]);
            }
        }
        // str_rate[gp][2][0]=rot_rate[gp][2][0]=0.0;                str_rate[gp][0][2]=rot_rate[gp][0][2]=0.0;        
        // str_rate[gp][2][2]=rot_rate[gp][2][2]=0.0;
    }
}

void calc_strain(double str_rate[m_gp_count][3][3], double dt, double strain[m_gp_count][3][3]) {

    for (int gp = 0; gp < m_gp_count; gp++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                strain[gp][i][j] = dt * str_rate[gp][i][j];

            }
        }
    }
}

void calc_pressure(double K_, double dstr[m_gp_count][3][3], double stress[m_gp_count][3][3], double pres[m_gp_count]) {
    double pi_ = 0.0;
    for (int gp = 0; gp < m_gp_count; gp++) {
        for (int i = 0; i < 3; i++) {
            pi_ += dstr[gp][i][i];
        }
    }
    pi_ = -pi_ / m_gp_count;
    for (int gp = 0; gp < m_gp_count; gp++) {
        pres[gp] = -1.0 / 3.0 * (stress[gp][0][0] + stress[gp][1][1] + stress[gp][2][2]) + K_ * pi_;
        //printf("pres %e ",pres[gp]);
    }
    
}

void dev(double t[3][3], double dev[3][3]) {
    for (int i = 0; i < 3; i++) 
        for (int j = 0; j < 3; j++) 
          dev[i][j]= t[i][j]- 1.0 / 3.0 * (t[0][0] + t[1][1] + t[2][2])*(i==j);
    
}

void calc_stress2(double str_rate[m_gp_count][3][3], double rot_rate[m_gp_count][3][3], double tau[m_gp_count][3][3], double p[m_gp_count], double dt, double stress[m_gp_count][3][3]) {
    double srt[m_gp_count][3][3];
    double rs[m_gp_count][3][3];
    double d[3][3];
    for (int gp = 0; gp < m_gp_count; gp++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                srt[gp][i][j] = rs[gp][i][j] = 0.0;
                for (int k=0;k<m_dim;k++){
                  srt[gp][i][j] += tau[gp][i][k] * rot_rate[gp][j][k];
                  rs[gp][i][j] += rot_rate[gp][i][k] * tau[gp][k][j];
                }
                dev(str_rate[gp],d);
                tau[gp][i][j] += dt * ((2.0 * mat_G *d[i][j]) + rs[gp][i][j] + srt[gp][i][j]);
                stress[gp][i][j] = tau[gp][i][j] - p[gp] * (i == j);
                //printf ("stress %e",stress[gp][i][j]);
            }
            printf("\n");
        }
    }
    // tensor3 Sigma = FromFlatSym(m_sigma,          0 );    

    // stress[0][0][0]=Sigma.xx;    stress[0][1][1]=Sigma.yy;    stress[0][2][2]=Sigma.zz;              
    // stress[0][0][1]=stress[0][1][0]=Sigma.xy;    
    // stress[0][1][2]=stress[0][2][1]=Sigma.yz;    
    // stress[0][0][2]=stress[0][2][0]=Sigma.xz; 
    // print(Sigma);    
}

void calc_forces(double stress[m_nodxelem][3][3], double dNdX[m_nodxelem][m_dim][m_nodxelem], double forces[m_nodxelem][m_dim]) {
    double B[m_dim][m_nodxelem];
    #ifdef BIDIMM
    double s2[2][2];
    #else
    double s2[3][3];      
    #endif
    
    for (int gp = 0; gp < m_gp_count; gp++) {
        for (int i = 0; i < m_dim; i++) 
            for (int j = 0; j < m_dim; j++) 
              s2[i][j]=stress[gp][i][j];

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
              forces[i][j] = 0.0;
              for (int k = 0; k < m_dim; k++) 
                //forces[i][j] += dNdX[gp][k][i] * s2[k][j]*m_detJ[gp] * gauss_weights[gp];
                forces[i][j] += getDerivative(0,gp,k,i) * s2[k][j]/**m_detJ[gp]*/ * gauss_weights[gp];
              
              //printf ("forces %e",forces[i][j]);
            }
        }
        
    }
}

void calc_hg_forces(double rho, double vol, double cs,double fhg[m_nodxelem][m_dim]){

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
        fhg[i][d] = 0.0;

  for (int j = 0; j < jmax; j++) 
    for (int i = 0; i < m_nodxelem; i++) 
      for (int d = 0; d < m_dim; d++) {
        hmod[d][j] +=v_[i][d]*Sig[j][i];
        //printf("hmod: %6e", hmod[d][j]);
      }

  double ch = 0.06 * pow (vol,0.6666666) * rho * 0.25 * cs;
  
      for (int j = 0; j < jmax; j++) 
        for (int i = 0; i < m_nodxelem; i++) 
          for (int d = 0; d < m_dim; d++) 
            fhg[i][d] -=ch*hmod[d][j]*Sig[j][i];
  

  
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
    mat_G = E / (2.0 * (1 + nu));
    K_mod = E / (3.0 * (1.0 - 2.0 * nu));

    double mat_cs = sqrt(K_mod/rho);
    
    printf("Imposing bc..\n");
    impose_bc(v_, a_);
    printf("Done");

    calcElemJAndDerivatives();
    //calc_jacobian(x_, J);
    printf("m_m_detJ %6e\n",m_detJ[0]);
    
    double vol_0 = calc_vol();
    cout << "vol 0 "<<vol_0<<endl;
    double nod_mass = vol_0 * rho / m_nodxelem;

    


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
                u_[i][j] = dt * (v_[i][j] + (0.5 - m_beta) * dt * prev_a_[i][j]);
                //NEW; global
                int ig = i*m_dim + j;
                u_dt[ig] = dt * (v[ig] + (0.5 - m_beta) * dt * prev_a[ig]);
            }
        }

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                v_[i][j] += (1.0 - m_gamma) * dt * prev_a_[i][j];
                
                v[m_dim*i+j] = v_[i][j]; // ASSUMING LOCAL     
                //v[m_dim*i+j] += (1.0 - m_gamma) * dt * prev_a[m_dim*i+j];                

            }
        }

        impose_bc(v_, a_);
        for (int i = 0; i < m_nodxelem; i++) 
            for (int j = 0; j < m_dim; j++) 
                v[m_dim*i+j] = v_[i][j]; // ASSUMING LOCAL     
        calcElemJAndDerivatives();

        // calcElemStrainRates();
        // calcElemPressure(); //CRASHES IN 2D
        // CalcStressStrain(dt);

        //calc_jacobian(x_, J);

        calc_str_rate(dNdX, v_, str_rate, rot_rate);
        double str_inc[m_nodxelem][3][3];
        calc_strain(str_rate, dt, str_inc);

        calc_pressure(K_mod, str_inc, stress, pres);
        
        calc_stress2(str_rate, rot_rate, tau, pres, dt, stress);

        calc_forces(stress, dNdX, a_);
        
        cout << "rho "<<rho<<"cs "<<mat_cs<<endl;
        cout << "mat cs " <<mat_cs<<endl;
        calc_hg_forces(rho, vol_0, mat_cs, f_hg);

        for (int i = 0; i < m_nodxelem; i++) 
          for (int j = 0; j < m_dim; j++) 
            m_f_elem[i*m_dim + j] = -a_[i][j] / nod_mass + f_hg[i][j]; //LOCAL
          
        //assemblyForces(); 

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                a_[i][j] = (-a_[i][j]+ f_hg[i][j]) / nod_mass - m_alpha * prev_a_[i][j];
                cout << "hg force "<<f_hg[i][j]<<endl;
                a_[i][j] /= (1.0 - m_alpha);
                v_[i][j] += m_gamma * dt * a_[i][j];
                v[m_dim*i+j] = v_[i][j]; // ASSUMING LOCAL     
                
                //GLOBAL: gives wrong results in 2D
                ////----
                // int ig = i*m_dim + j;
                // a[ig] = m_fi[ig] -m_alpha * prev_a[ig]; //GLOBAL
                // a[ig] /= (1.0-m_alpha);
                // v[ig] += m_gamma * dt * a[ig];
            }
        }

        impose_bc(v_, a_);

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                u_[i][j] += m_beta * dt * dt * a_[i][j];
                x_[i][j] += u_[i][j];
                //x[3*i+j] += u_[i][j];
                
                int ig = i*m_dim + j;
                u_dt[ig] += m_beta * dt * dt * a[ig];
                //x[ig] += u_dt[ig];
            }
        }

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                prev_a_[i][j] = a_[i][j];
                
                prev_a[m_dim*i+j] = a[m_dim*i+j];
            }
        }

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                u_tot_[i][j] += u_[i][j];
                u[m_dim*i+j] += u[m_dim*i+j];
            }
        }

        t += dt;
        cout << "Time "<<t<< ", dt" <<dt<<endl;
    }

     printf("\n -- DISP --\n");
    for (int i = 0; i < m_nodxelem; i++) {
        for (int j = 0; j < m_dim; j++) {
            printf("%.6e ", u_tot_[i][j]);
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
    
    for (int i = 0; i < m_nodxelem; i++) {
        for (int j = 0; j < m_dim; j++) {
            printf("%.6e ", a_[i][j]);
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
  
  
  dom_d->Solve();
}
