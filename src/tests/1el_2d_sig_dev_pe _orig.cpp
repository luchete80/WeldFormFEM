#include <stdio.h>
#include <math.h>
#include <iostream>

using namespace std;

#define m_dim 2
#define m_nodxelem 4
#define m_gp_count 1

double E = 206e9;  // Young's modulus in Pa
double nu = 0.3;   // Poisson's ratio
double rho = 7850.0;
double mat_G;
double K_mod;
int red_int = 0;
double element_length = 1.0; // Length of the element
int axi_symm = 0;            //FALSE: PLAIN STRAIN

double dt = 0.8e-5;
//double tf = 0.8e-5;
double tf = 1.0e-3;
double x[m_nodxelem][m_dim] = {{0.0,0.0},{0.1,0.0},{0.1,0.1},{0.0,0.1}};
double v[m_nodxelem][m_dim];
double a[m_nodxelem][m_dim];
double u[m_nodxelem][m_dim];
double u_tot[m_nodxelem][m_dim] ={{0.,0.},{0.,0.},{0.,0.},{0.,0.}};

double f_hg[m_nodxelem][m_dim];

// double gauss_points[m_nodxelem][2]={{-0.577350269, -0.577350269},
                                    // {0.577350269, -0.577350269},
                                    // { 0.577350269,  0.577350269},
                                    // {-0.577350269,  0.577350269}};
                                    
// double gauss_weights[m_gp_count] = {1.,1.,1.,1.};

double gauss_points[1][2] = {{0,0}};
double gauss_weights[m_gp_count] = {1.};

double detJ[m_gp_count];
double invJ[m_gp_count][2][2];
double dNdX[m_gp_count][m_dim][m_nodxelem];
double str_rate[m_gp_count][3][3];
double rot_rate[m_gp_count][3][3];
double strain[m_gp_count][3][3];
double tau[m_gp_count][3][3];
double pres[m_gp_count];
double stress[m_gp_count][3][3];
double radius[m_gp_count];

void impose_bc(double vel[m_nodxelem][m_dim], double accel[m_nodxelem][m_dim]) {
    vel[2][1] = vel[3][1] = -1.0;
    vel[0][0] = vel[0][1] = vel[1][1] = 0.0;

    accel[2][1] = accel[3][1] = 0.0;
    accel[0][0] = accel[0][1] = accel[1][1] = 0.0;
}

void shape_functions(double xi, double eta, double N[m_nodxelem], double dNdX_[m_dim][m_nodxelem]) {
    N[0] = (1 - xi) * (1 - eta) / 4;
    N[1] = (1 + xi) * (1 - eta) / 4;
    N[2] = (1 + xi) * (1 + eta) / 4;
    N[3] = (1 - xi) * (1 + eta) / 4;

    dNdX_[0][0] = -(1 - eta) / 4;
    dNdX_[0][1] = (1 - eta) / 4;
    dNdX_[0][2] = (1 + eta) / 4;
    dNdX_[0][3] = -(1 + eta) / 4;

    dNdX_[1][0] = -(1 - xi) / 4;
    dNdX_[1][1] = -(1 + xi) / 4;
    dNdX_[1][2] = (1 + xi) / 4;
    dNdX_[1][3] = (1 - xi) / 4;
    
}

void calc_jacobian(double pos[m_nodxelem][m_dim], double J[m_gp_count][2][2]) {
    double N[m_nodxelem];
    double dNdX_[m_dim][m_nodxelem];
    double xi, eta;
    for (int gp = 0; gp < m_gp_count; gp++) {
        xi = gauss_points[gp][0];
        eta = gauss_points[gp][1];
        shape_functions(xi, eta, N, dNdX_);        
        for (int i = 0; i < m_dim; i++) {
            for (int j = 0; j < m_dim; j++) {
                J[gp][i][j] = 0.0;
                for (int k = 0; k < m_nodxelem; k++) {
                    J[gp][i][j] += dNdX_[i][k] * pos[k][j];
                      
        // elem%jacob(e,gp,1,:) = -x2(1,:)+x2(2,:)+x2(3,:)-x2(4,:)
        // elem%jacob(e,gp,2,:) = -x2(1,:)-x2(2,:)+x2(3,:)+x2(4,:)                
                
  
                    // printf("pos %.6e", pos[k][j]);
                    // printf ("J %.6e", J[gp][i][j]);
                }
            }
        }
        //1gp
        // for (int i = 0; i < m_dim; i++){
          // J[gp][0][i] = 0.25*(-pos[0][i]+pos[1][i]+pos[2][i]-pos[3][i]);
          // J[gp][1][i] = 0.25*(-pos[0][i]-pos[1][i]+pos[2][i]+pos[3][i]);                     
        // }
        // printf ("J %.6e %.6e \n %.6e %.6e\n", J[gp][0][0], J[gp][0][1], J[gp][1][0], J[gp][1][1] );
        double adJ[2][2]; 
        adJ[0][0]= J[gp][1][1];adJ[1][1]= J[gp][0][0];
        adJ[0][1]=-J[gp][0][1];adJ[1][0]=-J[gp][1][0];
        
        detJ[gp] = J[gp][0][0] * J[gp][1][1] - J[gp][0][1] * J[gp][1][0];
        for (int i = 0; i < m_dim; i++) {
            for (int j = 0; j < m_dim; j++) {
                  invJ[gp][i][j] = 1.0/detJ[gp]*adJ[i][j];
            }
        }

        for (int i = 0; i < m_dim; i++) {
            for (int j = 0; j < m_nodxelem; j++) {
              dNdX[gp][i][j] = 0.0;
              for (int k = 0; k < m_dim; k++) 
                dNdX[gp][i][j] += invJ[gp][i][k]*dNdX_[k][j];
            }
        }                
        
        
        //printf ("detJ %.6e\n", detJ[gp]);
    }
}

double calc_vol(double detJ[m_gp_count]) {
    double vol = 0.0;
    for (int gp = 0; gp < m_gp_count; gp++) {
        vol += detJ[gp] * gauss_weights[gp];
    }
    return vol;
}

void velocity_gradient_tensor(double dNdX[m_gp_count][m_dim][m_nodxelem], double vel[m_nodxelem][m_dim], double grad_v[m_nodxelem][m_dim][m_dim]) {
    for (int gp = 0; gp < m_gp_count; gp++) {
        for (int I = 0; I < m_dim; I++) {
            for (int J = 0; J < m_dim; J++){ 
                grad_v[gp][I][J] = 0.0;
                for (int k = 0; k < m_nodxelem; k++) {
                    grad_v[gp][I][J] += dNdX[gp][J][k] * vel[k][I];
                }
                // printf ("grad v %e " , grad_v[gp][I][J]);
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
                 printf ("str rat %.6e ",str_rate[gp][i][j]);

            }
        }
        str_rate[gp][2][0]=rot_rate[gp][2][0]=0.0;                str_rate[gp][0][2]=rot_rate[gp][0][2]=0.0;        
        str_rate[gp][2][2]=rot_rate[gp][2][2]=0.0;
    }
}

void calc_strain(double str_rate[m_gp_count][3][3], double dt, double strain[m_gp_count][3][3]) {
    printf ("strain\n");
    for (int gp = 0; gp < m_gp_count; gp++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                strain[gp][i][j] = dt * str_rate[gp][i][j];
                printf ("strain %.6e ",strain[gp][i][j]);
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
        printf("pres %e ",pres[gp]);
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
                printf ("stress %e",stress[gp][i][j]);
            }
            printf("\n");
        }
    }
}

void calc_forces(double stress[m_nodxelem][3][3], double dNdX[m_nodxelem][m_dim][m_nodxelem], double detJ[m_nodxelem], double forces[m_nodxelem][m_dim]) {
    double B[m_dim][m_nodxelem];
    double s2[2][2];
    
    for (int gp = 0; gp < m_gp_count; gp++) {
        for (int i = 0; i < m_dim; i++) 
            for (int j = 0; j < m_dim; j++) 
              s2[i][j]=stress[gp][i][j];

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
              forces[i][j] = 0.0;
              for (int k = 0; k < m_dim; k++) 
                forces[i][j] += dNdX[gp][k][i] * s2[k][j]*detJ[gp] * gauss_weights[gp];
              printf ("forces %e",forces[i][j]);
            }
            printf("\n");
        }
        
    }
}

void calc_hg_forces(double rho, double vol, double cs,double fhg[m_nodxelem][m_dim]){

  double Sig [4][4] = {{1.,-1.,1.,-1.}, {1.,-1.,1.,-1.},{1.,-1.,1.,-1.},{1.,-1.,1.,-1.}};
  double hmod[m_dim][4]={{0.0,0.0,0.0,0.0},{0.0,0.0,0.0,0.0}};
  int jmax = 4;

  for (int j = 0; j < jmax; j++) 
    for (int i = 0; i < m_nodxelem; i++) 
      for (int d = 0; d < m_dim; d++) {
        hmod[d][j] +=v[i][d]*Sig[j][i];
        printf("hmod: %6e", hmod[d][j]);
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

int main() {
    printf("Begin..\n");
    mat_G = E / (2.0 * (1 + nu));
    K_mod = E / (3.0 * (1.0 - 2.0 * nu));
    
    printf("Imposing bc..\n");
    impose_bc(v, a);
    printf("Done");

    double prev_a[m_nodxelem][m_dim];

    double J[m_gp_count][2][2];
    calc_jacobian(x, J);

    double vol_0 = calc_vol(detJ);
    cout << "vol 0 "<<vol_0<<endl;
    double nod_mass = vol_0 * rho / m_nodxelem;

    double rho_b = 0.8182;
    double alpha = (2.0 * rho_b - 1.0) / (1.0 + rho_b);
    double beta = (5.0 - 3.0 * rho_b) / ((1.0 + rho_b) * (1.0 + rho_b) * (2.0 - rho_b));
    double gamma = 1.5 - alpha;
    
    double mat_cs = sqrt(K_mod/rho);

    double t = 0.0;
    while (t < tf) {
      printf ("Time: %.6e\n", t);

        // PREDICTION PHASE
        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                u[i][j] = dt * (v[i][j] + (0.5 - beta) * dt * prev_a[i][j]);
            }
        }

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                v[i][j] += (1.0 - gamma) * dt * prev_a[i][j];
            }
        }

        impose_bc(v, a);

        calc_jacobian(x, J);

        calc_str_rate(dNdX, v, str_rate, rot_rate);
        double str_inc[m_nodxelem][3][3];
        calc_strain(str_rate, dt, str_inc);

        calc_pressure(K_mod, str_inc, stress, pres);
        
        calc_stress2(str_rate, rot_rate, tau, pres, dt, stress);

        calc_forces(stress, dNdX, detJ, a);
        
        calc_hg_forces(rho, vol_0, mat_cs, f_hg);

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                a[i][j] = -a[i][j] / nod_mass + f_hg[i][j] - alpha * prev_a[i][j];
                a[i][j] /= (1.0 - alpha);
                v[i][j] += gamma * dt * a[i][j];
            }
        }

        impose_bc(v, a);

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                u[i][j] += beta * dt * dt * a[i][j];
                x[i][j] += u[i][j];
            }
        }

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                prev_a[i][j] = a[i][j];
            }
        }

        for (int i = 0; i < m_nodxelem; i++) {
            for (int j = 0; j < m_dim; j++) {
                u_tot[i][j] += u[i][j];
            }
        }

        t += dt;
    }

     printf("\n -- DISP --\n");
    for (int i = 0; i < m_nodxelem; i++) {
        for (int j = 0; j < m_dim; j++) {
            printf("%.6e ", u_tot[i][j]);
        }
      printf("\n");
    }

    for (int i = 0; i < m_nodxelem; i++) {
        for (int j = 0; j < m_dim; j++) {
            printf("%.6e ", a[i][j]);
        }
      printf("\n");
    }
    return 0;
}
