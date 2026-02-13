// ================================
// stokes_martins_axisymmetric.cpp
// Implementación del esquema de Martins para axisimetría
// Usando tu clase Matrix - Esquema de punto de silla directo
// ================================

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <memory>
#include "Solver_Eigen.h"
#include "Mesh.h"
#include "WallTimer.h"
#include "Matrix.h"
#include "Domain_d.h"

namespace MetFEM{

Matrix SolveLU_(Matrix& A, Matrix& b) {
    if (A.m_row != A.m_col || b.m_col != 1 || A.m_row != b.m_row) {
        Matrix empty;
        return empty;
    }

    const int n = A.m_row;
    Matrix LU = A;
    Matrix x(n, 1);
    int* pivot = new int[n];
    Matrix pb = b;

    // Inicializar pivotes
    for (int i = 0; i < n; ++i) pivot[i] = i;

    // --- Factorización LU con pivoteo parcial ---
    for (int i = 0; i < n; ++i) {
        // Pivoteo parcial - USANDO getVal()
        int max_row = i;
        double max_val = fabs(LU.getVal(i, i));
        for (int k = i + 1; k < n; ++k) {
            if (fabs(LU.getVal(k, i)) > max_val) {
                max_val = fabs(LU.getVal(k, i));
                max_row = k;
            }
        }
        
        // Intercambiar filas - USANDO getVal/setVal
        if (max_row != i) {
            for (int j = 0; j < n; ++j) {
                double temp = LU.getVal(i, j);
                LU.Set(i, j, LU.getVal(max_row, j));
                LU.Set(max_row, j, temp);
            }
            int temp_pivot = pivot[i];
            pivot[i] = pivot[max_row];
            pivot[max_row] = temp_pivot;
        }

        double diag = LU.getVal(i, i);
        if (fabs(diag) < 1e-12) {
            delete[] pivot;
            Matrix empty;
            return empty;
        }

        // Eliminación - USANDO getVal/setVal
        for (int k = i + 1; k < n; ++k) {
            double factor = LU.getVal(k, i) / diag;
            LU.Set(k, i, factor);
            for (int j = i + 1; j < n; ++j) {
                double val = LU.getVal(k, j) - factor * LU.getVal(i, j);
                LU.Set(k, j, val);
            }
        }
    }

    // --- Aplicar permutaciones a b - USANDO getVal/setVal ---
    Matrix b_permuted = b;
    for (int i = 0; i < n; ++i) {
        b_permuted.Set(i, 0, b.getVal(pivot[i], 0));
    }

    // --- Sustitución hacia adelante ---
    Matrix y(n, 1);
    for (int i = 0; i < n; ++i) {
        double sum = b_permuted.getVal(i, 0);
        for (int j = 0; j < i; ++j) {
            sum -= LU.getVal(i, j) * y.getVal(j, 0);
        }
        y.Set(i, 0, sum);
    }

    // --- Sustitución hacia atrás ---
    for (int i = n - 1; i >= 0; --i) {
        double sum = y.getVal(i, 0);
        for (int j = i + 1; j < n; ++j) {
            sum -= LU.getVal(i, j) * x.getVal(j, 0);
        }
        x.Set(i, 0, sum / LU.getVal(i, i));
    }

    delete[] pivot;
    return x;
}

// ================================
// PARÁMETROS FÍSICAMENTE CORRECTOS
// ================================

// MATERIAL: Acero caliente tipo Norton-Hoff
double Kmat = 50.0e6;          // Coeficiente de consistencia [Pa·s^n]
double nexp = 0.2;             // Exponente strain-rate sensitivity


// Parámetros numéricos
double omega_v = 0.4;          // Relajación para velocidades
double omega_p = 0.1;          // Relajación para presiones (MUY importante)


double tol = 1e-3;             // Tolerancia relajada para simulación práctica
int max_iter = 50;
int nsteps = 5;
double v_z_top = -1.0;         // Compresión en dirección Z

// ================================
// ESTRUCTURAS DE DATOS
// ================================
struct Point2D {
    double x, y;  // x = r (radial), y = z (axial)
    Point2D(double r_=0, double z_=0) : x(r_), y(z_) {}
};

// Variables globales
int ndof_v, ndof_p, ndof_total;

std::vector<Point2D> coords;
std::vector<std::vector<int>> elements;
std::vector<double> velocity;    // [vr0, vz0, vr1, vz1, ...]
std::vector<double> pressure;    // Presión por elemento (P0 discontinua)
std::vector<double> eps_bar;     // Deformación plástica acumulada
std::vector<std::vector<Point2D>> coords_history;

// ================================
// FUNCIONES AUXILIARES
// ================================
inline double max(double a, double b) { return a > b ? a : b; }
inline double min(double a, double b) { return a < b ? a : b; }
inline double sign(double x) { return (x > 0) ? 1.0 : ((x < 0) ? -1.0 : 0.0); }

// ================================
// FUNCIONES DE FORMA Q1 (AXISIMÉTRICO)
// ================================
void shape_functions(double xi, double eta, 
                     std::vector<double>& N,
                     Matrix& dNdxi) {
    // N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)]
    N.resize(4);
    N[0] = 0.25 * (1.0 - xi) * (1.0 - eta);
    N[1] = 0.25 * (1.0 + xi) * (1.0 - eta);
    N[2] = 0.25 * (1.0 + xi) * (1.0 + eta);
    N[3] = 0.25 * (1.0 - xi) * (1.0 + eta);
    
    // dNdxi [2x4]
    dNdxi = Matrix(2, 4);
    dNdxi.Set(0, 0, -0.25*(1.0-eta)); dNdxi.Set(0, 1,  0.25*(1.0-eta));
    dNdxi.Set(0, 2,  0.25*(1.0+eta)); dNdxi.Set(0, 3, -0.25*(1.0+eta));
    dNdxi.Set(1, 0, -0.25*(1.0-xi));  dNdxi.Set(1, 1, -0.25*(1.0+xi));
    dNdxi.Set(1, 2,  0.25*(1.0+xi));  dNdxi.Set(1, 3,  0.25*(1.0-xi));
}

struct JacobianResult {
    std::vector<double> N;
    Matrix dNdX;
    double detJ;
    Matrix J;
};

JacobianResult jacobian_and_gradients(const std::vector<Point2D>& pos, 
                                      double xi, double eta) {
    JacobianResult result;
    
    // Obtener funciones de forma
    Matrix dNdxi;
    shape_functions(xi, eta, result.N, dNdxi);
    
    // Construir matriz de posiciones [4x2]
    Matrix pos_mat(4, 2);
    for(int i = 0; i < 4; i++) {
        pos_mat.Set(i, 0, pos[i].x);  // r
        pos_mat.Set(i, 1, pos[i].y);  // z
    }
    
    // Jacobiano: J = dNdxi * pos_mat (2x4 * 4x2 = 2x2)
    result.J = Matrix(2, 2);
    for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 2; j++) {
            double sum = 0.0;
            for(int k = 0; k < 4; k++) {
                sum += dNdxi.getVal(i, k) * pos_mat.getVal(k, j);
            }
            result.J.Set(i, j, sum);
        }
    }
    
    // Determinante del Jacobiano
    double a = result.J.getVal(0, 0), b = result.J.getVal(0, 1);
    double c = result.J.getVal(1, 0), d = result.J.getVal(1, 1);
    result.detJ = a*d - b*c;
    
    if(std::abs(result.detJ) < 1e-12) {
        result.detJ = (result.detJ >= 0 ? 1.0 : -1.0) * 1e-12;
    }
    
    // Inversa del Jacobiano
    Matrix invJ(2, 2);
    invJ.Set(0, 0,  d/result.detJ); invJ.Set(0, 1, -b/result.detJ);
    invJ.Set(1, 0, -c/result.detJ); invJ.Set(1, 1,  a/result.detJ);
    
    // Gradientes espaciales: dN/dX = inv(J) * dN/dξ
    result.dNdX = Matrix(2, 4);
    for(int i = 0; i < 2; i++) {
        for(int j = 0; j < 4; j++) {
            double sum = 0.0;
            for(int k = 0; k < 2; k++) {
                sum += invJ.getVal(i, k) * dNdxi.getVal(k, j);
            }
            result.dNdX.Set(i, j, sum);
        }
    }
    
    return result;
}


//~ // ================================
//~ // FUNCIÓN DE RESOLUCIÓN LINEAL (necesita ser implementada)
//~ // ================================
std::vector<double> solve_linear_system(const Matrix& K, const std::vector<double>& F) {
    // Implementar con tu solver preferido
    // Este es un placeholder simple usando eliminación gaussiana
    
    //int n = K.getRows();
    int n = K.m_row;
    std::vector<double> solution(n, 0.0);
    
    // Hacer una copia para no modificar la original
    Matrix A = K;
    std::vector<double> b = F;
    
    // Eliminación gaussiana
    for(int i = 0; i < n; i++) {
        // Pivoteo
        int max_row = i;
        double max_val = fabs(A.getVal(i, i));
        for(int k = i+1; k < n; k++) {
            if(fabs(A.getVal(k, i)) > max_val) {
                max_val = fabs(A.getVal(k, i));
                max_row = k;
            }
        }
        
        if(max_row != i) {
            // Intercambiar filas
            for(int j = i; j < n; j++) {
                double temp = A.getVal(i, j);
                A.Set(i, j, A.getVal(max_row, j));
                A.Set(max_row, j, temp);
            }
            double temp = b[i];
            b[i] = b[max_row];
            b[max_row] = temp;
        }
        
        // Hacer ceros debajo de la diagonal
        for(int k = i+1; k < n; k++) {
            double factor = A.getVal(k, i) / A.getVal(i, i);
            for(int j = i; j < n; j++) {
                A.Set(k, j, A.getVal(k, j) - factor * A.getVal(i, j));
            }
            b[k] -= factor * b[i];
        }
    }
    
    // Sustitución hacia atrás
    for(int i = n-1; i >= 0; i--) {
        solution[i] = b[i];
        for(int j = i+1; j < n; j++) {
            solution[i] -= A.getVal(i, j) * solution[j];
        }
        solution[i] /= A.getVal(i, i);
    }
    
    return solution;
}

// AÑADIR esta función (ecuación 4.50 de Martins)
double equivalent_stress_norton(double eps_dot_eq) {
    double eps_dot_eq_safe = max(eps_dot_eq, 1e-12);
    double sqrt3_eps = sqrt(3.0) * eps_dot_eq_safe;
    return Kmat * pow(sqrt3_eps, nexp);
}

struct StrainRateResult {
    double eps_dot_eq;
    std::vector<double> strain_vec;  // [err, ezz, ett, erz] para axisimetría
    //std::vector<double> C_vec;       // Vector C de Martins para presión hidrostática
};

StrainRateResult calculate_strain_rate_martins(Matrix& dNdX,
                                              const std::vector<double>& vel_elem,
                                              double r_gp,
                                              const std::vector<double>& N) {
    StrainRateResult result;
    result.strain_vec.resize(4);  // 4 componentes para axisimetría
    
    // Inicializar gradientes y velocidad interpolada
    double dvr_dr = 0.0, dvr_dz = 0.0;
    double dvz_dr = 0.0, dvz_dz = 0.0;
    double vr_interp = 0.0;
    
    for(int a = 0; a < 4; a++) {
        double vr = vel_elem[2*a];      // vr del nodo a
        double vz = vel_elem[2*a + 1];  // vz del nodo a
        
        dvr_dr += dNdX.getVal(0, a) * vr;
        dvr_dz += dNdX.getVal(1, a) * vr;
        dvz_dr += dNdX.getVal(0, a) * vz;
        dvz_dz += dNdX.getVal(1, a) * vz;
        vr_interp += N[a] * vr;
    }
      
   // Tensor de deformación (COMPONENTES COMPLETAS, NO multiplicadas)
      double err = dvr_dr;                          // ε_rr
      double ezz = dvz_dz;                          // ε_zz
      double ett = (r_gp < 1e-8) ? dvr_dr : vr_interp / r_gp;  // ε_θθ
      double erz = 0.5 * (dvr_dz + dvz_dr);         // ε_rz (¡CON 0.5!)
      
      result.strain_vec[0] = err;
      result.strain_vec[1] = ezz;
      result.strain_vec[2] = ett;
      result.strain_vec[3] = erz;
      
      // --- CORREGIDO: Cálculo exacto de ε̇_eq según d_martins ---
      // d_martins = diag([2/3, 2/3, 2/3, 1/3])
      double eps_dot_eq_sq = (2.0/3.0) * (err*err + ezz*ezz + ett*ett) 
                           + (1.0/3.0) * (2.0*erz*erz);  // Nota: 2*erz² porque es cortante
      
      result.eps_dot_eq = sqrt(eps_dot_eq_sq);
      if(result.eps_dot_eq < 1e-12) result.eps_dot_eq = 1e-12;
      
      return result;
}

#include "VTKWriter_tiny.hpp"


// ================================
// ENSAMBLAJE MARTINS CON INTEGRACIÓN REDUCIDA
// ================================
// ================================
// ENSAMBLAJE MARTINS CON INTEGRACIÓN REDUCIDA Y PENALIZACIÓN
// ================================
  void Domain_d::assemble_martins_system(const std::vector<double>& vel_vec,
                            const std::vector<double>& press_vec,
                            Matrix& K_glob,
                            std::vector<double>& F_glob,
                            double& incompressibility_error,
                            double& max_div_v) {
    
    // Inicializar sistema global
    K_glob.SetZero();
    std::fill(F_glob.begin(), F_glob.end(), 0.0);
    incompressibility_error = 0.0;
    max_div_v = 0.0;
        
    // Puntos de Gauss: 4 puntos completos para P, 1 punto reducido para Q
    const double a = 1.0 / sqrt(3.0);
    std::vector<std::pair<double, double>> gp_full = {
        {-a, -a}, {a, -a}, {a, a}, {-a, a}
    };
    std::vector<double> w_full = {1.0, 1.0, 1.0, 1.0};
    
    std::vector<std::pair<double, double>> gp_reduced = {{0.0, 0.0}};
    std::vector<double> w_reduced = {4.0};  // peso para 1 punto
    
    for(int e_idx = 0; e_idx < m_elem_count; e_idx++) {
        const auto& conn = elements[e_idx];
        
        // Posiciones de los nodos del elemento
        std::vector<Point2D> pos(4);
        for(int i = 0; i < 4; i++) {//node

            int n =  getElemNode(e_idx,i);
            //pos[i].x = x[m_dim*n];
            //pos[i].y = x[m_dim*n+1];
            pos[i] = coords[conn[i]];
        }
        
        // DOFs de velocidades del elemento
        std::vector<int> dofs_v;
        for(int n : conn) {
            dofs_v.push_back(2*n);      // vr del nodo n
            dofs_v.push_back(2*n + 1);  // vz del nodo n
        }
        
        int dof_p = ndof_v + e_idx;  // Presión P0 en este elemento
        
        // Extraer velocidades del elemento
        std::vector<double> vel_elem(8);
        for(int i = 0; i < 8; i++) {
            vel_elem[i] = vel_vec[dofs_v[i]];
        }
        
        double press_elem = press_vec[e_idx];
        
        // Matrices elementales según esquema de Martins
        Matrix P_e(8, 8);   // σP en eq. 4.56
        Matrix Q_e(8, 1);   // Q en eq. 4.56
        P_e.SetZero();
        Q_e.SetZero();
        
        // ===========================================
        // 1. INTEGRACIÓN COMPLETA PARA P (4 puntos)
        // ===========================================
        for(int gp = 0; gp < 4; gp++) {
            double xi = gp_full[gp].first;
            double eta = gp_full[gp].second;
            double w = w_full[gp];
            
            auto jac_result = jacobian_and_gradients(pos, xi, eta);
            
            // Coordenada radial en el punto de Gauss
            double r_gp = 0.0;
            for(int i = 0; i < 4; i++) {
                r_gp += jac_result.N[i] * pos[i].x;
            }
            r_gp = max(r_gp, 1e-12);
            
            // Factor axisimétrico
            double dOmega = 2.0 * M_PI * r_gp * jac_result.detJ * w;
            
            // Calcular strain-rate
            auto strain_result = calculate_strain_rate_martins(
                jac_result.dNdX, vel_elem, r_gp, jac_result.N
            );
            
            // --- NUEVO: Matriz d de Martins (ecuación 4.51) ---
            Matrix d_martins(4, 4);
            d_martins.SetZero();
            d_martins.Set(0, 0, 2.0/3.0);
            d_martins.Set(1, 1, 2.0/3.0);
            d_martins.Set(2, 2, 2.0/3.0);
            d_martins.Set(3, 3, 1.0/3.0);  // Para componente cortante

            // Matriz B (4x8) para axisimetría según Martins
            Matrix B(4, 8);
            B.SetZero();
            
            // ε_rr = ∂vr/∂r
            for(int a = 0; a < 4; a++) {
                B.Set(0, 2*a, jac_result.dNdX.getVal(0, a));
            }
            
            // ε_zz = ∂vz/∂z
            for(int a = 0; a < 4; a++) {
                B.Set(1, 2*a + 1, jac_result.dNdX.getVal(1, a));
            }
            
            // ε_θθ = vr/r (o tratamiento especial para eje)
            if(r_gp < 1e-8) {
                // En el eje: ε_θθ = ∂vr/∂r
                for(int a = 0; a < 4; a++) {
                    B.Set(2, 2*a, jac_result.dNdX.getVal(0, a));
                }
            } else {
                // Fuera del eje: ε_θθ = vr/r
                for(int a = 0; a < 4; a++) {
                    B.Set(2, 2*a, jac_result.N[a] / r_gp);
                }
            }
            
            // ε_rz = ∂vr/∂z + ∂vz/∂r (SIN 0.5 - el tensor ya tiene 0.5 en la definición)
            for(int a = 0; a < 4; a++) {
                B.Set(3, 2*a, jac_result.dNdX.getVal(1, a));      // ∂vr/∂z
                B.Set(3, 2*a + 1, jac_result.dNdX.getVal(0, a)); // ∂vz/∂r
            }
     
            // --- NUEVO: Matriz k = B^T * d * B (ecuación 4.56) ---
            Matrix Bt = B.getTranspose();
            Matrix dB = Matrix(4, 8);
            
            // d * B
            for(int i = 0; i < 4; i++) {
                for(int j = 0; j < 8; j++) {
                    double sum = 0.0;
                    for(int k = 0; k < 4; k++) {
                        sum += d_martins.getVal(i, k) * B.getVal(k, j);
                    }
                    dB.Set(i, j, sum);
                }
            }
            
            // k_matrix = B^T * (d * B)
            Matrix k_matrix(8, 8);
            k_matrix.SetZero();
            for(int i = 0; i < 8; i++) {
                for(int j = 0; j < 8; j++) {
                    double sum = 0.0;
                    for(int k = 0; k < 4; k++) {
                        sum += Bt.getVal(i, k) * dB.getVal(k, j);
                    }
                    k_matrix.Set(i, j, sum);
                }
            }
            
            // --- NUEVO: Tensión equivalente de Norton-Hoff (ecuación 4.50) ---
            double sigma_eq = equivalent_stress_norton(strain_result.eps_dot_eq);
            
            // --- NUEVO: Contribución a P_e (ecuación 4.56) ---
            // P_e += (sigma_eq / eps_dot_eq) * k_matrix * dOmega
            double factor = sigma_eq / strain_result.eps_dot_eq;
            
            for(int i = 0; i < 8; i++) {
                for(int j = 0; j < 8; j++) {
                    double val = P_e.getVal(i, j) 
                               + k_matrix.getVal(i, j) * factor * dOmega;
                    P_e.Set(i, j, val);
                }
            }
        }
              
      // ===========================================
      // 2. INTEGRACIÓN REDUCIDA PARA Q (1 punto)
      //    Matriz Q = ∫ B^T C dV  (Martins eq. 4.55)
      //    Vector C = [1, 1, 1, 0]^T (Martins eq. 4.49)
      // ===========================================
      for(int gp = 0; gp < 1; gp++) {
          double xi = gp_reduced[gp].first;
          double eta = gp_reduced[gp].second;
          double w = w_reduced[gp];
          
          auto jac_result = jacobian_and_gradients(pos, xi, eta);
          
          // Coordenada radial en el punto de Gauss
          double r_gp = 0.0;
          for(int i = 0; i < 4; i++) {
              r_gp += jac_result.N[i] * pos[i].x;
          }
          r_gp = max(r_gp, 1e-12);
          
          // Factor axisimétrico: dΩ = 2π·r·|J|·w
          double dOmega = 2.0 * M_PI * r_gp * jac_result.detJ * w;
          
          // --- MATRIZ B COMPLETA (4x8) PARA AXISIMETRÍA ---
          Matrix B(4, 8);
          B.SetZero();
          
          // ε_rr = ∂vr/∂r
          for(int a = 0; a < 4; a++) {
              B.Set(0, 2*a, jac_result.dNdX.getVal(0, a));
          }
          
          // ε_zz = ∂vz/∂z
          for(int a = 0; a < 4; a++) {
              B.Set(1, 2*a + 1, jac_result.dNdX.getVal(1, a));
          }
          
          // ε_θθ = vr/r (o límite en el eje)
          if(r_gp < 1e-8) {
              // En el eje: ε_θθ = ∂vr/∂r
              for(int a = 0; a < 4; a++) {
                  B.Set(2, 2*a, jac_result.dNdX.getVal(0, a));
              }
          } else {
              // Fuera del eje: ε_θθ = vr/r
              for(int a = 0; a < 4; a++) {
                  B.Set(2, 2*a, jac_result.N[a] / r_gp);
              }
          }
                    
          // ε_rz = ∂vr/∂z + ∂vz/∂r (SIN 0.5 - el tensor ya tiene 0.5 en la definición)
          for(int a = 0; a < 4; a++) {
              B.Set(3, 2*a, jac_result.dNdX.getVal(1, a));      // ∂vr/∂z
              B.Set(3, 2*a + 1, jac_result.dNdX.getVal(0, a)); // ∂vz/∂r
          }    
          // --- VECTOR C DE MARTINS (ECUACIÓN 4.49) ---
          // C = [1, 1, 1, 0]^T  (rr, zz, θθ, rz)
          std::vector<double> C_vec = {1.0, 1.0, 1.0, 0.0};
          
          // --- CALCULAR B^T * C (VECTOR DE 8 COMPONENTES) ---
          std::vector<double> BT_C(8, 0.0);
          
          // B^T es de dimensión 8x4, C es 4x1, resultado 8x1
          // BT_C[i] = Σ_{k=0}^{3} B^T[i][k] * C[k] = Σ_{k=0}^{3} B[k][i] * C[k]
          for(int i = 0; i < 8; i++) {           // i = DOF local
              double sum = 0.0;
              for(int k = 0; k < 4; k++) {       // k = componente de deformación
                  sum += B.getVal(k, i) * C_vec[k];
              }
              BT_C[i] = sum;
          }
          
          // --- ENSAMBLAR MATRIZ Q ELEMENTAL ---
          // Q_e[:, 0] += ∫ (B^T C) dΩ
          for(int i = 0; i < 8; i++) {
              double current_val = Q_e.getVal(i, 0);
              Q_e.Set(i, 0, current_val + BT_C[i] * dOmega);
          }
          
          // --- DIAGNÓSTICO: CALCULAR DIVERGENCIA (OPCIONAL) ---
          // Solo para verificar incompresibilidad, no afecta el ensamblaje
          double dvr_dr = 0.0, dvz_dz = 0.0, vr_interp = 0.0;
          for(int a = 0; a < 4; a++) {
              dvr_dr += jac_result.dNdX.getVal(0, a) * vel_elem[2*a];
              dvz_dz += jac_result.dNdX.getVal(1, a) * vel_elem[2*a + 1];
              vr_interp += jac_result.N[a] * vel_elem[2*a];
          }
          
          double div_v;
          if(r_gp < 1e-8) {
              div_v = 2.0 * dvr_dr + dvz_dz;
          } else {
              div_v = dvr_dr + vr_interp/r_gp + dvz_dz;
          }
          
          incompressibility_error += std::abs(div_v) * dOmega;
          max_div_v = max(max_div_v, std::abs(div_v));
      }
        
        // ===========================================
        // ENSAMBLAJE GLOBAL CON PENALIZACIÓN VOLUMÉTRICA
        // ===========================================
        
        // 1. Bloque P (parte de velocidades)
        for(int i_local = 0; i_local < 8; i_local++) {
            int i_global = dofs_v[i_local];
            
            for(int j_local = 0; j_local < 8; j_local++) {
                int j_global = dofs_v[j_local];
                
                double val = K_glob.getVal(i_global, j_global)
                           + P_e.getVal(i_local, j_local);
                K_glob.Set(i_global, j_global, val);
            }
            
            // 2. Bloque Q (acoplamiento velocidad-presión)
            K_glob.Set(i_global, dof_p, 
                      K_glob.getVal(i_global, dof_p) + Q_e.getVal(i_local, 0));
            
            // 3. Bloque Q^T (transpuesta) 
            K_glob.Set(dof_p, i_global,
                      K_glob.getVal(dof_p, i_global) + Q_e.getVal(i_local, 0));
        }
        
        // 4. BLOQUE DE PENALIZACIÓN VOLUMÉTRICA K_pp (¡NUEVO!)
        // Agregar el término dΩ/κ en la diagonal de presión
        double dOmega_penalty = 0.0;
        
        // Calcular dΩ para la penalización (usando el punto reducido)
        auto jac_penalty = jacobian_and_gradients(pos, 0.0, 0.0);
        double r_penalty = 0.0;
        for(int i = 0; i < 4; i++) {
            r_penalty += jac_penalty.N[i] * pos[i].x;
        }
        r_penalty = max(r_penalty, 1e-12);
        dOmega_penalty = 2.0 * M_PI * r_penalty * jac_penalty.detJ * 4.0; // w_reduced[0] = 4.0
        
        // Agregar penalización volumétrica: ∫ (1/κ) dΩ
        //double current_val = K_glob.getVal(dof_p, dof_p);
        //K_glob.Set(dof_p, dof_p, current_val + dOmega_penalty / kappa);

        if(e_idx == 0) {  // Solo primer elemento
            std::cout << "\n=== P_e ELEMENTO 0 (ANTES DE ENSAMBLAR) ===" << std::endl;
            P_e.Print();
            std::cout << "=== FIN P_e ===\n" << std::endl;
            
            std::cout << "\n=== Q_e ELEMENTO 0 ===" << std::endl;
            Q_e.Print();
            std::cout << "=== FIN Q_e ===\n" << std::endl;
            
            std::cout << "dOmega_penalty = " << dOmega_penalty << std::endl;
        }


        
    } // Fin del bucle de elementos

    //~ for(int i = ndof_v; i < ndof_total; i++) {
        //~ double diag = K_glob.getVal(i, i);
        //~ if(fabs(diag) < 1e-12) {
            //~ K_glob.Set(i, i, 1e-8);  // Regularización
        //~ }
    //~ }
    
    //~ // Pequeña regularización por si acaso (mejor estabilidad numérica)
    //~ for(int i = ndof_v; i < ndof_total; ++i) {
        //~ double diag = K_glob.getVal(i, i);
        //~ if(std::abs(diag) < 1e-12) {
            //~ K_glob.Set(i, i, diag + 1e-12);
        //~ }
    //~ }

}


// ================================
// APLICAR CONDICIONES DE CONTORNO AXISIMÉTRICAS (FUNCIÓN LIBRE)
// ================================
void apply_boundary_conditions_axisym(Matrix& K_glob,
                                     std::vector<double>& F_glob,
                                     int ndof_v, 
                                     int ndof_total,
                                     int* bcx_nod, double* bcx_val, int n_bcx,
                                     int* bcy_nod, double* bcy_val, int n_bcy) {
    
    std::cout << "\n🔧 APLICANDO CONDICIONES AXISIMÉTRICAS:" << std::endl;
    std::cout << "   vr (X): " << n_bcx << " condiciones" << std::endl;
    std::cout << "   vz (Y): " << n_bcy << " condiciones" << std::endl;
    
    int count_applied = 0;
    double v_z_top_expected = -1.0; // De tu código original
    
    // Aplicar vr (dof = 2*nodo)
    for(int i = 0; i < n_bcx; i++) {
        int node = bcx_nod[i];
        double value = bcx_val[i];
        int dof = 2 * node;  // DOF de vr
        
        // Mostrar primeras 5 condiciones para debug
        if(i < 5) {
            std::cout << "   vr: nodo " << node << " -> dof " << dof 
                      << " = " << value << std::endl;
        }
        
        // Restar contribución de la columna
        for(int j = 0; j < ndof_total; j++) {
            F_glob[j] -= K_glob.getVal(j, dof) * value;
        }
        
        // Poner fila y columna en cero
        for(int j = 0; j < ndof_total; j++) {
            K_glob.Set(dof, j, 0.0);
            K_glob.Set(j, dof, 0.0);
        }
        
        // Poner 1 en la diagonal y asignar valor
        K_glob.Set(dof, dof, 1.0);
        F_glob[dof] = value;
        count_applied++;
    }
    
    // Aplicar vz (dof = 2*nodo + 1)
    int vz_top_count = 0;
    for(int i = 0; i < n_bcy; i++) {
        int node = bcy_nod[i];
        double value = bcy_val[i];
        int dof = 2 * node + 1;  // DOF de vz
        
        if(i < 5) {
            std::cout << "   vz: nodo " << node << " -> dof " << dof 
                      << " = " << value << std::endl;
        }
        
        if(std::abs(value - v_z_top_expected) < 1e-6) {
            vz_top_count++;
        }
        
        // Restar contribución de la columna
        for(int j = 0; j < ndof_total; j++) {
            F_glob[j] -= K_glob.getVal(j, dof) * value;
        }
        
        // Poner fila y columna en cero
        for(int j = 0; j < ndof_total; j++) {
            K_glob.Set(dof, j, 0.0);
            K_glob.Set(j, dof, 0.0);
        }
        
        // Poner 1 en la diagonal y asignar valor
        K_glob.Set(dof, dof, 1.0);
        F_glob[dof] = value;
        count_applied++;
    }
    
    std::cout << "   " << vz_top_count << " vz = " << v_z_top_expected 
              << " (COMPRESIÓN)" << std::endl;
    std::cout << "   TOTAL: " << count_applied << " condiciones aplicadas" << std::endl;
}



    SolveResult Domain_d::solve_step_martins(std::vector<double>& vel_guess,
                                  std::vector<double>& press_guess,
                                  Matrix& K_temp,
                                  std::vector<double>& F_temp) {
        
        std::vector<double> vel_prev = vel_guess;
        std::vector<double> press_prev = press_guess;
        
        SolveResult result;
        
    for(int iter = 0; iter < max_iter; iter++) {
        // Ensamblar sistema según Martins
        double max_mu, incompress_error, max_div;

        assemble_martins_system(vel_guess, press_guess, K_temp, F_temp, 
                               incompress_error, max_div);


          std::cout << "\n" << std::string(70, '=') << std::endl;
          std::cout << " DEBUG - MATRIZ GLOBAL C++ (primeros 10x10)" << std::endl;
          std::cout << std::string(70, '=') << std::endl;
          
          // Bloque velocidad-velocidad
          std::cout << "\n BLOQUE V-V (P_glob) - 5x5 superior izquierda:" << std::endl;
          for(int i = 0; i < std::min(5, ndof_v); i++) {
              for(int j = 0; j < std::min(5, ndof_v); j++) {
                  printf("%12.6e ", K_temp.getVal(i, j));
              }
              std::cout << std::endl;
          }
          
          // Bloque velocidad-presión
          std::cout << "\n BLOQUE V-P (Q_glob) - 5x5:" << std::endl;
          for(int i = 0; i < std::min(5, ndof_v); i++) {
              for(int j = 0; j < std::min(5, ndof_p); j++) {
                  printf("%12.6e ", K_temp.getVal(i, ndof_v + j));
              }
              std::cout << std::endl;
          }
          
          // Bloque presión-velocidad (Q^T)
          std::cout << "\n BLOQUE P-V (Q^T) - 5x5:" << std::endl;
          for(int i = 0; i < std::min(5, ndof_p); i++) {
              for(int j = 0; j < std::min(5, ndof_v); j++) {
                  printf("%12.6e ", K_temp.getVal(ndof_v + i, j));
              }
              std::cout << std::endl;
          }
          
          // Bloque presión-presión (K_pp)
          std::cout << "\n BLOQUE P-P (K_pp) - primeros 5 diag:" << std::endl;
          for(int i = 0; i < std::min(5, ndof_p); i++) {
              double val = K_temp.getVal(ndof_v + i, ndof_v + i);
              printf("  K_pp[%d,%d] = %12.6e\n", i, i, val);
          }

        
        // APLICAR CONDICIONES USANDO LOS MIEMBROS DE LA CLASE
        apply_boundary_conditions_axisym(K_temp, F_temp, 
                                        ndof_v, ndof_total,
                                        this->bcx_nod, this->bcx_val, this->bc_count[0],
                                        this->bcy_nod, this->bcy_val, this->bc_count[1]);
        
        try {
            // Resolver sistema de punto de silla
            std::vector<double> sol = solve_linear_system(K_temp, F_temp);
            
            // Extraer solución
            std::vector<double> vel_new(ndof_v);
            std::vector<double> press_new(ndof_p);
            
            for(int i = 0; i < ndof_v; i++) vel_new[i] = sol[i];
            for(int i = 0; i < ndof_p; i++) press_new[i] = sol[ndof_v + i];
            
            // Under-relaxation
            std::vector<double> vel_new_relax(ndof_v);
            std::vector<double> press_new_relax(ndof_p);
            
            for(int i = 0; i < ndof_v; i++) {
                vel_new_relax[i] = (1.0 - omega_v) * vel_guess[i] 
                                 + omega_v * vel_new[i];
            }
            for(int i = 0; i < ndof_p; i++) {
                press_new_relax[i] = (1.0 - omega_p) * press_guess[i] 
                                   + omega_p * press_new[i];
            }
            
            // Verificar convergencia
            double dv_norm = 0.0;
            double vel_norm = 0.0;
            for(int i = 0; i < ndof_v; i++) {
                double diff = vel_new_relax[i] - vel_guess[i];
                dv_norm += diff * diff;
                vel_norm += vel_new_relax[i] * vel_new_relax[i];
            }
            dv_norm = sqrt(dv_norm);
            vel_norm = sqrt(vel_norm);
            double dv_rel = dv_norm / (vel_norm + 1e-12);
            
            double dp_norm = 0.0;
            double press_norm = 0.0;
            for(int i = 0; i < ndof_p; i++) {
                double diff = press_new_relax[i] - press_guess[i];
                dp_norm += diff * diff;
                press_norm += press_new_relax[i] * press_new_relax[i];
            }
            dp_norm = sqrt(dp_norm);
            press_norm = sqrt(press_norm);
            double dp_rel = dp_norm / (press_norm + 1e-12);
            
            // Calcular error de incompresibilidad promedio
            //~ double domain_volume = M_PI * Lr * Lr * Lz;
            //~ double avg_div = incompress_error / (domain_volume + 1e-12);

            double current_volume = 0.0;
            for(int e_idx = 0; e_idx < m_elem_count; e_idx++) {
                const auto& conn = elements[e_idx];
                std::vector<Point2D> pos(4);
                for(int i = 0; i < 4; i++) pos[i] = coords[conn[i]];
                
                auto jac = jacobian_and_gradients(pos, 0.0, 0.0);
                double r = 0.0;
                for(int i = 0; i < 4; i++) r += jac.N[i] * pos[i].x;
                current_volume += 2.0 * M_PI * max(r, 1e-12) * jac.detJ * 4.0;
            }

            // Normalizar con volumen actual
            double avg_div = incompress_error / (current_volume + 1e-12);
            
            std::cout << "    Iter " << iter+1 
                      << ": Δv=" << dv_rel 
                      << ", Δp=" << dp_rel 
                      << ", div(v)=" << avg_div
                      << ", μ_max=" << max_mu/1e6 << " MPa·s" << std::endl;
            
            if(dv_rel < tol && dp_rel < tol && avg_div < 1e-2) {
                result.velocity = vel_new_relax;
                result.pressure = press_new_relax;
                result.converged = true;
                result.iterations = iter + 1;
                result.div_error = avg_div;
                result.max_div = max_div;
                return result;
            }
            
            // Actualizar para siguiente iteración
            vel_guess = vel_new_relax;
            press_guess = press_new_relax;
            
        } catch(const std::exception& e) {
            std::cout << "    Error en iteración " << iter+1 << ": " << e.what() << std::endl;
            result.velocity = vel_prev;
            result.pressure = press_prev;
            result.converged = false;
            result.iterations = iter + 1;
            return result;
        }
    }
    
    std::cout << "    No convergió en " << max_iter << " iteraciones" << std::endl;
    result.velocity = vel_guess;
    result.pressure = press_guess;
    result.converged = false;
    result.iterations = max_iter;
    return result;
}

//~ // ================================
//~ // VERIFICACIONES FÍSICAS (MARTINS)
//~ // ================================
//~ void perform_physical_checks(const std::vector<double>& vel,
                            //~ const std::vector<double>& press,
                            //~ const std::vector<double>& eps_bar_local) {
    
    //~ std::cout << "\n" << std::string(70, '=') << std::endl;
    //~ std::cout << "VERIFICACIONES FÍSICAS (Martins)" << std::endl;
    //~ std::cout << std::string(70, '=') << std::endl;
    
    //~ // 1. Volumen conservado
    //~ double initial_vol = M_PI * Lr * Lr * Lz;
    //~ double final_height = 0.0, final_radius = 0.0;
    //~ for(const auto& p : coords) {
        //~ final_height = max(final_height, p.y);
        //~ final_radius = max(final_radius, p.x);
    //~ }
    //~ double min_height = 1e100;
    //~ for(const auto& p : coords) {
        //~ min_height = min(min_height, p.y);
    //~ }
    //~ final_height -= min_height;
    
    //~ double final_vol = M_PI * final_radius * final_radius * final_height;
    //~ double vol_change = (final_vol - initial_vol) / initial_vol;
    
    //~ std::cout << "Volumen inicial: " << initial_vol << " m³" << std::endl;
    //~ std::cout << "Volumen final:   " << final_vol << " m³" << std::endl;
    //~ std::cout << "Cambio volumen:  " << vol_change*100 << "%" << std::endl;
    //~ if(std::abs(vol_change) < 0.01) {
        //~ std::cout << "✓ Conservación de volumen ACEPTABLE" << std::endl;
    //~ } else {
        //~ std::cout << "⚠ ADVERTENCIA: Pérdida significativa de volumen" << std::endl;
    //~ }
    
    //~ // 2. Incompresibilidad
    //~ std::cout << "\nINCOMPRESIBILIDAD:" << std::endl;
    //~ double div_v_max = 0.0;
    //~ double div_v_avg = 0.0;
    //~ double total_vol = 0.0;
    
    //~ for(int e_idx = 0; e_idx < m_elem_count; e_idx++) {
        //~ const auto& conn = elements[e_idx];
        //~ std::vector<Point2D> pos(4);
        //~ for(int i = 0; i < 4; i++) {
            //~ pos[i] = coords[conn[i]];
        //~ }
        
        //~ // Velocidades del elemento
        //~ std::vector<double> vel_elem(8);
        //~ for(int i = 0; i < 4; i++) {
            //~ int node = conn[i];
            //~ vel_elem[2*i] = velocity[2*node];
            //~ vel_elem[2*i + 1] = velocity[2*node + 1];
        //~ }
        
        //~ // Punto central
        //~ auto jac_result = jacobian_and_gradients(pos, 0.0, 0.0);
        //~ double r_center = 0.0;
        //~ for(int i = 0; i < 4; i++) {
            //~ r_center += jac_result.N[i] * pos[i].x;
        //~ }
        //~ r_center = max(r_center, 1e-12);
        
        //~ // Calcular divergencia
        //~ double dvr_dr = 0.0, dvz_dz = 0.0, vr_center = 0.0;
        //~ for(int a = 0; a < 4; a++) {
            //~ dvr_dr += jac_result.dNdX.getVal(0, a) * vel_elem[2*a];
            //~ dvz_dz += jac_result.dNdX.getVal(1, a) * vel_elem[2*a + 1];
            //~ vr_center += jac_result.N[a] * vel_elem[2*a];
        //~ }
        
        //~ double div_v;
        //~ if(r_center < 1e-8) {
            //~ div_v = 2.0 * dvr_dr + dvz_dz;
        //~ } else {
            //~ div_v = dvr_dr + vr_center/r_center + dvz_dz;
        //~ }
        
        //~ double elem_vol = jac_result.detJ * 2.0 * M_PI * r_center;
        //~ div_v_avg += std::abs(div_v) * elem_vol;
        //~ total_vol += elem_vol;
        //~ div_v_max = max(div_v_max, std::abs(div_v));
    //~ }
    
    //~ div_v_avg /= total_vol;
    //~ std::cout << "Divergencia promedio: " << div_v_avg << " 1/s" << std::endl;
    //~ std::cout << "Divergencia máxima:   " << div_v_max << " 1/s" << std::endl;
    //~ if(div_v_avg < 1e-2) {
        //~ std::cout << "✓ Incompresibilidad ACEPTABLE" << std::endl;
    //~ } else {
        //~ std::cout << "⚠ ADVERTENCIA: Incompresibilidad débil" << std::endl;
    //~ }
    
    //~ // 3. Verificar viscosidades
    //~ std::vector<double> eps_dot_values, mu_values;
    //~ for(int e_idx = 0; e_idx < m_elem_count; e_idx++) {
        //~ const auto& conn = elements[e_idx];
        //~ std::vector<Point2D> pos(4);
        //~ for(int i = 0; i < 4; i++) {
            //~ pos[i] = coords[conn[i]];
        //~ }
        
        //~ std::vector<double> vel_elem(8);
        //~ for(int i = 0; i < 4; i++) {
            //~ int node = conn[i];
            //~ vel_elem[2*i] = velocity[2*node];
            //~ vel_elem[2*i + 1] = velocity[2*node + 1];
        //~ }
        
        //~ auto jac_result = jacobian_and_gradients(pos, 0.0, 0.0);
        //~ double r_center = 0.0;
        //~ for(int i = 0; i < 4; i++) {
            //~ r_center += jac_result.N[i] * pos[i].x;
        //~ }
        
        //~ auto strain_result = calculate_strain_rate_martins(
            //~ jac_result.dNdX, vel_elem, r_center, jac_result.N
        //~ );
        
        //~ eps_dot_values.push_back(strain_result.eps_dot_eq);
        //~ //mu_values.push_back(effective_viscosity_norton(strain_result.eps_dot_eq));
    //~ }
    
    //~ double eps_min = *std::min_element(eps_dot_values.begin(), eps_dot_values.end());
    //~ double eps_max = *std::max_element(eps_dot_values.begin(), eps_dot_values.end());
    //~ double eps_mean = 0.0;
    //~ for(double val : eps_dot_values) eps_mean += val;
    //~ eps_mean /= eps_dot_values.size();
    
    //~ double mu_min = *std::min_element(mu_values.begin(), mu_values.end());
    //~ double mu_max = *std::max_element(mu_values.begin(), mu_values.end());

    
    //~ std::cout << "\nPROPIEDADES DEL MATERIAL:" << std::endl;
    //~ std::cout << "ε̇_eq: min=" << eps_min << ", max=" << eps_max 
              //~ << ", mean=" << eps_mean << " 1/s" << std::endl;
    
    //~ // 4. Fuerzas estimadas
    //~ double area_top = M_PI * final_radius * final_radius;
    //~ double press_mean = 0.0;
    //~ for(double p : press) press_mean += p;
    //~ press_mean /= press.size();
    //~ double force_top_est = press_mean * area_top;
    
    //~ std::cout << "\nFUERZAS ESTIMADAS:" << std::endl;
    //~ std::cout << "Área de contacto: " << area_top << " m²" << std::endl;
    //~ std::cout << "Presión promedio: " << press_mean/1e6 << " MPa" << std::endl;
    //~ std::cout << "Fuerza total estimada: " << force_top_est/1e3 << " kN" << std::endl;
//~ }




// ================================
// SIMULACIÓN PRINCIPAL (MARTINS)
// ================================

  void Domain_d::Solve_Martins_Picard(){
    
    dt = 0.001;
   
    // Inicialización
    //initialize_mesh();
    coords.reserve(m_node_count);
    elements.reserve(m_elem_count);
    
    // DOFs
    ndof_v = 2 * m_node_count;
    ndof_p = m_elem_count;
    ndof_total = ndof_v + ndof_p;
    
    // Inicializar variables
    velocity.resize(ndof_v, 0.0);
    pressure.resize(ndof_p, 0.0);
    eps_bar.resize(m_elem_count, 0.0);
    
    // Guardar malla inicial
    coords_history.push_back(coords);

    
    for (int n=0;n<m_node_count;n++){
      coords[n].x = x[m_dim*n  ];
      coords[n].y = x[m_dim*n+1];
    }
    for (int e=0;e<m_elem_count;e++){
      elements[e].resize(4);
      for (int i=0;i<4;i++){
          elements[e][i] = m_elnod[m_nodxelem*e+i];
      }
    }

    //auto fixed_dofs = setup_boundary_conditions();
       
    
    std::cout << std::string(70, '=') << std::endl;
    std::cout << "SIMULACIÓN AXISIMÉTRICA DE FORJA - ESQUEMA MARTINS" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
    std::cout << "Material: Norton-Hoff, K=" << Kmat << " Pa·s^" << nexp 
              << ", n=" << nexp << std::endl;
    std::cout << "Elementos: " << m_elem_count << ", Nodos: " << m_node_count << std::endl;
    std::cout << "DOFs: " << ndof_total << " (" << ndof_v 
              << " velocidades + " << ndof_p << " presiones P0)" << std::endl;
    std::cout << "Esquema: Martins (4.55-4.56)" << std::endl;
    std::cout << "Integración: P(4 puntos) + Q(1 punto reducido)" << std::endl;

    ///// SOLVER THINGS.
    Solver_Eigen *solver = new Solver_Eigen();
    m_solver = solver;
    m_solver->setDomain(this);
    
    m_solver->Allocate();                  
    


    // Solver_Eigen *solver = new Solver_Eigen();
    // m_solver = solver;
    // m_solver->setDomain(this);   // Le pasás la malla, DOFs, etc.
    // m_solver->Allocate();         // Reserva memoria
    
    // Matrices temporales para reuso
    Matrix K_temp(ndof_total, ndof_total);
    std::vector<double> F_temp(ndof_total, 0.0);
    
    // Bucle temporal
    for(int step = 0; step < nsteps; step++) {
        std::cout << "\n--- PASO " << step+1 << "/" << nsteps 
                  << " (dt=" << dt << " s) ---" << std::endl;
        
        auto result = solve_step_martins(velocity, pressure, K_temp, F_temp);
        
        velocity = result.velocity;
        pressure = result.pressure;
        
        if(!result.converged) {
            std::cout << "  ⚠ ADVERTENCIA: Paso no convergido completamente!" << std::endl;
        }
        
        // Actualizar deformación acumulada
        for(int e_idx = 0; e_idx < m_elem_count; e_idx++) {
            const auto& conn = elements[e_idx];
            std::vector<Point2D> pos(4);
            for(int i = 0; i < 4; i++) {
                pos[i] = coords[conn[i]];
            }
            
            // Velocidades del elemento
            std::vector<double> vel_elem(8);
            for(int i = 0; i < 4; i++) {
                int node = conn[i];
                vel_elem[2*i] = velocity[2*node];
                vel_elem[2*i + 1] = velocity[2*node + 1];
            }
            
            // Punto central (ξ=0, η=0)
            auto jac_result = jacobian_and_gradients(pos, 0.0, 0.0);
            double r_center = 0.0;
            for(int i = 0; i < 4; i++) {
                r_center += jac_result.N[i] * pos[i].x;
            }
            
            auto strain_result = calculate_strain_rate_martins(
                jac_result.dNdX, vel_elem, r_center, jac_result.N
            );
            
            eps_bar[e_idx] += dt * strain_result.eps_dot_eq;
        }
        
        // Updated Lagrangian: actualizar coordenadas
        for(int i = 0; i < m_node_count; i++) {
            coords[i].x += dt * velocity[2*i];      // r
            coords[i].y += dt * velocity[2*i + 1];  // z
        }
        coords_history.push_back(coords);
        
        // Estadísticas
        double vr_min = 1e100, vr_max = -1e100;
        double vz_min = 1e100, vz_max = -1e100;
        double p_min = 1e100, p_max = -1e100;
        
        for(int i = 0; i < m_node_count; i++) {
            vr_min = min(vr_min, velocity[2*i]);
            vr_max = max(vr_max, velocity[2*i]);
            vz_min = min(vz_min, velocity[2*i + 1]);
            vz_max = max(vz_max, velocity[2*i + 1]);
        }
        
        for(int i = 0; i < ndof_p; i++) {
            p_min = min(p_min, pressure[i]);
            p_max = max(p_max, pressure[i]);
        }
        
        std::cout << "  v_r: [" << vr_min << ", " << vr_max << "] m/s" << std::endl;
        std::cout << "  v_z: [" << vz_min << ", " << vz_max << "] m/s" << std::endl;
        std::cout << "  Presión: [" << p_min/1e6 << ", " << p_max/1e6 << "] MPa" << std::endl;
        std::cout << "  Presión media: " << (p_min + p_max)/(2.0*1e6) << " MPa" << std::endl;
    }
    
    cout << "Writing VTK"<<endl;
        // ESCRIBIR VTK DESPUÉS DE CADA PASO
        VTKWriter::writeVtkFile("forja", 1, coords, elements, 
                               velocity, pressure, eps_bar, 
                               m_node_count, m_elem_count);
    cout << "DONE."<<endl;
    // Verificaciones físicas finales
    //perform_physical_checks(velocity, pressure, eps_bar);
    
    std::cout << "\n" << std::string(70, '=') << std::endl;
    std::cout << "SIMULACIÓN COMPLETADA - ESQUEMA MARTINS OPTIMIZADO" << std::endl;
    std::cout << std::string(70, '=') << std::endl;
}



}//MetFEM
