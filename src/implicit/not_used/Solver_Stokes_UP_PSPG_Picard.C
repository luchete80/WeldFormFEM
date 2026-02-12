// ================================
// stokes_pspg_plainstrain.cpp
// Migración completa de Python a C++ con Plain Strain
// Usando tu clase Matrix
// ================================

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <memory>

// ================================
// PARÁMETROS FÍSICAMENTE CORRECTOS
// ================================
int nelr = 15, nelz = 15;
double Lx = 0.0127, Ly = 0.03;  // Dimensiones en X e Y (plain strain)

// MATERIAL: Acero caliente tipo Norton-Hoff
double Kmat = 50.0e6;          // Coeficiente de consistencia [Pa·s^n]
double nexp = 0.2;             // Exponente strain-rate sensitivity
double mu_ref = Kmat / 3.0;    // Viscosidad de referencia

// Parámetros numéricos
double beta_stab = 0.1;
double relax_factor = 0.6;
double tol = 1e-8;
int max_iter = 50;
double dt = 0.001;
int nsteps = 4;
double v_y_top = -1.0;  // Compresión en dirección Y

// ================================
// ESTRUCTURAS DE DATOS
// ================================
struct Point2D {
    double x, y;
    Point2D(double x_=0, double y_=0) : x(x_), y(y_) {}
};

// Variables globales
int nx_nodes, ny_nodes, nnodes, nelem;
int ndof_v, ndof_p, ndof_total;

std::vector<Point2D> coords;
std::vector<std::vector<int>> elements;
std::vector<double> velocity;    // [vx0, vy0, vx1, vy1, ...]
std::vector<double> pressure;    // Presión por elemento
std::vector<double> eps_bar;     // Deformación plástica acumulada
std::vector<std::vector<Point2D>> coords_history;

// ================================
// FUNCIONES AUXILIARES
// ================================
inline double max(double a, double b) { return a > b ? a : b; }
inline double min(double a, double b) { return a < b ? a : b; }

// ================================
// FUNCIONES DE FORMA Q1 (PLAIN STRAIN)
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
        pos_mat.Set(i, 0, pos[i].x);
        pos_mat.Set(i, 1, pos[i].y);
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

// ================================
// VISCOPLASTICIDAD NORTON-HOFF (PLAIN STRAIN)
// ================================
double effective_viscosity_norton(double eps_dot_eq) {
    // Evitar divisiones por cero
    double eps_dot_eq_safe = max(eps_dot_eq, 1e-12);
    return (Kmat / 3.0) * pow(eps_dot_eq_safe, nexp - 1.0);
}

struct StrainRateResult {
    double eps_dot_eq;
    std::vector<double> strain_vec;  // [exx, eyy, exy] para plain strain
    double eps_mean;
};

StrainRateResult calculate_strain_rate_plain_strain(const Matrix& dNdX,
                                                   const std::vector<double>& vel_elem,
                                                   const std::vector<double>& N) {
    StrainRateResult result;
    result.strain_vec.resize(3);
    
    // Inicializar gradientes
    double dvx_dx = 0.0, dvx_dy = 0.0;
    double dvy_dx = 0.0, dvy_dy = 0.0;
    
    for(int a = 0; a < 4; a++) {
        double vx = vel_elem[2*a];      // vx del nodo a
        double vy = vel_elem[2*a + 1];  // vy del nodo a
        
        dvx_dx += dNdX.getVal(0, a) * vx;
        dvx_dy += dNdX.getVal(1, a) * vx;
        dvy_dx += dNdX.getVal(0, a) * vy;
        dvy_dy += dNdX.getVal(1, a) * vy;
    }
    
    // Tensor de deformación para plain strain (2D)
    double exx = dvx_dx;
    double eyy = dvy_dy;
    double exy = 0.5 * (dvx_dy + dvy_dx);
    
    // En plain strain, ezz = 0 (deformación en z es cero)
    // Pero contribuye a la parte volumétrica
    double ezz = 0.0;
    
    // Parte esférica (volumétrica)
    result.eps_mean = (exx + eyy + ezz) / 3.0;
    
    // Parte desviadora
    double exx_dev = exx - result.eps_mean;
    double eyy_dev = eyy - result.eps_mean;
    double ezz_dev = ezz - result.eps_mean;
    double exy_dev = exy;  // La componente cortante ya es desviadora
    
    // Invariante J2 para plain strain
    // J2 = 1/2 * (exx_dev^2 + eyy_dev^2 + ezz_dev^2) + exy^2
    double J2 = 0.5 * (exx_dev*exx_dev + eyy_dev*eyy_dev + ezz_dev*ezz_dev) 
                + exy_dev*exy_dev;
    
    // Strain-rate equivalente: ε̇_eq = √(3J2)
    result.eps_dot_eq = sqrt(3.0 * J2);
    
    result.strain_vec[0] = exx;
    result.strain_vec[1] = eyy;
    result.strain_vec[2] = exy;
    
    return result;
}

// ================================
// ENSAMBLAJE PSPG (PLAIN STRAIN)
// ================================
void assemble_mixed_system_correct(const std::vector<double>& vel_vec,
                                  const std::vector<double>& press_vec,
                                  Matrix& K_glob,
                                  std::vector<double>& F_glob,
                                  double& max_mu_eff) {
    
    // Inicializar matriz global (sparse, pero usamos densa por simplicidad inicial)
    K_glob.SetZero();
    std::fill(F_glob.begin(), F_glob.end(), 0.0);
    max_mu_eff = 0.0;
    
    // Puntos de Gauss
    const double a = 1.0 / sqrt(3.0);
    std::vector<std::pair<double, double>> gp_full = {
        {-a, -a}, {a, -a}, {a, a}, {-a, a}
    };
    std::vector<double> w_full = {1.0, 1.0, 1.0, 1.0};
    
    for(int e_idx = 0; e_idx < nelem; e_idx++) {
        const auto& conn = elements[e_idx];
        
        // Posiciones de los nodos del elemento
        std::vector<Point2D> pos(4);
        for(int i = 0; i < 4; i++) {
            pos[i] = coords[conn[i]];
        }
        
        // DOFs de velocidades del elemento
        std::vector<int> dofs_v;
        for(int n : conn) {
            dofs_v.push_back(2*n);      // vx del nodo n
            dofs_v.push_back(2*n + 1);  // vy del nodo n
        }
        
        int dof_p = ndof_v + e_idx;
        
        // Extraer velocidades del elemento
        std::vector<double> vel_elem(8);
        for(int i = 0; i < 8; i++) {
            vel_elem[i] = vel_vec[dofs_v[i]];
        }
        
        double press_elem = press_vec[e_idx];
        
        // Matrices elementales
        Matrix K_e(8, 8);
        Matrix G_e(8, 1);
        Matrix Kstab_e(8, 8);
        K_e.SetZero();
        G_e.SetZero();
        Kstab_e.SetZero();
        
        // Tamaño característico del elemento
        double h_x = fabs(pos[1].x - pos[0].x) + fabs(pos[2].x - pos[3].x);
        double h_y = fabs(pos[3].y - pos[0].y) + fabs(pos[2].y - pos[1].y);
        double h_e = 0.5 * sqrt(h_x*h_x + h_y*h_y);
        
        double area_elem = 0.0;
        double mu_eff_max_elem = 0.0;
        
        // Integración en puntos de Gauss
        for(int gp = 0; gp < 4; gp++) {
            double xi = gp_full[gp].first;
            double eta = gp_full[gp].second;
            double w = w_full[gp];
            
            auto jac_result = jacobian_and_gradients(pos, xi, eta);
            
            // dOmega = detJ * w (sin factor 2πr de axisimetría)
            double dOmega = jac_result.detJ * w;
            area_elem += dOmega;
            
            // Calcular strain-rate y viscosidad
            auto strain_result = calculate_strain_rate_plain_strain(
                jac_result.dNdX, vel_elem, jac_result.N
            );
            
            double mu_eff = effective_viscosity_norton(strain_result.eps_dot_eq);
            mu_eff_max_elem = max(mu_eff_max_elem, mu_eff);
            
            // 1. Matriz B (3x8) para plain strain
            Matrix B(3, 8);
            B.SetZero();
            
            // ε_xx = ∂vx/∂x
            for(int a = 0; a < 4; a++) {
                B.Set(0, 2*a, jac_result.dNdX.getVal(0, a));
            }
            
            // ε_yy = ∂vy/∂y
            for(int a = 0; a < 4; a++) {
                B.Set(1, 2*a + 1, jac_result.dNdX.getVal(1, a));
            }
            
            // γ_xy = ∂vx/∂y + ∂vy/∂x
            for(int a = 0; a < 4; a++) {
                B.Set(2, 2*a, jac_result.dNdX.getVal(1, a));      // ∂vx/∂y
                B.Set(2, 2*a + 1, jac_result.dNdX.getVal(0, a));  // ∂vy/∂x
            }
            
            // 2. Tensor constitutivo D para plain strain (3x3)
            // Para material incompresible con penalización
            Matrix D(3, 3);
            D.SetZero();
            
            // Forma desviadora para viscosidad
            D.Set(0, 0, 2.0 * mu_eff);
            D.Set(1, 1, 2.0 * mu_eff);
            D.Set(2, 2, mu_eff);  // Componente cortante
            
            // 3. Contribución viscosa: B^T * D * B
            Matrix Bt = B.getTranspose();
            Matrix DxB = Matrix(3, 8);
            
            // Multiplicación D * B
            for(int i = 0; i < 3; i++) {
                for(int j = 0; j < 8; j++) {
                    double sum = 0.0;
                    for(int k = 0; k < 3; k++) {
                        sum += D.getVal(i, k) * B.getVal(k, j);
                    }
                    DxB.Set(i, j, sum);
                }
            }
            
            // B^T * (D * B)
            for(int i = 0; i < 8; i++) {
                for(int j = 0; j < 8; j++) {
                    double sum = 0.0;
                    for(int k = 0; k < 3; k++) {
                        sum += Bt.getVal(i, k) * DxB.getVal(k, j);
                    }
                    K_e.Set(i, j, K_e.getVal(i, j) + sum * dOmega);
                }
            }
            
            // 4. Matriz de acoplamiento presión-velocidad (B_vol)
            // Para plain strain: div(v) = ∂vx/∂x + ∂vy/∂y
            Matrix B_vol(1, 8);
            B_vol.SetZero();
            
            for(int a = 0; a < 4; a++) {
                B_vol.Set(0, 2*a, jac_result.dNdX.getVal(0, a));      // ∂vx/∂x
                B_vol.Set(0, 2*a + 1, jac_result.dNdX.getVal(1, a));  // ∂vy/∂y
            }
            
            // G = B_vol^T * N_p (N_p = 1 para P0)
            Matrix B_vol_t = B_vol.getTranspose();
            for(int i = 0; i < 8; i++) {
                G_e.Set(i, 0, G_e.getVal(i, 0) + B_vol_t.getVal(i, 0) * dOmega);
            }
            
            // 5. Estabilización PSPG
            double tau = beta_stab * (h_e*h_e) / (2.0 * mu_eff + 1e-12);
            
            // Kstab = tau * B_vol^T * B_vol
            for(int i = 0; i < 8; i++) {
                for(int j = 0; j < 8; j++) {
                    double sum = 0.0;
                    for(int k = 0; k < 1; k++) {
                        sum += B_vol_t.getVal(i, k) * B_vol.getVal(k, j);
                    }
                    Kstab_e.Set(i, j, Kstab_e.getVal(i, j) + tau * sum * dOmega);
                }
            }
        }
        
        // Penalización volumétrica
        double Kbulk_eff = 1000.0 * mu_eff_max_elem;
        double M_p_elem = area_elem / Kbulk_eff;
        
        max_mu_eff = max(max_mu_eff, mu_eff_max_elem);
        
        // ENSAMBLAJE GLOBAL
        // 1. Parte de velocidades (K_e + Kstab_e)
        for(int i_local = 0; i_local < 8; i_local++) {
            int i_global = dofs_v[i_local];
            
            for(int j_local = 0; j_local < 8; j_local++) {
                int j_global = dofs_v[j_local];
                
                double val = K_glob.getVal(i_global, j_global)
                           + K_e.getVal(i_local, j_local)
                           + Kstab_e.getVal(i_local, j_local);
                
                K_glob.Set(i_global, j_global, val);
            }
            
            // Acoplamiento velocidad-presión (G)
            K_glob.Set(i_global, dof_p, 
                      K_glob.getVal(i_global, dof_p) + G_e.getVal(i_local, 0));
            
            // Transpuesta de G
            K_glob.Set(dof_p, i_global,
                      K_glob.getVal(dof_p, i_global) + G_e.getVal(i_local, 0));
        }
        
        // Penalización de presión
        K_glob.Set(dof_p, dof_p,
                  K_glob.getVal(dof_p, dof_p) - M_p_elem);
    }
}

// ================================
// CONDICIONES DE CONTORNO (PLAIN STRAIN)
// ================================
std::unordered_map<int, double> setup_boundary_conditions() {
    std::unordered_map<int, double> fixed_dofs;
    
    // Base inferior fija (y = 0)
    for(int i = 0; i < nx_nodes; i++) {
        int node_id = i;  // Primera fila
        fixed_dofs[2*node_id] = 0.0;      // vx = 0
        fixed_dofs[2*node_id + 1] = 0.0;  // vy = 0
    }
    
    // Tapa superior: velocidad impuesta (y = Ly)
    for(int i = 0; i < nx_nodes; i++) {
        int node_id = (ny_nodes - 1) * nx_nodes + i;
        fixed_dofs[2*node_id] = 0.0;           // vx = 0
        fixed_dofs[2*node_id + 1] = v_y_top;   // vy = velocidad impuesta
    }
    
    // Lado izquierdo: simetría o condición de deslizamiento
    for(int j = 0; j < ny_nodes; j++) {
        int node_id = j * nx_nodes;
        fixed_dofs[2*node_id] = 0.0;  // vx = 0 (simetría)
        // vy libre
    }
    
    // Lado derecho: libre (o condición según necesites)
    // Para compresión pura, podrías dejar libre o imponer vx = 0
    
    return fixed_dofs;
}

void apply_boundary_conditions(Matrix& K_glob,
                              std::vector<double>& F_glob,
                              const std::unordered_map<int, double>& fixed_dofs) {
    
    for(const auto& [dof, value] : fixed_dofs) {
        // Restar contribución de la columna
        for(int i = 0; i < ndof_total; i++) {
            F_glob[i] -= K_glob.getVal(i, dof) * value;
        }
        
        // Poner fila y columna en cero
        for(int i = 0; i < ndof_total; i++) {
            K_glob.Set(dof, i, 0.0);
            K_glob.Set(i, dof, 0.0);
        }
        
        // Poner 1 en la diagonal
        K_glob.Set(dof, dof, 1.0);
        F_glob[dof] = value;
    }
}

// ================================
// SOLUCIÓN CON RELAJACIÓN (PICARD)
// ================================
struct SolveResult {
    std::vector<double> velocity;
    std::vector<double> pressure;
    bool converged;
    int iterations;
};

SolveResult solve_step_with_relaxation(std::vector<double>& vel_guess,
                                      std::vector<double>& press_guess,
                                      const std::unordered_map<int, double>& fixed_dofs,
                                      Matrix& K_temp,  // Matriz temporal para reuso
                                      std::vector<double>& F_temp) {  // Vector temporal
    
    std::vector<double> vel_prev = vel_guess;
    std::vector<double> press_prev = press_guess;
    
    SolveResult result;
    
    for(int iter = 0; iter < max_iter; iter++) {
        // Ensamblar sistema
        double max_mu;
        assemble_mixed_system_correct(vel_guess, press_guess, K_temp, F_temp, max_mu);
        
        // Aplicar condiciones de contorno
        apply_boundary_conditions(K_temp, F_temp, fixed_dofs);
        
        try {
            // Resolver sistema (usando tu solver)
            // Asumo que tienes una función solve_linear_system
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
                vel_new_relax[i] = (1.0 - relax_factor) * vel_guess[i] 
                                 + relax_factor * vel_new[i];
            }
            
            for(int i = 0; i < ndof_p; i++) {
                press_new_relax[i] = (1.0 - relax_factor) * press_guess[i] 
                                   + relax_factor * press_new[i];
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
            
            std::cout << "    Iter " << iter+1 
                      << ": Δv=" << dv_rel 
                      << ", Δp=" << dp_rel 
                      << ", μ_max=" << max_mu << " Pa·s" << std::endl;
            
            if(dv_rel < tol && dp_rel < tol) {
                result.velocity = vel_new_relax;
                result.pressure = press_new_relax;
                result.converged = true;
                result.iterations = iter + 1;
                return result;
            }
            
            // Actualizar para siguiente iteración
            vel_guess = vel_new_relax;
            press_guess = press_new_relax;
            vel_prev = vel_guess;
            press_prev = press_guess;
            
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

// ================================
// INICIALIZACIÓN DE MALLADO
// ================================
void initialize_mesh() {
    nx_nodes = nelr + 1;
    ny_nodes = nelz + 1;
    nnodes = nx_nodes * ny_nodes;
    nelem = nelr * nelz;
    
    // DOFs
    ndof_v = 2 * nnodes;
    ndof_p = nelem;
    ndof_total = ndof_v + ndof_p;
    
    // Coordenadas
    coords.reserve(nnodes);
    for(int j = 0; j < ny_nodes; j++) {
        for(int i = 0; i < nx_nodes; i++) {
            double x = i * Lx / nelr;
            double y = j * Ly / nelz;
            coords.push_back(Point2D(x, y));
        }
    }
    
    // Elementos Q4
    elements.reserve(nelem);
    for(int j = 0; j < nelz; j++) {
        for(int i = 0; i < nelr; i++) {
            int n1 = j * nx_nodes + i;
            int n2 = n1 + 1;
            int n3 = n2 + nx_nodes;
            int n4 = n1 + nx_nodes;
            elements.push_back({n1, n2, n3, n4});
        }
    }
    
    // Inicializar variables
    velocity.resize(ndof_v, 0.0);
    pressure.resize(ndof_p, 0.0);
    eps_bar.resize(nelem, 0.0);
    
    // Guardar malla inicial
    coords_history.push_back(coords);
}

// ================================
// SIMULACIÓN PRINCIPAL
// ================================
void run_simulation() {
    std::cout << "==============================================================" << std::endl;
    std::cout << "SIMULACIÓN PLAIN STRAIN DE FORJA - FÍSICAMENTE CORRECTA" << std::endl;
    std::cout << "==============================================================" << std::endl;
    std::cout << "Material: Norton-Hoff, K=" << Kmat << " Pa·s^" << nexp 
              << ", n=" << nexp << std::endl;
    std::cout << "Elementos: " << nelem << ", Nodos: " << nnodes << std::endl;
    std::cout << "DOFs: " << ndof_total << " (" << ndof_v 
              << " velocidades + " << ndof_p << " presiones)" << std::endl;
    
    // Inicialización
    initialize_mesh();
    auto fixed_dofs = setup_boundary_conditions();
    
    // Matrices temporales para reuso
    Matrix K_temp(ndof_total, ndof_total);
    std::vector<double> F_temp(ndof_total, 0.0);
    
    // Bucle temporal
    for(int step = 0; step < nsteps; step++) {
        std::cout << "\n--- PASO " << step+1 << "/" << nsteps 
                  << " (dt=" << dt << " s) ---" << std::endl;
        
        auto result = solve_step_with_relaxation(velocity, pressure, 
                                                fixed_dofs, K_temp, F_temp);
        
        velocity = result.velocity;
        pressure = result.pressure;
        
        if(!result.converged) {
            std::cout << "  ¡ADVERTENCIA: Paso no convergido completamente!" << std::endl;
        }
        
        // Actualizar deformación acumulada
        for(int e_idx = 0; e_idx < nelem; e_idx++) {
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
            auto strain_result = calculate_strain_rate_plain_strain(
                jac_result.dNdX, vel_elem, jac_result.N
            );
            
            eps_bar[e_idx] += dt * strain_result.eps_dot_eq;
        }
        
        // Updated Lagrangian: actualizar coordenadas
        for(int i = 0; i < nnodes; i++) {
            coords[i].x += dt * velocity[2*i];
            coords[i].y += dt * velocity[2*i + 1];
        }
        coords_history.push_back(coords);
        
        // Estadísticas
        double vx_min = 1e100, vx_max = -1e100;
        double vy_min = 1e100, vy_max = -1e100;
        double p_min = 1e100, p_max = -1e100;
        
        for(int i = 0; i < nnodes; i++) {
            vx_min = min(vx_min, velocity[2*i]);
            vx_max = max(vx_max, velocity[2*i]);
            vy_min = min(vy_min, velocity[2*i + 1]);
            vy_max = max(vy_max, velocity[2*i + 1]);
        }
        
        for(int i = 0; i < ndof_p; i++) {
            p_min = min(p_min, pressure[i]);
            p_max = max(p_max, pressure[i]);
        }
        
        std::cout << "  v_x: [" << vx_min << ", " << vx_max << "] m/s" << std::endl;
        std::cout << "  v_y: [" << vy_min << ", " << vy_max << "] m/s" << std::endl;
        std::cout << "  Presión: [" << p_min << ", " << p_max << "] Pa" << std::endl;
    }
    
    std::cout << "\n==============================================================" << std::endl;
    std::cout << "SIMULACIÓN COMPLETADA" << std::endl;
    std::cout << "==============================================================" << std::endl;
}

// ================================
// FUNCIÓN PRINCIPAL
// ================================
int main() {
    run_simulation();
    return 0;
}

// ================================
// FUNCIÓN DE RESOLUCIÓN LINEAL (EJEMPLO)
// Necesitas implementar según tu solver
// ================================
std::vector<double> solve_linear_system(const Matrix& K, const std::vector<double>& F) {
    // Esto es un placeholder - debes implementar con tu solver preferido
    
    std::vector<double> solution(F.size(), 0.0);
    
    // Ejemplo simple: resolver con eliminación gaussiana
    // En la práctica, usarías Eigen, PETSc, o tu propio solver
    int n = K.getRows();
    
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
        // ================================
// stokes_pspg_plainstrain.cpp
// Migración completa de Python a C++ con Plain Strain
// Usando tu clase Matrix
// ================================

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <memory>

// ================================
// PARÁMETROS FÍSICAMENTE CORRECTOS
// ================================
int nelr = 15, nelz = 15;
double Lx = 0.0127, Ly = 0.03;  // Dimensiones en X e Y (plain strain)

// MATERIAL: Acero caliente tipo Norton-Hoff
double Kmat = 50.0e6;          // Coeficiente de consistencia [Pa·s^n]
double nexp = 0.2;             // Exponente strain-rate sensitivity
double mu_ref = Kmat / 3.0;    // Viscosidad de referencia

// Parámetros numéricos
double beta_stab = 0.1;
double relax_factor = 0.6;
double tol = 1e-8;
int max_iter = 50;
double dt = 0.001;
int nsteps = 4;
double v_y_top = -1.0;  // Compresión en dirección Y

// ================================
// ESTRUCTURAS DE DATOS
// ================================
struct Point2D {
    double x, y;
    Point2D(double x_=0, double y_=0) : x(x_), y(y_) {}
};

// Variables globales
int nx_nodes, ny_nodes, nnodes, nelem;
int ndof_v, ndof_p, ndof_total;

std::vector<Point2D> coords;
std::vector<std::vector<int>> elements;
std::vector<double> velocity;    // [vx0, vy0, vx1, vy1, ...]
std::vector<double> pressure;    // Presión por elemento
std::vector<double> eps_bar;     // Deformación plástica acumulada
std::vector<std::vector<Point2D>> coords_history;

// ================================
// FUNCIONES AUXILIARES
// ================================
inline double max(double a, double b) { return a > b ? a : b; }
inline double min(double a, double b) { return a < b ? a : b; }

// ================================
// FUNCIONES DE FORMA Q1 (PLAIN STRAIN)
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
        pos_mat.Set(i, 0, pos[i].x);
        pos_mat.Set(i, 1, pos[i].y);
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

// ================================
// VISCOPLASTICIDAD NORTON-HOFF (PLAIN STRAIN)
// ================================
double effective_viscosity_norton(double eps_dot_eq) {
    // Evitar divisiones por cero
    double eps_dot_eq_safe = max(eps_dot_eq, 1e-12);
    return (Kmat / 3.0) * pow(eps_dot_eq_safe, nexp - 1.0);
}

struct StrainRateResult {
    double eps_dot_eq;
    std::vector<double> strain_vec;  // [exx, eyy, exy] para plain strain
    double eps_mean;
};

StrainRateResult calculate_strain_rate_plain_strain(const Matrix& dNdX,
                                                   const std::vector<double>& vel_elem,
                                                   const std::vector<double>& N) {
    StrainRateResult result;
    result.strain_vec.resize(3);
    
    // Inicializar gradientes
    double dvx_dx = 0.0, dvx_dy = 0.0;
    double dvy_dx = 0.0, dvy_dy = 0.0;
    
    for(int a = 0; a < 4; a++) {
        double vx = vel_elem[2*a];      // vx del nodo a
        double vy = vel_elem[2*a + 1];  // vy del nodo a
        
        dvx_dx += dNdX.getVal(0, a) * vx;
        dvx_dy += dNdX.getVal(1, a) * vx;
        dvy_dx += dNdX.getVal(0, a) * vy;
        dvy_dy += dNdX.getVal(1, a) * vy;
    }
    
    // Tensor de deformación para plain strain (2D)
    double exx = dvx_dx;
    double eyy = dvy_dy;
    double exy = 0.5 * (dvx_dy + dvy_dx);
    
    // En plain strain, ezz = 0 (deformación en z es cero)
    // Pero contribuye a la parte volumétrica
    double ezz = 0.0;
    
    // Parte esférica (volumétrica)
    result.eps_mean = (exx + eyy + ezz) / 3.0;
    
    // Parte desviadora
    double exx_dev = exx - result.eps_mean;
    double eyy_dev = eyy - result.eps_mean;
    double ezz_dev = ezz - result.eps_mean;
    double exy_dev = exy;  // La componente cortante ya es desviadora
    
    // Invariante J2 para plain strain
    // J2 = 1/2 * (exx_dev^2 + eyy_dev^2 + ezz_dev^2) + exy^2
    double J2 = 0.5 * (exx_dev*exx_dev + eyy_dev*eyy_dev + ezz_dev*ezz_dev) 
                + exy_dev*exy_dev;
    
    // Strain-rate equivalente: ε̇_eq = √(3J2)
    result.eps_dot_eq = sqrt(3.0 * J2);
    
    result.strain_vec[0] = exx;
    result.strain_vec[1] = eyy;
    result.strain_vec[2] = exy;
    
    return result;
}

// ================================
// ENSAMBLAJE PSPG (PLAIN STRAIN)
// ================================
void assemble_mixed_system_correct(const std::vector<double>& vel_vec,
                                  const std::vector<double>& press_vec,
                                  Matrix& K_glob,
                                  std::vector<double>& F_glob,
                                  double& max_mu_eff) {
    
    // Inicializar matriz global (sparse, pero usamos densa por simplicidad inicial)
    K_glob.SetZero();
    std::fill(F_glob.begin(), F_glob.end(), 0.0);
    max_mu_eff = 0.0;
    
    // Puntos de Gauss
    const double a = 1.0 / sqrt(3.0);
    std::vector<std::pair<double, double>> gp_full = {
        {-a, -a}, {a, -a}, {a, a}, {-a, a}
    };
    std::vector<double> w_full = {1.0, 1.0, 1.0, 1.0};
    
    for(int e_idx = 0; e_idx < nelem; e_idx++) {
        const auto& conn = elements[e_idx];
        
        // Posiciones de los nodos del elemento
        std::vector<Point2D> pos(4);
        for(int i = 0; i < 4; i++) {
            pos[i] = coords[conn[i]];
        }
        
        // DOFs de velocidades del elemento
        std::vector<int> dofs_v;
        for(int n : conn) {
            dofs_v.push_back(2*n);      // vx del nodo n
            dofs_v.push_back(2*n + 1);  // vy del nodo n
        }
        
        int dof_p = ndof_v + e_idx;
        
        // Extraer velocidades del elemento
        std::vector<double> vel_elem(8);
        for(int i = 0; i < 8; i++) {
            vel_elem[i] = vel_vec[dofs_v[i]];
        }
        
        double press_elem = press_vec[e_idx];
        
        // Matrices elementales
        Matrix K_e(8, 8);
        Matrix G_e(8, 1);
        Matrix Kstab_e(8, 8);
        K_e.SetZero();
        G_e.SetZero();
        Kstab_e.SetZero();
        
        // Tamaño característico del elemento
        double h_x = fabs(pos[1].x - pos[0].x) + fabs(pos[2].x - pos[3].x);
        double h_y = fabs(pos[3].y - pos[0].y) + fabs(pos[2].y - pos[1].y);
        double h_e = 0.5 * sqrt(h_x*h_x + h_y*h_y);
        
        double area_elem = 0.0;
        double mu_eff_max_elem = 0.0;
        
        // Integración en puntos de Gauss
        for(int gp = 0; gp < 4; gp++) {
            double xi = gp_full[gp].first;
            double eta = gp_full[gp].second;
            double w = w_full[gp];
            
            auto jac_result = jacobian_and_gradients(pos, xi, eta);
            
            // dOmega = detJ * w (sin factor 2πr de axisimetría)
            double dOmega = jac_result.detJ * w;
            area_elem += dOmega;
            
            // Calcular strain-rate y viscosidad
            auto strain_result = calculate_strain_rate_plain_strain(
                jac_result.dNdX, vel_elem, jac_result.N
            );
            
            double mu_eff = effective_viscosity_norton(strain_result.eps_dot_eq);
            mu_eff_max_elem = max(mu_eff_max_elem, mu_eff);
            
            // 1. Matriz B (3x8) para plain strain
            Matrix B(3, 8);
            B.SetZero();
            
            // ε_xx = ∂vx/∂x
            for(int a = 0; a < 4; a++) {
                B.Set(0, 2*a, jac_result.dNdX.getVal(0, a));
            }
            
            // ε_yy = ∂vy/∂y
            for(int a = 0; a < 4; a++) {
                B.Set(1, 2*a + 1, jac_result.dNdX.getVal(1, a));
            }
            
            // γ_xy = ∂vx/∂y + ∂vy/∂x
            for(int a = 0; a < 4; a++) {
                B.Set(2, 2*a, jac_result.dNdX.getVal(1, a));      // ∂vx/∂y
                B.Set(2, 2*a + 1, jac_result.dNdX.getVal(0, a));  // ∂vy/∂x
            }
            
            // 2. Tensor constitutivo D para plain strain (3x3)
            // Para material incompresible con penalización
            Matrix D(3, 3);
            D.SetZero();
            
            // Forma desviadora para viscosidad
            D.Set(0, 0, 2.0 * mu_eff);
            D.Set(1, 1, 2.0 * mu_eff);
            D.Set(2, 2, mu_eff);  // Componente cortante
            
            // 3. Contribución viscosa: B^T * D * B
            Matrix Bt = B.getTranspose();
            Matrix DxB = Matrix(3, 8);
            
            // Multiplicación D * B
            for(int i = 0; i < 3; i++) {
                for(int j = 0; j < 8; j++) {
                    double sum = 0.0;
                    for(int k = 0; k < 3; k++) {
                        sum += D.getVal(i, k) * B.getVal(k, j);
                    }
                    DxB.Set(i, j, sum);
                }
            }
            
            // B^T * (D * B)
            for(int i = 0; i < 8; i++) {
                for(int j = 0; j < 8; j++) {
                    double sum = 0.0;
                    for(int k = 0; k < 3; k++) {
                        sum += Bt.getVal(i, k) * DxB.getVal(k, j);
                    }
                    K_e.Set(i, j, K_e.getVal(i, j) + sum * dOmega);
                }
            }
            
            // 4. Matriz de acoplamiento presión-velocidad (B_vol)
            // Para plain strain: div(v) = ∂vx/∂x + ∂vy/∂y
            Matrix B_vol(1, 8);
            B_vol.SetZero();
            
            for(int a = 0; a < 4; a++) {
                B_vol.Set(0, 2*a, jac_result.dNdX.getVal(0, a));      // ∂vx/∂x
                B_vol.Set(0, 2*a + 1, jac_result.dNdX.getVal(1, a));  // ∂vy/∂y
            }
            
            // G = B_vol^T * N_p (N_p = 1 para P0)
            Matrix B_vol_t = B_vol.getTranspose();
            for(int i = 0; i < 8; i++) {
                G_e.Set(i, 0, G_e.getVal(i, 0) + B_vol_t.getVal(i, 0) * dOmega);
            }
            
            // 5. Estabilización PSPG
            double tau = beta_stab * (h_e*h_e) / (2.0 * mu_eff + 1e-12);
            
            // Kstab = tau * B_vol^T * B_vol
            for(int i = 0; i < 8; i++) {
                for(int j = 0; j < 8; j++) {
                    double sum = 0.0;
                    for(int k = 0; k < 1; k++) {
                        sum += B_vol_t.getVal(i, k) * B_vol.getVal(k, j);
                    }
                    Kstab_e.Set(i, j, Kstab_e.getVal(i, j) + tau * sum * dOmega);
                }
            }
        }
        
        // Penalización volumétrica
        double Kbulk_eff = 1000.0 * mu_eff_max_elem;
        double M_p_elem = area_elem / Kbulk_eff;
        
        max_mu_eff = max(max_mu_eff, mu_eff_max_elem);
        
        // ENSAMBLAJE GLOBAL
        // 1. Parte de velocidades (K_e + Kstab_e)
        for(int i_local = 0; i_local < 8; i_local++) {
            int i_global = dofs_v[i_local];
            
            for(int j_local = 0; j_local < 8; j_local++) {
                int j_global = dofs_v[j_local];
                
                double val = K_glob.getVal(i_global, j_global)
                           + K_e.getVal(i_local, j_local)
                           + Kstab_e.getVal(i_local, j_local);
                
                K_glob.Set(i_global, j_global, val);
            }
            
            // Acoplamiento velocidad-presión (G)
            K_glob.Set(i_global, dof_p, 
                      K_glob.getVal(i_global, dof_p) + G_e.getVal(i_local, 0));
            
            // Transpuesta de G
            K_glob.Set(dof_p, i_global,
                      K_glob.getVal(dof_p, i_global) + G_e.getVal(i_local, 0));
        }
        
        // Penalización de presión
        K_glob.Set(dof_p, dof_p,
                  K_glob.getVal(dof_p, dof_p) - M_p_elem);
    }
}

// ================================
// CONDICIONES DE CONTORNO (PLAIN STRAIN)
// ================================
std::unordered_map<int, double> setup_boundary_conditions() {
    std::unordered_map<int, double> fixed_dofs;
    
    // Base inferior fija (y = 0)
    for(int i = 0; i < nx_nodes; i++) {
        int node_id = i;  // Primera fila
        fixed_dofs[2*node_id] = 0.0;      // vx = 0
        fixed_dofs[2*node_id + 1] = 0.0;  // vy = 0
    }
    
    // Tapa superior: velocidad impuesta (y = Ly)
    for(int i = 0; i < nx_nodes; i++) {
        int node_id = (ny_nodes - 1) * nx_nodes + i;
        fixed_dofs[2*node_id] = 0.0;           // vx = 0
        fixed_dofs[2*node_id + 1] = v_y_top;   // vy = velocidad impuesta
    }
    
    // Lado izquierdo: simetría o condición de deslizamiento
    for(int j = 0; j < ny_nodes; j++) {
        int node_id = j * nx_nodes;
        fixed_dofs[2*node_id] = 0.0;  // vx = 0 (simetría)
        // vy libre
    }
    
    // Lado derecho: libre (o condición según necesites)
    // Para compresión pura, podrías dejar libre o imponer vx = 0
    
    return fixed_dofs;
}

void apply_boundary_conditions(Matrix& K_glob,
                              std::vector<double>& F_glob,
                              const std::unordered_map<int, double>& fixed_dofs) {
    
    for(const auto& [dof, value] : fixed_dofs) {
        // Restar contribución de la columna
        for(int i = 0; i < ndof_total; i++) {
            F_glob[i] -= K_glob.getVal(i, dof) * value;
        }
        
        // Poner fila y columna en cero
        for(int i = 0; i < ndof_total; i++) {
            K_glob.Set(dof, i, 0.0);
            K_glob.Set(i, dof, 0.0);
        }
        
        // Poner 1 en la diagonal
        K_glob.Set(dof, dof, 1.0);
        F_glob[dof] = value;
    }
}

// ================================
// SOLUCIÓN CON RELAJACIÓN (PICARD)
// ================================
struct SolveResult {
    std::vector<double> velocity;
    std::vector<double> pressure;
    bool converged;
    int iterations;
};

SolveResult solve_step_with_relaxation(std::vector<double>& vel_guess,
                                      std::vector<double>& press_guess,
                                      const std::unordered_map<int, double>& fixed_dofs,
                                      Matrix& K_temp,  // Matriz temporal para reuso
                                      std::vector<double>& F_temp) {  // Vector temporal
    
    std::vector<double> vel_prev = vel_guess;
    std::vector<double> press_prev = press_guess;
    
    SolveResult result;
    
    for(int iter = 0; iter < max_iter; iter++) {
        // Ensamblar sistema
        double max_mu;
        assemble_mixed_system_correct(vel_guess, press_guess, K_temp, F_temp, max_mu);
        
        // Aplicar condiciones de contorno
        apply_boundary_conditions(K_temp, F_temp, fixed_dofs);
        
        try {
            // Resolver sistema (usando tu solver)
            // Asumo que tienes una función solve_linear_system
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
                vel_new_relax[i] = (1.0 - relax_factor) * vel_guess[i] 
                                 + relax_factor * vel_new[i];
            }
            
            for(int i = 0; i < ndof_p; i++) {
                press_new_relax[i] = (1.0 - relax_factor) * press_guess[i] 
                                   + relax_factor * press_new[i];
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
            
            std::cout << "    Iter " << iter+1 
                      << ": Δv=" << dv_rel 
                      << ", Δp=" << dp_rel 
                      << ", μ_max=" << max_mu << " Pa·s" << std::endl;
            
            if(dv_rel < tol && dp_rel < tol) {
                result.velocity = vel_new_relax;
                result.pressure = press_new_relax;
                result.converged = true;
                result.iterations = iter + 1;
                return result;
            }
            
            // Actualizar para siguiente iteración
            vel_guess = vel_new_relax;
            press_guess = press_new_relax;
            vel_prev = vel_guess;
            press_prev = press_guess;
            
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

// ================================
// INICIALIZACIÓN DE MALLADO
// ================================
void initialize_mesh() {
    nx_nodes = nelr + 1;
    ny_nodes = nelz + 1;
    nnodes = nx_nodes * ny_nodes;
    nelem = nelr * nelz;
    
    // DOFs
    ndof_v = 2 * nnodes;
    ndof_p = nelem;
    ndof_total = ndof_v + ndof_p;
    
    // Coordenadas
    coords.reserve(nnodes);
    for(int j = 0; j < ny_nodes; j++) {
        for(int i = 0; i < nx_nodes; i++) {
            double x = i * Lx / nelr;
            double y = j * Ly / nelz;
            coords.push_back(Point2D(x, y));
        }
    }
    
    // Elementos Q4
    elements.reserve(nelem);
    for(int j = 0; j < nelz; j++) {
        for(int i = 0; i < nelr; i++) {
            int n1 = j * nx_nodes + i;
            int n2 = n1 + 1;
            int n3 = n2 + nx_nodes;
            int n4 = n1 + nx_nodes;
            elements.push_back({n1, n2, n3, n4});
        }
    }
    
    // Inicializar variables
    velocity.resize(ndof_v, 0.0);
    pressure.resize(ndof_p, 0.0);
    eps_bar.resize(nelem, 0.0);
    
    // Guardar malla inicial
    coords_history.push_back(coords);
}

// ================================
// SIMULACIÓN PRINCIPAL
// ================================
void run_simulation() {
    std::cout << "==============================================================" << std::endl;
    std::cout << "SIMULACIÓN PLAIN STRAIN DE FORJA - FÍSICAMENTE CORRECTA" << std::endl;
    std::cout << "==============================================================" << std::endl;
    std::cout << "Material: Norton-Hoff, K=" << Kmat << " Pa·s^" << nexp 
              << ", n=" << nexp << std::endl;
    std::cout << "Elementos: " << nelem << ", Nodos: " << nnodes << std::endl;
    std::cout << "DOFs: " << ndof_total << " (" << ndof_v 
              << " velocidades + " << ndof_p << " presiones)" << std::endl;
    
    // Inicialización
    initialize_mesh();
    auto fixed_dofs = setup_boundary_conditions();
    
    // Matrices temporales para reuso
    Matrix K_temp(ndof_total, ndof_total);
    std::vector<double> F_temp(ndof_total, 0.0);
    
    // Bucle temporal
    for(int step = 0; step < nsteps; step++) {
        std::cout << "\n--- PASO " << step+1 << "/" << nsteps 
                  << " (dt=" << dt << " s) ---" << std::endl;
        
        auto result = solve_step_with_relaxation(velocity, pressure, 
                                                fixed_dofs, K_temp, F_temp);
        
        velocity = result.velocity;
        pressure = result.pressure;
        
        if(!result.converged) {
            std::cout << "  ¡ADVERTENCIA: Paso no convergido completamente!" << std::endl;
        }
        
        // Actualizar deformación acumulada
        for(int e_idx = 0; e_idx < nelem; e_idx++) {
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
            auto strain_result = calculate_strain_rate_plain_strain(
                jac_result.dNdX, vel_elem, jac_result.N
            );
            
            eps_bar[e_idx] += dt * strain_result.eps_dot_eq;
        }
        
        // Updated Lagrangian: actualizar coordenadas
        for(int i = 0; i < nnodes; i++) {
            coords[i].x += dt * velocity[2*i];
            coords[i].y += dt * velocity[2*i + 1];
        }
        coords_history.push_back(coords);
        
        // Estadísticas
        double vx_min = 1e100, vx_max = -1e100;
        double vy_min = 1e100, vy_max = -1e100;
        double p_min = 1e100, p_max = -1e100;
        
        for(int i = 0; i < nnodes; i++) {
            vx_min = min(vx_min, velocity[2*i]);
            vx_max = max(vx_max, velocity[2*i]);
            vy_min = min(vy_min, velocity[2*i + 1]);
            vy_max = max(vy_max, velocity[2*i + 1]);
        }
        
        for(int i = 0; i < ndof_p; i++) {
            p_min = min(p_min, pressure[i]);
            p_max = max(p_max, pressure[i]);
        }
        
        std::cout << "  v_x: [" << vx_min << ", " << vx_max << "] m/s" << std::endl;
        std::cout << "  v_y: [" << vy_min << ", " << vy_max << "] m/s" << std::endl;
        std::cout << "  Presión: [" << p_min << ", " << p_max << "] Pa" << std::endl;
    }
    
    std::cout << "\n==============================================================" << std::endl;
    std::cout << "SIMULACIÓN COMPLETADA" << std::endl;
    std::cout << "==============================================================" << std::endl;
}

// ================================
// FUNCIÓN PRINCIPAL
// ================================
int main() {
    run_simulation();
    return 0;
}

// ================================
// FUNCIÓN DE RESOLUCIÓN LINEAL (EJEMPLO)
// Necesitas implementar según tu solver
// ================================
std::vector<double> solve_linear_system(const Matrix& K, const std::vector<double>& F) {
    // Esto es un placeholder - debes implementar con tu solver preferido
    
    std::vector<double> solution(F.size(), 0.0);
    
    // Ejemplo simple: resolver con eliminación gaussiana
    // En la práctica, usarías Eigen, PETSc, o tu propio solver
    int n = K.getRows();
    
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