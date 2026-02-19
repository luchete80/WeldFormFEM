// Solver_Eigen_UP.cpp
#include "Solver_Eigen_UP.h"
#include "Matrix.h"
#include "Domain_d.h"
#include <iomanip>

namespace MetFEM {

void Solver_Eigen_UP::Allocate() {
    // Inicializar DOFs
    m_dof_v = m_dom->m_node_count * m_dom->m_dim;
    m_dof_p = m_dom->m_elem_count;
    m_dof = m_dof_v + m_dof_p;
    
    std::cout << "Allocate UP: " << m_dof_v << " velocity DOFs, "
              << m_dof_p << " pressure DOFs (total: " << m_dof << ")" << std::endl;
    
    K.resize(m_dof, m_dof);
    U = Eigen::VectorXd::Zero(m_dof);
    R = Eigen::VectorXd::Zero(m_dof);
    
    // Reservar espacio para tripletas
    size_t est_nnz_vv = m_dom->m_elem_count * 64;  // 8x8 = 64
    size_t est_nnz_vp = m_dom->m_elem_count * 8;   // 8x1 = 8
    size_t est_nnz_pp = m_dom->m_elem_count;       // 1x1 por elemento
    
    m_triplets_vv.reserve(est_nnz_vv);
    m_triplets_vp.reserve(est_nnz_vp);
    m_triplets_pp.reserve(est_nnz_pp);
}

void Solver_Eigen_UP::assemblyGlobalMatrix() {
    // Limpiar tripletas
    m_triplets_vv.clear();
    m_triplets_vp.clear();
    m_triplets_pp.clear();
    
    // Inicializar RHS
    R = Eigen::VectorXd::Zero(m_dof);
    
    double incompressibility_error = 0.0;
    double max_div_v = 0.0;
    
    for (int e = 0; e < m_dom->m_elem_count; ++e) {
        // Calcular matrices elementales
        Matrix Ke_vv(8, 8);
        Matrix Ke_vp(8, 1);
        double volume = 0.0;
        
        calculateElementMatrices(e, Ke_vv, Ke_vp, volume, 
                                 incompressibility_error, max_div_v);
        
        // Ensamblar usando el método auxiliar
        assembleUPElement(e, Ke_vv, Ke_vp, volume);
    }
    
    // Combinar todas las tripletas
    std::vector<T> all_triplets;
    size_t total_size = m_triplets_vv.size() + m_triplets_vp.size() + m_triplets_pp.size();
    all_triplets.reserve(total_size);
    
    all_triplets.insert(all_triplets.end(), m_triplets_vv.begin(), m_triplets_vv.end());
    all_triplets.insert(all_triplets.end(), m_triplets_vp.begin(), m_triplets_vp.end());
    all_triplets.insert(all_triplets.end(), m_triplets_pp.begin(), m_triplets_pp.end());
    
    // Construir matriz global
    K.setFromTriplets(all_triplets.begin(), all_triplets.end());
    K.makeCompressed();
    
    std::cout << "Matriz ensamblada: " << K.nonZeros() << " elementos no cero" << std::endl;
    std::cout << "Error incompresibilidad: " << incompressibility_error 
              << ", max div(v): " << max_div_v << std::endl;
}

// IMPLEMENTACIÓN DEL MÉTODO FALTANTE
void Solver_Eigen_UP::assembleElementBlock(
    const std::vector<int>& row_dofs,
    const std::vector<int>& col_dofs,
    const Matrix& Ke,
    std::vector<T>& triplets)
{
    int nr = row_dofs.size();
    int nc = col_dofs.size();
    
    for (int i = 0; i < nr; ++i) {
        int I = row_dofs[i];
        for (int j = 0; j < nc; ++j) {
            int J = col_dofs[j];
            double val = Ke.at(i, j);
            if (std::abs(val) > 1e-15) {
                triplets.emplace_back(I, J, val);
            }
        }
    }
}

void Solver_Eigen_UP::assembleUPElement(int e, const Matrix& Ke_vv, 
                                        const Matrix& Ke_vp, double volume) {
    // Obtener DOFs globales
    std::vector<int> dofs_v = getVelocityDOFs(e);
    int dof_p = m_dof_v + e;
    
    // Ensamblar bloque velocidad-velocidad
    assembleElementBlock(dofs_v, dofs_v, Ke_vv, m_triplets_vv);
    
    // Ensamblar bloque velocidad-presión (Q)
    std::vector<int> col_p = {dof_p};
    assembleElementBlock(dofs_v, col_p, Ke_vp, m_triplets_vp);
    
    // Añadir transpuesta (Q^T) - mantener simetría
    // Nota: Esto asume que Ke_vp es columna, su transpuesta será fila
    for (int i = 0; i < 8; ++i) {
        double val = Ke_vp.at(i, 0);
        if (std::abs(val) > 1e-15) {
            m_triplets_vp.emplace_back(dof_p, dofs_v[i], val);
        }
    }
    
    // Penalización volumétrica
    if (m_use_penalty) {
        double penalty_val = volume / m_penalty_factor;
        if (std::abs(penalty_val) > 1e-15) {
            m_triplets_pp.emplace_back(dof_p, dof_p, penalty_val);
        }
    }
}

std::vector<int> Solver_Eigen_UP::getVelocityDOFs(int e) {
    std::vector<int> dofs;
    int nnodes = m_dom->m_nodxelem;
    int dim = m_dom->m_dim;
    dofs.reserve(nnodes * dim);
    
    for (int a = 0; a < nnodes; ++a) {
        int node = m_dom->getElemNode(e, a);
        for (int d = 0; d < dim; ++d) {
            dofs.push_back(node * dim + d);
        }
    }
    return dofs;
}

void Solver_Eigen_UP::calculateElementMatrices(int e, Matrix& Ke_vv, Matrix& Ke_vp, 
                                              double& volume, double& incomp_error, 
                                              double& max_div) {
                                                /*
    // IMPLEMENTACIÓN DE MARTINS (adaptada de tu assemble_martins_system)
    // Puntos de Gauss
    const double a = 1.0 / sqrt(3.0);
    std::vector<std::pair<double, double>> gp_full = {
        {-a, -a}, {a, -a}, {a, a}, {-a, a}
    };
    std::vector<double> w_full = {1.0, 1.0, 1.0, 1.0};
    
    std::vector<std::pair<double, double>> gp_reduced = {{0.0, 0.0}};
    std::vector<double> w_reduced = {4.0};
    
    // Obtener coordenadas de los nodos del elemento
    std::vector<Point2D> pos(4);
    for (int i = 0; i < 4; ++i) {
        int node = m_dom->getElemNode(e, i);
        pos[i].x = m_dom->x[node * m_dom->m_dim];
        pos[i].y = m_dom->x[node * m_dom->m_dim + 1];
    }
    
    // Extraer velocidades actuales (para diagnóstico)
    std::vector<double> vel_elem(8);
    std::vector<int> dofs_v = getVelocityDOFs(e);
    for (int i = 0; i < 8; ++i) {
        vel_elem[i] = U[dofs_v[i]];  // U tiene la solución actual
    }
    
    // Inicializar matrices
    Ke_vv.SetZero();
    Ke_vp.SetZero();
    
    // ========== INTEGRACIÓN COMPLETA PARA P (4 puntos) ==========
    for (int gp = 0; gp < 4; ++gp) {
        double xi = gp_full[gp].first;
        double eta = gp_full[gp].second;
        double w = w_full[gp];
        
        auto jac = m_dom->jacobian_and_gradients(e, xi, eta);
        
        // Coordenada radial
        double r_gp = 0.0;
        for (int i = 0; i < 4; ++i) {
            r_gp += jac.N[i] * pos[i].x;
        }
        r_gp = std::max(r_gp, 1e-12);
        
        double dOmega = 2.0 * M_PI * r_gp * jac.detJ * w;
        
        // Matriz d de Martins (ecuación 4.51)
        Matrix d_martins(4, 4);
        d_martins.SetZero();
        d_martins.Set(0, 0, 2.0/3.0);
        d_martins.Set(1, 1, 2.0/3.0);
        d_martins.Set(2, 2, 2.0/3.0);
        d_martins.Set(3, 3, 1.0/3.0);
        
        // Matriz B (4x8)
        Matrix B(4, 8);
        B.SetZero();
        
        for (int a = 0; a < 4; ++a) {
            // ε_rr
            B.Set(0, 2*a, jac.dNdX.getVal(0, a));
            // ε_zz
            B.Set(1, 2*a+1, jac.dNdX.getVal(1, a));
            
            if (r_gp < 1e-8) {
                // En el eje: ε_θθ = ∂vr/∂r
                B.Set(2, 2*a, jac.dNdX.getVal(0, a));
            } else {
                // Fuera del eje: ε_θθ = vr/r
                B.Set(2, 2*a, jac.N[a] / r_gp);
            }
            
            // ε_rz
            B.Set(3, 2*a, jac.dNdX.getVal(1, a));      // ∂vr/∂z
            B.Set(3, 2*a+1, jac.dNdX.getVal(0, a));    // ∂vz/∂r
        }
        
        // Calcular k_matrix = B^T * d * B
        Matrix Bt = B.getTranspose();
        Matrix dB(4, 8);
        
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 8; ++j) {
                double sum = 0.0;
                for (int k = 0; k < 4; ++k) {
                    sum += d_martins.getVal(i, k) * B.getVal(k, j);
                }
                dB.Set(i, j, sum);
            }
        }
        
        Matrix k_matrix(8, 8);
        k_matrix.SetZero();
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                double sum = 0.0;
                for (int k = 0; k < 4; ++k) {
                    sum += Bt.getVal(i, k) * dB.getVal(k, j);
                }
                k_matrix.Set(i, j, sum);
            }
        }
        
        // Aquí iría el cálculo de sigma_eq (Norton-Hoff)
        double sigma_eq = 1.0;  // Placeholder
        double eps_dot_eq = 1.0; // Placeholder
        
        double factor = sigma_eq / eps_dot_eq;
        
        // Ensamblar Ke_vv
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                double val = Ke_vv.getVal(i, j) + k_matrix.getVal(i, j) * factor * dOmega;
                Ke_vv.Set(i, j, val);
            }
        }
    }
    
    // ========== INTEGRACIÓN REDUCIDA PARA Q (1 punto) ==========
    auto jac_red = m_dom->jacobian_and_gradients(e, 0.0, 0.0);
    double r_red = 0.0;
    for (int i = 0; i < 4; ++i) {
        r_red += jac_red.N[i] * pos[i].x;
    }
    r_red = std::max(r_red, 1e-12);
    volume = 2.0 * M_PI * r_red * jac_red.detJ * 4.0;
    
    // Matriz B reducida
    Matrix B_red(4, 8);
    B_red.SetZero();
    
    for (int a = 0; a < 4; ++a) {
        B_red.Set(0, 2*a, jac_red.dNdX.getVal(0, a));
        B_red.Set(1, 2*a+1, jac_red.dNdX.getVal(1, a));
        
        if (r_red < 1e-8) {
            B_red.Set(2, 2*a, jac_red.dNdX.getVal(0, a));
        } else {
            B_red.Set(2, 2*a, jac_red.N[a] / r_red);
        }
        
        B_red.Set(3, 2*a, jac_red.dNdX.getVal(1, a));
        B_red.Set(3, 2*a+1, jac_red.dNdX.getVal(0, a));
    }
    
    // Vector C de Martins [1,1,1,0]^T
    std::vector<double> C_vec = {1.0, 1.0, 1.0, 0.0};
    
    // Calcular B^T * C
    for (int i = 0; i < 8; ++i) {
        double sum = 0.0;
        for (int k = 0; k < 4; ++k) {
            sum += B_red.getVal(k, i) * C_vec[k];
        }
        double val = Ke_vp.getVal(i, 0) + sum * volume;
        Ke_vp.Set(i, 0, val);
    }
    
    // Diagnóstico de divergencia
    double dvr_dr = 0.0, dvz_dz = 0.0, vr_interp = 0.0;
    for (int a = 0; a < 4; ++a) {
        dvr_dr += jac_red.dNdX.getVal(0, a) * vel_elem[2*a];
        dvz_dz += jac_red.dNdX.getVal(1, a) * vel_elem[2*a+1];
        vr_interp += jac_red.N[a] * vel_elem[2*a];
    }
    
    double div_v;
    if (r_red < 1e-8) {
        div_v = 2.0 * dvr_dr + dvz_dz;
    } else {
        div_v = dvr_dr + vr_interp/r_red + dvz_dz;
    }
    
    incomp_error += std::abs(div_v) * volume;
    max_div = std::max(max_div, std::abs(div_v));

  */
}

int Solver_Eigen_UP::Solve() {
    // Aplicar condiciones de contorno
    applyDirichletBCs();
    
    // Configurar solver
    solver.analyzePattern(K);
    solver.setPivotThreshold(1e-8);
    
    // Factorizar
    //~ if (solver.factorize(K) != Eigen::Success) {
        //~ std::cerr << "Error en factorización" << std::endl;
        //~ return -1;
    //~ }
    
    // Resolver
    U = solver.solve(R);
    
    if (solver.info() != Eigen::Success) {
        std::cerr << "Error en resolución" << std::endl;
        return -1;
    }
    
    // Actualizar velocidades en Domain_d
    for (int i = 0; i < m_dof_v; ++i) {
        m_dom->v[i] = U[i];
    }
    
    std::cout << "Solución UP completada. ||U|| = " << U.norm() << std::endl;
    return 0;
}

void Solver_Eigen_UP::applyDirichletBCs() {
    // Versión simplificada - puedes expandir según necesidades
    std::vector<std::pair<int, double>> bc_dofs;
    
    // BCs de velocidad
    for (int dim = 0; dim < m_dom->m_dim; ++dim) {
        int count = m_dom->bc_count[dim];
        int* nodes = nullptr;
        double* values = nullptr;
        
        if (dim == 0) { nodes = m_dom->bcx_nod; values = m_dom->bcx_val; }
        else if (dim == 1) { nodes = m_dom->bcy_nod; values = m_dom->bcy_val; }
        else if (dim == 2) { nodes = m_dom->bcz_nod; values = m_dom->bcz_val; }
        
        for (int i = 0; i < count; ++i) {
            int dof = nodes[i] * m_dom->m_dim + dim;
            bc_dofs.emplace_back(dof, values[i]);
        }
    }
    
    // Aplicar BCs (método simplificado)
    for (const auto& [dof, val] : bc_dofs) {
        // Limpiar fila/columna
        for (int k = K.outerIndexPtr()[dof]; k < K.outerIndexPtr()[dof+1]; ++k) {
            if (K.innerIndexPtr()[k] != dof) K.valuePtr()[k] = 0.0;
        }
        K.coeffRef(dof, dof) = 1.0;
        R[dof] = val;
    }
    
    K.makeCompressed();
}

void Solver_Eigen_UP::applyPressureBCs() {
    // Implementar si hay BCs de presión
}

} // namespace MetFEM
