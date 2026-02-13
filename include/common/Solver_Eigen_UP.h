// Solver_Eigen_UP.h
#ifndef _SOLVER_EIGEN_UP_H
#define _SOLVER_EIGEN_UP_H

#include "Solver_Eigen.h"

namespace MetFEM {

class Solver_Eigen_UP : public Solver_Eigen {
public:
    Solver_Eigen_UP() : m_penalty_factor(1e6), m_use_penalty(true) {}
    virtual ~Solver_Eigen_UP() {}

    // Override métodos clave
    virtual void Allocate() override;
    virtual void assemblyGlobalMatrix() override;
    virtual int Solve() override;
    virtual void applyDirichletBCs() override;
    
    // Métodos específicos para UP
    void setPenaltyFactor(double factor) { m_penalty_factor = factor; }
    void usePenaltyMethod(bool use) { m_use_penalty = use; }
    
    // Getters para componentes
    Eigen::VectorXd getVelocity() { return U.head(m_dof_v); }
    Eigen::VectorXd getPressure() { return U.tail(m_dof_p); }

private:
    double m_penalty_factor;
    bool m_use_penalty;
    
    // Estructuras para almacenamiento separado
    std::vector<T> m_triplets_vv;  // Bloque velocidad-velocidad
    std::vector<T> m_triplets_vp;  // Bloque velocidad-presión
    std::vector<T> m_triplets_pp;  // Bloque presión-presión
    
    // ¡¡¡DECLARAR TODOS LOS MÉTODOS AUXILIARES!!!
    void assembleElementBlock(const std::vector<int>& row_dofs,
                              const std::vector<int>& col_dofs,
                              const Matrix& Ke,
                              std::vector<T>& triplets);
    
    std::vector<int> getVelocityDOFs(int e);
    
    void calculateElementMatrices(int e, Matrix& Ke_vv, Matrix& Ke_vp, 
                                  double& volume, double& incomp_error, 
                                  double& max_div);
    
    void assembleUPElement(int e, const Matrix& Ke_vv, 
                           const Matrix& Ke_vp, double volume);
    
    void applyPressureBCs();
};

} // namespace MetFEM

#endif
