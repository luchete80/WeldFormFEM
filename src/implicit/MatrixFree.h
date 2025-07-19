#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>  // Para GMRES/BiCGSTAB

typedef Eigen::VectorXd Vector;
typedef Eigen::SparseMatrix<double> SparseMat;

// ===================== OPERADOR MATRIX-FREE =====================
class MatrixFreeOperator {
public:
    MatrixFreeOperator(YourSolverClass* solver) : m_solver(solver) {}

    // Función que Eigen usará para calcular K * v (matrix-free)
    void multiply(const Eigen::VectorXd& v, Eigen::VectorXd& Kv) const {
        // 1. Convertir Eigen::VectorXd a tu formato de matriz
        Matrix v_internal(m_node_count * m_dim, 1);
        for (int i = 0; i < v.size(); ++i) {
            v_internal.Set(i, 0, v[i]);
        }

        // 2. Calcular K * v usando tu lógica existente (paralelizable)
        Matrix Kv_internal = ComputeKvAction(v_internal);  // <-- Tu función matrix-free

        // 3. Convertir de vuelta a Eigen::VectorXd
        Kv.resize(Kv_internal.Rows());
        for (int i = 0; i < Kv_internal.Rows(); ++i) {
            Kv[i] = Kv_internal.Get(i, 0);
        }
    }

private:
    YourSolverClass* m_solver;
};

// ===================== PRECONDICIONADOR (EJ: MASA DIAGONAL) =====================
class MassPreconditioner {
public:
    MassPreconditioner(YourSolverClass* solver) : m_solver(solver) {
        // Precomputar M⁻¹ (diagonal)
        m_inv_mass.resize(m_solver->m_node_count * m_solver->m_dim);
        for (int n = 0; n < m_solver->m_node_count; ++n) {
            for (int d = 0; d < m_solver->m_dim; ++d) {
                int idx = n * m_solver->m_dim + d;
                m_inv_mass[idx] = 1.0 / (m_solver->m_mdiag[n] + 1e-12);  // Evitar división por cero
            }
        }
    }

    // Aplicar el precondicionador: z = M⁻¹ * r
    void solve(const Vector& r, Vector& z) const {
        z = r.cwiseProduct(m_inv_mass);  // Multiplicación elemento a elemento
    }

private:
    YourSolverClass* m_solver;
    Vector m_inv_mass;
};

// ===================== SOLVER IMPLÍCITO =====================
void YourSolverClass::SolveImplicitStep() {
    // 1. Configurar el solver de Eigen (GMRES o BiCGSTAB)
    Eigen::GMRES<MatrixFreeOperator, MassPreconditioner> gmres;
    gmres.setTolerance(1e-6);
    gmres.setMaxIterations(1000);

    // 2. Inicializar operador matrix-free y precondicionador
    MatrixFreeOperator op(this);
    MassPreconditioner prec(this);

    // 3. Convertir el residual a Eigen::VectorXd
    Vector r_eigen(m_node_count * m_dim);
    Matrix r_global = ComputeGlobalResidual();  // <-- Tu función existente
    for (int i = 0; i < r_global.Rows(); ++i) {
        r_eigen[i] = r_global.Get(i, 0);
    }

    // 4. Resolver el sistema: K Δu = -r
    Vector delta_u_eigen(m_node_count * m_dim);
    gmres.compute(op).solveWithGuess(-r_eigen, delta_u_eigen);

    // 5. Actualizar desplazamientos (convertir delta_u_eigen a tu formato)
    for (int i = 0; i < m_node_count * m_dim; ++i) {
        u[i] += delta_u_eigen[i];
    }
}

Matrix YourSolverClass::ComputeKvAction(const Matrix& v) {
    Matrix Kv(m_node_count * m_dim, 1);
    Kv.SetZero();

    par_loop(e, m_elem_count) {
        // 1. Obtener desplazamientos locales del elemento
        Matrix v_e = gatherElementDisplacements(e, v);

        // 2. Calcular K_e * v_e (sin ensamblar K_e explícitamente)
        Matrix B = getElemBMatrix(e);
        Matrix D = getConstitutiveMatrix(e);
        Matrix K_e_v = MatMul(B.Transpose(), MatMul(D, B)) * vol[e] * v_e;

        // 3. Ensamblar en el vector global
        scatterAdd(Kv, K_e_v, e);
    }

    return Kv;
}