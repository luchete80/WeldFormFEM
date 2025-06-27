import numpy as np

# === Material properties ===
E = 200e9       # Young's modulus [Pa]
nu = 0.3        # Poisson's ratio
dim = 3
n_nodes = 4

# === Nodal coordinates of tetrahedron (4 nodes) ===
coords_0 = np.array([
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0]
])

def constitutive_matrix(E, nu):
    lmbda = (E * nu) / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    C = np.zeros((6, 6))
    C[:3, :3] = lmbda
    np.fill_diagonal(C[:3, :3], lmbda + 2 * mu)
    C[3:, 3:] = np.eye(3) * mu
    return C

def compute_B(coords):
    x0, x1, x2, x3 = coords
    J = np.column_stack((x1 - x0, x2 - x0, x3 - x0))
    detJ = np.linalg.det(J)
    volume = abs(detJ) / 6.0
    invJ = np.linalg.inv(J)

    dN = np.array([
        [-1, -1, -1],
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]
    ])
    grads = dN @ invJ

    B = np.zeros((6, 12))
    for i in range(4):
        base = i * 3
        dN_dx, dN_dy, dN_dz = grads[i]
        B[0, base + 0] = dN_dx
        B[1, base + 1] = dN_dy
        B[2, base + 2] = dN_dz
        B[3, base + 0] = dN_dy
        B[3, base + 1] = dN_dx
        B[4, base + 1] = dN_dz
        B[4, base + 2] = dN_dy
        B[5, base + 2] = dN_dx
        B[5, base + 0] = dN_dz

    return B, grads, volume

def F_from_u(grads, u):
    F = np.eye(3)
    for a in range(4):
        u_vec = u[a*3:(a+1)*3]
        F += np.outer(u_vec, grads[a])
    return F

def almansi_strain(F):
    b = F @ F.T
    b_inv = np.linalg.inv(b)
    e = 0.5 * (np.eye(3) - b_inv)
    return e

def strain_to_voigt(e):
    return np.array([
        e[0, 0], e[1, 1], e[2, 2],
        2 * e[0, 1], 2 * e[1, 2], 2 * e[0, 2]
    ])

def apply_dirichlet(K, f, fixed_dofs):
    for dof in fixed_dofs:
        K[dof, :] = 0
        K[:, dof] = 0
        K[dof, dof] = 1
        f[dof] = 0
    return K, f

def compute_Kgeo(grads, sigma, volume):
    Kgeo = np.zeros((12, 12))  # 4 nodes × 3 dof each

    for a in range(4):
        grad_a = grads[a].reshape(3, 1)  # 3x1
        for b in range(4):
            grad_b = grads[b].reshape(3, 1)  # 3x1

            # Compute 3x3 block: grad_a^T * sigma * grad_b
            block = grad_a.T @ sigma @ grad_b  # scalar

            # Outer product to expand to 3x3 block
            outer = np.outer(grad_a.flatten(), grad_b.flatten())  # 3x3

            # Multiply by stress to get the contribution
            kab = (sigma @ grad_b) @ grad_a.T  # 3x3

            # Place in the correct block in Kgeo
            Kgeo[3*a:3*a+3, 3*b:3*b+3] += volume * kab

    return Kgeo

##ALTERNATIVE
def compute_Kgeo_voigt(B, stress, volume):
    S = np.zeros((6, 6))
    S[:3, :3] = np.diag(stress[:3])
    S[3:, 3:] = np.diag(stress[3:])
    Kgeo = B.T @ S @ B * volume
    return Kgeo

def voigt_to_tensor(voigt):
    """
    Convert 6x1 Voigt vector to 3x3 symmetric tensor.
    Voigt notation:
    [σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_zx]
    """
    return np.array([
        [voigt[0], voigt[3], voigt[5]],
        [voigt[3], voigt[1], voigt[4]],
        [voigt[5], voigt[4], voigt[2]]
    ])   



# === External force (final target) ===
f_ext_total = np.zeros(dim * n_nodes)
f_ext_total[11] = -1000.0  # Apply force at Node 3, z-direction
# === Fixed DOFs (Dirichlet BCs) ===
fixed_dofs = [0, 1, 2, 4, 5, 8]


# === Initialize ===
C = constitutive_matrix(E, nu)
u = np.zeros(12)  # Total displacement

n_steps = 10

for step in range(1, n_steps + 1):
    print(f"\n--- Load Step {step}/{n_steps} ---")
    f_ext = f_ext_total * (step / n_steps)

    for iteration in range(10):  # Newton-Raphson loop
        coords_current = coords_0 + u.reshape((4, 3))
        B, grads, volume = compute_B(coords_current)

        F = F_from_u(grads, u)
        e_tensor = almansi_strain(F)
        e_voigt = strain_to_voigt(e_tensor)
        stress = C @ e_voigt

        f_int = B.T @ stress * volume
        Kmat = B.T @ C @ B * volume
        sigma_tensor = voigt_to_tensor(stress)
        Kgeo = compute_Kgeo(grads, sigma_tensor, volume)
        K_total = Kmat + Kgeo
        
        #print ("Kgeo", Kgeo)

        # Residual
        r = f_ext - f_int

        # Apply BCs
        K_mod, r_mod = apply_dirichlet(K_total.copy(), r.copy(), fixed_dofs)

        # Solve for du
        du = np.linalg.solve(K_mod, r_mod)
        u += du

        norm_du = np.linalg.norm(du)
        print(f"  Iter {iteration+1}: ||du|| = {norm_du:.3e}")
        if norm_du < 1e-6:
            break

# === Final Output ===
print("\nFinal Displacements (Updated Lagrangian):")
for i in range(n_nodes):
    ux, uy, uz = u[3*i:3*i+3]
    print(f"Node {i}: Ux = {ux:.4e}, Uy = {uy:.4e}, Uz = {uz:.4e}")