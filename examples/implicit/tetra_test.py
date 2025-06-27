import numpy as np

# === Material properties ===
E = 200e9       # Young's modulus [Pa]
nu = 0.3        # Poisson's ratio
dim = 3
n_nodes = 4

# === Nodal coordinates of tetrahedron (4 nodes) ===
coords = np.array([
    [0.0, 0.0, 0.0],  # Node 0
    [1.0, 0.0, 0.0],  # Node 1
    [0.0, 1.0, 0.0],  # Node 2
    [0.0, 0.0, 1.0]   # Node 3
])

# === External forces (global vector) ===
f_ext = np.zeros(dim * n_nodes)
f_ext[1] = -1000.0  # Apply -1000 N in Y direction at node 0

# === Construct constitutive matrix for isotropic linear elastic material ===
def constitutive_matrix(E, nu):
    lmbda = (E * nu) / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))
    C = np.zeros((6, 6))
    C[:3, :3] = lmbda
    np.fill_diagonal(C[:3, :3], lmbda + 2 * mu)
    C[3:, 3:] = np.eye(3) * mu
    print ("C matrix ",C)
    return C

# === Compute B matrix and element stiffness ===
def compute_B_and_Ke(coords, C):
    # Build the volume and shape function derivatives using Jacobian
    x0, x1, x2, x3 = coords
    J = np.column_stack((x1 - x0, x2 - x0, x3 - x0))
    detJ = np.linalg.det(J)
    volume = abs(detJ) / 6.0

    invJ = np.linalg.inv(J)

    # Shape function gradients in reference tetrahedron (constant)
    dN = np.array([
        [-1, -1, -1],  # N0
        [1,  0,  0],   # N1
        [0,  1,  0],   # N2
        [0,  0,  1]    # N3
    ])
    grads = dN @ invJ  # shape function gradients in physical coords

    # Build B matrix (6 x 12)
    B = np.zeros((6, 12))
    for i in range(4):
        base = i * 3
        dN_dx, dN_dy, dN_dz = grads[i]
        # Normal strains
        B[0, base + 0] = dN_dx
        B[1, base + 1] = dN_dy
        B[2, base + 2] = dN_dz
        # Shear strains
        B[3, base + 0] = dN_dy
        B[3, base + 1] = dN_dx
        B[4, base + 1] = dN_dz
        B[4, base + 2] = dN_dy
        B[5, base + 2] = dN_dx
        B[5, base + 0] = dN_dz

    # Element stiffness matrix: Ke = Bᵀ * C * B * volume
    Ke = B.T @ C @ B * volume
    return Ke, B, volume

# === Solve KU = F with simple Dirichlet BCs ===
def apply_dirichlet(K, f, fixed_dofs):
    for dof in fixed_dofs:
        K[dof, :] = 0
        K[:, dof] = 0
        K[dof, dof] = 1
        f[dof] = 0
    return K, f

# === Main ===
C = constitutive_matrix(E, nu)
Ke, B, vol = compute_B_and_Ke(coords, C)

# Apply fixed boundary conditions: fix all DOFs of node 0
fixed_dofs = [0,1,2,  4,5 ,8]  # Node 0: x, y, z
K_mod, f_mod = apply_dirichlet(Ke.copy(), f_ext.copy(), fixed_dofs)

#f_mod[ 9] = -1.0e3
#f_mod[10] = -1.0e3
f_mod[11] = -1.0e3
# Solve system
u = np.linalg.solve(K_mod, f_mod)

# === Output results ===
print("Displacements (U):")
print ("K_mod",K_mod)
for i in range(n_nodes):
    ux, uy, uz = u[3*i:3*i+3]
    print(f"Node {i}: Ux = {ux:.4e} m, Uy = {uy:.4e} m, Uz = {uz:.4e} m")

