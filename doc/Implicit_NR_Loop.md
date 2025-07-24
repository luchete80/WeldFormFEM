# Implicit Dynamics Solver Algorithm

```python
while Time < end_time:
    # 1. Initialize time step
    prev_v = v
    delta_v = 0
    x_initial = x
    
    # 2. Newton-Raphson iterations
    while not converged:
        # Velocity update
        v = prev_v + delta_v
        
        # Apply boundary conditions
        ImposeBCV(v)  
        
        # Displacement and position
        u = dt * v
        x = x_initial + u
        
        # Acceleration (Newmark-β)
        a = (v - prev_v)/(γ*dt) - (1-γ)/γ * prev_a
        
        # Solve linear system
        K, R = AssembleSystem()
        ApplyBCsToMatrix(K, R)
        dv = Solve(K, R)
        
        # Update correction
        delta_v += dv
        
        # Check convergence
        converged = (norm(dv) < tol) and (norm(R) < ftol)
    
    # 3. Finalize step
    prev_v = v
    prev_a = a
    Time += dt
```

## Key Components

### 1. Time Step Initialization
- `prev_v`: Stores velocities from previous step
- `delta_v`: Initialized to zero (NR correction accumulator)
- `x_initial`: Reference configuration for displacement

### 2. Newton-Raphson Loop
| Operation | Mathematical Form | Purpose |
|-----------|-------------------|---------|
| Velocity Update | `v = prev_v + delta_v` | Current velocity estimate |
| BC Application | `v[BC_nodes] = v_BC` | Enforces constraints |
| Displacement | `u = dt * v` | Integrated motion |
| Acceleration | `a = (v-prev_v)/(γΔt) - (1-γ)/γ * prev_a` | Newmark time integration |
| Linear Solve | `K·dv = R` | Computes velocity correction |

### 3. Newmark Parameters
```math
\begin{cases}
\beta = 0.25 & \text{(Unconditional stability)} \\
\gamma = 0.5 & \text{(Second-order accuracy)}
\end{cases}
```

### Boundary Condition Handling
- Applied at **two levels**:
  1. Direct velocity overwrite (`ImposeBCV`)
  2. Linear system constraints (`ApplyBCsToMatrix`)

## Visualization
![Implicit Solver Animation](implicit_solver_anim.gif)

> Generated using Manim with the provided Python code