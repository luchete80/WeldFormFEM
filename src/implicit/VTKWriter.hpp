// ================================
// SIMULACIÓN PRINCIPAL (MARTINS)
// ================================
void run_simulation_martins() {
    // Inicialización
    initialize_mesh();
    auto fixed_dofs = setup_boundary_conditions();
    
    // ... (código existente) ...
    
    // Escribir VTK del paso inicial
    VTKWriter::writeVtkFile("forja", 0, coords, elements, 
                           velocity, pressure, eps_bar, 
                           nnodes, nelem);
    
    // Bucle temporal
    for(int step = 0; step < nsteps; step++) {
        std::cout << "\n--- PASO " << step+1 << "/" << nsteps 
                  << " (dt=" << dt << " s) ---" << std::endl;
        
        auto result = solve_step_martins(velocity, pressure, 
                                        fixed_dofs, K_temp, F_temp);
        
        velocity = result.velocity;
        pressure = result.pressure;
        
        // ... (actualizar deformación y coordenadas) ...
        
        // ESCRIBIR VTK DESPUÉS DE CADA PASO
        VTKWriter::writeVtkFile("forja", step+1, coords, elements, 
                               velocity, pressure, eps_bar, 
                               nnodes, nelem);
    }
    
    // Escribir colección PVD para animación
    VTKWriter::writePVDCollection("forja", nsteps);
    
    // Verificaciones físicas finales
    perform_physical_checks(velocity, pressure, eps_bar);
    
    // ...
}
