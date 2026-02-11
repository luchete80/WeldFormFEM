// ================================
// VTK_Writer.h - Escritor de archivos VTK para visualización
// ================================

#ifndef VTK_WRITER_TINY_H
#define VTK_WRITER_TINY_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

struct Point2D;  // Forward declaration

class VTKWriter {
public:
    // Escribe un archivo VTK en formato legacy
    static void writeVtkFile(const std::string& filename,
                            int step,
                            const std::vector<Point2D>& coords,
                            const std::vector<std::vector<int>>& elements,
                            const std::vector<double>& velocity,
                            const std::vector<double>& pressure_elem,  // P0 por elemento
                            const std::vector<double>& eps_bar_elem,  // Deformación por elemento
                            int nnodes, int nelem) {
        
        // Crear nombre de archivo con step
        std::string fullname = filename + "_step" + std::to_string(step) + ".vtk";
        std::ofstream vtk(fullname);
        
        if(!vtk.is_open()) {
            std::cerr << "Error: No se pudo crear el archivo " << fullname << std::endl;
            return;
        }
        
        // ========================================
        // CABECERA VTK
        // ========================================
        vtk << "# vtk DataFile Version 3.0\n";
        vtk << "Resultados de simulacion de forja - Paso " << step << "\n";
        vtk << "ASCII\n";
        vtk << "DATASET UNSTRUCTURED_GRID\n";
        
        // ========================================
        // COORDENADAS DE NODOS
        // ========================================
        vtk << "POINTS " << nnodes << " float\n";
        vtk << std::scientific << std::setprecision(6);
        
        for(int i = 0; i < nnodes; i++) {
            // En axisimetría, escribimos como 3D (r, 0, z) para visualización
            vtk << coords[i].x << " 0.0 " << coords[i].y << "\n";
        }
        
        // ========================================
        // CONECTIVIDAD DE ELEMENTOS
        // ========================================
        int nodes_per_elem = 4;
        int size_per_elem = nodes_per_elem + 1;  // +1 para el tipo de celda
        vtk << "CELLS " << nelem << " " << nelem * (nodes_per_elem + 1) << "\n";
        
        for(int e = 0; e < nelem; e++) {
            vtk << nodes_per_elem;
            for(int i = 0; i < nodes_per_elem; i++) {
                vtk << " " << elements[e][i];
            }
            vtk << "\n";
        }
        
        // ========================================
        // TIPOS DE CELDA (VTK_QUAD = 9)
        // ========================================
        vtk << "CELL_TYPES " << nelem << "\n";
        for(int e = 0; e < nelem; e++) {
            vtk << "9\n";  // VTK_QUAD
        }
        
        // ========================================
        // DATOS EN NODOS (PUNTO)
        // ========================================
        vtk << "POINT_DATA " << nnodes << "\n";
        
        // ---- VELOCIDAD (VECTOR) ----
        vtk << "VECTORS velocity float\n";
        for(int i = 0; i < nnodes; i++) {
            vtk << velocity[2*i] << " 0.0 " << velocity[2*i + 1] << "\n";
        }
        
        // ---- MAGNITUD DE VELOCIDAD (SCALAR) ----
        vtk << "SCALARS velocity_magnitude float 1\n";
        vtk << "LOOKUP_TABLE default\n";
        for(int i = 0; i < nnodes; i++) {
            double vr = velocity[2*i];
            double vz = velocity[2*i + 1];
            double mag = sqrt(vr*vr + vz*vz);
            vtk << mag << "\n";
        }
        
        // ========================================
        // DATOS EN ELEMENTOS (CELDA)
        // ========================================
        vtk << "CELL_DATA " << nelem << "\n";
        
        // ---- PRESIÓN POR ELEMENTO (P0) ----
        vtk << "SCALARS pressure_P0 float 1\n";
        vtk << "LOOKUP_TABLE default\n";
        for(int e = 0; e < nelem; e++) {
            vtk << pressure_elem[e] / 1e6 << "\n";  // Convertir a MPa
        }
        
        // ---- DEFORMACIÓN PLÁSTICA ACUMULADA ----
        vtk << "SCALARS plastic_strain float 1\n";
        vtk << "LOOKUP_TABLE default\n";
        for(int e = 0; e < nelem; e++) {
            vtk << eps_bar_elem[e] << "\n";
        }
        
        // ---- TASA DE DEFORMACIÓN EQUIVALENTE ----
        vtk << "SCALARS strain_rate float 1\n";
        vtk << "LOOKUP_TABLE default\n";
        
        // Calcular ε̇_eq para cada elemento
        for(int e = 0; e < nelem; e++) {
            const auto& conn = elements[e];
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
            
            // Punto central
            auto jac_result = jacobian_and_gradients(pos, 0.0, 0.0);
            double r_center = 0.0;
            for(int i = 0; i < 4; i++) {
                r_center += jac_result.N[i] * pos[i].x;
            }
            r_center = max(r_center, 1e-12);
            
            auto strain_result = calculate_strain_rate_martins(
                jac_result.dNdX, vel_elem, r_center, jac_result.N
            );
            
            vtk << strain_result.eps_dot_eq << "\n";
        }
        
        // ---- VISCOSIDAD EFECTIVA ----
        vtk << "SCALARS effective_viscosity float 1\n";
        vtk << "LOOKUP_TABLE default\n";
        for(int e = 0; e < nelem; e++) {
            const auto& conn = elements[e];
            std::vector<Point2D> pos(4);
            for(int i = 0; i < 4; i++) {
                pos[i] = coords[conn[i]];
            }
            
            std::vector<double> vel_elem(8);
            for(int i = 0; i < 4; i++) {
                int node = conn[i];
                vel_elem[2*i] = velocity[2*node];
                vel_elem[2*i + 1] = velocity[2*node + 1];
            }
            
            auto jac_result = jacobian_and_gradients(pos, 0.0, 0.0);
            double r_center = 0.0;
            for(int i = 0; i < 4; i++) {
                r_center += jac_result.N[i] * pos[i].x;
            }
            r_center = max(r_center, 1e-12);
            
            auto strain_result = calculate_strain_rate_martins(
                jac_result.dNdX, vel_elem, r_center, jac_result.N
            );
            
            double mu = effective_viscosity_norton(strain_result.eps_dot_eq);
            vtk << mu / 1e6 << "\n";  // Convertir a MPa·s
        }
        
        // ---- DIVERGENCIA (incompresibilidad) ----
        vtk << "SCALARS divergence float 1\n";
        vtk << "LOOKUP_TABLE default\n";
        for(int e = 0; e < nelem; e++) {
            const auto& conn = elements[e];
            std::vector<Point2D> pos(4);
            for(int i = 0; i < 4; i++) {
                pos[i] = coords[conn[i]];
            }
            
            std::vector<double> vel_elem(8);
            for(int i = 0; i < 4; i++) {
                int node = conn[i];
                vel_elem[2*i] = velocity[2*node];
                vel_elem[2*i + 1] = velocity[2*node + 1];
            }
            
            auto jac_result = jacobian_and_gradients(pos, 0.0, 0.0);
            double r_center = 0.0;
            for(int i = 0; i < 4; i++) {
                r_center += jac_result.N[i] * pos[i].x;
            }
            r_center = max(r_center, 1e-12);
            
            double dvr_dr = 0.0, dvz_dz = 0.0, vr_center = 0.0;
            for(int a = 0; a < 4; a++) {
                dvr_dr += jac_result.dNdX.getVal(0, a) * vel_elem[2*a];
                dvz_dz += jac_result.dNdX.getVal(1, a) * vel_elem[2*a + 1];
                vr_center += jac_result.N[a] * vel_elem[2*a];
            }
            
            double div_v;
            if(r_center < 1e-8) {
                div_v = 2.0 * dvr_dr + dvz_dz;
            } else {
                div_v = dvr_dr + vr_center/r_center + dvz_dz;
            }
            
            vtk << div_v << "\n";
        }
        
        vtk.close();
        std::cout << "  VTK escrito: " << fullname << std::endl;
    }
    
    // Versión simplificada para un solo paso
    static void writeVtkFile(const std::string& filename,
                            const std::vector<Point2D>& coords,
                            const std::vector<std::vector<int>>& elements,
                            const std::vector<double>& velocity,
                            const std::vector<double>& pressure_elem,
                            const std::vector<double>& eps_bar_elem,
                            int nnodes, int nelem) {
        writeVtkFile(filename, 0, coords, elements, velocity, 
                    pressure_elem, eps_bar_elem, nnodes, nelem);
    }
    
    // Escribe un archivo .pvd que agrupa múltiples VTK para animación en Paraview
    static void writePVDCollection(const std::string& basename, int nsteps) {
        std::string filename = basename + ".pvd";
        std::ofstream pvd(filename);
        
        pvd << "<?xml version=\"1.0\"?>\n";
        pvd << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
        pvd << "  <Collection>\n";
        
        for(int step = 0; step <= nsteps; step++) {
            pvd << "    <DataSet timestep=\"" << step * 0.001  // dt
                << "\" file=\"" << basename << "_step" << step << ".vtk\"/>\n";
        }
        
        pvd << "  </Collection>\n";
        pvd << "</VTKFile>\n";
        
        pvd.close();
        std::cout << "  PVD collection escrito: " << filename << std::endl;
    }
};

#endif // VTK_WRITER_H
