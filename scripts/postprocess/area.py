import numpy as np

def parse_vtk(file_path):
    """Parsea un archivo VTK de forma robusta."""
    with open(file_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    # --- Extraer PUNTOS ---
    points = []
    points_start = None
    num_points = 0

    for i, line in enumerate(lines):
        if line.upper().startswith("POINTS"):
            parts = line.split()
            num_points = int(parts[1])  # Ej: "POINTS 123 float" → 123
            points_start = i + 1
            break

    if points_start is None:
        raise ValueError("Sección 'POINTS' no encontrada.")

    # Leer todas las coordenadas hasta alcanzar 3*num_points valores
    current_line = points_start
    while len(points) < 3 * num_points and current_line < len(lines):
        parts = lines[current_line].split()
        for part in parts:
            try:
                points.append(float(part))
            except ValueError:
                break  # Si falla, asumimos que terminaron los puntos
        current_line += 1

    points = np.array(points).reshape(-1, 3)
    
    print("POINTS: ",len(points))
    # --- Extraer CELDAS (Tetraedros) ---
    tetrahedra = []
    cells_start = None

    for i, line in enumerate(lines):
        if line.upper().startswith("CELLS"):
            cells_start = i + 1
            break

    if cells_start is not None:
        current_line = cells_start
        while current_line < len(lines):
            line_content = lines[current_line]
            # Detenerse si se encuentra una nueva sección (ej: "CELL_TYPES")
            if line_content.upper().startswith(("CELL_TYPES", "POINT_DATA", "CELL_DATA")):
                break
            parts = line_content.split()
            # Solo procesar líneas que comienzan con un número (cantidad de nodos)
            if parts and parts[0].isdigit():
                num_nodes = int(parts[0])
                if num_nodes == 4:  # Tetraedro (4 nodos)
                    tetrahedra.append(list(map(int, parts[1:5])))
            current_line += 1
    
    print ("Cells ",len(tetrahedra))
    # --- Extraer CELL_DATA (elem_area) ---
    elem_area = None
    cell_data_start = None

    for i, line in enumerate(lines):
        if line.upper().startswith("CELL_DATA"):
            cell_data_start = i
            break

    if cell_data_start is not None:
        for i in range(cell_data_start, min(cell_data_start + 20, len(lines))):
            if "ele_area" in lines[i].lower():
                # Los datos comienzan después de "LOOKUP_TABLE default"
                data_start = i + 2  # Saltar "SCALARS..." y "LOOKUP_TABLE..."
                elem_area = []
                while data_start < len(lines):
                    line = lines[data_start].strip()
                    # Detener si encontramos una nueva sección o línea vacía
                    if not line or line.startswith(("SCALARS", "LOOKUP", "POINT_DATA", "CELL_TYPES")):
                        break
                    try:
                        # Convertir cada línea a float (una valor por línea)
                        elem_area.append(float(line))
                    except ValueError:
                        break  # Si falla la conversión, terminar
                    data_start += 1
                break
    return points, tetrahedra, elem_area
def compute_face_area(points, face):
    """Calcula el área de una cara triangular."""
    p0, p1, p2 = points[face[0]], points[face[1]], points[face[2]]
    cross = np.cross(p1 - p0, p2 - p0)
    return 0.5 * np.linalg.norm(cross)

def main():
    file_path = "out_0.000000.vtk"  # Cambia esto
    points, tetrahedra, elem_area = parse_vtk(file_path)
    
    if elem_area is None:
        print("Error: 'ele_area' no encontrado en CELL_DATA.")
        return
    
    z_threshold = 0.03  # Valor central
    z_tolerance = 1e-4  # ± tolerancia para Z (ajústala)
    area_tolerance = 0.01  # Tolerancia del 1% para el área
    
    elem_area_calculated = np.zeros(len(tetrahedra))
    total_area_calculated = 0.0
    #total_area_vtk = sum(elem_area)
    total_area_vtk = 0.0

    elcount = 0
                
    print ("Points size: ",len(points))
    for elem_id, tetra in enumerate(tetrahedra):
        faces = [
            [tetra[0], tetra[1], tetra[2]],
            [tetra[0], tetra[1], tetra[3]],
            [tetra[1], tetra[2], tetra[3]],
            [tetra[0], tetra[2], tetra[3]]
        ]
        for face in faces:

            #print("Points array shape:", points.shape if points is not None else "None")

            z_values = [points[node][2] for node in face]
            if all(abs(z - z_threshold) <= z_tolerance for z in z_values):
                area = compute_face_area(points, face)
                elem_area_calculated[elem_id] += area
                total_area_calculated += area
                total_area_vtk += elem_area[elem_id]
                elcount +=1
                print (elem_id,",",elem_area[elem_id])
                
    print ("Elements length: ",len(elem_area))
    # Comparación por elemento y total
    discrepancies = []
    for elem_id in range(len(tetrahedra)):
        diff = abs(elem_area_calculated[elem_id] - elem_area[elem_id])
        rel_diff = diff / elem_area[elem_id] if elem_area[elem_id] != 0 else 0
        if rel_diff > area_tolerance:
            discrepancies.append((elem_id, elem_area[elem_id], elem_area_calculated[elem_id], rel_diff))
    
    # Resultados
    print(f"\n→ Área total calculada: {total_area_calculated:.6e}")
    print(f"→ Área total en VTK:    {total_area_vtk:.6e}")
    print(f"→ Diferencia total:     {abs(total_area_calculated - total_area_vtk):.6e} ({abs(total_area_calculated - total_area_vtk)/total_area_vtk*100:.2f}%)\n")
    print("Elementos en zona: ", elcount)
    print(f"→ Elementos analizados: {len(elem_area)}")
    print(f"→ Elementos con error > {area_tolerance*100}%: {len(discrepancies)}")
    
    if discrepancies:
        print("\nDetalle de discrepancias (ID, elem_area VTK, elem_area calculada, % error):")
        for entry in discrepancies[:15]:  # Muestra las 10 primeras
            print(f"  {entry[0]:4d}: {entry[1]:.6e}  vs  {entry[2]:.6e}  ({entry[3]*100:.2f}%)")
    else:
        print("¡Todas las áreas coinciden dentro de la tolerancia!")

if __name__ == "__main__":
    main()
