import pandas as pd

# Cargar datos
df = pd.read_csv('area.csv')

# 1. Suma total
total_area = df['ele_area'].sum()

# 2. Estadísticas básicas
stats = df['ele_area'].describe()

# 3. Filtrado de valores no cero
nonzero_area = df[df['ele_area'] > 0]['ele_area']

# Resultados
print(f"→ Suma total: {total_area:.6e}")
print(f"→ Elementos con área > 0: {len(nonzero_area)}")
print(f"→ Área promedio (no cero): {nonzero_area.mean():.6e}")
print("\nEstadísticas completas:")
print(stats.to_string())

# Opcional: Guardar resultados
results = pd.DataFrame({
    'Metrica': ['Suma total', 'Elementos con área >0', 'Área promedio (no cero)'],
    'Valor': [total_area, len(nonzero_area), nonzero_area.mean()]
})
 #results.to_csv('resultados_area.csv', index=False)
