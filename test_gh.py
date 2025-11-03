import os
import subprocess
from openai import OpenAI

# Ruta del archivo a mejorar

file_path = "src/common/ReMesher.C"

# Crear cliente con tu API Key

client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))

# Leer contenido original

# Obtener lista de archivos del repo
files = subprocess.check_output(["git", "ls-files"]).decode().splitlines()

for fpath in files:
    if fpath.endswith(".py"):
        with open(fpath, "r", encoding="utf-8") as f:
            code = f.read()  # ← estaba mal indentado
        print(f"\nAnalizando {fpath}...")

        prompt = f"Explicá brevemente qué hace este código:\n\n{code}"

        response = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=messages
        )   
        print(response.choices[0].message["content"])

# Nombre de la nueva rama

branch_name = "gpt-improvement-remesh"

# Crear rama nueva

subprocess.run(["git", "checkout", "-b", branch_name], check=False)

# Pedir sugerencias al modelo

messages = [
{"role": "system", "content": "Sos un experto en C++ y CUDA, especializado en algoritmos de remallado y optimización de código numérico."},
{"role": "user", "content": f"Analizá este archivo y mejorá la claridad, eficiencia y seguridad del código. Mantené la lógica original pero optimizá si es posible:\n\n{code}"}
]

response = client.chat.completions.create(
model="gpt-5",
messages=messages
)

improved_code = response.choices[0].message.content

# Guardar el código mejorado

with open(file_path, "w", encoding="utf-8") as f:
  f.write(improved_code)

# Commit y creación del Pull Request

subprocess.run(["git", "add", file_path])
subprocess.run(["git", "commit", "-m", "Mejoras automáticas de GPT en remesh.C"])
subprocess.run(["git", "push", "-u", "origin", branch_name])
subprocess.run(["gh", "pr", "create", "--fill", "--base", "main", "--head", branch_name])

print("Pull Request creado con las mejoras propuestas por GPT.")
