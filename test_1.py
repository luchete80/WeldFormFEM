import os
import subprocess
from openai import OpenAI

# Inicializar cliente
client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))

# Listar archivos del repo
files = subprocess.check_output(["git", "ls-files"]).decode().splitlines()

for fpath in files:
    if fpath.endswith(".py"):  # Podés cambiar la extensión
        with open(fpath, "r", encoding="utf-8") as f:
            code = f.read()

        print(f"\nAnalizando {fpath}...")

        # <-- Definimos 'messages' antes de usarlo
        messages = [
            {"role": "system", "content": "Sos un experto en Python y CUDA."},
            {"role": "user", "content": f"Explicá brevemente qué hace este código:\n\n{code}"}
        ]

        # Llamada a la API usando la nueva clase OpenAI
        response = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=messages
        )

        print(response.choices[0].message.content)

