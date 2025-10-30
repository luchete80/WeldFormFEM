#!/usr/bin/env python3
import os
import re
import sys

# Patrón: elimina desde la primera línea /****** hasta la última línea igual
HEADER_PATTERN = re.compile(
    r"(?s)^/\*{5,}.*?\n\s*/\*{5,}\*/\n*",  # Captura el bloque entero
    re.MULTILINE
)

def remove_header(filepath):
    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()

    new_content, num_subs = HEADER_PATTERN.subn("", content, count=1)

    if num_subs > 0:
        lines_removed = content.count("\n") - new_content.count("\n")
        with open(filepath, "w", encoding="utf-8") as f:
            f.write(new_content)
        print(f"🧹 Removed header from: {filepath} ({lines_removed} lines)")
    else:
        print(f"— No header found in: {filepath}")

def process_path(target_path):
    if os.path.isfile(target_path):
        # Si es un archivo directo
        if target_path.endswith((".cpp", ".c", ".h", ".hpp")):
            remove_header(target_path)
        else:
            print(f"⚠️  Skipping non-C/C++ file: {target_path}")
    elif os.path.isdir(target_path):
        # Si es un directorio, recorrer recursivamente
        for root, _, files in os.walk(target_path):
            if any(skip in root for skip in ["build", "external", "third_party"]):
                continue
            for file in files:
                if file.endswith((".cpp", ".c", ".h", ".hpp")):
                    filepath = os.path.join(root, file)
                    remove_header(filepath)
    else:
        print(f"❌ Path not found: {target_path}")

def main():
    # Si no se pasa argumento → usa la raíz actual
    target_path = sys.argv[1] if len(sys.argv) > 1 else "."
    process_path(target_path)

if __name__ == "__main__":
    main()
