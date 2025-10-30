#!/usr/bin/env python3
import os
import sys
import datetime
import re

PROJECT_NAME = "WeldformFEM"
AUTHOR = "Luciano Buglioni"
LICENSE = "GNU General Public License v3.0 or later"
EMAIL = "weldform.sph@gmail.com"
WEB = "https://www.opensourcemech.com"
CURRENT_YEAR = datetime.datetime.now().year

HEADER_TEMPLATE = """\
/*************************************************************************/
/*  {filename:<60} */
/*  {project_name} - High-Performance Explicit & Implicit FEM Solvers     */
/*  (CPU/GPU, C++/CUDA)                                                  */
/*                                                                       */
/*  {email}                                                              */
/*  {web}                                                                */
/*                                                                       */
/*  Copyright (c) {start_year}-{end_year} {author:<25} */
/*                                                                       */
/*  This file is part of the {project_name} project.                     */
/*  Licensed under the {license}. See the LICENSE file in the project    */
/*  root for full license information.                                   */
/*************************************************************************/
"""

def generate_header(filename, start_year):
    return HEADER_TEMPLATE.format(
        filename=os.path.basename(filename),
        project_name=PROJECT_NAME,
        email=EMAIL,
        web=WEB,
        start_year=start_year,
        end_year=CURRENT_YEAR,
        author=AUTHOR,
        license=LICENSE
    )

def has_header(content):
    return "This file is part of the" in content and "Copyright" in content

def extract_start_year(content):
    m = re.search(r"Copyright \(c\) (\d{4})", content)
    if m:
        return m.group(1)
    return str(CURRENT_YEAR)

def remove_header(content):
    """Elimina el bloque de header que empieza con /******** y termina con */"""
    return re.sub(r"(?s)^/\*{9,}.*?\*/\s*", "", content, count=1)

def process_file(filepath, action):
    with open(filepath, "r", encoding="utf-8") as f:
        content = f.read()

    if action == "add":
        if has_header(content):
            # Actualiza año si ya tiene header
            start_year = extract_start_year(content)
            updated_header = generate_header(filepath, start_year)
            content = re.sub(r"(?s)^/\*{9,}.*?\*/", updated_header.strip(), content, count=1)
            print(f"Updated header in: {filepath}")
        else:
            start_year = str(CURRENT_YEAR)
            new_header = generate_header(filepath, start_year)
            content = new_header + "\n\n" + content
            print(f"Added header to: {filepath}")

    elif action == "remove":
        if has_header(content):
            content = remove_header(content)
            print(f"Removed header from: {filepath}")

    with open(filepath, "w", encoding="utf-8") as f:
        f.write(content)

def main():
    if len(sys.argv) < 2:
        print("Usage: python update_headers.py [--add | --remove] [optional_path]")
        sys.exit(1)

    action = sys.argv[1].replace("-", "")
    if action not in ("add", "remove"):
        print("Error: action must be --add or --remove")
        sys.exit(1)

    path = sys.argv[2] if len(sys.argv) > 2 else "."

    for root, _, files in os.walk(path):
        for file in files:
            if file.endswith((".cpp", ".c", ".h", ".hpp", "CMakeLists.txt", "README.md")):
                filepath = os.path.join(root, file)
                process_file(filepath, action)

if __name__ == "__main__":
    main()
