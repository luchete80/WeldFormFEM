file.write("<?xml version=\"1.0\"?>\n")
file.write("<VTKFile type=\"UnstructuredGrid\">\n")
file.write("<UnstructuredGrid>\n")
file.write("<Piece NumberOfPoints=\"" << grid.Num_Verts() << "\" NumberOfCells=\"" << grid.Num_Cells() << "\">\n")
file.write("<Points>\n")
file.write("<DataArray NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\" >\n")
for (int n=0;n<grid.Num_Verts();++n) {
    for (int i=0; i<3; ++i) file<< setw(16) << setprecision(8) << scientific << grid.Vertex(n).Comp()[i] << endl;
}
file.write("</DataArray>\n")
file.write("</Points>\n")
file.write("<Cells>\n")

file.write("<DataArray Name=\"connectivity\" type=\"Int32\" format=\"ascii\" >\n")
for (int c=0;c<grid.Num_Cells();++c) {
    for (int n=0;n<grid.Cell(c).Num_Vertex();++n) {
        file.write(grid.Cell(c).Vert(n) << "\t";
    }
    file.write(endl;
}

file.write("</DataArray>\n")
file.write("<DataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\" >\n")
int offset=0;
for (int c=0;c<grid.Num_Cells();++c) {
    offset+=grid.Cell(c).Num_Vertex();
    file.write(offset << endl;
}
file.write("</DataArray>\n")

file.write("<DataArray Name=\"types\" type=\"UInt8\" format=\"ascii\" >\n")
for (int c=0;c<grid.Num_Cells();++c) {
    if (grid.Cell(c).Num_Vertex()==4) file.write("10\n") // Tetra
    if (grid.Cell(c).Num_Vertex()==8) file.write("12\n") // Hexa
    if (grid.Cell(c).Num_Vertex()==6) file.write("13\n") // Prism
    if (grid.Cell(c).Num_Vertex()==5) file.write("14\n") // Pyramid (Wedge)
}
file.write(endl;
file.write("</DataArray>\n");

file.write("</Cells>\n")

file.write("<CellData Scalars=\"scalars\" format=\"ascii\">\n")

//Begin data field output
file.write("<DataArray Name=\"";

file.write("Var";
file.write("\" type=\"Float32\" format=\"ascii\" >\n")

//If scalars
for (int n=0;n<grid.Num_Cells();n++)
    file.write(field.Val(n).Comp()[0] << endl;

file.write("</DataArray>\n")

// End of data output
file.write("</CellData>\n")

file.write("</Piece>\n")
file.write("</UnstructuredGrid>\n")
file.write("</VTKFile>\n")
file.close();