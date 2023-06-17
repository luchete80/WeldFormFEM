void readstep_(char* inputFile, void* outShape){

    STEPControl_Reader reader;
    reader = STEPControl_Reader();

    int succeed = reader.ReadFile(inputFile);

    if(!succeed){
        std::cout << "There was an error with the input file" << std::endl;
        return;
    }

    reader.NbRootsForTransfer();
    reader.TransferRoots();

    TopoDS_Shape myShape = reader.OneShape();
    TopoDS_Shape* myShapePtr = new TopoDS_Shape();
    (*myShapePtr) = myShape;

    outShape = myShapePtr;

    return;
}