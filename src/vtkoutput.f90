Module VTKOutput
  use ModPrecision, only : fp_kind
use NodeData
use ElementData
use Domain 
! <VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">
  ! <UnstructuredGrid>
    ! <Piece NumberOfPoints="8" NumberOfCells="6">
      ! <PointData Scalars="scalars">
        ! <DataArray type="Float32" Name="pressure" Format="ascii">
        
contains

subroutine WriteMeshVTU 
  implicit none
  integer :: e,dof,d,n
  open (1,file='output.vtu')!, position='APPEND')  
  write (1,*)  "<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""BigEndian"">"
  write (1,*)  "  <UnstructuredGrid>"
  !!!! FOR NOT RETURN CARRIAGE
  write (1, *) "<Piece NumberOfPoints=""" ,node_count, """ NumberOfCells=""",elem_count,""""  !!!!Note that an explicit format descriptor is needed when using
  !write (1, '(A,2x,I5)') '<Piece NumberOfPoints="'
  write (1,*) "<PointData Scalars=""scalars"">"
  write (1,*) "<DataArray type=""Float32"" Name=""Position"" NumberOfComponents=""3"" Format=""ascii"">"

  do n =1, node_count
    write (1,*) nod%x(n,1), nod%x(n,2), nod%x(n,3)
  end do
  
  write (1,*) "  </DataArray>"
  write (1,*) "</PointData>"
  
  write (1,*) "<DataArray type=""Int32"" Name=""connectivity"" Format=""ascii"">"
  do e=1, elem_count
    do n=1,nodxelem
      !write (1,*) elem%elnod(e,n)
      write(1,"(I,2x,I5)",advance="no") elem%elnod(e,n)
    end do
  end do
  
  close(1)
end subroutine WriteMeshVTU

End Module VTKOutput
 