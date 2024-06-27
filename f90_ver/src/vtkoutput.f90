!!! https://github.com/bhagath555/mesh_VTU
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

subroutine WriteMeshVTU (fname)
  implicit none
  character(len=*), intent(in) :: fname
  integer :: e,dof,d,n, offs
  real (fp_kind) :: z
  
  open (1,file=fname, action="write")!, position='APPEND')  
  write (1,*)  "<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""BigEndian"">"
  write (1,*)  "  <UnstructuredGrid>"
  !!!! FOR NOT RETURN CARRIAGE
  write (1,*)  "    <Piece NumberOfPoints=""" ,node_count, """ NumberOfCells=""",elem_count,""">"  !!!!Note that an explicit format descriptor is needed when using
  !write (1, '(A,2x,I5)') '<Piece NumberOfPoints="'
  write (1,*) "      <Points>"
  write (1,*) "        <DataArray type=""Float32"" Name=""Position"" NumberOfComponents=""3"" Format=""ascii"">"

  z = 0.0d0
  do n =1, node_count
    if (dim .eq. 3) then
          z = nod%x(n,3)
    end if
    write (1,*) nod%x(n,1), nod%x(n,2), z !!! No problem with putting the three component at point coords
    

  end do 
  
  write (1,*) "         </DataArray>"
  write (1,*) "       </Points>"
  
  write (1,*) "       <Cells>" 
  write (1,*) "         <DataArray type=""Int32"" Name=""connectivity"" Format=""ascii"">"
  do e=1, elem_count
    do n=1,nodxelem
      write(1,"(I10,1x)",advance="no") elem%elnod(e,n) - 1
    end do
  end do
  write (1,*) "" !!!! END OF LINE
  write (1,*) "        </DataArray>"
  write (1,*) "        <DataArray type=""Int32"" Name=""offsets"" Format=""ascii"">"
  
  offs = nodxelem
  do e=1, elem_count
      write(1,"(I10,1x)",advance="no") offs
      offs = offs + nodxelem
  end do
  write (1,*) "" !!!! END OF LINE
  write (1,*) "        </DataArray>"

  write (1,*) "        <DataArray type=""Int32"" Name=""types"" Format=""ascii"">"  
  do e=1, elem_count
      if (dim .eq. 2) then
        write(1,"(I3,1x)",advance="no") 9
      else 
        write(1,"(I3,1x)",advance="no") 12
      end if
  end do
  write (1,*) "" !!!! END OF LINE
  write (1,*) "        </DataArray>"
  write (1,*) "      </Cells>"
  
    !!!!!!! POINT DATA
  write (1,*) "      <PointData Scalars=""scalars"">"
  !!!!!!!!!!!_--------------------------------------
  
  write (1,*) "        <DataArray type=""Float32"" Name=""u"" NumberOfComponents=""",dim,""" Format=""ascii"">"
  do n =1, node_count
		do d =1, dim
    write (1,"(e13.5)") nod%u(n,d)
		end do
	end do   
  write (1,*) "        </DataArray>"  

  write (1,*) "        <DataArray type=""Float32"" Name=""v"" NumberOfComponents=""", dim ,""" Format=""ascii"">"
  do n =1, node_count
		do d =1, dim
    write (1,*) nod%v(n,d)
		end do
	end do   
  write (1,*) "        </DataArray>"  

  write (1,*) "        <DataArray type=""Float32"" Name=""a"" NumberOfComponents=""",dim,""" Format=""ascii"">"
  do n =1, node_count
		do d =1, dim
			write (1, "(e13.4)") nod%a(n,d)
		end do
	end do   
  write (1,*) "        </DataArray>" 

  ! write (1,*) "        <DataArray type=""Float32"" Name=""f_hourg"" NumberOfComponents=""3"" Format=""ascii"">"
  ! do n =1, node_count
    ! write (1,*) nod%f_hour(n,1), nod%f_hour(n,2), nod%f_hour(n,3)
  ! end do   
  ! write (1,*) "        </DataArray>" 
  
  !write (1,*) "      </PointData>"

  !write (1,*) "      <PointData Scalars=""scalars"">"
  write (1,*) "        <DataArray type=""Float32"" Name=""density"" NumberOfComponents=""1"" Format=""ascii"">"
  do n =1, node_count
    write (1,*) nod%rho(n)
  end do   
  write (1,*) "        </DataArray>"  

  write (1,*) "        <DataArray type=""Float32"" Name=""mass"" NumberOfComponents=""1"" Format=""ascii"">"
  do n =1, node_count
    write (1,*) nod%m(n)
  end do   
  write (1,*) "        </DataArray>"  
  
  write (1,*) "        <DataArray type=""Float32"" Name=""sigma_eq"" NumberOfComponents=""1"" Format=""ascii"">"
  do n =1, node_count
    write (1,*) nod%sigma_eq(n)
  end do   
  write (1,*) "        </DataArray>"    

  write (1,*) "        <DataArray type=""Float32"" Name=""pl_strain"" NumberOfComponents=""1"" Format=""ascii"">"
  do n =1, node_count
    write (1,*) nod%pl_strain(n)
  end do   
  write (1,*) "        </DataArray>"     

  ! write (1,*) "        <DataArray type=""Float32"" Name=""sigma"" NumberOfComponents=""6"" Format=""ascii"">"
  ! do n =1, node_count
    ! write (1,*) nod%sigma(n,1,1), nod%sigma(n,2,2), nod%sigma(n,3,3), &
                ! nod%sigma(n,1,2), nod%sigma(n,2,3), nod%sigma(n,3,1)!, &
                ! !nod%sigma(n,3,1), nod%sigma(n,3,2), nod%sigma(n,3,3)
  ! end do   
  ! write (1,*) "        </DataArray>"    
   
  
  
  write (1,*) "      </PointData>"
  
  ! <DataArray type="Int32" Name="offsets" Format="ascii">
    ! 3  6  9  12  15  18
  ! </DataArray>
  ! <DataArray type="Int32" Name="types" Format="ascii">
    ! 5  5  5  5  5  5
  ! </DataArray>
 
  write (1,*) "    </Piece>"
  write (1,*) "  </UnstructuredGrid>"
  write (1,*) "</VTKFile>"
  
  close(1)
end subroutine WriteMeshVTU

End Module VTKOutput
 