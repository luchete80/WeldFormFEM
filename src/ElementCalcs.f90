subroutine calculate_strain ()
  integer :: e
  do while (e < elem_count)
    !!Assemble element displacement
    !elem%strain(e,:) = matmul(elem%bl(e,:),elem)
    e = e + 1 
  end do   
end subroutine

subroutine calculate_strain_rate ()

end subroutine