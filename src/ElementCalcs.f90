subroutine calculate_strains ()
  integer :: e
  do while (e < elem_count)
    !!Assemble element displacement
    !elem%strain(e,:) = matmul(elem%bl(e,:),elem)
    e = e + 1 
  end do   
end subroutine