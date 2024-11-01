


Axisymmetric

'''
  if (bind_dom_type .eq. 3) then
    call calculate_element_shapeMat() !ONLY FOR VOLUMETRIC CALCS
    print *, "shape mat", elem%math(1,1,:,:)
    print *, "calc radis"
    call calculate_element_radius()   #THIS IS USED!!!!
    call calculate_element_MassMat ()
  end if

'''


'''
subroutine calc_nodal_masses ()
...
      if (bind_dom_type .eq. 3) then
        if (axisymm_vol_weight .eqv. .true.) then 
        nod%m(n) = nod%m(n) +  elem%vol(nod%nodel(n,ne)) * elem%rho(nod%nodel(n,ne),gp)/4.0d0 * elem%radius(nod%nodel(n,ne),gp)!!WEIGHT
        else
          nod%m(n) = nod%m(n) +  elem%vol(nod%nodel(n,ne)) * elem%rho(nod%nodel(n,ne),gp)/ &
                                                           (4.0d0 * elem%radius(nod%nodel(n,ne),gp) ) !! AREA WEIGHT
                                                           
'''