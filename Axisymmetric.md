



Axisymmetric

1. 
Addad bind_dom_type and axisymm_vol_weight

- Added radius calc via N functions

2. 
- Added nodal mass calc (calc_nodal_masses) by using radius
  - If are weight: sum (rho x V x r / 4)
  _ If vol weight: 
  
3. Calc Vol (used?), only in vol weight

3.
- Internal Forces: (cal_elem_forces () inside mechanical)


4. Strain Rate calc 
- 
        !!! er hoop = vr/r
        if (dim .eq. 2 .and. bind_dom_type .eq. 3) then 
          ! if (elem%gausspc(e) .eq. 1) then
            elem%str_rate(e,gp, 3,3) = elem%str_rate(e,gp, 3,3) + 0.25d0*elem%vele (e,dim*(n-1)+1,1) / elem%radius(e,gp) !! 0.25 is shapemat
          ! print *, "hoop er", elem%str_rate(e,gp, 3,3) 
          elem%rot_rate(e,gp, 3,3) = 0.0d0
          ! end if
        end if 


Used:
flags/type:
bind_dom_type, 

Vars
vol
m_radius



CODE

'''
  if (
Addad bind_dom_type and  .eq. 3) then
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
-----------------------------------------
cal_elem_forces () in Mechanical (in C version is dev_t void Domain_d::calcElemForces(){)


Elem Forces:

      if (axisymm_vol_weight .eqv. .true.) then
        fc = elem%radius(e,gp)
      end if      
      
      
                  if (axisymm_vol_weight .eqv. .true.) then      
            !!! RADIUS IS CANCELLED IN THE SECOND TERMS             
            elem%f_int(e,n,1) = elem%f_int(e,n,1) + (elem%dHxy_detJ(e,gp,2,n) * elem%sigma (e,gp, 1,2)) &
                                                    * elem%radius(e,gp) &
                                                    + 0.25d0*(elem%sigma (e,gp, 1,1) - &
                                                    elem%sigma (e,gp, 3,3) ) * elem%detJ(e,gp)
              ! print *, "term 1 ", (elem%dHxy_detJ(e,gp,2,n) * elem%sigma (e,gp, 1,2)) &
                                                    ! * elem%radius(e,gp) &
              ! print *, "term 2 ", 0.25d0*(elem%sigma (e,gp, 1,1) - &
                                                    ! elem%sigma (e,gp, 3,3) ) * elem%detJ(e,gp)
            elem%f_int(e,n,2) = elem%f_int(e,n,2) + (elem%dHxy_detJ(e,gp,1,n) * elem%sigma (e,gp, 1,2)) & 
                                                    * elem%radius(e,gp) &
                                                    + 0.25d0*elem%sigma (e,gp, 1,2) * elem%detJ(e,gp)
            else
              fa = 0.25d0/elem%radius(e,gp) * elem%detJ(e,gp) !!! THEN IS WEIGHTED BY 4 in case of gauss point =1
              !!! AREA WEIGHTED, BENSON EQN 2.4.3.2
              !!! 2.4.3.2 remains sig * Area/(4 r0), which is (4detJ)/(4r0) = detJ /r0
              !!! LATER IS MULTIPLIED BY WEIGHT WICH GIVES THE AREA

              elem%f_int(e,n,1) = elem%f_int(e,n,1) + elem%dHxy_detJ(e,gp,2,n) * elem%sigma (e,gp, 1,2) - &
                                                     (elem%sigma (e,gp, 1,1) - elem%sigma (e,gp, 3,3) ) * fa
                                                     
              elem%f_int(e,n,2) = elem%f_int(e,n,2) + elem%dHxy_detJ(e,gp,1,n) * elem%sigma (e,gp, 1,2) - &
                                                     elem%sigma (e,gp, 1,2) * fa   
