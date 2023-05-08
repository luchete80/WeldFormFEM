!NODAL CALCULATION
!Calculate Lumped matrix

!Predictor 
!uest_n+1 = un + dt v_n + dt2/2 a_n
!Estimate u and vel from previous steps
!Solve eqns of motion at t_n+1 = tn +dt
!Calculate a from  M dacc = fext (tn+1) - fint(uest, vest) -fcont
!Solve motion at t_n+1
!Update vel with CURRENT ACCELERATION