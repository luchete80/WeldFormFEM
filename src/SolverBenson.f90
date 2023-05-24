!!!!BENSON 1992

! (1) Knowing the stress, pressure, hourglass forces and shock viscosity at t” in each zone or
! element, the forces at the nodes are calculated. The accelerations of the nodes are
! calculated by dividing the nodal forces by the nodal masses.
! (2) The acceleration is integrated to give the velocity at t”+l”.
! (3) The velocity is integrated to give the displacement at t”+‘.
! (4) The constitutive model for the strength of the material is integrated from t” to t”+’ now
! that the motion of the material is known.
! (5) The artificial shock viscosity and hourglass viscosity are calculated from u”+“~.
! (6) The internal energy is updated based on the work done between t” and t”+‘.
! (7) Based on the density and energy at t”+l, the pressure is calculated from the equation of
! state.
! (8) A new time step size is calculated based on the speed of sound through each of the
! elements and their geometry.
! (9) Advance the time and return to step (1)