  ! // Compute the Jacobian​

  ! model->computeUnderJacobian(true);​

! ​

  ! // Compute the Mass Matrix if not already computed​

  ! model->computeMassMatrix();​

! ​

  ! computeTimeStep(true);​

  ! model->computeStrains();​

  ! model->computePressure();​ !!!FROM STRAIN INCREMENT!!

  ! model->computeStress(timeStep);​

  ! model->computeFinalRotation();​

  ! model->computeInternalForces();
  
  
  !!! AND DENSITY WHERE IS COMPUTED????
  