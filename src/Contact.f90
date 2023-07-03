module Contact

use ModPrecision, only : fp_kind
use NodeData
use ElementData
use Domain
use class_ContMesh
 
contains

subroutine CalcCoontactForces(trimesh, dom) !!!! TODO: ADD DOMAIN
    implicit none
    type (Mesh), intent(in) :: trimesh
    type (dom_type), intent (in):: dom
    integer :: sn, me ! master node, slave element
		! Vec3_t xij, qj;
		! double h,K;
    !!!! LOOP FROM MASTER 
		! for (size_t a = 0; a < ContPairs[k].Size();a++) {

    real(fp_kind) :: vr(3)
      ! bool is_first = false;
      ! for (int m=0;m<meshcount;m++){
        ! if (Particles[ContPairs[k][a].first]->ID == contact_surf_id[m])
          ! is_first = true;
      ! }
      

      ! if (is_first) 	{ 	//Cont Surf is partcicles from FEM
        ! P1 = ContPairs[k][a].second; P2 = ContPairs[k][a].first; 	}
      ! else {
        ! P1 = ContPairs[k][a].first; P2 = ContPairs[k][a].second; } 
    
      ! vr = Particles[P1]->v - Particles[P2]->v;		//Fraser 3-137
      !Calculate velocity at element face centroid
      do sn=1, size(dom%slavenod)
        do me=1, size(dom%slavenod)
          !vr = dom%x(dom%slavenod(mn))- trimesh%v(e) 
			
      ! delta_ = - dot( Particles[P2]->normal , vr);	//Penetration rate, Fraser 3-138
      

        ! m = Particles[P2]->mesh;

				! e = trimesh[m]-> element[Particles[P2]->element];
        
        ! x_pred = Particles[P1]->x + Particles[P1]->v * deltat + Particles[P1]->a * deltat * deltat/2.0;
        ! vr_pred = Particles[P1]->v + Particles[P1]->a * deltat - Particles[P2]->v;
        
        ! dist =  dot (Particles[P2]->normal, x_pred ) - trimesh[m]-> element[Particles[P2]->element] -> pplane;

        ! if( dist  < Particles[P1]->h) {

          ! qj = Particles[P1]->x - dist * Particles[P2]->normal;
                                 ! //Check if it is inside triangular element
					! //Find a vector 
					! //Fraser 3-147
					! inside = true;
          ! if (trimesh[m]->dimension == 3){
            ! i=0;		
            ! while (i<3 && inside){
              ! j = i+1;	if (j>2) j = 0;
              ! crit = dot (cross ( *trimesh[m]->node[e -> node[j]] 
                                          ! - *trimesh[m]->node[e -> node[i]],
                                          ! qj  - *trimesh[m]->node[e -> node[i]]),
                                ! Particles[P2]->normal);
              ! if (crit < 0.0) inside = false;
              ! i++;
            ! }
          ! } else { //MESH DIMENSION = 2
            ! i=0;
            ! while (i<2 && inside){
              ! j = i+1;	if (j>1) j = 0;
              ! crit = dot ( *trimesh[m]->node[e -> node[j]] 
                                          ! - *trimesh[m]->node[e -> node[i]],
                                          ! qj  - *trimesh[m]->node[e -> node[i]]);
              ! if (crit < 0.0) inside = false;
              ! i++;
            ! }
          ! }
					
					! if (inside ) { //Contact point inside element, contact proceeds
            ! inside_geom++;
            ! end=true;

            ! delta = Particles[P1]->h - dist;

              ! kij = Particles[P1]->Mass / (deltat * deltat);
						

						! omega = sqrt (kij/Particles[P1]->Mass);
						! psi_cont = 2. * Particles[P1]->Mass * omega * DFAC; // Fraser Eqn 3-158
            

						! omp_set_lock(&Particles[P1]->my_lock);
                ! Particles[P1] -> contforce = (kij * delta - psi_cont * delta_) * Particles[P2]->normal; // NORMAL DIRECTION, Fraser 3-159    
              ! Particles[P1] -> delta_cont = delta;
						! omp_unset_lock(&Particles[P1]->my_lock);
            

            ! tgvr = vr + delta_ * Particles[P2]->normal;  // -dot(vr,normal) * normal, FRASER 3-168
            ! norm_tgvr = norm(tgvr);  
            ! tgdir = tgvr / norm_tgvr;              
            ! atg = Particles[P1] -> a - dot (Particles[P1] -> a,Particles[P2]->normal)*Particles[P2]->normal;


						! dt_fext = contact_force_factor * (Particles[P1]->Mass * 2. * norm(Particles[P1]->v) / norm (Particles[P1] -> contforce));


            ! Particles[P1] -> a += Particles[P1] -> contforce / Particles[P1] -> Mass; 
    end do !master element
  end do ! Mesh Slave Nodes
end subroutine CalcCoontactForces

end Module Contact