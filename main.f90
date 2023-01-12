program main

   use parameters
   use cartesian_mesh
   use flux
   use solver
   use high_order
   use initial_condition
   use boundary_conditions
   use time_scheme
   use computing_error
   use vtk

   implicit none

   ! Déclaration des variables
   integer :: i, j, k, kms, n, iteration
   real(PR), dimension(:), allocatable :: x, y, xm, ym
   real(PR), dimension(:,:), allocatable :: Rho, u, v, kappa, p
   real(PR), dimension(:,:), allocatable :: Rhop1, up1, vp1, sigma
   real(PR), dimension(:,:), allocatable :: Rho_ex
   real(PR), dimension(:,:,:), allocatable :: Un, Unp1, Uex, Sn, phi
   real(PR) :: t0, t, dt
   real(PR) :: PD, SPEED
   real(PR) :: l1error, l2error
   character(100) :: numero, num

   ! Creating the file in which we are gonna store the L1 and L2 errors for different meshes
   open(25, file='Convergence_errors.dat', access = 'append')


   ! Initialisation
   k = 0
   kms = 640
   iteration = 0
   t0 = 0._PR
   t = 0._PR
   PD = 0._PR
   SPEED = 0._PR
   l1error = 0._PR
   l2error = 0._PR



   ! Allocation des vecteurs 1D pour le maillage
   allocate(x(1:imax+1))
   allocate(y(1:jmax+1))
   allocate(xm(1:imax+1))
   allocate(ym(1:jmax+1))

   ! Allocation des vecteurs pour la solution
   allocate(Rho(1:imax+1,1:jmax+1))
   allocate(u(1:imax+1,1:jmax+1))
   allocate(v(1:imax+1,1:jmax+1))
   allocate(p(1:imax+1,1:jmax+1))
   allocate(sigma(1:imax+1,1:jmax+1))
   allocate(Rhop1(1:imax+1,1:jmax+1))
   allocate(up1(1:imax+1,1:jmax+1))
   allocate(vp1(1:imax+1,1:jmax+1))
   allocate(kappa(1:imax+1,1:jmax+1))
   allocate(Un(1:imax+1,1:jmax+1,1:3))
   allocate(Sn(1:imax+1,1:jmax+1,1:3))
   allocate(Unp1(1:imax+1,1:jmax+1,1:3))
   allocate(phi(1:imax+1,1:jmax+1,1:3))

   ! Chargement du maillage 2D cartésien
   call load_mesh(x,y,xm,ym)

   ! Initialisation
   Rho = 0._PR; u = 0._PR; v = 0._PR; p = 0._PR
   Rhop1 = 0._PR; up1 = 0._PR; vp1 = 0._PR; sigma = 0._PR
   Un = 0._PR; Unp1 = 0._PR; phi = 0._PR

   ! Conditions initiales
   call init(Rho,u,v,p,sigma,xm,ym,kms,t0,Un)

   ! Densité exacte
   Rho_ex = Rho

   !num = 1000
   write(num,*) 123456789
   call sortie_vtk(num,imax,jmax,x,y,Rho,u,v)

   ! Boucle en temps
   do while (t <= tf)

      ! Computing time step
      !CFL <= min(dx, dy) / max(|u| + c, |v| + c) où c = sqrt(p_prime)
      PD = pressure_prime(maxval(Rho))
      dt = dx*dy/(8._pr * maxval((abs(u) + sqrt(PD))/dx + (abs(v) + sqrt(PD))/dy) * (dx + dy) &
         & + maxval(sigma)*dx*dy)

      print*, "The time step is equal to:", dt, "."

      ! Paramètre thete_e permettant la limitation
      theta_e(1:3) = (/1._PR/(1._PR + maxval(sigma)*t), 1._PR, 1._PR/)

      ! Boundary conditions
      call Neumann(Rho,u,v)

      ! Solution using an explicit Euler solver or RK3 scheme
      !call RK3(Un,t,dt,kms,xm,ym,Unp1)
      call explicit_euler(Un,t,dt,kms,xm,ym,Unp1)

      do j = 2, jmax

         do i = 2, imax

            ! Swapping solutions
            call conservative_to_non_conservative(Unp1(i,j,1:3),Rhop1(i,j),up1(i,j),vp1(i,j))
            call swapping(Rhop1(i,j),up1(i,j),vp1(i,j),Rho(i,j),u(i,j),v(i,j))
            call non_conservative_to_conservative(Rho(i,j),u(i,j),v(i,j),Un(i,j,1:3))
            
         end do
         
      end do

      if (mod(iteration,pas_affichage) == 0) then ! On stocke la solution par fréquence de 
                                          ! 100 pas de temps  
         k=k+1

         write(numero,*) k

         if (k <= 9) then

            numero = '000'//trim(adjustl(numero))

         elseif (9 < k .and. k <= 99) then

            numero = '00'//trim(adjustl(numero))

         elseif (99 < k .and. k <= 999) then

            numero = '0'//trim(adjustl(numero))

         end if
         
         call sortie_vtk(numero,imax,jmax,x,y,Rho,u,v)

      end if

      iteration = iteration + 1
      t = t + dt
      print*, "On est à la ", iteration, "itération."
      print*, "Le temps est égal à :", t, "."

   end do

   ! Calcul des erreurs L1 et L2
   l1error = L1_Error(Rho_ex,Rho)
   l2error = L2_Error(Rho_ex,Rho)
   write(25,*) l1error, l2error
  
end program main