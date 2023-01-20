module parameters

    implicit none
    
    character(3) :: case    ! one := First test case : convergence to a manufactured solution
                            ! two := Second test case : convergence to the diffusive limit
                            ! thr := Blast with friction

    integer, parameter :: PR = selected_real_kind(15,70)
    real(PR), parameter :: PI = 4._PR*atan(1._PR)
    real(PR) :: tf
    real(PR), public :: gamma
    real(PR), public :: theta_e
    integer :: imax, jmax, pas_affichage, kms
    real(PR) :: Lx, Ly, dx, dy, CFL

contains

    subroutine read_params()

        open(1,file = 'params.dat')
        read(1,*) case
        print*, case
        read(1,*) tf
        print*, tf
        read(1,*) gamma
        print*, gamma
        read(1,*) imax
        print*, imax
        read(1,*) jmax
        print*, jmax
        read(1,*) pas_affichage
        print*, pas_affichage
        read(1,*) Lx
        print*, Lx
        read(1,*) Ly
        print*, Ly
        read(1,*) CFL
        print*, CFL
        read(1,*) kms
        print*, kms
        dx = Lx/imax
        print*, dx
        dy = Ly/jmax
        print*, dy
        
    end subroutine read_params

    
end module parameters