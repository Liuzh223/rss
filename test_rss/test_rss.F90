 
PROGRAM main
    
USE MOD_Soilsurface_resistance
USE MOD_Precision

    IMPLICIT NONE
    
    integer :: i 

    INTEGER ::           nl_soil                       ! upper bound of array


    REAL(r8) ::          porsl        ! soil porosity [-]
    REAL(R8) ::          bsw          ! Clapp-Hornberger "B"
    REAL(R8) ::          psi0         ! saturated soil suction (mm) (NEGATIVE)
    REAL(R8) ::          dz_soisno    ! layer thickness (m)
    REAL(R8) ::          t_soisno     ! soil/snow skin temperature (K)
    REAL(R8) ::          wliq_soisno  ! liquid water (kg/m2)
    REAL(R8) ::          wice_soisno  ! ice lens [kg/m2]
    REAL(R8) ::          hksati       ! hydraulic conductivity at saturation (mm h2o/s)
    REAL(R8) ::          qg           !
    REAL(R8) ::          forc_rhoair                   ! density air [kg/m**3]
    REAL(R8) ::          wfc          !
    REAL(R8) ::          fsno
    REAL(R8) ::          rss                           ! soil surface resistance(m/s)
    REAL(R8) ::          dsl    
        nl_soil         = 10
        forc_rhoair     = 1.29_r8
        hksati          = 1.25e-6_r8 
        porsl           = 0.45_r8
        bsw             = 7.12_r8
        psi0            = -1.21e-3_r8 
        dz_soisno       = 1.75e-2_r8
        t_soisno        = 297.15_r8
        qg              = 1.5e-2_r8
        wfc             = 0.31
        fsno            = 0.
        do i = 1,1000
           wliq_soisno = (i/1000.)*porsl*1000.*dz_soisno
           wice_soisno = 1.1e-4_r8  
           call soilsurface_resistance (nl_soil,forc_rhoair,hksati,porsl,bsw,psi0,&
                  dz_soisno,t_soisno,wliq_soisno,wice_soisno,wfc,fsno,qg,dsl,rss)
          write(*,*) rss
        enddo  
end program main
