#undef RSS_SL14
#undef RSS_SZ09
#define RSS_TR13
#undef Soilbeta
#undef BCC
#undef P_WLR
#undef MI_WLR
#undef MA_WLR
#define M_Q

MODULE MOD_SoilSurfaceResistance
  ! -----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate the soil surface resistance by using three parameterization schemes
  !
  ! ORIGINAL:
  ! Zhuo Liu, June, 2023
  !
  ! -----------------------------------------------------------------------
  ! !USE

   USE MOD_Precision
   !USE MOD_Vars_Global
   !USE MOD_Const_Physical, only : denice, denh2o  ! physical constant
   IMPLICIT NONE
   SAVE

   PUBLIC :: SoilSurfaceResistance


CONTAINS
!-----------------------------------------------------------------------

 SUBROUTINE SoilSurfaceResistance (nl_soil,forc_rhoair,hksati,porsl,bsw,psi0,&
                   dz_soisno,t_soisno,wliq_soisno,wice_soisno,wfc,fsno,qg,dsl,rss)
  !=======================================================================
  ! !DESCRIPTION:
  !
  ! REFERENCES:
  !
  !
  !=======================================================================

   !USE MOD_Precision
   !USE MOD_Const_Physical, only : denice, denh2o  ! physical constant
   IMPLICIT NONE


!-----------------------Argument-----------------------------------------

   integer, intent(in) :: &
        nl_soil            ! upper bound of array

   real(r8), intent(in) :: &
        porsl           , &! soil porosity [-]
        bsw             , &! Clapp-Hornberger "B"
        psi0            , &! saturated soil suction (mm) (NEGATIVE)
        dz_soisno       , &! layer thickness (m)
        t_soisno        , &! soil/snow skin temperature (K)
        wliq_soisno     , &! liquid water (kg/m2)
        wice_soisno     , &! ice lens [kg/m2]
        hksati          , &! hydraulic conductivity at saturation (mm h2o/s)
        qg              , &! ground specific humidity [kg/kg]
        wfc             , &!
        fsno            , &!
        forc_rhoair        ! density air [kg/m**3]

   real(r8), intent(out) :: &
        rss             , &! soil surface resistance(m/s)
        dsl                ! dry soil layer thickness


!-----------------------Local Variables------------------------------

   real(r8) :: &
        vol_liq,          & ! vol_liq
        smp_node,         & ! matrix potential
        eff_porosity,     & ! effective porosity = porosity - vol_ice
        aird,             & ! air free pore space
        d0,               & ! molecular diffusivity of water vapor in air (m2/s)
        eps,              & ! air filled pore space
        dg,               & ! d0*(tortuosity of the vapor flow paths through the soil matrix)
        dw,               & ! aqueous diffusivity (m2/s)
        hk,               & ! hydraulic conductivity
        d1,               & !
        beta,             & !
        tao,              & !
        B                   ! liquid water density / water vapor density

!-----------------------End Variables list---------------------------

   ! calculate the top soil volumetric water content (m3/m3), soil matrix potential and soil hydraulic conductivity
   !NOTE: max(wliq_soisno,1e-6)/...?
   !vol_liq  = max(wliq_soisno/(1000.*dz_soisno),0.001)
   vol_liq  = max(wliq_soisno,1e-6)/(1000.*dz_soisno)
   smp_node = psi0*(vol_liq/porsl)**(-bsw)
   hk       = hksati*(vol_liq/porsl)**(2.*bsw+3.)

   ! eff_porosity not calculated til SoilHydrolog
   eff_porosity = max(0.01_r8,porsl-min(porsl, wice_soisno/(dz_soisno*917.)))

   ! calculate diffusivity (dg, dw) and air free pore space
   aird        = porsl*(psi0/-1.e7_r8)**(1./bsw)
   d0          = 2.12e-5*(t_soisno/273.15)**1.75 ![Bitelli et al., JH, 08]
   eps         = porsl - aird

#ifdef BCC
   tao         = eps*eps*(eps/porsl)**(3._r8/max(3._r8,bsw))
#endif

#ifdef P_WLR
   tao         = 0.66*eps*(eps/porsl)
#endif

#ifdef MI_WLR
   tao         = eps**(4._r8/3._r8)*(eps/porsl)
#endif

#ifdef MA_WLR
   tao         = eps**(3./2.)*(eps/porsl)
#endif

#ifdef M_Q
   tao         = eps**(4._r8/3._r8)*(eps/porsl)**(2.0_r8)
#endif

   dg          = d0*tao
   dw          = -hk*bsw*smp_node/vol_liq


#ifdef RSS_SL14
   ! calculate dsl by SL14
   dsl         = dz_soisno*max(0.001_r8,(0.8*eff_porosity - vol_liq)) &
               / max(0.001_r8,(0.8*porsl- aird))

   dsl         = max(dsl,0._r8)
   dsl         = min(dsl,0.2_r8)
   rss         = dsl/dg + 20._r8
#endif

#ifdef RSS_SZ09
   ! calculate dsl by SZ09
   dsl         = dz_soisno*((exp((1._r8 - vol_liq/eff_porosity)**5) - 1._r8)/ (2.71828 - 1._r8))
   dsl         = min(dsl,0.2_r8)
   dsl         = max(dsl,0._r8)
   rss         = dsl/dg
#endif

#ifdef RSS_TR13
   ! calculate dsl by TR13
   B           = 1000._r8/(qg*forc_rhoair)
   d1          = (B*vol_liq*dw + eps*dg)/(B*vol_liq+eps)
   dsl         = dz_soisno*(1._r8/(2._r8*(B*vol_liq+eps)))
   dsl         = min(dsl,0.2_r8)
   dsl         = max(dsl,0._r8)
   rss         = dsl/d1
#endif


#ifdef Soilbeta
   wx   = (wliq_soisno/1000.+wice_soisno/917.)/dz_soisno

   !fac  = min(1._r8, wx/porsl)
   !fac  = max( fac, 0.01_r8 )
   !! Lee and Pielke 1992 beta, added by K.Sakaguchi
   IF (vol_liq < wfc ) THEN  !when water content of ths top layer is less than that at F.C.
      beta = (1._r8-fsno)*0.25_r8*(1._r8 - cos(vol_liq*3.1415926/wfc))**2._r8
   ELSE   !when water content of ths top layer is more than that at F.C.
      beta = 1._r8
   ENDIF
   ! raw = 50, NOTE: raw = 50. is just one specific case
   rss         = 50._r8 * (1._r8/beta-1._r8)
#endif

   rss         = min(1.e6_r8,rss)

 END Subroutine SoilSurfaceResistance

END MODULE MOD_SoilsurfaceResistance
