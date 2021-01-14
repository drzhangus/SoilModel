! SSCFM (Surface Soil Contaminant Fate Model)
! Define system of differential equations 
!===========================================================================================================================     
  subroutine derivative(T, PLt, Rst, Srt, Rnst, DT, YWRK, Y, YDOT)
    !
    use Global 
    use HydroVariable 
    implicit none
    !
    real(kind=8), Intent( INout ) :: T, PLt, Rst, Srt, Rnst, DT   
    real(kind=8), Intent( INout ) :: YWRK(*), Y(*)
    real(kind=8), Intent( OUT)    :: YDOT(*)
    real(kind=8) :: Fr1, Fe1, Fl1, Fdecay1, Fvol1, Zcheck  
    integer i, j, k

    ! DE 1: di - particle size
    if ( MMC == 1) then
      Y(1)  = 0.D0
      alpha = 0.D0
    else  
      ! not miscible
      if (sshape == 1) then
        x = 3.D0         
      else
        x = 2.D0
      endif

      if (Y(2) < Mmin) Y(2) = Mmin 
      if (abs(YWRK(2)) > nullerror) then
        Y(1) = YWRK(1) * ((Y(2)/YWRK(2))**(1.D0/x))           
      else
        Y(1) = Dmin
      endif          
      if (Y(1) < Dmin) Y(1) = Dmin           
      if (Y(1) > D0)   Y(1) = D0      
          
      if (sshape ==1) then        
        x = 3.D0  ! sphere
        if (abs(Y(1)) > nullerror) then
          alpha = 6.D0 / (Rous*Y(1)*1000000.D0)     
        else
          alpha = 0.D0
        endif
      else      
        x = 2.D0   ! cylinder
        if ( abs(Y(1)) > nullerror ) then  
          alpha = 1.D0 / (Rous*1000000.D0) * (2.D0/sl + 4.D0/Y(1))
        else
          alpha = 0.D0
        endif
      endif 
      !
    endif
    
    !! No DE for diameter
    YDOT(1) = 0.D0
    
    ! DE 2: Ms - solid mass (g)
    Cdp = Y(3) / (Rtmp*thetwtmp)
    if (MMC == 1) then
      Fprep = 0.D0
    else
      ! solubility limits
      ! A multiplier of 1000 is used for force equilibrium without rediculousl small time steps, Dortch 2015    
      if (Cdp > Cs) then                                    
        Fprep = 1000.D0 * PA * Zb * thetwtmp * (Cdp - Cs) / DT         
      else
        Fprep = 0.D0
      endif
    endif
    
    ! Dissolution flux (g/yr or g/day))
    if (MMC == 1) then
      Fdis = PLt
    else
      Fdis = RItmp * alpha * Cs * Y(2)          ! use rainfall rather than precip for dissolution 
    endif     
    if (dataset == 2) Fdis = Fdis * Zb / Zmax   
    
    if (Esolid == 1) then   
      if (dataset /= 2) then
       Fes = Y(2) * Etmp / Zb    
      elseif (dataset == 2) then
       Fes = Y(2) * Etmp / Zmax
      endif
    else
      Fes = 0.D0
    endif 
    
    !! YDOT(2) = DMs/dt
    YDOT(2) = PLt - Fes - Fdis + Fprep - Rst*Y(2) - Srt  
    
    ! DE 3: Ctt - total concentration per total volume (g/m3)
    if (dataset == 0) then
      if (Rotmp > 0) then
        pk = (soila*thets*RItmp*Fdptmp) / (de*Roub*thetwtmp*PNItmp)
        Fr1 = de * (1.D0 - exp(-pk)) * PNItmp * PA  ! Surface runoff flux (m3/yr)
      else
        Fr1 = 0.D0
      endif
      !
    else    
      Rer = 0.D0
      do i = 1,24
        pk  = (soila*PI24tmp(i)*thets*Fdptmp) / (de*Roub*thetwtmp)
        Rer = Rer + (de*(1.D0 - Exp(-pk)))
      enddo
      !
      if (Rotmp > 0) then
        Fr1 = Rer * PA                            ! Surface runoff flux (m3/day)
      else
        Fr1 = 0.D0
      endif
    endif
       
    Fe1     = Etmp * PA                           ! Erosion flux (m3/yr or m3/day)    
    Fl1     = qwtmp * Fdptmp / thetwtmp * PA      ! Leaching flux (m3/yr or m3/day)   
    Fdecay1 = (PKl*Fdptmp + PKa*Fpptmp) * (PA*Zb) ! Degradation flux (m3/yr or m3/day)      
    Fvol1   = PKv * Faptmp * PA                   ! Volatilization rate (m3/yr or m3/day)
 
    if (dataset == 2) then   
      Zcheck = Zmin + 0.001
      if (Zb < Zcheck) then
        Fdecay1 = 0.D0
        Fvol1 = 0.D0
      endif
    endif 
    
    !! YDOT(3) = DCtt/dt
    if (dataset /= 2) then
      YDOT(3) = ( Fdis - Fprep - (Fr1 + Fe1 + Fl1 + Fdecay1 + Fvol1 + (Rnst*PA*Zb))*Y(3) ) / (PA*Zb)       
    elseif (dataset == 2) then   
      YDOT(3) = ( Fdis - Fprep - (Fr1 + Fe1 + Fl1 + Fdecay1 + Fvol1 + (Rnst*PA*Zmax))*Y(3) ) / (PA*Zmax)      
    endif
    
end subroutine