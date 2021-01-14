! Fourth-order Runge-Kutta integration subroutine (RK4)
!===========================================================================================================================   
  subroutine RK4(NEQ, QPTime, DT, Y) 
    ! 
    use Global 
    use HydroVariable 
    implicit none
    !
    integer,      Intent( INout ) :: NEQ
    real(kind=8), Intent( INout ) :: QPTime, DT, Y(*) 
    real(kind=8) :: YT(NEQ), YDOT1(NEQ), YDOT2(NEQ), YDOT3(NEQ), YDOT4(NEQ), YWRK(NEQ), YDOT(NEQ)
    real(kind=8) :: TSAVED, PLt, Rst, Srt, Rnst, outQPTime, Zcheck  
    real(kind=8), save   :: TI22, TI33     
    integer i, j, k
                   
    do i = 1,NEQ              
      YWRK(i) = Y(i)            
    enddo
        
    TSAVED = QPTime
    if (dataset == 1 .or. dataset == 2) then      
      do while (.not.((QPTime >= Fracyeardaily(Idays)) .AND. (QPTime < Fracyeardaily(Idays + 1))))
        Idays = Idays + 1
      enddo 
      !
      do while (.not.((QPTime >= Fracyearhourly(Ihours)) .AND. (QPTime < Fracyearhourly(Ihours + 1))))
        Ihours = Ihours + 1
      enddo
      !    
      do i=1, 24
        PI24tmp(i) = PI(((Idays-1)*24) + i)
      enddo 
      !
      if (dataset == 2) then
        Zb = thickness(Idays)  
      endif
      
      PItmp    = PI(Ihours)
      thetwtmp = thetw(Idays)
      !PNItmp   = PNI(Idays)        
      PIttmp   = PIt(Idays) 
      RItmp    = RI(Idays)
      Rotmp    = Ro(Idays)      
      Etmp     = E(Idays)       
      qwtmp    = qw(Idays)             
      Fiftmp   = Fif(Idays)               
    endif
         
    if (dataset == 0) then      
      thetwtmp = thetw(1)
      PNItmp   = PNI(1)        
      PIttmp   = PIt(1)      
      !PItmp    = PI(1)
      RItmp    = RI(1)
      Rotmp    = Ro(1)      
      Etmp     = E(1)
      qwtmp    = qw(1)
      Fiftmp   = Fif(1)   
    endif
 
    Rtmp   = R(Idays)
    Fdptmp = Fdp(Idays)
    Fpptmp = Fpp(Idays)
    Faptmp = Fap(Idays)
    
    ! set loadings
    if (dataset == 0) then 
      PresTime = QPTime
    else      
      PresTime = QPTime/365.25D0 
    endif
    
    if (PresTime < yr(1)) PLt = 0.        
    if ((PresTime >= yr(1)) .AND. (PresTime < yr(nyrs))) then
      do i = 1,nyrs-1
        if ( (PresTime >= yr(i)) .AND. (PresTime < yr(i + 1))) PLt = PL(i)
      enddo       
    endif   
    if (PresTime >= yr(nyrs)) PLt = PL(nyrs)  
    
    ! BMP
    if (BMPflag == 1) then
      if (PresTime < BMPyr(1)) then
        Rst  = 0.
        Srt  = 0.
        Rnst = 0.
      endif
                        
      if ((PresTime >= BMPyr(1)) .AND. (PresTime < BMPyr(BMPnyrs))) then
        do i = 1,BMPnyrs-1
          if ((PresTime >= BMPyr(i)) .AND. (PresTime < BMPyr(i + 1))) then
            Rst  = Rs(i)
            Srt  = Sr(i)
            Rnst = Rns(i)
          endif  
        enddo
      endif
      
      if (PresTime >= BMPyr(BMPnyrs)) then
        Rst  = Rs(BMPnyrs)
        Srt  = Sr(BMPnyrs)
        Rnst = Rns(BMPnyrs)
      endif      
      !
    else
      Rst  = 0.
      Srt  = 0.
      Rnst = 0.
    endif
    
    !---------------------------------------------------------------------------------------------------  
    call derivative(QPTime, PLt, Rst, Srt, Rnst, DT, YWRK, Y, YDOT1)      
    do i=1,NEQ
      YT(i) = Y(i) + 0.5*DT*YDOT1(i)
    enddo 
    
    QPTime = TSAVED + 0.5*DT
    
    call derivative(QPTime, PLt, Rst, Srt, Rnst, DT, YWRK, YT, YDOT2)   
    do i = 1,NEQ
      YT(i) = Y(i) + 0.5*DT*YDOT2(i)
    enddo   
    
    QPTime = TSAVED + 0.5*DT
    
    call derivative(QPTime, PLt, Rst, Srt, Rnst, DT, YWRK, YT, YDOT3)         
    do i = 1,NEQ
      YT(i) = Y(i) + DT*YDOT3(i)
    enddo
    
    QPTime = TSAVED + DT
    
    call derivative(QPTime, PLt, Rst, Srt, Rnst, DT, YWRK, YT, YDOT4)             
    do i = 1,NEQ
      Y(i) = Y(i) + DT * (YDOT1(i) + 2.*YDOT2(i) + 2.*YDOT3(i) + YDOT4(i)) / 6.
    enddo
           
    if (Y(3) < 0) Y(3) = D0
        
    if (MMC == 1) then
      Y(1) = D0
    else    
      if (Y(2) < Mmin) Y(2) = Mmin
      
      if (abs(YWRK(2)) > nullerror) then
        Y(1) = YWRK(1) * ((Y(2) / YWRK(2))**(1./x))      
      else
        Y(1) = Dmin
      endIf
      !
      if (Y(1) > D0)   Y(1) = D0                                 
      if (Y(1) < Dmin) Y(1) = Dmin
    endif
    
    if (MMC == 1) then
      alpha = D0
    else
      if (sshape == 1) then      
        x = 3.0D0   ! sphere
        if (abs(Y(1)) > nullerror) then
          alpha = 6.0D0/(Rous*Y(1)*1000000.D0)
        else
          alpha = D0
        endif
      else       
        x = 2.0D0   ! cylinder sl = m
        if (abs(Y(1)) > nullerror) then
          alpha = 1.0D0 / (Rous*1000000.D0)*(2.0D0/sl + 4.0D0/Y(1))
        else
          alpha = D0
        endif
      endif  
    endif
    
    if (MMC == 1) then
      Fdis = PLt
    else
      Fdis = RItmp * alpha * Cs * Y(2)
    endif  
    if (dataset == 2) Fdis = Fdis * Zb / Zmax   
    
    if (Esolid == 1) then
      if (dataset /= 2) then
       Fes = Y(2) * Etmp / Zb    
      elseif (dataset == 2) then
       Fes = Y(2) * Etmp / Zmax
      endif
    else
      Fes = 0.0     
    endif  
    
    if (dataset == 0) then
      if (Rotmp > 0) then              
        Fr = de * (1.0D0 - exp(-pk)) * PNItmp * PA * Y(3)
      else
        Fr = 0.0
      endif
    else
      if (Rotmp > 0) then
        Fr = Rer * PA * Y(3)
      else
        Fr = 0.0  
      endif
    endif
             
    Fe     = Etmp * PA * Y(3)
    Fl     = qwtmp * Fdptmp / thetwtmp * PA * Y(3)
    Fdecay = (PKl*Fdptmp + PKa*Fpptmp) * (PA*Zb) * Y(3)
    Fvol   = PKv * Faptmp * PA * Y(3)
    Fiw    = Fiftmp * qwtmp * Fdptmp / thetwtmp * PA * Y(3)
    Fgw    = Fl - Fiw
    Cdp    = Y(3) / (Rtmp*thetwtmp)
    
    if (dataset == 2) then   
      Zcheck = Zmin + 0.001
      if (Zb < Zcheck) then
        Fdecay = 0.D0
        Fvol = 0.D0
      endif
      CTs    = (Y(2) + PA*Zmax*Y(3)) / (Roub*PA*Zmax)
    elseif (dataset /= 2) then
      CTs    = (Y(2) + PA*Zb*Y(3)) / (Roub*PA*Zb)
    endif 
            
    if (MMC == 1) then
      Fprep = 0.D0
    else
      if (Cdp > Cs) then
        Fprep = 1000.0D0 * PA * Zb * thetwtmp * (Cdp - Cs) / DT
      else
        Fprep = 0.D0
      endif
    endif
    
    if (dataset == 0) then  
      if (BMPflag == 1) then
        Qro  = Rotmp*PA / PNItmp                              !Qro is runoff per rainfall day, m3/day
        Qinf = qwtmp*PA / 365.25                              !Qinf = infiltration flow per day, m3/day
        Fro  = (Fes + Fe + Fr) / PNItmp    
        if ( Qro > 0.D0) then
          TSSi = 1000000.D0 * Roub * PA * Etmp / (Rotmp*PA)   !TSSi concentration (mg/L)
        else
          TSSi = 0.D0
        endif
      endif
      !
    else             
      if (BMPflag == 1) then  
        Qro  = Rotmp*PA                                       !runoff flow Qro(m^3/d)
        Qinf = qwtmp*PA                                       !infilt. Qinf(m^3/d)
        Fro  = Fes + Fe + Fr    
        if ( Qro > 0.D0) then
          TSSi = 1000000.D0 * Roub * PA * Etmp / (Rotmp*PA) 
        else
          TSSi = 0.D0
        endif
      endif 
    endif
        
    if (dataset /= 0) then
      if (QPTime == DT) then     
        TI22 = QPTime       
      else
        if (QPTime >= (TI22+1.D0)) then
          TI22 = QPTime     
          if (dataset /= 2) then   
            write(21,'(F14.4, 17(ES11.2))'), QPTime/365.25, Y(1),Y(2),PA*Zb*Y(3),Y(3), &
                                             Cdp,CTs,Fdis,Fes,Fe,Fl,Fr,Fdecay,Fvol,Fgw,Fes+Fe,Fr+Fiw,Fprep
          elseif (dataset == 2) then
            write(21,'(F14.4, 17(ES11.2))'), QPTime/365.25, Y(1),Y(2),PA*Zmax*Y(3),Y(3), &
                                             Cdp,CTs,Fdis,Fes,Fe,Fl,Fr,Fdecay,Fvol,Fgw,Fes+Fe,Fr+Fiw,Fprep  
          endif
          if (BMPflag == 1) write(26,'(F14.4, 5(ES17.2))'), QPTime/365.25, Qro,Qinf,Fro,Fl,TSSi          
          !
        elseif (QPTime >= (TF-DTMIN)) then
          if (dataset /= 2) then   
            write(21,'(F14.4, 17(ES11.2))'), QPTime/365.25, Y(1),Y(2),PA*Zb*Y(3),Y(3), &
                                             Cdp,CTs,Fdis,Fes,Fe,Fl,Fr,Fdecay,Fvol,Fgw,Fes+Fe,Fr+Fiw,Fprep
          elseif (dataset == 2) then
            write(21,'(F14.4, 17(ES11.2))'), QPTime/365.25, Y(1),Y(2),PA*Zmax*Y(3),Y(3), &
                                             Cdp,CTs,Fdis,Fes,Fe,Fl,Fr,Fdecay,Fvol,Fgw,Fes+Fe,Fr+Fiw,Fprep  
          endif
          if (BMPflag == 1) write(26,'(F14.4, 5(ES17.2))'), QPTime/365.25, Qro,Qinf,Fro,Fl,TSSi        
        endif
      endif
    endif
    
    if (dataset == 0) then
      outQPTime = QPTime    
      if (QPTime == DT) then   
        TI33 = QPTime         
        write(21,'(F14.4, 17(ES11.2))'), outQPTime, Y(1),Y(2),PA*Zb*Y(3),Y(3), &
                                         Cdp,CTs,Fdis,Fes,Fe,Fl,Fr,Fdecay,Fvol,Fgw,Fes+Fe,Fr+Fiw,Fprep
        if (BMPflag == 1) write(26,'(F14.4, 5(ES17.2))'), QPTime/365.25, Qro,Qinf,Fro,Fl,TSSi   
        !
      elseif (QPTime > DT) then      
        if (QPTime .ge. (TI33+0.1D0)) then      
          TI33 = QPTime     
          write(21,'(F14.4, 17(ES11.2))'), outQPTime, Y(1),Y(2),PA*Zb*Y(3),Y(3), &
                                           Cdp,CTs,Fdis,Fes,Fe,Fl,Fr,Fdecay,Fvol,Fgw,Fes+Fe,Fr+Fiw,Fprep
          if (BMPflag == 1) write(26,'(F14.4, 5(ES17.2))'), QPTime/365.25, Qro,Qinf,Fro,Fl,TSSi   
        else
          if (QPTime >= (TF-DTMIN)) then
            write(21,'(F14.4, 17(ES11.2))'), outQPTime, Y(1),Y(2),PA*Zb*Y(3),Y(3), &
                                             Cdp,CTs,Fdis,Fes,Fe,Fl,Fr,Fdecay,Fvol,Fgw,Fes+Fe,Fr+Fiw,Fprep
            if (BMPflag == 1) write(26,'(F14.4, 5(ES17.2))'), QPTime/365.25, Qro,Qinf,Fro,Fl,TSSi   
          endif
        endif
      endif
    endif
        
end subroutine