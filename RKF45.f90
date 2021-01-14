! Runge-Kutta-Fehlberg 4(5) subroutine 
! Fourth-order solution with fifth-order error control (RK45)
!===========================================================================================================================         
  subroutine RKF45(NEQ, TI, TFS, H, TETOL, Y)
    !
    use Global 
    use HydroVariable 
    implicit none
    !
    integer, Intent( INout )      :: NEQ
    real(kind=8), Intent( INout ) :: TI, TFS, H, TETOL, Y(*)
    integer ICOEFF, kcount, i, j, k, kk, IFLAG, IFLAG2, L, ii 
    !
    real(kind=8) :: SDT, DT, DT1, EDS, TWRK, DTopt  
    real(kind=8) :: ALPH(6), CH(6), CE(5), betaA(6,5)
    real(kind=8) :: YDOT(NEQ), YWRK(NEQ), F(NEQ, 6), Z(NEQ), TOL(NEQ), ER(NEQ), Yfirst(NEQ)
    real(kind=8) :: PLt, Rst, Srt, Rnst, Temp, Cpore, DT2x, DTnext, TI1, Zcheck   
    real(kind=8), save   :: TI2, outTI     

    ICOEFF = 0
    if ( ICOEFF == 0) then
      ! dimension and define integration coefficients          
      ALPH(2) = real(1./4.,  8)
      ALPH(3) = real(3./8.,  8)
      ALPH(4) = real(12./13.,8)
      ALPH(5) = 1.D0
      ALPH(6) = 0.5D0
      
      betaA(2, 1) = real(1./4., 8)
      betaA(3, 1) = real(3./32.,8)
      betaA(3, 2) = real(9./32.,8)
      betaA(4, 1) = real(1932./2197., 8)
      betaA(4, 2) = real(-7200./2197.,8)
      betaA(4, 3) = real(7296./2197., 8)
      betaA(5, 1) = real(439./216.,   8)
      betaA(5, 2) = -8.D0
      betaA(5, 3) = real(3680./513., 8)
      betaA(5, 4) = real(-845./4104.,8)
      
      betaA(6, 1) = real(-8./27.,8)
      betaA(6, 2) = 2.D0
      betaA(6, 3) = real(-3544./2565.,8)
      betaA(6, 4) = real(1859./4104., 8)
      betaA(6, 5) = real(-11./40.,    8)

      CH(1) = real(16./135.,8)
      CH(2) = 0.D0
      CH(3) = real(6656./12825., 8)
      CH(4) = real(28561./56430.,8)
      CH(5) = real(-9./50.,8)
      CH(6) = real(2./55., 8)

      CE(1) = real(25./216.,8)
      CE(2) = 0.D0
      CE(3) = real(1408./2565.,8)
      CE(4) = real(2197./4101.,8)
      CE(5) = real(-1./5.,8)
       
      ICOEFF = 1       
    endif
    SDT = abs((TFS-TI)/(TFS-TI))      
    DT  = Abs(H) * SDT
    DT1 = DT
    EDS = TETOL/20.D0
    kcount = 1
    IFLAG  = 2      ! = 2 not converged, = 1 converged
    IFLAG2 = 1
   
    TWRK = TI
    do i = 1,NEQ
      YWRK(i) = Y(i)
    enddo

    do     
      if (IFLAG == 1) EXIT   
                 
      if (dataset == 1 .or. dataset == 2) then              
        if (abs(TI-10.D0) < nullerror) write(25,*) TI, Fracyearhourly(Ihours), Fracyearhourly(Ihours + 1)
        
        do while (.not.((TI >= Fracyeardaily(Idays)) .AND. (TI < Fracyeardaily(Idays + 1))))
          Idays = Idays + 1
        enddo
        !
        do while (.not.((TI >= Fracyearhourly(Ihours)) .AND. (TI < Fracyearhourly(Ihours + 1))))
          Ihours = Ihours + 1
        enddo 
        !
        do i=1,24
          PI24tmp(i) = PI(((Idays-1)*24) + i)
        enddo  
        
        if (dataset == 2) then
          Zb = thickness(Idays)  
        endif
        
        PItmp    = PI(Ihours)
        RItmp    = RI(Idays)
        Rotmp    = Ro(Idays)       
        Etmp     = E(Idays)        
        qwtmp    = qw(Idays)       
        !PNItmp   = PNI(Idays)        
        PIttmp   = PIt(Idays)      
        thetwtmp = thetw(Idays)      
        Fiftmp   = Fif(Idays)      
      endif
         
      if (dataset == 0) then
        !PItmp    = PI(1)
        RItmp    = RI(1)   
        Rotmp    = Ro(1)
        Etmp     = E(1)
        qwtmp    = qw(1)
        PNItmp   = PNI(1)        
        PIttmp   = PIt(1) 
        thetwtmp = thetw(1)       
        Fiftmp   = Fif(1)     
      endif
          
      Rtmp   = R(Idays)
      Fdptmp = Fdp(Idays)
      Fpptmp = Fpp(Idays)
      Faptmp = Fap(Idays)
      
      ! set loadings
      if (dataset == 0) then 
        PresTime = TI
      else      
        PresTime = TI/365.25D0 
      endif
           
      if (PresTime < yr(1)) PLt = 0.             
      if ((PresTime >= yr(1)) .AND. (PresTime < yr(nyrs))) then
        do i = 1,nyrs-1
          if ( (PresTime >= yr(i)) .AND. (PresTime < yr(i + 1))) PLt = PL(i)
        enddo       
      endif    
      if (PresTime >= yr(nyrs)) PLt = PL(nyrs)
        
      if (BMPflag == 1) then
        if (PresTime < BMPyr(1)) then
          Rst  = 0.D0
          Srt  = 0.D0
          Rnst = 0.D0
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
        Rst  = 0.D0
        Srt  = 0.D0
        Rnst = 0.D0
      endif
      
      ! evaluate system of differential equations
      call derivative(TI, PLt, Rst, Srt, Rnst, DT, YWRK, Y, YDOT)    
      do i = 1,NEQ
        F(i,1) = YDOT(i)
      enddo
           
      do k = 2,6
        kk = k-1
        do i=1,NEQ
          Temp = 0.D0
          do j=1,kk
            Temp = Temp + betaA(k,j)*F(i,j)
          enddo
          Y(i) = YWRK(i) + DT*Temp
        enddo 
        
        TI = TWRK + ALPH(k)*DT       
        call derivative(TI, PLt, Rst, Srt, Rnst, DT, YWRK, Y, YDOT)
        do j = 1,NEQ
          F(j,k) = YDOT(j)
        enddo
      enddo
      !
      ! Compute 5th order estimate, Z(n+1)
      do i = 1,NEQ
        Temp = 0.D0
        do L = 1,6
          Temp = Temp + CH(L)*F(i,L)
        enddo
        Z(i) = YWRK(i) + DT*Temp
      enddo
         
      ! Compute 4th order estimate, Y(n+1)
      do i = 1,NEQ
        Temp = 0.D0
        do L = 1,5
          Temp = Temp + CE(L)*F(i,L)
        enddo
        Y(i) = YWRK(i) + DT*Temp
      enddo
      IFLAG = 1
      DT1   = DT
         
      do i = 1,NEQ
        ER(i) = abs(Z(i)-Y(i))
        if (kcount==1) Yfirst(i) = Y(i)
        TOL(i) = 0.001D0 * abs(YWRK(i)-Yfirst(i)) + (1.0D-20)         
        if (TOL(i) < TETOL) TOL(i) = TETOL
        if (ER(i) > TOL(i)) then
          IFLAG = 2
          DT    = 0.5D0 * DT1
          DTopt = DT1 * 0.84D0 * (TOL(i)*DT1/(ER(i) + 0.0000000001D0))**(1.D0/4.D0)
          if (DTopt < DT) DT = DTopt
        endif
      enddo
      
      if (IFLAG2 < 2) then
        Cpore = Y(3) / (thetwtmp*Rtmp)
        if (Cpore > Cs) then
          IFLAG  = 2
          DT     = DTMIN
          IFLAG2 = 2
        endif
      endif
    
      if (kcount > 4) IFLAG = 1
  
      if (IFLAG == 2) then
        kcount = kcount + 1
        TI = TWRK
        do i = 1, NEQ
          Y(i) = YWRK(i)
        enddo
              
        if (DT < DTMIN) DT = DTMIN        
        if (DT > DTMAX) DT = DTMAX              
        if (Y(3) < 0.D0) Y(3) = 0.D0    
      else
        TI     = TWRK + DT1
        DT2x   = 2.0D0*DT1
        DTnext = DT2x
        
        if (DTnext < DTMIN) DTnext = DTMIN        
        if (DTnext > DTMAX) DTnext = DTMAX          
        if (abs(DTnext) > abs(TFS-TI)) DTnext = TFS-TI   
        H   = DTnext
        TI1 = TWRK
        if ( Y(3) < 0.D0 ) Y(3) = 0.D0  
            
        if (MMC == 1) then
          Y(1) = 0.D0
        else
          if (Y(2) < Mmin) Y(2) = Mmin               
          if (abs(YWRK(2)) > nullerror) then    
            Y(1) = YWRK(1) * ((Y(2)/YWRK(2))**(1.0D0/x))
          else
            Y(1) = Dmin
          endif
               
          if ( Y(1) > D0 )  Y(1) = D0           
          if ( Y(1) < Dmin) Y(1) = Dmin  
        endif
           
        if (MMC == 1) then
          alpha = 0.D0
        else
          if (sshape == 1) then
            x = 3.0D0
            if (abs(Y(1)) > nullerror) then   
              alpha = 6.D0/(Rous*Y(1)*1000000.D0)
            else
              alpha = 0.D0
            endif
          else
            x = 2.0D0
            if (abs(Y(1)) > nullerror) then    
              alpha = 1.D0/(Rous*1000000.D0)*(2.D0/sl+4.D0/Y(1))
            else
              alpha = 0.D0
            endif
          endif 
          !
        endif   
     
        if (MMC == 1) then
          Fdis = PLt
        else
          Fdis = RItmp * alpha * Cs * Y(2)
        end if  
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
        
        if (dataset == 0) then
          if (Rotmp > 0) then
            Fr = de * (1.D0 - exp(-pk)) * PNItmp * PA * Y(3)
          else
            Fr = 0.D0
          end if
          !
        else    
         if ( Rotmp > 0) then
            Fr = Rer * PA*Y(3)
         else
            Fr = 0.D0   
          endif
        end if
           
        Fe     = Etmp * PA * Y(3)
        Fl     = qwtmp * Fdptmp / thetwtmp * PA * Y(3)
        Fdecay = (PKl * Fdptmp + PKa * Fpptmp) * (PA * Zb) * Y(3)
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
          CTs    = (Y(2) + PA * Zb * Y(3))/(Roub * PA * Zb)
        endif
            
        if (MMC == 1) then
          Fprep = 0.D0
        else
          if (Cdp > Cs) then
            Fprep = 1000.D0 * PA * Zb * thetwtmp * (Cdp - Cs)/DT1
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
            Qro  = Rotmp*PA  
            Qinf = qwtmp*PA   
            Fro  = Fes + Fe + Fr  
            if ( Qro > 0.D0) then
              TSSi = 1000000.D0 * Roub * PA * Etmp / (Rotmp*PA)
            else
              TSSi = 0.D0
            endif
          endif 
        endif
        
        NoTime = NoTime + 1
               
        if (dataset /= 0) then
          outTI = TI 
          write(100,*) TI, TI2, TI/365.25
          if (NoTime == 1) then
            TI2 = TI
            if (dataset /= 2) then       
              write(21,'(F14.4, 17(ES11.2))'), outTI/365.25, Y(1),Y(2),PA*Zb*Y(3),Y(3), &
                                               Cdp,CTs,Fdis,Fes,Fe,Fl,Fr,Fdecay,Fvol,Fgw,Fes+Fe,Fr+Fiw,Fprep
            elseif (dataset == 2) then
              write(21,'(F14.4, 17(ES11.2))'), outTI/365.25, Y(1),Y(2),PA*Zmax*Y(3),Y(3), &
                                               Cdp,CTs,Fdis,Fes,Fe,Fl,Fr,Fdecay,Fvol,Fgw,Fes+Fe,Fr+Fiw,Fprep  
            endif
            if (BMPflag == 1) write(26,'(F14.4, 5(ES17.2))'), outTI/365.25, Qro,Qinf,Fro,Fl,TSSi            
          elseif (Notime > 1) then
              if (TI .ge. (TI2+1.0) ) then   
                TI2 = TI
                if (dataset /= 2) then   
                  write(21,'(F14.4, 17(ES11.2))'), outTI/365.25, Y(1),Y(2),PA*Zb*Y(3),Y(3), &
                                                   Cdp,CTs,Fdis,Fes,Fe,Fl,Fr,Fdecay,Fvol,Fgw,Fes+Fe,Fr+Fiw,Fprep
                elseif (dataset == 2) then
                  write(21,'(F14.4, 17(ES11.2))'), outTI/365.25, Y(1),Y(2),PA*Zmax*Y(3),Y(3), &
                                                   Cdp,CTs,Fdis,Fes,Fe,Fl,Fr,Fdecay,Fvol,Fgw,Fes+Fe,Fr+Fiw,Fprep  
                endif
                if (BMPflag == 1) write(26,'(F14.4, 5(ES17.2))'), outTI/365.25, Qro,Qinf,Fro,Fl,TSSi
                !
              elseif (TI >= (TFS-DTMIN)) then
                if (dataset /= 2) then   
                  write(21,'(F14.4, 17(ES11.2))'), outTI/365.25, Y(1),Y(2),PA*Zb*Y(3),Y(3), &
                                                   Cdp,CTs,Fdis,Fes,Fe,Fl,Fr,Fdecay,Fvol,Fgw,Fes+Fe,Fr+Fiw,Fprep
                elseif (dataset == 2) then
                  write(21,'(F14.4, 17(ES11.2))'), outTI/365.25, Y(1),Y(2),PA*Zmax*Y(3),Y(3), &
                                                   Cdp,CTs,Fdis,Fes,Fe,Fl,Fr,Fdecay,Fvol,Fgw,Fes+Fe,Fr+Fiw,Fprep  
                endif
                if (BMPflag == 1) write(26,'(F14.4, 5(ES17.2))'), outTI/365.25, Qro,Qinf,Fro,Fl,TSSi
              endif
          endif
        endif
               
        if (dataset == 0) then
          outTI = TI 
          if (NoTime == 1) then
            TI2 = TI
            write(21,'(F14.4, 17(ES11.2))'), outTI, Y(1),Y(2),PA*Zb*Y(3),Y(3), &
                                             Cdp,CTs,Fdis,Fes,Fe,Fl,Fr,Fdecay,Fvol,Fgw,Fes+Fe,Fr+Fiw,Fprep
            if (BMPflag == 1) write(26,'(F14.4, 5(ES17.2))'), outTI/365.25, Qro,Qinf,Fro,Fl,TSSi
            !
          else
            if (Notime > 1) then
              if (TI >= (TI2+0.1D0)) then
                TI2 = TI
                write(21,'(F14.4, 17(ES11.2))'), outTI, Y(1),Y(2),PA*Zb*Y(3),Y(3), &
                                                 Cdp,CTs,Fdis,Fes,Fe,Fl,Fr,Fdecay,Fvol,Fgw,Fes+Fe,Fr+Fiw,Fprep
                if (BMPflag == 1) write(26,'(F14.4, 5(ES17.2))'), outTI/365.25, Qro,Qinf,Fro,Fl,TSSi
              else
                if (TI >= (TFS-DTMIN)) then
                write(21,'(F14.4, 17(ES11.2))'), outTI, Y(1),Y(2),PA*Zb*Y(3),Y(3), &
                                                 Cdp,CTs,Fdis,Fes,Fe,Fl,Fr,Fdecay,Fvol,Fgw,Fes+Fe,Fr+Fiw,Fprep
                if (BMPflag == 1) write(26,'(F14.4, 5(ES17.2))'), outTI/365.25, Qro,Qinf,Fro,Fl,TSSi
                endif
              endif
            endif
          endif
        endif 
        !
      endif     
    enddo
     
  end subroutine