! SSCFM (Surface Soil Contaminant Fate Model) Version 2.0
! simulates user-defined multi-phase constituents (solid, aqueous, gas, soil sorded) in the AOI
! developed by Zhonglong Zhang and Mark Dortch
!
! This verion was developed from three verions of VB6 codes 
! 1) AA verion reads average annual hydrology inputs and the computation time step is based on fraction of "year"
! 2) TV version reads daily hydrology and hourly rainfall inputs and the computation time step is based on fraction of "day"
!
! Updated January 2021
!===========================================================================================================================  
  program Main
    use Global 
    implicit none
    
    integer iPosition, TFF, rk45, i, j, k, kk, flag, firstyear, nyrs1
    integer NEQ, NumHourPerDay	   
    integer istartyear, iyear, imonth, ihour, jhours	
    integer(kind=4) :: ntotdays, itotdays, iday, jdays
    logical SkipLoop
    !
    real(kind=8) :: H, HH, DT, TI, TETOL         
    real(kind=8) :: bmptemp2, bmptemp3, valtemp, loading 
    real(kind=8), allocatable :: Y(:), YE(:), yr1(:), PL1(:)
    !
    character(len=255) :: Temp, xTemp, Temp1, Temp2, Temp3    
    character(len=255) :: inputfile, BMPfile, cmd, arg
    character(len=255) :: astrSplitFiles(5)
    character(len=255), allocatable :: sArgs(:)       	 
       
    NEQ = 3
    NumHourPerDay = 24
     
    Zmin = 0.0508     
    
    allocate(Y(NEQ))
    allocate(YE(NEQ))
 
    i = command_argument_count()              
    if (i == 0) STOP 'Enter the inp file'  
    allocate(sArgs(i))
       
    i = 1
    do
      call get_command_argument(i, arg)
      if (LEN_TRIM(arg) == 0) EXIT          
      sArgs(i) = TRIM(arg)
      i = i+1
    enddo
       
    inputfile = sArgs(1)
   
    if (INDEX(inputfile,'~senso') > 0 ) then
      call  StrReplace(inputfile,'~senso','ofile',inputfile)
    endif
    
    if (LEN_TRIM(inputfile) /= 0 ) then
      open(22, file=inputfile)
      print *, ' Input file is :',trim(inputfile)
    else
      STOP 'NO inp file'
		endif

    !---------------------------------------------------------------------------------------------------  
    !21 inputfile.out 
    !22 inputfile.inp 
    !23 dailyRainfall.out
    !24 hourlyRainfall.out    
    !20 ofileDaily.inp
		!
	  !Read header
	  SkipLoop = .FALSE.
	  do while(.NOT. SkipLoop)
		read(22,'(a)') Temp
		if(index(Temp, "#") == 0) SkipLoop = .TRUE.
	  end do
	  backspace(22)    
       
    read(22,*)  dataset
    write(*,10) 'Hydro input', dataset  
            
    read(22,*) istartyear
    write(*,10) 'Start year=', istartyear   
    read(22,*) TF     ! total simulation duration (years for AA or days for TV)
    TFF = floor(TF)     
        
    read(22,*) H
    write(*,9) 'H=', H 
    if (abs(H) < nullerror)  STOP 'Time step is 0'  
    if (dataset /= 0 ) then     
      if (H < 0.036525) H = 0.036525D0
    else
      if (H < 0.0001) H = 0.0001D0
    endif
      
    HH    = H
    DTMIN = HH
    DTMAX = 1.D0
    DT    = H
           
    read(22,*) rk45   !! (1) for RKF45, (0) for RK4                            
    write(*,10) 'rk45=', rk45       
    if ((rk45 == 1) .AND. (H < nullerror)) then
      H = 0.01D0
      print *, 'H has been changed to 0.01, due to rk45=1 and H < nullerror'
    endif
                    
    read(22,*) PA
    write(*,9) 'PA=', PA 
    
    if (dataset == 0 .or. dataset == 1)  then   
      read(22,*) Zb
      write(*,9) 'Zb=', Zb  
    elseif (dataset == 2) then     
      read(22,*) Zmax
      write(*,9) 'Zmax=', Zmax
    endif
    
    read(22,*) de
    write(*,9) 'de=', de  
    
    !read(22,*) thets
    !write(*,9) 'thets=', thets         
    !read(22,*) Roub
    !write(*,9) 'Roub=', Roub         
    !read(22,*) soila
    !write(*,9) 'soila=', soila 

    allocate(Fracyeardaily(TFF+2))
    allocate(R(TFF+2))
    allocate(Fdp(TFF+2))
    allocate(Fpp(TFF+2))
    allocate(Fap(TFF+2))          
    if (dataset == 1 .or. dataset == 2) then   
      read(22,*) thets
      write(*,9) 'thets=', thets         
      read(22,*) Roub
      write(*,9) 'Roub=', Roub         
      read(22,*) soila
      write(*,9) 'soila=', soila 
    
      Temp = inputfile       
      call StringSplit(Temp,'.',astrSplitFiles,i)    
      Temp1 = trim(astrSplitFiles(1))//'Daily.inp'  ! Daily thet, precip, Ro, E, gw, Fif                
      open(20, file = Temp1)
      read(20, '(a)') Temp1
    
      Temp2 = 'dailyRainfall.out'                   ! Daily rainfall (m/d)     
      Temp3 = 'hourlyRainfall.out'                  ! Hourly rainfall (m/hr)    
      open(23, file = Temp2)  
      read(23, '(a)') Temp2  
      open(24, file = Temp3)                 
      read(24, '(a)') Temp3                
                
      allocate (PI((TFF*NumHourPerDay) + 2))                
      allocate (Fracyearhourly((TFF*NumHourPerDay) + 2))
      allocate (thetw(TFF+2))	
      !allocate (PNI(TFF+2))     ! not used 
      allocate (PIt(TFF+2))   
      allocate (RI(TFF+2))
      allocate (Ro(TFF+2))   
      allocate (E(TFF+2))   
      allocate (qw(TFF+2)) 
      allocate (Fif(TFF+2))
      if (dataset == 2) allocate (thickness(TFF+2))   
      !
    elseif (dataset == 0) then                    
      i = 1 
      allocate(thetw(i))
      allocate(PNI(i))
      !allocate(PI(i))           ! not used 
      allocate(PIt(i))                     
      allocate(RI(i))     
      allocate(Ro(i))   
      allocate(E(i))   
      allocate(qw(i))
      allocate(Fif(i))      
      !     
      read(22,*) thetw(1)
      write(*,9) 'thetw=', thetw(1)
      
      read(22,*) thets
      write(*,9) 'thets=', thets         
      read(22,*) Roub
      write(*,9) 'Roub=', Roub         
      read(22,*) soila
      write(*,9) 'soila=', soila 
      
      read(22,*) PNI(1)
      write(*,9) 'PNI=', PNI(1)          
      read(22,*) RI(1)
      write(*,9) 'RI=', RI(1)             
      read(22,*) PIt(1) 
      write(*,9) 'PIt=', PIt(1)                
      read(22,*) RO(1) 
      write(*,9) 'RO=', RO(1)             
      read(22,*) E(1)
      write(*,9) 'E=', E(1)             
      read(22,*) qw(1)
      write(*,9) 'qw=', qw(1)                       
      read(22,*) Fif(1)
      write(*,9) 'Fif=', Fif(1)
                       
    else    
      write(*,10) 'Hydro input', dataset
      STOP 'Invalid input for rainfall dataset'                  
    endif
    
    read(22,*) BMPflag
    write(*,10) '  BMP flag=', BMPflag  
    !if (BMPflag == 1) then
    !endif
    
    !---------------------------------------------------------------------------------------------------  
    read(22,*)  NMC
    write(*,10) 'NMC=',NMC  
         
    if (dataset /=0 ) write(*,*) 'Reading the daily data for thetw,PIt,Ro,E,qw,Fif,RI, daily and hourly rainfall.' 

    if (dataset == 2) then    
      Temp = 'soillayerThickness.txt'
      open(19,file = Temp)
      read(19,'(a)') Temp
      read(19,'(a)') Temp
      read(19,'(a)') Temp
      read(19,'(a)') Temp
      read(19,'(a)') Temp
    endif    
    
    Temp = inputfile   
    call StringSplit(Temp,'.',astrSplitFiles,i) 
    Temp = trim(astrSplitFiles(1))//'.out'    
    open(21,file=Temp)
    
    !---------------------------------------------------------------------------------------------------  
    Idays  = 0
    Ihours = 0   
    do jdays = 1,TFF
      Idays = Idays+1
      !     
      if (dataset == 1 .or. dataset == 2) then     
        read(20,*) iyear,imonth,iday,thetw(Idays),PIt(Idays),Ro(Idays),E(Idays),qw(Idays),Fif(Idays)
        read(23,*) iyear,imonth,iday,RI(Idays) 
        if (dataset == 2)  then                   
          read(19,*) iyear, imonth, iday, thickness(Idays),flag
          if (thickness(Idays) .lt. Zmin) then
            thickness(Idays) = Zmin
          elseif (thickness(Idays) .gt. Zmax) then
            thickness(Idays) = Zmax
          endif
        endif
            
        Fracyeardaily(Idays) = Idays - 1           
        do jhours = 1,24
          Ihours = Ihours+1
          read(24,*) iyear,imonth,iday,ihour,PI(Ihours) 
          Fracyearhourly(Ihours) = Fracyeardaily(Idays) + ((jhours-1)/24.0D0)     
        enddo 
      endif    
    enddo
    
    if (dataset == 1 .or. dataset == 2) then   
      Fracyeardaily(Idays+1)   = Idays       
      Fracyearhourly(Ihours+1) = Idays
    endif
            
    write(*,'(A30)') 'Complete reading general inputs.'
    
    ntotdays = Idays
		
    !if (dataset == 0) then
      write(21,'(A14,14(A11),2(A13),A11)'), &
                'Time(yr)','di(m)','Ms(g)','Mns(g)','Ctt(g/m3)','Cdp(g/m3)','CTs(mg/kg)', &
                'Fdis(g/yr)','Fes(g/yr)','Fe(g/yr)','Fl(g/yr)','Fr(g/yr)','Fdcay(g/yr)', &
                'Fvol(g/yr)','Fgw(g/yr)','Fes+Fe(g/yr)','Fr+Fiw(g/yr)','Fprep(g/yr)'    
    !else
    !  write(21,'(A14,14(A11),2(A12),A11)'), &
    !            'Time(yr)','di(m)','Ms(g)','Mns(g)','Ctt(g/m3)','Cdp(g/m3)','CTs(mg/kg)', &
    !            'Fdis(g/d)','Fes(g/d)','Fe(g/d)','Fl(g/d)','Fr(g/d)','Fdcay(g/d)',       &
    !            'Fvol(g/d)','Fgw(g/d)','Fes+Fe(g/d)','Fr+Fiw(g/d)','Fprep(g/d)'      
    !endif		
		
    TI    = 0.D0
    TETOL = 0.0000000001D0
    Dmin  = 0.0000001D0
    Mmin  = 0.D0    

    !---------------------------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------------------------  
    do j = 1,NMC
      Ihours = 1
      Idays  = 1
          
      read(22,*)  MC
      write(*,10) 'MC=',MC
      !        
      read(22,*)  MMC
      write(*,10) 'MMC=', MMC
      
      ! Miscible MC 
      if (MMC == 1) then
        Esolid = 2
        sshape = 1
        Y(1) = 0.D0
        Y(2) = 0.D0            
      else
        read(22,*) Esolid  
        read(22,*) sshape       
        if (sshape == 2 ) read(22,*) sl
        
        read(22,*) Y(1)
        D0 = Y(1)
        
        read(22,*) Y(2)
        if (Y(2) < 0.000001D0) Y(2) = 0.000001D0                
        M0 = Y(2)          
      endif
      
      read(22,*) Y(3)
            
      write(*,10) 'Esolid=',  Esolid
      write(*,10) 'sshape=',  sshape
      write(*,9)  'di/Y(1)=', Y(1)
      write(*,9)  'Ms/Y(2)=', Y(2)
      write(*,9)  'Ctt/Y(3)=',Y(3)           
      if (( MMC /=1) .AND. (sshape ==2)) write(*,9) 'sl=',sl
                 
      read(22,*) Rous
      write(*,9) 'Rous=', Rous    
      read(22,*) PKH
      write(*,9) 'PKH=', PKH 
      read(22,*) PKd
      write(*,9) 'PKd=', PKd 
      read(22,*) PKl
      write(*,9) 'PKl=', PKl 
      read(22,*) PKa
      write(*,9) 'PKa=', PKa 
      read(22,*) PKv
      write(*,9) 'PKv=', PKv 
      read(22,*) Cs
      write(*,9) 'Cs=', Cs 
      
      !---------------------------------------------------------------------------------------------------  
      read(22,*) nyrs1
      write(*,10) 'nyrs=', nyrs1     
      read (22, *) firstyear,loading
      write(*,10) 'Beginning year of MC loading', firstyear
      if(firstyear < istartyear) STOP 'Beginning year of loading is smaller than the start year of simulation.'
      nyrs = nyrs1 + firstyear-istartyear      
      
      if (allocated(yr1)) deallocate(yr1, yr, PL, PL1)
      allocate (yr1(nyrs1),yr(nyrs),PL(nyrs),PL1(nyrs1))      
      if (firstyear == istartyear) then
        yr1(1) = firstyear
        PL(1) = loading
        if(dataset .ne. 0) PL(1) = PL(1)/365.25D0
        yr(1) = 0.0
      elseif (firstyear > istartyear) then
        do k = 1, firstyear-istartyear
          yr(k) = k-1
          PL(k) = 0.0
        enddo 
        yr(1+firstyear-istartyear) = firstyear-istartyear
        PL(1+firstyear-istartyear) = loading
        if(dataset .ne. 0) PL(1 + firstyear-istartyear) = PL(1 + firstyear-istartyear)/365.25D0
      endif
      !
      do i = 2,nyrs1
        read(22,*) yr1(i), PL1(i)
        if (firstyear == istartyear) then
          yr(i) = yr1(i) - firstyear
          PL(i) = PL1(i)
          if(dataset .ne. 0) PL(i) = PL(i)/365.25D0
        elseif (firstyear > istartyear) then
          kk = i + firstyear-istartyear
          yr(kk) = yr1(i) - istartyear
          PL(kk) = PL1(i)
          if(dataset .ne. 0) PL(kk) = PL(kk)/365.25D0
        endif   
      enddo 
      !---------------------------------------------------------------------------------------------------  
          
      write(*,11) 'Complete reading MC', MC
      !---------------------------------------------------------------------------------------------------  
                
      do itotdays = 1,ntotdays                
        if (dataset == 0) then       
          R(itotdays)   = 1.D0 + ((thets-thetw(1))*PKH + Roub*PKd) / thetw(1)
          Fdp(itotdays) = 1.D0 / R(itotdays)
          Fpp(itotdays) = Roub*PKd / (R(itotdays)*thetw(1))
          Fap(itotdays) = (thets-thetw(1)) * PKH/(R(itotdays)*thetw(1))
        else
          R(itotdays)   = 1.D0 + ((thets-thetw(itotdays))*PKH + Roub*PKd) / thetw(itotdays)
          Fdp(itotdays) = 1.D0 / R(itotdays)
          Fpp(itotdays) = Roub*PKd / (R(itotdays)*thetw(itotdays))
          Fap(itotdays) = (thets-thetw(itotdays)) * PKH/(R(itotdays)*thetw(itotdays))          
        endif             
      enddo
     
      write(21,*)   
      write(21,20), 'MC constituent', MC  
      
      if (BMPflag == 1) write(26,20), 'MC constituent', MC  

      do while (abs(TF-TI) > DT) 
        valtemp = TF-TI  
        if (abs(valtemp) < DT) DT = valtemp*(valtemp/abs(valtemp)) 
        if (rk45 == 1) then  
          call RKF45(NEQ, TI, TF, H, TETOL, Y)
        else           
          call RK4(NEQ, TI, H, Y)
        endif 
      enddo
      
      TI = 0.D0    
      TF = real(TFF,8)            
      H  = HH
      NoTime = 0   
              
      print *, 'MC constituent', MC, 'Finished.'           
    enddo
    !---------------------------------------------------------------------------------------------------  
    !---------------------------------------------------------------------------------------------------  
    
    close(19) 
		close(20)
    close(21)   
    close(22)   
    close(23)   
    close(24)         
         
    deallocate(Y, YE, sArgs)      
    deallocate(thetw, PIt, RI, Ro, E, qw, Fif)    
    deallocate(R, Fdp, Fpp, Fap)
    if (allocated(PI))  deallocate(PI)
    if (allocated(PNI)) deallocate(PNI)
    if (allocated(Fracyeardaily))  deallocate(Fracyeardaily)
    if (allocated(Fracyearhourly)) deallocate(Fracyearhourly)
  
9   format(A15,G20.6)
10  format(A14,I10)
11  format(A19,I4)   
20  format(A15,I10)    
    
    print *, 'All Finished.'
    
end program 
    