! SSCFM (Surface Soil Contaminant Fate Model)
! Global variables
!===========================================================================================================================  
  module Global
    integer dataset, NMC, MMC, MC, Esolid, sshape, nyrs, RFyears
    integer BMPflag, BMPnyrs   
    integer(kind=4) :: NoTime, Ihours, Idays      
    !
    real(kind=8) :: nullerror = 1.D-10  
    real(kind=8) :: yr0, TF
    real(kind=8) :: D0, Dmin, DTMIN, DTMAX, M0, Mmin, sl, x
    real(kind=8) :: PA, Zb, de, Zmin, Zmax   
    real(kind=8) :: thets, Roub, soila
    real(kind=8) :: PKH, PKd, Rous
    real(kind=8) :: PKl, PKa, PKv
    real(kind=8) :: Tw, Cs0, Cs, CTs
    real(kind=8) :: Fdptmp, Fpptmp, Faptmp, Rtmp, Cdp
    real(kind=8) :: alpha, pk
    real(kind=8) :: Fdis, Fes, Fe, Fl, Fr, Fdecay, Fvol            
    !
    real(kind=8), allocatable :: Fsw, Fgw, Fiw, Fprep 
    real(kind=8), allocatable :: PNI(:), thetw(:), PI(:), PIt(:), RI(:), Ro(:), E(:), qw(:)    
    real(kind=8), allocatable :: Fdp(:), Fpp(:), Fap(:), R(:), Fif(:), thickness(:) 
    real(kind=8), allocatable :: yr(:), PL(:)
    real(kind=8), allocatable :: BMPyr(:), BMPyr1(:), Rs(:), Rns(:), Sr(:)   
    real(kind=8), allocatable :: Fracyeardaily(:), Fracyearhourly(:)  ! when yearly dataset used, Fracyeardaily is used for yearly
  end module  
  
!===========================================================================================================================     
  module HydroVariable    
    real(kind=8) :: PNItmp 
    real(kind=8) :: thetwtmp, PItmp, PIttmp, RItmp, Rotmp, Etmp, qwtmp, PI24tmp(24)
    real(kind=8) :: PresTime
    real(kind=8) :: Fiftmp, Qro, Qinf, Rer, Fro, TSSi   
  end module