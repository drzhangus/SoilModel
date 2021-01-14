!=========================================================================================================================== 
  subroutine StrReplace(InStr,OldChar,NewChar,OutStr)
    implicit none
    !
    character(len = *) , Intent( IN ) :: InStr
    character(len = *) , Intent( IN ) :: OldChar
    character(len = LEN(OldChar)) , Intent( IN )  :: NewChar
    character(len = LEN(InStr)) , Intent( INOUT ) :: OutStr
    integer :: i  ! loop variable
    !
    OutStr=InStr
    i=INDEX(OutStr,OldChar)
    do while(i>0)
	    OutStr(i:i+LEN(OldChar)-1)=NewChar
	    i=INDEX(OutStr,OldChar)
    end do
  end subroutine
  
!===========================================================================================================================     
  subroutine StringSplit(InStr,delimiter,StrArray,nsize)
    implicit none
    !
    character(len = *) ,  Intent( IN ) :: InStr
    character(len = *)  , Intent( IN ) :: delimiter
    character(len = LEN(InStr)),dimension(LEN(InStr)),Intent( OUT ) :: StrArray
    integer, Intent( OUT ) :: nsize ! effective size of StrArray
    integer:: i,j    ! loop variable
    integer:: istart ! split index for Start Position
    !
    nsize=0
    istart=1
    do i=1,LEN(InStr)
	    do j=1,LEN(delimiter)
		    if (InStr(i:i) == delimiter(j:j)) then
			    if (istart == i) then
			      istart=i+1 
			    endif
			    if (istart<i) then
				    nsize=nsize+1
				    StrArray(nsize)=InStr(istart:i-1)
				    istart=i+1
			    endif
		    endif
	    enddo
    enddo
    !
    if (nsize>0) then
	    if (istart<LEN(InStr)) then
		    nsize=nsize+1
		    StrArray(nsize)=InStr(istart:LEN(InStr))
	    endif
    endif
    ! 
    if ( (nsize<1) .AND. (LEN(TRIM(InStr)) > 0 )) then
		  nsize=1
		  StrArray(1)=InStr
    endif
  end subroutine
 
!===========================================================================================================================     
  real*8 function sgn(d)
    implicit none
    real*8::d
    !
    if (d > 0.D0) then
      sgn = 1.
    elseif (d == 0.D0) then
      sgn = 0.
    elseif (d < 0.D0) then
      sgn = -1.
    else
      write(*,*) "false"
    endif
  end function