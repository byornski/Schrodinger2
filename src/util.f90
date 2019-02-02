module util

contains


  
  
  subroutine assert(cond,msg)
    !Checks the assertion is true or STOPs and gives error message
    logical, intent(in) :: cond
    character(len=*), intent(in) :: msg
    if (.not. cond) then
       write(*,*) msg
       STOP 'Failed Assertion'
    end if
  end subroutine assert


  

end module util
