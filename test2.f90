program slash_replace
  implicit none
  character(len=20) :: str = "Hello/World"
  print *, "Original string:", str
  print *, "Replaced string:", replace_slash(str)
contains
  function replace_slash(s) result(r)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: r
    integer :: i
    do i = 1, len(s)
      if (s(i:i) == "/") then
        r(i:i) = "\"
      else
        r(i:i) = s(i:i)
      end if
    end do
  end function replace_slash
end program slash_replace
