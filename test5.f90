program TypeArrayExample
    implicit none
  
    ! 定义一个包含type的数组
    type :: MyType
      real :: value
      integer :: count
    end type
  
    type(MyType), dimension(5) :: myArray
    integer :: i
  
    ! 初始化数组
    do i = 1, 5
      myArray(i)%value = 0.5 * i
      myArray(i)%count = i
    end do
  
    ! 遍历数组
    do i = 1, 5
      print *, "Element", i, ": Value =", myArray(i)%value, ", Count =", myArray(i)%count
    end do
  
  end program TypeArrayExample
  