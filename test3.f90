program Example
    implicit none
  
    integer, parameter :: m = 10  ! 数组的大小
    real :: arr1(m), arr2(m), arr3(m), arr4(m), arr5(m)
    integer :: i
  
    ! 初始化数组，这里只是一个例子，实际情况应该根据需求进行赋值
    arr1 = 1.0
    arr2 = 2.0
    arr3 = 3.0
    arr4 = 4.0
    arr5 = 5.0
  
    ! 写入输出，将每个数组按行排列
    ! do i = 1, m
    !    write(*,*) arr1(i), arr2(i), arr3(i), arr4(i), arr5(i)
    ! end do
    write(*,*) (arr1(i), arr2(i), arr3(i), arr4(i), arr5(i), i = 1, m)

  end program Example
  