program DynamicMemoryAllocation
    integer, dimension(:), allocatable :: my_array
  
    ! 调用子程序来分配内存并填充数组
    call allocateAndFillArray(my_array)
  
    ! 打印数组内容
    print *, my_array
  
    ! 释放内存
    deallocate(my_array)
  
  contains
  
    subroutine allocateAndFillArray(arr)
      integer, dimension(:), allocatable, intent(out) :: arr
      integer :: i
  
      allocate(arr(3))
  
      do i = 1, 3
        arr(i) = i
      end do
    end subroutine allocateAndFillArray
  
  end program DynamicMemoryAllocation
  