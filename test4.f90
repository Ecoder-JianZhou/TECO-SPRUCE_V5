program DynamicFormat
    implicit none
    integer, parameter :: n = 5
    integer :: i
    ! 显式声明格式字符串
    character(len=50) :: format_string

    ! 定义整型变量
    integer :: data(n)
    data = [1, 2, 3, 4, 5]
  
    
  
    ! 构建格式字符串
    format_string = '(' // repeat('(i5, ",")', n-1) // 'i5)'
    print*, format_string
  
    ! 使用动态格式写入数据
    write(*, format_string) (data(i), i=1,n)
  
  end program DynamicFormat
  