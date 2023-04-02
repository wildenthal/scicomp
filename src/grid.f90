! Author: Gaston Barboza
module grid
      implicit none
      private
      
      integer, parameter :: K = kind(0.d0)
      integer, parameter :: file_number = 11
      real(K) :: xL, yL
      integer :: n
      public initialize_grid, save_grid, xpos, xint

      contains
      
      ! general grid operations
      subroutine initialize_grid(grid,lengthX,lengthY,dots)
            ! allocate grid of water heights and set it to 0
            real(K), allocatable :: grid(:,:)
            real(K) :: lengthX, lengthY
            integer :: i,j, dots
            n = dots
            allocate(grid(n,n))
            do j = 1, n
                  do i = 1, n
                        grid(i,j) = 0.0_K
                  end do
            end do
            xL = lengthX
            yL = lengthY
      end

      subroutine save_grid(grid)
            integer :: i,j
            real(K) :: axis_value(2)

            real(K), dimension(:,:) :: grid
            
            ! write height of water to grid.dat file
            open(unit=file_number, file='grid.dat', status='replace')
            do j = 1, n
                  do i = 1, n
                        write(file_number, '(f10.3)',advance='no') grid(i,j)
                        write(file_number, '(a)', advance='no') ' '
                  end do
                  write(file_number,*)
            end do
            close(file_number)
            
            ! write coordinates of x axis to file
            open(unit=file_number,file='x_axis.dat', status='replace')
            do i = 1, n
                  axis_value = xpos((/i,1/))
                  write(file_number, *) axis_value(1)
            end do
            close(file_number)
            
            ! write coordinates of y axis to file
            open(unit=file_number,file='y_axis.dat', status='replace')
            do j = 1, n
                  axis_value = xpos((/1,j/))
                  write(file_number, *) axis_value(2)
            end do
            close(file_number)
      end

      function xpos(integer_position)
            ! grid is indexed by integer numbers; this converts indices
            !   to coordinate positions
            integer, dimension(2) :: integer_position
            real(K), dimension(2) :: xpos
            
            xpos(1) = -xL/2.0_K + xL*(integer_position(1)-1)/(n-1)
            xpos(2) = -yL/2.0_K + yL*(integer_position(2)-1)/(n-1)
      end

      function xint(x)
            ! this converts coordinate positions to grid indices
            real(K), dimension(2) :: x
            integer, dimension(2) :: xint

            xint(1) = (x(1)/xL + 0.5_K)*(n-1)+1
            xint(2) = (x(1)/xL + 0.5_K)*(n-1)+1
      end 
end
