! Author: Gaston Barboza
program cheerios
      use grid 
      use shapes
      use forces
      implicit none
      
      !general parameters
      integer, parameter   :: K = kind(0.d0)
      integer, parameter   :: file_number = 11
      real(K), parameter   :: pi = 4*ATAN(1.0_K)
      real(K), allocatable :: history(:,:) ! keeps track of positions
      integer              :: step, max_steps 
      real(K)              :: time_step ! time delta for integrator

      ! grid parameters
      integer              :: n  ! number of gridpoints (n by n)
      real(K)              :: xL ! x length of box
      real(K)              :: yL ! y lenght of box
      
      ! grid variables
      real(K), allocatable :: grid(:,:)

      ! initializations
      type(sphere)         :: sphere1
      
      call read_params()
      
      allocate(history(2,max_steps))
      
      call initialize_grid(grid,xL,yL,n)
      
      print *
      write (*,'(A)',advance='no') "Initialized a grid of size"
      write (*,'(f10.3)',advance='no') xL
      write (*,'(A)',advance='no') " by"
      write (*,'(f10.3)',advance='no') yL
      print *, "meters"

      sphere1 = read_shape('params/sphere1.param')
      
      write (*,'(A)',advance='no') "Added sphere of radius"
      write (*,'(f10.3)',advance='no') sphere1%radius
      write (*,'(A)',advance='no') " meters at position"
      write (*,'(f10.3)') sphere1%coordinates(1), sphere1%coordinates(2)
      print *

      history(:,1) = sphere1%coordinates
      history(:,2) = sphere1%coordinates
      
      ! evolve system
      call balance(sphere1)

      write (*,'(A)',advance='no') "Sphere is of density"
      write (*,'(f10.3)',advance='no') sphere1%density
      write (*,'(A)',advance='no') " kg/m^3 and is floating at a height of"
      write (*,'(f10.3)',advance='no') sphere1%height
      print *, "meters"
      print *

      call tilt_water(sphere1)
      print *, "Tilted water to simulate a capillary effect"
      call make_dip(sphere1)
      do step = 3,max_steps
            call integrate_normals(sphere1)
            call integrate_time(sphere1,step,time_step,history)
      end do

      call save_results()
      print *, "Trajectory of sphere on water surface saved to file."
      print *

      contains

      subroutine read_params()
            integer            :: line = 0
            integer            :: ios = 0
            integer            :: data_position
            character(len=100) :: buffer, label

            open(unit=11,file="params/system.param")

            do while (ios == 0)
                  read(11, '(A)', iostat=ios) buffer
                  if (ios == 0) then
                        line = line + 1
                        data_position = scan(buffer, ' ')
                        label = buffer(1:data_position)
                        buffer = buffer(data_position+1:)

                        select case(label)
                        case('x_length')
                              read(buffer, *, iostat=ios) xL
                        case('y_length')
                              read(buffer, *, iostat=ios) yL
                        case('grid_res')
                              read(buffer, *, iostat=ios) n
                        case('max_steps')
                              read(buffer, *, iostat=ios) max_steps
                        case('time_step')
                              read(buffer, *, iostat=ios) time_step
                        case default
                              print *, "invalid label at line", line
                        end select
                  end if
            end do
            
            close(11)
      end
      
      subroutine save_results()
            integer :: i
            open(unit=12, file="trajectory.dat", status='replace')
                  do i = 1, max_steps
                        write(12, *) history(:,i)
                  enddo
            close(12)
      end

end program
