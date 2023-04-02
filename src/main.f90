program cheerios
      use grid 
      use shapes
      use forces
      implicit none
      
      !general parameters
      integer, parameter   :: K = kind(0.d0)
      integer, parameter   :: file_number = 11
      real(K), parameter   :: pi = 4*ATAN(1.0_K)
      real(K), allocatable :: history(:,:)
      integer              :: step, max_steps
      real(K)              :: time_step

      ! grid parameters
      integer              :: n
      real(K)              :: xL
      real(K)              :: yL
      
      ! grid variables
      real(K), allocatable :: grid(:,:)

      ! object parameters
      integer              :: resolution

      ! initializations
      type(sphere)         :: sphere1
      call read_params()
      allocate(history(2,max_steps))
      call initialize_grid(grid,xL,yL,n)
      sphere1 = read_shape('params/sphere1.param')
      history(:,1) = sphere1%coordinates
      history(:,2) = sphere1%coordinates
      
      ! evolve system
      call tilt_water(sphere1)
      call make_dip(sphere1)
      do step = 3,max_steps
            call integrate_normals(sphere1)
            call integrate_time(sphere1,step,time_step,history)
      end do
!      print *, history

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
                        case('local_res')
                              read(buffer, *, iostat=ios) resolution
                        case default
                              print *, "invalid label at line", line
                        end select
                  end if
            end do
      end

end program
