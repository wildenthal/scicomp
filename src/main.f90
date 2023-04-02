program cheerios
      use grid 
      use shapes
      use forces
      implicit none
      
      !general parameters
      integer, parameter :: K = kind(0.d0)
      integer, parameter :: file_number = 11
      real(K), parameter :: pi = 4*ATAN(1.0_K)
      integer, parameter :: max_steps = 100
      real(K)            :: history(2,max_steps)
      integer            :: time_step

      ! grid parameters
      integer, parameter :: n = 200
      real(K), parameter :: xL = 10.0_K
      real(K), parameter :: yL = 10.0_K
      
      ! grid variables
      real(K), allocatable :: grid(:,:)

      ! object parameters
      integer, parameter :: resolution = 100

      ! initializations
      type(sphere), save :: sphere1
      call initialize_grid(grid,xL,yL,n)
      sphere1 = read_shape('params/sphere1.param')
      history(:,1) = sphere1%coordinates
      history(:,2) = sphere1%coordinates
      call make_dip(sphere1)
      
      ! evolve system
      do time_step = 3,max_steps
            call integrate_normals(sphere1)
            call integrate_time(sphere1,time_step,0.01_K,history)
      end do
end program
