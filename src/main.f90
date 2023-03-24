program cheerios
      use grid_ops 
      use shapes
      implicit none
      
      !general parameters
      integer, parameter :: K = kind(0.d0)
      integer, parameter :: file_number = 11
      real(K), parameter :: pi = 4*ATAN(1.0_K)

      ! grid parameters
      integer, parameter :: n = 200
      real(K), parameter :: xL = 10.0_K
      real(K), parameter :: yL = 10.0_K
      
      ! grid variables
      real(K), allocatable :: grid(:,:)

      ! object parameters
      integer, parameter :: resolution = 100

      ! initializations
      type(sphere) :: sphere1
      type(sphere) :: sphere2
      call initialize_grid(grid,xL,yL,n)
      sphere1 = read_sphere('params/sphere1.param')
      call make_dip_sphere(sphere1)
      call integrate_normals(sphere1)

      !print *, sphere1%radius, sphere1%coordinates, sphere1%height, &
      !      &sphere1%contact_angle, sphere1%resolution
      
      ! operations
      !call make_dip(sphere1,grid)
      !call make_dip(sphere2)
      !call save_grid(grid)
      !call integrate_normals(sphere1)

      print *, 'Program end'
end program
