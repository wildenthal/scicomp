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

      ! initializations
      type(sphere) :: sphere1
      type(sphere) :: sphere2
      call initialize_grid(grid,xL,yL,n)
      sphere1%coordinates = (/2.0_K,2.0_K/)
      sphere1%intcoord = xint(sphere1%coordinates,grid)
      sphere1%height = 0.25_K
      
      sphere2%coordinates = (/2.5_K,2.5_K/)
      sphere2%intcoord = xint(sphere2%coordinates,grid)
      sphere2%height = 0.25_K

      
      ! operations
      call make_dip(sphere1,grid)
      !call make_dip(sphere2)
      call save_grid(grid)
      !call integrate_normals(sphere1)

      print *, 'Program end'
end program
