module shape_forces
      use shapes
      implicit none
      private
      
      integer, parameter :: K = kind(1d0)
      real(K), parameter :: pi = 4*ATAN(1.0_K)

      public balance, integrate_normals
      
      interface balance
            module procedure balance_sphere
      end interface
      
      interface integrate_normals
            module procedure integrate_normals_sphere
      end interface

      contains

      subroutine integrate_normals_sphere(obj)
            type(sphere) :: obj
            integer      :: i, j, n, normal(2), added, skipped 
            real(K)      :: h, pos(2)

           n = size(obj%reference_frame,1)

            if (obj%height < 0) then
                        stop "not implemented: &
                                         & center of mass below water level"
            end if
            
            normal = 0
            skipped = 0
            added = 0

            do j = -n/2,n/2 
                  do i = -n/2,n/2
                        pos = local_pos((/i,j/),obj)
                        h = obj%radius**2 - pos(1)**2 - pos(2)**2
                        h = obj%height-sqrt(h)
                        if (abs(obj%reference_frame(i,j)-h) < 1d-9) then
                              normal = normal + (/i,j/)
                              added = added + 1
                        else
                              skipped = skipped + 1
                        end if
                  end do
            end do
            obj%acceleration = -normal
      end

      subroutine balance_sphere(obj)
            type(sphere) :: obj
      end

end
