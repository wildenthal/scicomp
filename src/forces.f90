module forces
      use shapes
      implicit none
      private
      
      integer, parameter :: K   = kind(1d0)
      real(K), parameter :: pi  = 4*ATAN(1.0_K)
      real(K), parameter :: g   = 9.8_K          ! gravity
      real(K), parameter :: rho = 1000_K         ! density of water
      real(K), parameter :: gam = 72.8e-3        ! air-water surface tension

      public balance, integrate_normals, integrate_time, tilt_water
      
      interface balance
            module procedure balance_sphere
      end interface
      
      interface integrate_normals
            module procedure integrate_normals_sphere
      end interface

      contains

      subroutine integrate_time(obj,step,time_step, history)
            type(sphere) :: obj
            integer      :: step
            real(K)      :: time_step, history(:,:)
            obj%coordinates = 2*obj%coordinates-history(:,step-2) + &
                                                obj%acceleration*time_step**2
            history(:,step) = obj%coordinates
      end

      subroutine integrate_normals_sphere(obj)
            type(sphere) :: obj
            integer      :: i, j, n
            real(K)      :: h, norm, water, pos(2), normal(3)

           n = size(obj%reference_frame,1)

            if (obj%height < 0) then
                        stop "not implemented: &
                                         & center of mass below water level"
            end if
            
            normal = 0

            do j = -n/2,n/2 
                  do i = -n/2,n/2
                        pos = local_pos((/i,j/),obj)
                        h = obj%radius**2 - pos(1)**2 - pos(2)**2
                        h = obj%height-sqrt(h)
                        water = obj%reference_frame(i,j)
                        if (abs(water-h) < 1d-9) then
                              normal(1) = normal(1) - pos(1)
                              normal(2) = normal(2) - pos(2)
                              normal(3) = normal(3) - water
                        end if
                  end do
            end do
            norm = sqrt(normal(1)**2+normal(2)**2+normal(3)**2)
            normal = normal/norm
            obj%acceleration = -normal(1:2)*normal(3)*obj%mass*g
      end

      subroutine balance_sphere(obj)
            type(sphere) :: obj
            real(K)      :: weight, bouyancy, force, eps
            weight = obj%mass*g

            do
                  bouyancy = submerged_volume(obj)*rho*g
                  force = bouyancy - weight
            
                  eps = abs(force)/weight

                  if (eps < 1e-7) exit
                  if (force > 0) then
                        obj%height = obj%height + obj%radius*eps
                  else
                        obj%height = obj%height - obj%radius*eps
                  end if
            end do
      end

      subroutine tilt_water(obj)
            type(sphere) :: obj
            integer      :: i, j, n
            real(K)      :: slope
            n = size(obj%reference_frame,1)
            do j = -n/2, n/2
                  do i = -n/2, n/2
                        slope = real(i)/n*obj%radius
                        obj%reference_frame(i,j) = obj%reference_frame(i,j)+slope
                  end do
            end do
      end
end
