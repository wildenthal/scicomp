module forces
      use shapes
      implicit none
      private
      
      integer, parameter :: K   = kind(1d0)
      real(K), parameter :: pi  = 4*ATAN(1.0_K)
      real(K), parameter :: g   = 9.8_K          ! gravity
      real(K), parameter :: rho = 1000_K         ! density of water
      real(K), parameter :: gam = 72.8e-3        ! air-water surface tension

      public balance, integrate_normals, integrate_time
      
      interface balance
            module procedure balance_sphere
      end interface
      
      interface integrate_normals
            module procedure integrate_normals_sphere
      end interface

      contains

      subroutine integrate_time(obj,time_step,time_delta,history)
            type(sphere) :: obj
            integer      :: time_step
            real(K)      :: time_delta, history(:,:)

            obj%coordinates = 2*obj%coordinates-history(time_step-2,:) + &
                                                obj%acceleration*time_delta**2
            history(time_step,:) = obj%coordinates
      end

      subroutine integrate_normals_sphere(obj)
            type(sphere) :: obj
            integer      :: i, j, n, added, skipped 
            real(K)      :: h, norm, water, pos(2), normal(3)

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
                        water = obj%reference_frame(i,j)
                        if (abs(water-h) < 1d-9) then
                              normal(1) = normal(1) + pos(1)
                              normal(2) = normal(2) + pos(2)
                              normal(3) = normal(3) - water
                        else
                              skipped = skipped + 1
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
            integer      :: counter
            counter = 0
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
                  counter = counter + 1
            end do
      end
end
