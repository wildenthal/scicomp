module forces
      use shapes
      implicit none
      private
      
      integer, parameter :: K   = kind(1d0)
      real(K), parameter :: pi  = 4*ATAN(1.0_K)
      real(K), parameter :: g   = 9.8_K          ! gravity
      real(K), parameter :: rho = 1000_K         ! density of water
      real(K), parameter :: gam = 72.8e-3        ! air-water surface tension

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
            print *, obj%height
            print *, counter
            call water_height_sphere(obj)
            print *, obj%water_height
      end
      
      subroutine water_height_sphere(obj)
            type(sphere) :: obj
            real(K)      :: z, a
            integer      :: i
            
            z = obj%water_height
            a = sqrt(rho*g/gam)
            
            print *, obj%contact_angle
            print *, z/obj%radius
            print *, asin(z/obj%radius)
            print *, tan(obj%contact_angle)
            print *, a

            do i = 1, 100
                  z = a*tan(obj%contact_angle - asin(z/obj%radius))
                  print *, z
            end do

            obj%water_height = z
      end

end
