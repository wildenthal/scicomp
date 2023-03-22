program cheerios
      use grid_ops 
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

      ! sphere variables
      type sphere
            real(K) :: coordinates(2)
            integer :: intcoord(2)
            real(K) :: height
            real(K) :: radius = 1.3_K
            real(K) :: contact_angle = pi/2
      end type

      ! initializations
      type(sphere) :: sphere1
      type(sphere) :: sphere2
      call initialize_grid(grid,xL,yL,n)
      print *, "grid done!"
      print *, grid(1,1)
      sphere1%coordinates = (/2.0_K,2.0_K/)
      sphere1%intcoord = xint(sphere1%coordinates,grid)
      sphere1%height = 0.25_K
      sphere2%coordinates = (/2.5_K,2.5_K/)
      sphere2%intcoord = xint(sphere2%coordinates,grid)
      sphere2%height = 0.25_K

      
      ! operations
      call make_dip(sphere1)
      !call make_dip(sphere2)
      !call save_grid(grid)
      call integrate_normals(sphere1)


      contains
      
      ! sphere grid operations
      subroutine make_dip(obj)
            type(sphere) :: obj
            integer      :: radius_int
            integer      :: i, j, iL, iH, jL, jH
            real(K)      :: h
            real(K)      :: pos(2)

            if (obj%height < 0) then
                        stop "not implemented: center of mass below water level"
            end if

            radius_int = obj%radius/xL*(n-1)

            iL = obj%intcoord(1)-radius_int
            iH = obj%intcoord(1)+radius_int
            jL = obj%intcoord(2)-radius_int
            jH = obj%intcoord(2)+radius_int
            
            do j = jL, jH
                  do i = iL, iH
                        pos = xpos((/i,j/),grid)
                        h = obj%radius**2 - (pos(1)-obj%coordinates(1))**2 &
                                          - (pos(2)-obj%coordinates(2))**2
                        h = min(0.0_K,obj%height - sqrt(h))
                        grid(i,j) = h
                  end do
            end do
            print *, 'dip made'
      end
      
      real(K) function capillary_dip(obj,pos)
            type(sphere)          :: obj
            real(K), dimension(2) :: pos
            capillary_dip = 0.0_K
      end

      subroutine integrate_normals(obj)
            type(sphere) :: obj
            integer      :: radius_int
            integer      :: i, j, iL, iH, jL, jH, l
            real(K)      :: h
            real(K)      :: pos(2)
            integer      :: contact
            integer      :: normal(2)
            integer      :: counter 

            normal = 0
            counter = 1

            if (obj%height < 0) then
                        stop "not implemented: center of mass below water level"
            end if

            radius_int = obj%radius/xL*(n-1)
            print *, radius_int

            iL = obj%intcoord(1)-radius_int
            iH = obj%intcoord(1)+radius_int
            jL = obj%intcoord(2)-radius_int
            jH = obj%intcoord(2)+radius_int

            open(unit=file_number,file='normal.dat',status='replace')
            do j = jL, jH
                  do i = iL, iH
                        pos = xpos((/i,j/),grid)
                        h = obj%radius**2 - (pos(1)-obj%coordinates(1))**2 &
                                          - (pos(2)-obj%coordinates(2))**2
                        h = obj%height - sqrt(h)
                        if (abs(grid(i,j)-h) < 1d-9) then
                              normal = normal + (/i,j/) - obj%intcoord 
                              if (counter == 4) then
                                    normal = normal - (/1,1/)
                                    counter = 1
                              else
                                    counter = counter + 1
                              end if
                              write (file_number,*) (/i,j/)-obj%intcoord
                        end if
                  end do
            end do
            close(file_number)
            print *, 'integral over normal vectors:'
            print *, normal
      end
      
      subroutine contact_position(obj)
            type(sphere) :: obj 
      ! z * sin phi = B*Sigma(D,theta)
      ! phi = pi - theta + arctan z
      ! sin( pi - theta + arctan z ) = - sin( arctan(z) - theta)
      ! - sin( arctan(z) - theta) ) = 1/sqrt(1+z**2) * ( - z cos(theta) + sin(theta))
      
      ! solve for z: 
      ! z/sqrt(1+z**2)*( - z cos(theta) + sin(theta)) = B*Sigma(D,theta)

      
      end
end program
