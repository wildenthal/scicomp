module shapes
      use grid_ops
      implicit none
      private
      
      integer, parameter :: K = kind(1d0)
      real(K), parameter :: pi = 4*ATAN(1.0_K)
      integer, parameter :: file_unit=11

      public sphere, make_dip
      
      type sphere
            real(K)              :: radius = 1.0_K
            real(K)              :: coordinates(2) = (/0.0_K,0.0_K/)
            real(K)              :: height = 1.0_K
            real(K)              :: contact_angle = pi/2
            integer              :: intcoord(2) = (/0,0/)
            integer              :: resolution
            real(K), allocatable :: reference_frame
      end type

      contains

      subroutine initialize_sphere(obj,file_name)
            type(sphere)       :: obj
            character(len=*)   :: file_name
            integer            :: line = 0
            integer            :: ios = 0
            integer            :: data_position
            character(len=100) :: buffer, label

            open(unit=file_unit,file=file_name)
            
            do while (ios == 0)
                  read(file_unit, '(A)', iostat=ios) buffer
                  if (ios == 0) then
                        line = line + 1
                        data_position = scan(buffer, ' ')
                        label = buffer(1:data_position)
                        buffer = buffer(data_position+1:)

                        select case(label)
                        case('radius')
                              read(buffer, *, iostat=ios) obj%radius
                        case('x_coord')
                              read(buffer, *, iostat=ios) obj%coordinates(1)
                        case('y_coord')
                              read(buffer, *, iostat=ios) obj%coordinates(2)
                        case('height')
                              read(buffer, *, iostat=ios) obj%height
                        case('contact_angle')
                              read(buffer, *, iostat=ios) obj%contact_angle
                        case('ref_frame_resolution')
                              read(buffer, *, iostat=ios) obj%resolution
                        end select
                  end if
            end do
      end

      subroutine make_dip(obj,grid)
            type(sphere) :: obj
            integer      :: radius_int
            integer      :: i, j, iL, iH, jL, jH, n
            real(K)      :: h, xL
            real(K)      :: pos(2)
            real(K), dimension(:,:) :: grid

            n = 200
            xL = 10.0_K

            if (obj%height < 0) then
                        stop "not implemented: &
                                         & center of mass below water level"
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
      end
      
      real(K) function capillary_dip(obj,pos)
            type(sphere)          :: obj
            real(K), dimension(2) :: pos
            capillary_dip = 0.0_K
      end

      subroutine integrate_normals(obj,grid)
            type(sphere) :: obj
            integer      :: radius_int
            integer      :: i, j, iL, iH, jL, jH, l
            real(K)      :: h,xL
            real(K)      :: pos(2)
            integer      :: contact
            integer      :: normal(2)
            integer      :: counter, n
            real(K), dimension(:,:) :: grid

            normal = 0
            counter = 1
            n = 200
            xL=10.0_K

            if (obj%height < 0) then
                        stop "not implemented: center of mass below water level"
            end if

            radius_int = obj%radius/xL*(n-1)
            print *, radius_int

            iL = obj%intcoord(1)-radius_int
            iH = obj%intcoord(1)+radius_int
            jL = obj%intcoord(2)-radius_int
            jH = obj%intcoord(2)+radius_int

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
                        end if
                  end do
            end do
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

end
