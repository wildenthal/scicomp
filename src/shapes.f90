module shapes
      implicit none
      private
      
      integer, parameter :: K = kind(1d0)
      real(K), parameter :: pi = 4*ATAN(1.0_K)
      integer, parameter :: file_unit=11
      integer, parameter :: n = 100

      public sphere, make_dip_sphere, read_sphere, integrate_normals, &
                                                                  apply_force
      
      type sphere
            real(K)              :: radius, height, contact_angle
            real(K)              :: coordinates(2), velocity(2), acceleration(2)
            integer              :: intcoord(2)
            real(K)              :: reference_frame(-n/2:n/2,-n/2:n/2)
      end type

      contains

      type(sphere) function read_sphere(file_name)
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
                        case default
                              print *, "invalid label at line", line
                        end select
                  end if
            end do
            
            obj%contact_angle = obj%contact_angle*pi/180

            read_sphere = obj
      end

      function local_pos_sphere(integer_position,obj)
            type(sphere)          :: obj
            integer, dimension(2) :: integer_position
            real(K), dimension(2) :: local_pos_sphere

            local_pos_sphere(1) = obj%radius*integer_position(1)/(n/2)
            local_pos_sphere(2) = obj%radius*integer_position(2)/(n/2)
      end

      subroutine make_dip_sphere(obj)
            type(sphere) :: obj
            integer      :: i, j, counter
            real(K)      :: h, pos(2)
            

            counter = 0
            if (obj%height < 0) then
                        stop "not implemented: &
                                         & center of mass below water level"
            end if
            
            do j = -n/2,n/2 
                  do i = -n/2,n/2
                        pos = local_pos_sphere((/i,j/),obj)
                        h = obj%radius**2 - pos(1)**2 - pos(2)**2
                        h = min(0.0_K,obj%height-sqrt(h))
                        obj%reference_frame(i,j) = h
                        if (h < 0.0_K) then
                              counter = counter + 1
                        end if
                  end do
            end do
            !print *, counter, 'dips'
      end

      subroutine integrate_normals(obj)
            type(sphere) :: obj
            integer      :: i, j, normal(2), added, skipped 
            real(K)      :: h, pos(2)

            if (obj%height < 0) then
                        stop "not implemented: &
                                         & center of mass below water level"
            end if
            
            normal = 0
            skipped = 0
            added = 0

            do j = -n/2,n/2 
                  do i = -n/2,n/2
                        pos = local_pos_sphere((/i,j/),obj)
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
      
      subroutine apply_force(obj1,obj2)
            type(sphere) :: obj1, obj2
            real(K)      :: displacement(2), pos(2)
            integer      :: i,j

            displacement = obj1%coordinates-obj2%coordinates
            do j = -n/2,n/2
                  do i = -n/2,n/2
                        pos = local_pos_sphere((/i,j/),obj1)
                        if (abs(pos(1)**2+pos(2)**2-&
                                          0.25*obj1%radius**2)<1e-6) then
                              displacement = displacement + pos
                              !print *, pos
                              !print *, obj1%reference_frame(i,j)
                              obj1%reference_frame(i,j) = &
                                             obj1%reference_frame(i,j) - &
                                             &capillary_dip(pos,obj2)
                              !print *, obj1%reference_frame(i,j)
                              !print *
                        end if
                  end do
            end do
      end

      real(K) function capillary_dip(pos,obj)
            real(K), dimension(2) :: pos
            type(sphere) :: obj
            capillary_dip =0.0_K
            if (pos(1) > 0.0_K) then
                  capillary_dip = 0.1_K
            end if
      end


      subroutine contact_position(obj)
            type(sphere) :: obj 
            ! z * sin phi = B*Sigma(D,theta)
            ! phi = pi - theta + arctan z
            ! sin( pi - theta + arctan z ) = - sin( arctan(z) - theta)
            ! - sin( arctan(z) - theta) ) =
            !                 1/sqrt(1+z**2) * ( - z cos(theta) + sin(theta))
            ! solve for z: 
            ! z/sqrt(1+z**2)*( - z cos(theta) + sin(theta)) = B*Sigma(D,theta)
            
      end

end
