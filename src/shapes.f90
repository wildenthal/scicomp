! Author: Gaston Barboza
module shapes
      implicit none
      private
      
      ! general parameters
      integer, parameter :: K = kind(1d0)
      real(K), parameter :: pi = 4*ATAN(1.0_K)
      integer, parameter :: file_unit=11
      integer, parameter :: n = 100 ! local grid resolution

      public sphere, make_dip, read_shape, local_pos
      public submerged_volume
      
      ! global interfaces: allows extending procedures to different
      !   shape types
      interface make_dip
            module procedure make_dip_sphere
      end interface

      interface read_shape
            module procedure read_sphere
      end interface
      
      interface local_pos
             module procedure local_pos_sphere
      end interface

      interface submerged_volume
            module procedure submerged_volume_sphere
      end interface
      
      ! shapes go here
      type sphere
            real(K)              :: radius, height, contact_angle 
            real(K)              :: coordinates(2), acceleration(2)
            integer              :: intcoord(2)
            real(K)              :: reference_frame(-n/2:n/2,-n/2:n/2)
            real(K)              :: volume, mass, density
      end type

      contains

      type(sphere) function read_sphere(file_name)
            ! loads all shape characteristics from a file
            ! taken from jblevins.org/log/control-file

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
                        case('density')
                              read(buffer, *, iostat=ios) obj%density
                        case default
                              print *, "invalid label at line", line
                        end select
                  end if
            end do
            
            obj%contact_angle = obj%contact_angle*pi/180
            
            obj%volume = 4*pi*obj%radius**3/3

            obj%mass = obj%volume*obj%density

            read_sphere = obj
      end

      function local_pos_sphere(integer_position,obj)
            ! local reference frame is indexed by integers,
            !   this converts them to coordinate positions

            type(sphere)          :: obj
            integer, dimension(2) :: integer_position
            real(K), dimension(2) :: local_pos_sphere

            local_pos_sphere(1) = obj%radius*integer_position(1)/(n/2)
            local_pos_sphere(2) = obj%radius*integer_position(2)/(n/2)
      end

      subroutine make_dip_sphere(obj)
            ! makes a dip in the water where the sphere is floating

            type(sphere) :: obj
            integer      :: i, j
            real(K)      :: h, pos(2)
            

            if (obj%height < 0) then
                        stop "not implemented: &
                                         & center of mass below water level"
            end if
            
            do j = -n/2,n/2 
                  do i = -n/2,n/2
                        pos = local_pos_sphere((/i,j/),obj)
                        h = obj%radius**2 - pos(1)**2 - pos(2)**2
                        h = min(obj%reference_frame(i,j),obj%height-sqrt(h))
                        obj%reference_frame(i,j) = h
                  end do
            end do
      end 

      real(K) function submerged_volume_sphere(obj)
            ! calculates the volume of the submerged part of the sphere
            !   using the formula for the volume of a spherical cap
            !   of height h: vol = pi*h**2*(radius-h/3)

            type(sphere) :: obj
            real(K)      :: svs, h
            
            if (obj%height>obj%radius) then
                  stop "Sphere ejected from water"
            else if (obj%height<-obj%radius) then
                  stop "Sphere sunk below water"
            end if

            if (obj%height>0) then
                  h = obj%radius-obj%height
                  svs = pi*h**2*(obj%radius-h/3)
            else
                  svs = 4*pi*obj%radius**3/3 
                  h = obj%radius+obj%height
                  svs = svs - pi*h**2*(obj%radius - h/3)
            end if

            submerged_volume_sphere = svs
      end
end
