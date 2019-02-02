module well_parameters
  USE ISO_FORTRAN_ENV, only : REAL64
  implicit none

  private

  !Type to hold parameters about a potential well
  type well_params
     !Height of potential walls
     character(len=20) :: filename
     character(len=1)  :: type  !=S for square well
     real(REAL64) :: depth
     real(REAL64) :: width

     !A function to calculate the potential at some position
     procedure(calculate_potential), pointer :: V => null()

  end type well_params


  !Declare an interface for a well potential function (F2003)
  abstract interface
     real(REAL64) function calculate_potential(well_data,x)
       use ISO_FORTRAN_ENV, only : REAL64
       import well_params
       implicit none
       class(well_params), intent(in) :: well_data
       real(REAL64),       intent(in) :: x
     end function calculate_potential
  end interface

  public :: well_params
  public :: read_params
  
contains

  real(REAL64) function well_calculate_potential_square(well_data,x) result(V)
    !Returns the potential V at position x
    class(well_params), intent(in) :: well_data
    real(REAL64),       intent(in) :: x
    
    logical :: inside_well
    
    !First calculate V(y)
    inside_well = abs(x) < (well_data%width / 2.0_REAL64)

    !Potential is zero inside the well
    if (inside_well) then
       V = 0
    else !and some other value outside
       V = well_data%depth
    end if
   
  end function well_calculate_potential_square
  
  type(well_params) function read_params(filename)
    !Reads a set of well parameters from a file
    character(len=*), intent(in) :: filename

    integer :: stat
    integer, parameter :: file_unit = 101

    !Open the file and check for errors
    open(unit=file_unit,file=filename,status='OLD',action='READ',iostat=stat)
    if (stat/=0) then
       write(*,*) 'Failed to open file: ', filename
       STOP 'Failed to read file'
    end if

    !Use the input filename for output seed
    read_params%filename = filename

    !Read the well type
    read(file_unit,'(A)') read_params%type

    !Now read the rest of the data from the file
    select case(read_params%type)
    case('S')

       !Set the potential function pointer
       read_params%V => well_calculate_potential_square
       
       !Square well has only depth and width
       read(file_unit,*,iostat=stat) read_params%depth, read_params%width
       if (stat/=0) STOP "Failed to parse well file"
       
    case default

       STOP 'Invalid well type'
       
    end select
       
    !Finally check the function V has been assigned
    if (.not. associated(read_params%V)) STOP 'Error in developer: Potential function not set!'

    
    !Close the file
    close(unit=file_unit)
    
  end function read_params




  
end module well_parameters

