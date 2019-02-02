module potential
  use ISO_FORTRAN_ENV, only : REAL64
  use well_parameters, only : well_params
  use basis,           only : basis_type

  implicit none
  
  type pot_type
     real(REAL64), allocatable :: coeffs(:)
     type(basis_type)          :: basis
     real(REAL64)              :: offset
   contains
     procedure :: apply_offset => pot_apply_offset
  end type pot_type



  

contains

  subroutine pot_apply_offset(pot,offset)
    class(pot_type), intent(inout) :: pot
    real(REAL64),    intent(in)    :: offset

    real(REAL64) :: relative_offset

    !Find out how much offset we need to apply
    relative_offset = offset - pot%offset

    !Now modify the potential
    pot%coeffs = pot%coeffs - relative_offset

    !And save the total offset
    pot%offset = offset

  end subroutine pot_apply_offset
  
  type(pot_type) function pot_from_file(basis,filename) result(pot)
    use well_parameters, only : well_params, read_params
    type(basis_type),  intent(in) :: basis
    character(len=*),  intent(in) :: filename
    type(well_params) :: wp

    !Read the file
    wp = read_params(filename)

    !Convert to pot on grid
    pot = pot_from_well(basis,wp)

  end function pot_from_file
  
  type(pot_type) function pot_from_well(basis, wp) result(pot)
    type(basis_type),  intent(in) :: basis
    type(well_params), intent(in) :: wp
    
    real(REAL64) :: x
    integer      :: ns

    !Allocate output space
    allocate(pot%coeffs(basis%num_samples))

    !Save the basis
    pot%basis = basis
   
    !Now calculate V at each point
    do ns=1,basis%num_samples
       x = basis%pos(ns)
       pot%coeffs(ns) = wp%V(x)
    end do

    !Set offset to zero
    pot%offset = 0.0_REAL64
    
  end function pot_from_well

end module potential
