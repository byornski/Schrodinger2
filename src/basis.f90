module basis
  use ISO_FORTRAN_ENV, only : REAL64
  use util,            only : assert
  implicit none
  
  type basis_type
     integer :: num_samples
     real(REAL64) :: x_start
     real(REAL64) :: x_end
     real(REAL64) :: step_size

   contains
     !Translation to position
     procedure :: pos => basis_pos

     !Integration routines
     procedure :: integrate => basis_integrate
     procedure :: dot => basis_dot

     !Derivatives
     procedure :: second_deriv => basis_2nd_deriv

     !Equality check
     procedure :: basis_equals
     generic :: operator(.eq.) => basis_equals

     !Assignment (ie copy)
     procedure :: basis_copy
     generic :: assignment(=) => basis_copy
     
  end type basis_type

contains

  function basis_2nd_deriv(basis,vals) result(derivs)
    class(basis_type), intent(in) :: basis
    real(REAL64),      intent(in) :: vals(basis%num_samples)
    real(REAL64)                  :: derivs(basis%num_samples)

    integer :: last, ns
    
    last = size(derivs)

    !Do first and last values assuming outside values are zero
    derivs(1) = vals(2) - 2.0_REAL64 * vals(1)
    derivs(last) = vals(last-1) - 2.0_REAL64 * vals(last)
    
    !Now do the rest of the values properly
    do ns=2,last-1
       derivs(ns) = vals(ns+1) + vals(ns-1) - 2.0_REAL64 * vals(ns)        
    end do

    !Finally normalise
    derivs(:) = derivs(:) / basis%step_size**2
    
  end function basis_2nd_deriv

  function basis_multiply(basis,vals1,vals2) result(vals_mul)
    class(basis_type), intent(in) :: basis
    real(REAL64), dimension(basis%num_samples), intent(in) :: vals1, vals2
    real(REAL64), dimension(basis%num_samples) :: vals_mul
    
    vals_mul(:) = vals1(:) * vals2(:)
    
  end function basis_multiply
  
  real(REAL64) function basis_pos(basis,n) result(x)
    !Translates a basis index n into a x position
    class(basis_type), intent(in) :: basis
    integer,           intent(in) :: n

    x = basis%x_start + basis%step_size * (n-1)

  end function basis_pos

  
  subroutine basis_copy(b1,b2)
    class(basis_type), intent(out) :: b1
    class(basis_type), intent(in)  :: b2

    b1%num_samples = b2%num_samples
    b1%x_start     = b2%x_start
    b1%x_end       = b2%x_end
    b1%step_size   = b2%step_size
    
  end subroutine basis_copy
  

  type(basis_type) function basis_new(x_start, x_end, num_steps) result(basis)
    real(REAL64), intent(in) :: x_start, x_end
    integer,      intent(in) :: num_steps
    
    !Initialises a new basis
    basis%x_start = x_start
    basis%x_end = x_end
    basis%num_samples = num_steps

    !Set up the step size
    basis%step_size = (x_end - x_start) / real(num_steps-1,REAL64)
    write(*,*) basis%step_size
    
  end function basis_new


  real(REAL64) function basis_integrate(basis,values)
    class(basis_type), intent(in) :: basis
    real(REAL64),     intent(in) :: values(:)

    integer :: ns
    
    !Check the number of values is the same as num_samples
    call assert(size(values).eq.basis%num_samples,"basis_integrate: Number of samples not consistent")

    !Do the actual integration
    basis_integrate = 0.0_REAL64
    do ns=1,basis%num_samples
       basis_integrate = basis_integrate + values(ns)
    end do
    
    !Normalise
    basis_integrate = basis_integrate * basis%step_size
    
  end function basis_integrate

  logical function basis_equals(b1,b2)
    class(basis_type), intent(in) :: b1,b2

    !Check basis domains are the same
    basis_equals = (abs(b1%x_start - b2%x_start) .lt. 1d-6) &
         .and. (abs(b1%x_start - b2%x_start) .lt. 1d-6) &
         .and. (b1%num_samples .eq. b2%num_samples)

  end function basis_equals


  
  real(REAL64) function basis_dot(basis,vals1,vals2) 
    class(basis_type), intent(in) :: basis
    real(REAL64),     intent(in) :: vals1(:), vals2(:)

    integer :: ns
    
    !Check the number of values is the same as num_samples
    call assert(size(vals1).eq.basis%num_samples,'wvfn 1 is not of the right size')
    call assert(size(vals2).eq.basis%num_samples,'wvfn 2 is not of the right size')

    !Do the actual integration
    basis_dot = 0.0_REAL64
    do ns=1,basis%num_samples
       basis_dot = basis_dot + vals1(ns) * vals2(ns)
    end do
    
    !Normalise
    basis_dot = basis_dot * basis%step_size

  end function basis_dot


  
  
end module basis
