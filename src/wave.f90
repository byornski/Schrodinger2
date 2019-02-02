module wave
  use ISO_FORTRAN_ENV, only : REAL64
  use basis,           only : basis_type
  use potential,       only : pot_type
  use util,            only : assert
  
  implicit none


  !Type for a single band (has a pointer to a basis)
  type band_type
     real(REAL64), allocatable :: coeffs(:)
     real(REAL64)              :: eigenvalue
     type(basis_type)          :: basis 
   contains
     procedure :: allocate => band_allocate
     procedure :: randomize => band_randomize

     !Addition
     procedure :: add => band_add
     generic   :: operator(+) => add
     
     !Dot product
     procedure :: dot => band_dot
     generic   :: operator(*) => dot

     !Normalisation
     procedure :: norm => band_norm
     procedure :: normalize => band_normalize
     procedure :: scale => band_scale

     !Applying potentials
     procedure :: apply_H => band_apply_H
     procedure :: apply_T => band_apply_T
     procedure :: apply_pot => band_apply_pot

     !Copying/assignment
     procedure :: copy => band_copy
     generic   :: assignment(=) => copy
     
  end type band_type
  

  
  type wvfn_type
     integer :: num_bands
     type(basis_type) :: basis
     type(band_type), allocatable :: bands(:)
     
   contains

     procedure :: wvfn_add
     generic   :: operator(+) => wvfn_add

     procedure :: apply_H => wvfn_apply_H
     procedure :: apply_T => wvfn_apply_T
     procedure :: apply_pot => wvfn_apply_pot

     
     procedure :: randomize => wvfn_randomize
     procedure :: normalize => wvfn_normalize
     procedure :: norm => wvfn_norm

     procedure :: scale => wvfn_scale
     
     procedure :: write_file => wvfn_write_file

     procedure :: wvfn_dot
     
  end type wvfn_type



contains


  subroutine band_copy(bcopy, band)
    class(band_type), intent(in) :: band
    class(band_type), intent(out) :: bcopy
    call bcopy%allocate(band%basis)
    bcopy%coeffs = band%coeffs
  end subroutine band_copy

    

  
  subroutine wvfn_write_file(wvfn,filename)
    class(wvfn_type), intent(in) :: wvfn
    character(len=*), intent(in) :: filename

    integer :: nb,ns
    
    open(unit=100,file=filename)

    do nb=1,wvfn%num_bands
       do ns=1,wvfn%basis%num_samples
          write(100,*) wvfn%basis%pos(ns), wvfn%bands(nb)%coeffs(ns)
       end do
       write(100,*)
    end do

    close(unit=100)
    
  end subroutine wvfn_write_file
  
  type(wvfn_type) function wvfn_allocate_like(other_wvfn) result(wvfn)
    !Allocate a wvfn the same size as another
    type(wvfn_type), intent(in) :: other_wvfn
    wvfn = wvfn_allocate(other_wvfn%num_bands,other_wvfn%basis)
  end function wvfn_allocate_like
  
  type(wvfn_type) function wvfn_allocate(num_bands,basis) result(wvfn)
    integer,          intent(in) :: num_bands
    type(basis_type), intent(in) :: basis

    integer :: nb
    
    !Save params
    wvfn%num_bands = num_bands

    !Integration region parameters
    wvfn%basis = basis
    
    !Allocate space for the bands
    allocate(wvfn%bands(num_bands))

    !Now allocate each
    do nb=1,num_bands
       call wvfn%bands(nb)%allocate(basis)
    end do
    
  end function wvfn_allocate


  subroutine band_allocate(band,basis)
    !Allocates the space for a band
    class(band_type), intent(inout) :: band
    type(basis_type), target, intent(in)    :: basis
    band%eigenvalue = 0.0_REAL64
    allocate(band%coeffs(basis%num_samples))
    band%basis = basis
  end subroutine band_allocate



  
  
  subroutine wvfn_randomize(wvfn)
    !Initialise a wvfn with random values
    class(wvfn_type), intent(inout) :: wvfn
    integer :: nb
    do nb=1,wvfn%num_bands
       call wvfn%bands(nb)%randomize()
       !call random_number(wvfn%coeffs(ns,nb))
    end do

  end subroutine wvfn_randomize

  subroutine band_randomize(band)
    class(band_type), intent(inout) :: band
    integer :: ns
    real(REAL64) :: rn
    do ns=1,band%basis%num_samples
       call random_number(rn)
       band%coeffs(ns) = 2_REAL64*rn - 1.0_REAL64
    end do
  end subroutine band_randomize


  real(REAL64) function band_dot(band1,band2)
    class(band_type), intent(in) :: band1, band2
    band_dot = band1%basis%dot(band1%coeffs,band2%coeffs)
  end function band_dot
  
  real(REAL64) function wvfn_dot(wvfn1,wvfn2,nb)
    !Calculates <wvfn1|wvfn2>
    class(wvfn_type), intent(in) :: wvfn1, wvfn2
    integer, intent(in) :: nb
    
    !First check the start and end points are the same
    call assert(wvfn1%basis.eq.wvfn2%basis,"Wvfns dont have the same domain")

    !Now do the integration
    wvfn_dot = wvfn1%basis%dot(wvfn1%bands(nb)%coeffs,wvfn2%bands(nb)%coeffs)
    
  end function wvfn_dot

  real(REAL64) function band_norm(band)
    class(band_type), intent(in) :: band
    band_norm = band_dot(band,band)
  end function band_norm

  real(REAL64) function wvfn_norm(wvfn,nb)
    class(wvfn_type), intent(in) :: wvfn
    integer, intent(in) :: nb
    wvfn_norm = wvfn%bands(nb)%norm()
  end function wvfn_norm
  
  subroutine wvfn_normalize(wvfn)
    !Normalise each band in a wavefunction
    class(wvfn_type), intent(inout) :: wvfn
    integer :: nb

    do nb=1,wvfn%num_bands
       call wvfn%bands(nb)%normalize()
    end do
       
  end subroutine wvfn_normalize


  subroutine band_normalize(band)
    class(band_type), intent(inout) :: band
    real(REAL64) :: norm
    norm = band%norm()
    band%coeffs = band%coeffs / sqrt(norm)
  end subroutine band_normalize

  
  type(wvfn_type) function wvfn_apply_H(wvfn,pot) result(H_wvfn)
    !Applies the hamiltonian to the wvfn H = (T + V)
    class(wvfn_type), intent(in) :: wvfn
    type(pot_type), intent(in) :: pot

    H_wvfn = wvfn_apply_T(wvfn) + wvfn_apply_pot(wvfn,pot) 
    
  end function wvfn_apply_H


  type(band_type) function band_apply_H(band,pot) result(H_band)
    class(band_type), intent(in) :: band
    type(pot_type),   intent(in) :: pot
    H_band = band_apply_T(band) + band_apply_pot(band,pot)
  end function band_apply_H
  
  
  type(wvfn_type) function wvfn_apply_pot(wvfn,pot) result(pot_wvfn)
    !Apply the potential to a wvfn
    use potential, only : pot_type
    
    class(wvfn_type), intent(in) :: wvfn
    type(pot_type),  intent(in) :: pot

    integer :: nb
    
    !Check basis is ok
    call assert(wvfn%basis .eq. pot%basis, 'Pot and wvfn basis not the same!')

    !Allocate output space
    pot_wvfn = wvfn_allocate_like(wvfn)

    !Apply potential
    do nb=1,wvfn%num_bands
       pot_wvfn%bands(nb) = band_apply_pot(wvfn%bands(nb),pot)
    end do
    
  end function wvfn_apply_pot


  type(band_type) function band_apply_pot(band,pot) result(pot_band)
    class(band_type), intent(in) :: band
    type(pot_type),   intent(in) :: pot

    !Check basis is ok
    call assert(band%basis .eq. pot%basis, 'band_apply_pot: Pot and wvfn basis not the same!')

    !Allocate output space
    call pot_band%allocate(band%basis)

    !Apply potential
    pot_band%coeffs = pot%coeffs * band%coeffs
    
  end function band_apply_pot


  
  type(wvfn_type) function wvfn_apply_T(wvfn) result(T_wvfn)
    class(wvfn_type), intent(in) :: wvfn
    integer :: nb
    
    !Allocate output space
    T_wvfn = wvfn_allocate_like(wvfn)

    !Apply T
    do nb=1,wvfn%num_bands
       T_wvfn%bands(nb) = band_apply_T(wvfn%bands(nb))
    end do
    
  end function wvfn_apply_T


  type(band_type) function band_apply_T(band) result(T_band)
    class(band_type), intent(in)  :: band

    !Allocate output space
    call T_band%allocate(band%basis)

    !Apply T
    T_band%coeffs = -band%basis%second_deriv(band%coeffs)

  end function band_apply_T


  !ADDING
  type(wvfn_type) function wvfn_add(wvfn1,wvfn2)
    !Adds two wavefunctions
    class(wvfn_type), intent(in) :: wvfn1, wvfn2
    integer :: nb
    
    call assert(wvfn1%num_bands.eq.wvfn2%num_bands,"Wvfn_add: Num bands not the same")
    
    !Allocate output space
    wvfn_add = wvfn_allocate_like(wvfn1)
    
    !Add
    do nb=1,wvfn1%num_bands
       wvfn_add%bands(nb) = wvfn1%bands(nb) + wvfn2%bands(nb)
    end do
  end function wvfn_add

  type(band_type) function band_add(band1,band2)
    !Adds two bands
    class(band_type), intent(in) :: band1, band2
    
    !Check basis
    call assert(band1%basis.eq.band2%basis,"Wvfn_add: Basis not the same")

    !Allocate
    call band_add%allocate(band1%basis)
    
    !Add
    band_add%coeffs = band1%coeffs + band2%coeffs
    
  end function band_add
  


  
  !SCALING
  subroutine wvfn_scale(wvfn,scale)
    !Scales a wvfn
    class(wvfn_type), intent(inout) :: wvfn
    real(REAL64),     intent(in)    :: scale
    integer :: nb
    do nb=1,wvfn%num_bands
       call wvfn%bands(nb)%scale(scale)
    end do
  end subroutine wvfn_scale
  
  subroutine band_scale(band,scale)
    !Scale a band
    class(band_type), intent(inout) :: band
    real(REAL64),     intent(in)    :: scale
    band%coeffs = band%coeffs * scale
  end subroutine band_scale
    
  
end module wave
