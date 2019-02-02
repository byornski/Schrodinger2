program test
  
  use ISO_FORTRAN_ENV, only : REAL64
  use wave,  only : wvfn_type, wvfn_allocate
  use basis, only : basis_type, basis_new
  use potential, only : pot_type, pot_from_file
  use hamiltonian, only : diagonalise
  
  implicit none
  
  type(wvfn_type) :: wvfn
  type(basis_type) :: main_basis
  type(pot_type)  :: well_pot

  integer :: ns

  call init_random_seed()
  
  !Set up the basis
  main_basis = basis_new(-100._REAL64, 100._REAL64, 2000)

  !Create the potential
  well_pot = pot_from_file(main_basis,"test.inp")
  call well_pot%apply_offset(50.0_REAL64)
  
  !Now set up the wvfn
  wvfn = wvfn_allocate(1,main_basis)

  !Randomize it
  call wvfn%randomize()

  !Try and diag
  call diagonalise(wvfn,well_pot)
  
  !Write it
  call wvfn%write_file("test.wvfn")

  !Write the pot
  open(file="pot.dat",unit=106)
  do ns=1,well_pot%basis%num_samples
     write(106,*) well_pot%basis%pos(ns), well_pot%coeffs(ns)
  end do
  close(106)

  
contains
  SUBROUTINE init_random_seed()
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
            
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    
    CALL SYSTEM_CLOCK(COUNT=clock)
    
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)
          
    DEALLOCATE(seed)
  END SUBROUTINE init_random_seed
  
end program test
