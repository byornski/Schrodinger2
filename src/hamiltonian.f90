module hamiltonian
  use ISO_FORTRAN_ENV, only : REAL64

  use potential, only : pot_type
  use wave,      only : wvfn_type, band_type
  
  implicit none

  



contains
  
  subroutine diagonalise(wvfn,pot)
    type(wvfn_type), intent(inout) :: wvfn
    type(pot_type),  intent(in)    :: pot

    type(band_type) :: h_band, h_band_step, step_band
    
    real(REAL64), parameter :: step_size = -1e-4_REAL64
    !real(REAL64) ,parameter :: fake_step = -1e-11_REAL64
    real(REAL64) :: step

    integer, parameter :: convergence_window = 3
    real(REAL64), parameter :: convergence_tol = 1d-6

    integer :: convergence_num
    
    
    real(REAL64) :: E
    integer :: nb, ni

    real(REAL64) :: delta_e, last_e

    logical :: first_energy

    call step_band%allocate(wvfn%basis)
    
    do nb=1,wvfn%num_bands

       last_e = 0.0_REAL64
       convergence_num = 0
       step = step_size


       first_energy = .true.
       
       do ni=0,200000

          !Copy step band
          step_band = wvfn%bands(nb)

          !Calculate H|psi>
          h_band = wvfn%bands(nb)%apply_H(pot)
          h_band_step = h_band

          !Calculate E = <psi|H|psi>
          E = wvfn%bands(nb) * h_band

          delta_E = E - last_e
          last_e = e

          if (first_energy) then
             write(*,*) 'Step',ni,'Energy',E + pot%offset
             first_energy = .false.
          else
             write(*,*) 'Step',ni,'Energy',E + pot%offset, 'dE:', delta_E
          end if
          
          if (abs(delta_e) < convergence_tol) then
             convergence_num = convergence_num + 1
          else
             convergence_num = 0
          end if
          
          if (convergence_num .ge. convergence_window) exit

          !Line min
          call line_min(wvfn%bands(nb),h_band,pot,step_size)

       end do

       write(*,*) 'finished band',nb
       
    end do

    
  end subroutine diagonalise

  subroutine line_min(band,h_band,pot,trial_step_size,debug_output)
    type(band_type), intent(inout) :: band
    type(band_type), intent(in)    :: h_band
    type(pot_type), intent(in) :: pot
    real(REAL64), intent(in) :: trial_step_size

    logical, intent(in), optional :: debug_output
    
    type(band_type) :: step_band
    real(REAL64) :: E, last_e

    real(REAL64) :: enew, de

    logical :: debug

    debug = .false.
    if (present(debug_output)) debug = debug_output
    
    !Calculate E
    E = band * h_band
    last_e = E

    !Copy h_band
    call step_band%allocate(h_band%basis)
    step_band = h_band

    !Scale small
    call step_band%scale(trial_step_size)

    !Calculate e
    step_band = band + step_band
    call step_band%normalize()
    enew = calc_energy(step_band,pot)
    de = enew - last_e
    
    if (de > 0) then
       block
         logical :: found_step_size
         found_step_size = .false.

         do while (.not. found_step_size)

            !Shrink step
            call step_band%scale(0.1_REAL64)
            
            !Calculate e
            step_band = band + step_band
            call step_band%normalize()
            enew = calc_energy(step_band,pot)
            de = enew - last_e

            found_step_size = de < 0
            
         end do
       end block
    end if

    
    if (debug) write(*,*) 'lilstep', enew,de
    
    last_e = enew
          
    !Now walk in this direction while the energy goes down
    do while(de < 0)
       step_band = step_band + step_band
       call step_band%normalize()
       enew = calc_energy(step_band,pot)
       de = enew - last_e
       last_e = enew
       if (debug) write(*,*) 'lilstep', enew, de
    end do
    
    band = step_band
    

  end subroutine line_min

  

  real(REAL64) function calc_energy(band,pot)
    type(band_type), intent(in) :: band
    type(pot_type), intent(in) :: pot
    type(band_type) :: h_band


    h_band = band%apply_H(pot)

    calc_energy = band * h_band

    

  end function calc_energy

end module hamiltonian
