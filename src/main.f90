PROGRAM Bloch
  !$ use omp_lib
  use CONSTANTS
  use WAVE_FUNC
  use TD_CALC
  use COEFF

  implicit none
  integer :: ib
  !$ double precision st, en
  !real(8) :: zuH0zu, zuEtzu
  !complex(kind(0d0)),dimension(:) :: zu_ltov(0:Nx-1,1:Nk)
  !real(8) :: coef
  !integer :: iexp

  select case(igauge)
  case(1)
    write(*,*) "--In Velocity Gauge."
  case(2)
    write(*,*) "--In Length Gauge."
  case(3)
    write(*,*) "--In Length Gauge Expanding."
  end select

  call Initialize()

  write(*,*) "dx : ",dx
  write(*,*) "dt : ",dt
  write(*,*) "dk : ",dk
  write(*,*) "Nx : ",Nx
  write(*,*) "Nk : " , Nk

  !Grand state calculation
  write(*,*)
  write(*,*)"Calculate Base Wave Function"
  write(*,*)
  write(*,*)"V0 : ", V0
  write(*,*)"Width : ", width
  write(*,*)"a0 : ",a0
  !write(*,*) "Nexp_base : " , Nexp_base
  write(*,*) "Niter : ", Niter
  write(*,*) "Nb : ",Nb

  !$ st = omp_get_wtime()
  !call Imaginary_time_method()
    call ground_state
  !$ en = omp_get_wtime()
  !$ print *, "Elapsed time [sec]:", en-st

  !When Length gauge, phase matching
  select case(igauge)
  case(1)
    call td_init()
  case(2)
    call Phase_matching()
    call td_init()
    call fd_coef_k()
  case(3)
    call Phase_matching()
    call fd_coef_k()
    call Coef_init()
  end select

  call base_plot()


  !$ st = omp_get_wtime()
  write(*,*)
  write(*,*) "Time Development"
  write(*,*)
  write(*,*)"Nt : ",Nt
  write(*,*)"Tau : ",tau
  write(*,*)"E0 : ", E0
  write(*,*)"Omega : ",omega
  !write(*,*)"Tad : ",Tad
  write(*,*) "Nexp : ", Nexp
  write(*,*) "Nbt : ", Nbt

  norm(:)=0.d0
  cur(:,:)=0.d0
  hav(:,:)=0.d0

  open(9,file='./output/tdphi.dat')
  open(8,file='./output/td.dat')
  do it = 0, (Nt), 1

    select case(igauge)
    case(1)
      At_in=At(it)
      Et_in=0
      call Gauge_Velocity()
    case(2)
      At_in=0
      Et_in=Et(it)
      call Gauge_Length()
    case(3)
      At_in=0
      Et_in=Et(it)
      call Coef_Length()
    end select

    !$ if(it==0) write(*,*)"Nthreads : ",nthreads

    if(mod(it,5000)==0 .OR. it==(Nt) )then
       write(*,'(<1>i,<1>e)') it,sum(norm(:))
       !do ix = 0, Nx-1, 1
        ! write(9,'(<1>i,<4>e)') ix*dx,(sum(abs(zu(ix,:,1))**2))/Nk,(sum(abs(zu(ix,:,1))**2)+sum(abs(zu(ix,:,2))**2))/(2.d0*Nk),pot(ix),sum(abs(zu_ltov(ix,:)))/Nk
       !end do
       !write(9,*)''
       !write(9,*)''
    end if
    write(8,'(<1>i,<2>e)') it,norm(1),sum(norm(:))
  end do
  close(8)
  close(9)
  !$ en = omp_get_wtime()
  !$ print *, "Elapsed time [sec]:", en-st
  write(*,*) "Excitation Energy :"
  do ib = 1, Nbt, 1
    write(*,'(<1>i,<1>e)') ib ,(hav(ib,(Nt))-hav_base(ib))
  end do
  open(8,file="./output/td_V.dat")
  do ib = 1, Nbt, 1
    do it = 0, (Nt), 1
      write(8,'(<1>i,<8>e)') it, it*dt,(cur(ib,it)),sum(cur(1:ib,it)),hav(ib,it),sum(hav(1:ib,it)),Et(it),At(it)
    end do
    write(8,*)""
    write(8,*)""
  end do
  close(8)


  !Fourie transformation
  !$ st = omp_get_wtime()
  call FourieT()
  !$ en = omp_get_wtime()
  !$ print *, "Elapsed time [sec]:", en-st

  call gnuplot()

  write(*,*)"finish!"

END PROGRAM
