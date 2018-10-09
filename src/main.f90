PROGRAM Bloch
  !$ use omp_lib
  use CONSTANTS
  use WAVE_FUNC
  use TD_CALC
  use COEFF

  implicit none
  integer :: ib,it
  !$ double precision st, en
  !real(8) :: zuH0zu, zuEtzu
  !complex(kind(0d0)),dimension(:) :: zu_ltov(0:Nx-1,1:Nk)
  !real(8) :: coef
  !integer :: iexp

  call read_input_file()

  select case(igauge)
  case(1)
    write(*,*) "--In Velocity Gauge."
  case(2)
    write(*,*) "--In Length Gauge."
  case(3)
    write(*,*) "--In Length Gauge Expanding."
  end select

  call Initialize()

  !$ write(*,*)"Nthreads : ",NUMBER_THREADS
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
    !call Coef_init()
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

  norm(:,:)=0.d0
   cur(:,:)=0.d0
   hav(:,:)=0.d0

  do it = 0, (Nt), 1

    select case(igauge)
    case(1)
      call Gauge_Velocity(it,At(it))
    case(2)
      call Gauge_Length(it,Et(it))
    case(3)
      !call Coef_Length(it,Et(it))
    end select


    if(mod(it,5000)==0 .OR. it==(Nt) )then
       write(*,'(<1>i,<1>e)') it,sum(norm(:,it))
    end if
  end do
  !$ en = omp_get_wtime()
  !$ print *, "Elapsed time [sec]:", en-st

  write(*,*) "Excitation Energy :"
  do ib = 1, Nbt, 1
    write(*,'(<1>i,<1>e)') ib ,(hav(ib,(Nt))-hav_base(ib))
  end do
  open(8,file="./output/td.dat")
  do it=0,Nt
    write(8,'(<1>i,<4>e)') it,real(it)*dt,norm(:,it)
  end do
  close(8)
  open(8,file="./output/td_V.dat")
  write(8,*) "it, it*dt,(cur(ib,it)),sum(cur(1:ib,it)),hav(ib,it),sum(hav(1:ib,it)),Et(it),At(it)"
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
