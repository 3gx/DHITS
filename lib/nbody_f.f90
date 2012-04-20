program nbody_program
  implicit none

  integer, parameter :: NBODY_MAX = 128

  integer :: nbody, nbody_out
  real*8  :: time, dt, Tscale;
  real*8  :: Tend, dt_log, dt_out, dt_snap
  integer :: iteration, di_iter

  real*8  :: Mcentre
  real*8  :: mass_scale, pos_scale, vel_scale
  real*8, dimension(NBODY_MAX) :: id
  real*8, dimension(NBODY_MAX) :: mass
  real*8, dimension(NBODY_MAX) :: x
  real*8, dimension(NBODY_MAX) :: y
  real*8, dimension(NBODY_MAX) :: z
  real*8, dimension(NBODY_MAX) :: vx
  real*8, dimension(NBODY_MAX) :: vy
  real*8, dimension(NBODY_MAX) :: vz
  
  real*8, dimension(NBODY_MAX) :: x1
  real*8, dimension(NBODY_MAX) :: y1
  real*8, dimension(NBODY_MAX) :: z1
  real*8, dimension(NBODY_MAX) :: vx1
  real*8, dimension(NBODY_MAX) :: vy1
  real*8, dimension(NBODY_MAX) :: vz1

  integer :: i

  real*8 :: t_log, t_out, t_snap
  real*8 :: Ekin, Epot
  real*8 :: E0, Ep, E1, dE, ddE, dEmax

  real*8 :: wt_beg
  real*8 :: wt0, wt1

  integer :: nstep, iteration1, iter1



  read (*, *) nbody, nbody_out
  write (*, '(a,i4,a,i4)') 'nbody=', nbody, '  nbody_out=', nbody_out
  if (nbody < 0) stop 'nbody is not positive. Please check your input file'

  read (*, *) time, Tscale
  print *, 'time=   ', time, '  Tscale=', Tscale

  read (*, *) Tend, dt_log, dt_out, dt_snap, di_iter;
  print *, 'Tend=   ', Tend;
  print *, 'dTout=  ', dt_out;
  print *, 'dt_log= ', dt_log
  print *, 'dt_snap=', dt_snap
  print *, 'di_iter=', di_iter
  if (di_iter .eq. -1) di_iter = 1000000000

  read (*, *) dt
  print *, 'dt=     ', dt
  dt = dt * Tscale;

  read (*, *) mass_scale, pos_scale, vel_scale;
  print *, " scaling mass     by", mass_scale
  print *, " scaling position by", pos_scale
  print *, " scaling velocity by", vel_scale
  print *, " "

  read (*, *) Mcentre
  Mcentre = Mcentre * mass_scale
  print *, ' Mcentre=', Mcentre

  do i = 1, nbody
    read (*, *) id(i), mass(i), x(i), y(i), z(i), vx(i), vy(i), vz(i)

    mass(i) = mass(i) * mass_scale
    x   (i) = x   (i) *  pos_scale
    y   (i) = y   (i) *  pos_scale
    z   (i) = z   (i) *  pos_scale
    vx  (i) = vx  (i) *  vel_scale
    vy  (i) = vy  (i) *  vel_scale
    vz  (i) = vz  (i) *  vel_scale
    print *, '  i=', i, ' id=', id(i), 'mass=', mass(i)
  end do

  ! here one needs to make sure that the data has
  !   heliocentric coordinates &
  !   barycentric positions
  ! the library does not check this, so this is user's responsibility

  call dh_open(Mcentre);

  call dh_ekin(Ekin, nbody, mass, x,y,z, vx,vy,vz)
  call dh_epot(Epot, nbody, mass, x,y,z, vx,vy,vz)
  E0 = Ekin + Epot
  print *, 'E0=', E0

  Ep = E0
  
  t_log  = time/Tscale
  t_out  = time/Tscale
  t_snap = time/Tscale

  dEmax = 0.0;

  call dh_cvt2symp(nbody, dt, mass, x,y,z, vx,vy,vz)

  call get_wtime(wt0);
  wt_beg = wt0;

  iteration  = 0
  iteration1 = 0
  iter1      = 0

  do while (time/Tscale < Tend) 

    nstep = 32
    iter1 = 0;
    call dh_iterate(nbody, dt, nstep, mass, x,y,z, vx,vy,vz, time, iter1)
    iteration1 = iteration1 + iter1
    iteration  = iteration  + 1

    call get_wtime(wt1)
    wt1 = wt1 + 1.0e-16     ! adding a small number in case wt1 = wt0

    if (time/Tscale > t_log .or. &
        time/Tscale >= Tend .or.  &
        mod(iteration, di_iter) .eq. 0) then
        do i = 1, nbody
          x1 (i) = x (i)
          y1 (i) = y (i)
          z1 (i) = z (i)
          vx1(i) = vx(i)
          vy1(i) = vy(i)
          vz1(i) = vz(i)
        end do
        call dh_cvt2phys(nbody, dt, mass, x1,y1,z1, vx1,vy1,vz1)
        call dh_ekin(Ekin, nbody, mass, x1,y1,z1, vx1,vy1,vz1)
        call dh_epot(Epot, nbody, mass, x1,y1,z1, vx1,vy1,vz1)
        E1  = Ekin + Epot
        dE  = E1 - E0
        ddE = E1 - Ep
        Ep  = E1
        dEmax = max(dEmax, abs(dE/E0))

        write (*, '(a,i9,a,g15.7,a,f8.4,a,3g15.7,a,f7.3,a,a,f4.2,a)') &
               'iter= ', iteration1,  &
               ' Time= ', time, &
               ' dt= ', dt, &
               ' dE= ', dE, ddE, dEmax, &
               ' :: tW= ', (wt1 - wt_beg)/3600.0, ' hr', &
               ' <T>=   ', wt1 - wt0, ' sec'
        
        wt0 = wt1
        t_log = t_log + dt_log
    end if


  end do 




end
