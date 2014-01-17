!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Кластер иридия (Sutton-Chen potensial)
!    Квазидинамический режим (атомные единицы)
! 
!  Input data:
!    na - number of atoms
!    FCC.dat - coordinates of the atoms
!
!  Literature:
!  [1] http://www.wag.caltech.edu/home-pages/tahir/psfiles/51.ps

  Program SC_QDM_N

  implicit NONE

  integer(4), parameter::itM =10000 !  максимальное число итераций
  integer(4), parameter::nmax=  itM !  максимальное число точек

  integer(4), parameter::nf1 =   17
  integer(4), parameter::nf3 =   23
  integer(4), parameter::nf4 =   24
  integer(4), parameter::nf5 =   25

  integer(4), parameter::nat = 2123 !  R00(3,nat)

  integer(4) npar,mpar
  real(8)    epsilon,c,a0,massa

  real(8)    dt,Epot(itM),Ekin(itM),s

  real(8),   allocatable::P(:,:),P0(:,:),P1(:,:),   &
                          R0(:,:),R(:,:),dWtot(:,:),E(:,:)

  real(8)    R00(3*nat), dR(3),tm0,Tm1,t0,t,TimMin,eps,Fm,    &
             Wtot,g0,Cenergy,Clength,Cmass, R$(3)


  logical(4) Kmin,GlobIt

  integer(4) jpr,nnn,iter,itfin,itest,i,j,mode,     &
             m,mm,ii,it,ion,k,nt,na,ierr

  integer(4) frequency,duration
  character  inf*79,CharTrans2*4,name1*7/'Ir_QDM_'/,name*11

  character  text1*14/' Ir   77.0    '/,  &
             text2*14/'      SBK     '/

  character(20) :: filename = 'tetrahedron.mat'//CHAR(0)

  common/blmin/jpr,nnn,iter
  common/blpar/epsilon,c,a0,massa,npar,mpar

  data inf/' Sutton-Chen potential, a.u.'/

!       -------------- Атомные единицы -------------
!       1 a.е.энергии =        27.2116  эВ
!       1 a.е.длины   =        0.529177 A
!       1 a.е.заряда  =        4.80324d-10 е.з. СГСЭ
!       1 a.е.дип.мом.=        2.541764133 Дебая
!       1 a.e.поляризуемости = 0.529177**3 A^3
!       1 а.е.времени =        2.419d-17 c
!       1 а.е.массы =          9.1095d-28 г
!       Mproton/Melectron =    1836.15
!       --------------------------------------------

  call read_input(filename, R00, na)

  name=name1//CharTrans2(na)

  Open(nf1,FILE=name//'.res',FORM='formatted',STATUS='unknown')
  Open(nf3,FILE=name//'.xyz',FORM='formatted',STATUS='unknown')
  Open(nf4,FILE=name//'.ang',FORM='formatted',STATUS='unknown')
  Open(nf5,FILE=name//'.ene',FORM='formatted',STATUS='unknown')

  allocate(P(3,na),P0(3,na),P1(3,na),        &
           R0(3,na),R(3,na),dWtot(3,na),E(3,na),stat=ierr)
  if(ierr /= 0)then
   STOP ' SC_QDM_N: P(3,na),P0(3,na),... allocation failed!'
  endif


  tm0=TimMin()

  Cenergy=27.2116d0    ! http://en.wikipedia.org/wiki/Hartree
  Clength=0.529177d0   ! http://en.wikipedia.org/wiki/Bohr_radius
  Cmass  =1836.15d0    ! Mproton/Melectron

  epsilon=epsilon/Cenergy ! перевод в а.е. энергии
  a0=a0/Clength           ! перевод в а.е. длины
  massa=massa*Cmass

! ====================================================================
  dt   =10.d0        !   временной шаг (а.е.)
  eps  =1.d-8        !   нижняя граница пикового значения Екин
  mode =1            !   0, прогон без релаксации
  itest=0            !   тест для произв. полн. энергии ( 0 - нет )
! ====================================================================
  write(nf1,1111)
  write(*  ,1111)
  1111 format(/1x,79(1h=)/23x,'I R I D I U M ( quasidynamical mode )' &
                                                         /1x,79(1h=))

  write(nf1,1112)na,eps
  1112 format(' number of atoms = ',i12/ &
              ' eps             = ',1pd12.5/)

! Генерация координат ( в нулевом прибл. ) и эфф. зарядов ионов МК;
! задание степеней свободы ионов при релаксации

  do i=1,na

   ! Print coordinates
   write(nf1,11)i,R00((i-1)*3+1),R00((i-1)*3+2),R00((i-1)*3+3),massa
   11 format(1x,i3,')',5x,'( ',f8.2,3x,f8.2,3x,f8.2,' )',5x,'mass=',f12.5)

   do m=1,3
    R0(m,i)=a0*R00((i-1)*3+m)
   enddo

  enddo

  do m=1,3
   do j=1,na
    R(m,j)=R0(m,j)
   enddo
  enddo

! write(nf1,22)((R(mm,ii),mm=1,3),ii=1,2)
! 22 format(/' R',2x,1p3d12.5,2x,3d12.5)

! Вычисление потенциальной энергии и ее частных производных
! для стартовой пространственной конфигурации
  Call FgSC(na,R,Wtot,dWtot,g0)
  write(nf1,8200)Cenergy*Wtot
  write(nf1,8201)(dWtot(mm,1),mm=1,3)
  write(nf1,8202)g0

  do i=1,na
   P1(1,i)=0.d0
   P1(2,i)=0.d0
   P1(3,i)=0.d0
  enddo

  do it=1,itM
   Epot(i)=+1.d0
   Ekin(i)=-1.d0
  enddo

  it=1
  Ekin(1)=0.d0
  Epot(1)=Wtot
  Kmin=.FALSE.
  GlobIt=.FALSE.

  t0=TimMin()
  100 continue  !  Hачало итераций /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

  if(Kmin)then ! Kmin=.TRUE. - локальный минимум //////////////////
! print *,char(7)

   t=TimMin()-t0
   write(nf1,10)it,Ekin(it),Epot(it),Ekin(it)+Epot(it)
   write(*  ,10)it,Ekin(it),Epot(it),Ekin(it)+Epot(it)
   10  format(5x,i5,')',5x,'Ekin=',1pd12.5,3x,  &
                           'Epot=',d12.5,3x,'Etot=',d12.5)
   t0=TimMin()
   if(Ekin(it).lt.eps)then
    itfin=it
    goto 200  ! выход из итерационного процесса @@@@@@@@
   endif
   if(it.ge.itM)goto 202
   do i=1,na   ! обнуление импульсов в локальном минимуме
    P1(1,i)=0.d0
    P1(2,i)=0.d0
    P1(3,i)=0.d0
   enddo
   Kmin=.FALSE.
   GlobIt=.FALSE.
  endif ! if(Kmin)...

  it=it+1

  if(it.gt.1)then
   Call FgSC(na,R,Wtot,dWtot,g0)
  endif
  do i=1,na
   P0(1,i)=P1(1,i)
   P0(2,i)=P1(2,i)
   P0(3,i)=P1(3,i)
   P1(1,i)=P0(1,i)-dt*dWtot(1,i)
   R(1,i)=R(1,i)+dt*P1(1,i)/massa
   P1(2,i)=P0(2,i)-dt*dWtot(2,i)
   R(2,i)=R(2,i)+dt*P1(2,i)/massa
   P1(3,i)=P0(3,i)-dt*dWtot(3,i)
   R(3,i)=R(3,i)+dt*P1(3,i)/massa
   P(1,i)=0.5d0*(P1(1,i)+P0(1,i))
   P(2,i)=0.5d0*(P1(2,i)+P0(2,i))
   P(3,i)=0.5d0*(P1(3,i)+P0(3,i))
  enddo

  s=0.d0
  do i=1,na
   s=s+(P(1,i)**2+P(2,i)**2+P(3,i)**2)/massa
  enddo

  Ekin(it)=0.5d0*s
! Вычисление потенциальной энергии и ее частных производных
! для текущей пространственной конфигурации
  Call FgSC(na,R,Wtot,dWtot,g0)
  Epot(it)=Wtot

  if(it.gt.2 .AND. GlobIt)then
   if(Ekin(it).lt.Ekin(it-1))Kmin=.TRUE.
  endif

  GlobIt=.TRUE.

  if(it.lt.itM)goto 100

  202 continue
  write(nf1,201)
  write(*  ,201)
  201 format(//' Self-consistency is not reached...')
  stop

  200 continue ! /////////////////////////////////////////////

  500 continue

  Wtot=Cenergy*Wtot
  write(nf1,8200)Wtot
  8200 format(/' Wtot=',1pd15.8,' eV')
  write(nf1,8201)(dWtot(mm,1),mm=1,3)
  8201 format(/' dWtot(1,1)=',1pd16.9/' dWtot(2,1)=',d16.9/ &
                                      ' dWtot(3,1)=',d16.9)
  write(nf1,8202)g0
  8202 format(/' gmax=',1pd15.8)

  write(nf1,8300)
  8300 format(/1x,79(1h-)/t24,'Coordinates and Displacements (a0) '/)

  do i=1,na
   do m=1,3
    dR(m)=(R(m,i)-R0(m,i))/a0
   enddo
   if(mode.eq.0)goto 300
   !write(nf1,8307)i,(R(m,i)/a0,m=1,3),(dR(m),m=1,3)
   write(nf1,8307)i,(R(m,i)/a0,m=1,3),(R(m,i),m=1,3)
   8307 format(1x,i3,')',3f12.5,3x,3f12.5)
  enddo

  300 continue

  if(mode.eq.0)goto 777

! перевод кинетической и потенциальной энергий в eV
  do i=1,itfin
   Ekin(i)=Cenergy*Ekin(i)
   Epot(i)=Cenergy*Epot(i)
  enddo

! Файл для сохранения координат кластера в ед. a0
   write(nf3,773)((R(m,i)/a0,m=1,3),i=1,na)
   773 format(1x,3f18.7)

! Файл для сохранения координат кластера в ангстремах (для GAMESS)
   do i=1,na
    write(nf4,774)text1,(R(m,i)*Clength,m=1,3)
    774 format(1x,a14,3f18.7)
    write(nf4,775)text2
    775 format(1x,a14/)
   enddo

! Файл для сохранения Ekin и Epot
  write(nf5,776)(Ekin(i),Epot(i),i=1,itfin)
  776 format(1x,f15.5,3x,f15.5)

  777 tm1=TimMin()-tm0

  write(nf1,700)tm1
  print 700,tm1
  700 format(/' <<< Run is over, Tcalc = ',1pd9.3,' min >>>'/  &
                                                  1x,79(1h=))
  print *,char(7)

  deallocate(P,P0,P1,R0,R,dWtot,E,stat=ierr)
  if(ierr /= 0)then
   STOP ' SC_QDM_N: P,P0,... deallocation failed!'
  endif

  stop
  END

!///////////////////////////////////////////////////////////////////

  real*8 function Fi1(r)
  implicit NONE
  integer(4) n,m
  real(8)    r,epsilon,c,a,massa
  common/blpar/epsilon,c,a,massa,n,m
  Fi1=-(a*m*(a/r)**(m-1))/(r*r)
  return
  end

  real*8 function V1(r)
  implicit NONE
  integer(4) n,m
  real(8)    r,epsilon,c,a,massa
  common/blpar/epsilon,c,a,massa,n,m
  V1=-(a*n*(a/r)**(n-1))/(r*r)
  return
  end

  subroutine SubFr(r1,r2,Fr)
  implicit NONE
  real(8)    r1(3),r2(3),Fr(3),d
  integer(4) m
  d=dsqrt((r1(1)-r2(1))**2+(r1(2)-r2(2))**2+(r1(3)-r2(3))**2)
  do m=1,3
   Fr(m)=(r1(m)-r2(m))/d
  enddo
  return
  end

  real*8 function Ro(i,R,na)
  implicit NONE
  integer(4) n,m,i,j,na
  real(8)    epsilon,c,a,massa,sum,R(3,na),d
  common/blpar/epsilon,c,a,massa,n,m
  sum=0.d0
  do j=1,na
   if(j==i)cycle
   d=dsqrt((R(1,i)-R(1,j))**2+    &
           (R(2,i)-R(2,j))**2+    &
           (R(3,i)-R(3,j))**2)
   sum=sum+(a/d)**m
  enddo
  Ro=sum
  return
  end

  integer(4) function del(i,j)
  if(i==j)then
   del=1
  else
   del=0
  endif
  return
  end

  Subroutine FgSC(na,R,Wtot,dWtot,g0)
  implicit NONE

  integer*4  na,i,j,k,m,l,nf1,iter,nnn,jpr
  integer(4) npar,mpar,del,np
  real(8)    epsilon,c,a0,massa
  real(8)    R(3,na),Rmin,V,Vij,W,Wtot,dWtot(3,na),g0,gg, &
         s1(3),s2(3),d,ri(3),rj(3),rk(3),dr(3), &
         den,sden,dum,Ro,V1,Fi1

  !common blocks : http://www.obliquity.com/computer/fortran/common.html

  common/bl1/nf1
  common/blmin/jpr,nnn,iter
  common/blpar/epsilon,c,a0,massa,npar,mpar

  data Rmin/0.1d0/

  W=0.d0

  do m=1,3
   do j=1,na
    dWtot(m,j)=0.d0
   enddo
  enddo

  do i=1,na
   den=Ro(i,R,na)
   sden=dsqrt(den)
   dum=0.5d0*c/sden ! TODO: Why multiplication on 0.5?
   do m=1,3
    ri(m)=R(m,i)
   enddo
   V=0.d0
   do j=1,na
    if(j==i)cycle
    do m=1,3
     rj(m)=R(m,j)
    enddo
    d=dsqrt((ri(1)-rj(1))**2+   &
            (ri(2)-rj(2))**2+   &
            (ri(3)-rj(3))**2)
    if(d <  Rmin)goto 13
    Vij=(a0/d)**npar
    V=V+Vij
    Call SubFr(ri,rj,dr)
    do m=1,3
      s1(m)=0.5d0*V1(d)*dr(m)
    enddo
    do l=1,na
     do m=1,3
      dWtot(m,l)=dWtot(m,l)+s1(m)*(del(i,l)-del(j,l))
     enddo
    enddo ! l
   enddo ! j
   do k=1,na
    if(k==i)cycle
    do m=1,3
     rk(m)=R(m,k)
    enddo
    d=dsqrt((ri(1)-rk(1))**2+   &
            (ri(2)-rk(2))**2+   &
              (ri(3)-rk(3))**2)
    Call SubFr(ri,rk,dr)
    do m=1,3
      s2(m)=dum*Fi1(d)*dr(m)
    enddo
    do l=1,na
     do m=1,3
      dWtot(m,l)=dWtot(m,l)-s2(m)*(del(i,l)-del(k,l))
     enddo
    enddo ! l
   enddo ! k
   W=W+0.5d0*V-c*sden
  enddo ! i

  Wtot=epsilon*W

  g0=0.d0
  do l=1,na
   do m=1,3
    dWtot(m,l)=epsilon*dWtot(m,l)
    g0=g0+dWtot(m,l)**2
   enddo
  enddo
  g0=dsqrt(g0)
  iter=iter+1
  nnn=nnn+1
  if(nnn.lt.jpr)return
  nnn=0

  return

  13 write(nf1,113)i,j,d
  113 format(/' FgSC :       d very small'/    &
              '              i=',i3/           &
              '              j=',i3/           &
              '              d=',1pd12.6/)
  1313 write(nf1,1113)i,k,d
  1113 format(/' FgSC :      d very small'/    &
              '              i=',i3/           &
              '              k=',i3/           &
              '              d=',1pd12.6/)
  Stop
  END

  character(4) function CharTrans2(x)
  ! Преобразование целого положительного
  ! числа < 10000 в строку
  implicit NONE
  integer(4), parameter::m=4
  integer(4) x,xx,i,j,n1(m),k
  character*(m) ss
  character*1 s(10)
  data s/'1','2','3','4','5','6','7','8','9','0'/

  xx = x

  do i=m-1,0,-1
   j=10**i
   n1(m-i) = int(xx/j)
   xx = xx - j * n1(m-i)
!  print*,i,j,n1(m-i),xx
  enddo

  k=n1(1)
  if(k.eq.0)k=10
  ss=s(k)
  do i=2,m
   j=n1(i)
   if(j.eq.0)j=10
   ss=ss(1:i-1)//s(j)
  enddo

  CharTrans2=ss
  return
  END

  Real*8 Function TimMin()
!
!      <<<<< FPS4 >>>>>
!
  implicit None
  real*4 SECNDS,t
! t is number of hundredths of second since previous midnight
  t=SECNDS(0.0)
  TimMin=DBLE(t)/60.d0
  Return
  END

  BLOCK DATA ParIr ! Sutton-Chen potential, a.u.

  integer(4) n,m
  real(8)    epsilon,c,a,massa

  common/blpar/epsilon,c,a,massa,n,m

  data n       /13/
  data m       / 6/

  ! data from table 1 ([1])
  data epsilon /3.7674d-3  / ! в эВ
  data c       /224.8150d0/
  data a       /3.8344d0  / ! в ангстремах
  data massa   /192.2170d0/ ! атомный вес

! Иридий состоит из двух стабильных изотопов:
! Ir191 (37.3%) и Ir193 (62.7%);
! атомный вес 192.217
    end
