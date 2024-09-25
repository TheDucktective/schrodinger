module parametres
    real    , parameter :: L=1.5 , xstep=0.001 , Tmax = 0.005 
    real    , parameter :: dt = 0.00000025 , pi= acos(-1.0) 
    real    , parameter :: kx = 100,ky=0 , sigma = 0.1 
    real    , parameter :: x0 = 0.3, y0 = L/2
    integer , parameter :: disc = int(L/xstep)+1
    integer , parameter :: steps = Tmax/dt 
    integer , parameter :: dimension = 1 , n = 2*dimension
    integer , parameter :: set = steps/100
    complex , parameter :: i0=complex(0,1)
    character(len=20)   :: potential_shape = 'young'
    real , dimension(disc,disc) :: V
  end module parametres

  module square_well_parameters
    use parametres
    real , parameter :: potential_position = L/2 , V0 = 1E6
    integer :: potential_length = 10, potential_gap = 30
    integer ::  potential_distance = 60
    
  end module square_well_parameters

  module func
    use parametres
    use square_well_parameters
    contains
    subroutine write_in_file(i,psi)
      use parametres
      implicit none                        
      integer , intent(in)                    :: i             
      complex(4), dimension(disc,disc)        :: psi           
      integer                                 :: j,p           
      integer                                 :: file_unit     
      character(len=20)                       :: file_name     
      file_unit = i+10
      !call wave_function(norm,x,y)
      if (file_unit>=10)then
        open(file_unit,file=file_name(i),status='replace')
        do j=1,disc,1
          if (modulo(j,10)==0)then
            do p=1,disc,1
              if (modulo(p,10)==0) then
                write(file_unit,*)(p-1)*xstep,(j-1)*xstep,cabs(psi(p,j))
              endif
            end do
            write(file_unit,*)' '
          endif 
        
        enddo
      endif
      close(file_unit)
    end subroutine write_in_file

    ! Clearing data files from previous run
    subroutine clear_file()
      open(9,file='Tcoef',status='replace')
      open(11,file='data\file_1.dat',status='replace')
      close(9)
      close(11)
    end subroutine clear_file
  end module func

  program main
    use parametres
    use func
    implicit none
    complex(4) , dimension(disc,disc) :: psi
    real                           :: t=0
    call clear_file()
    call init(psi)
    call propagate(t,psi)    
  end program main

  ! Propagation subroutine 
  subroutine propagate(t,psi)
    use parametres
    use func
    use square_well_parameters
    implicit none
    real                        :: t   
    complex(4) , dimension(disc,n) :: x,y             
    complex(4), dimension(disc,disc) :: psi
    integer           :: j,i=1      
    integer           :: pot_x1 , pot_x2
    external          :: deriv           
    i=1
    j=2
    pot_x1 = int(potential_position/xstep)
    pot_x2 = int(potential_position/xstep+potential_length)
    call pot()   
    do while (j <= steps)                           
      call rk4(t,psi,dt,deriv)              
      psi(1,:)=0 ; psi(disc,:)=0                    
      psi(:,1)=0 ; psi(:,disc)=0
      if (modulo(j,set)==0)then                      
        write(*,*)'Writing in file :',i
        i=i+1
        call write_in_file(i,psi)    
        call normal_check(psi)
      endif
      t=t+dt                                                   
      j=j+1
    enddo
  end subroutine propagate


  ! Normalizing subroutine
  subroutine normalize(psi)
    use parametres
    implicit none
    complex(4), dimension(disc,disc), intent(inout) :: psi  
    real :: sum                                                  
    integer :: i,j                                            
    sum=0
    do j=1,disc,1
        do i=1,disc,1
        sum=sum+cabs(psi(i,j))                 
        enddo
    enddo
    psi=psi/(sum)                                                
  end subroutine normalize

  subroutine normal_check(psi)
    use parametres
    implicit none
    complex, dimension(disc,disc), intent(in) :: psi
    real    :: sum                                          
    integer :: i,j                                           
    sum=0
    do j=1,disc,1
        do i=1,disc,1
        sum=sum+cabs(psi(i,j))                                     
        enddo
    enddo
    write(*,*)'Norm = ',sum
  end subroutine normal_check



   ! Potential defining subroutine for the V(x) vector
  subroutine pot()
    use parametres
    use square_well_parameters
    implicit none
    integer      :: i , j , p, c, a 
    if (potential_shape=='square')then        
      i = int(potential_position/xstep)                    
      j = potential_length
      V=0
      V(i:i+j,:)=V0
    else if (potential_shape == 'hole')then
        i = int(potential_position/xstep)                  
        j = potential_length
        p = potential_gap/2
        c = int(disc/2)
        V=0
        V(i:i+j,:)=V0
        V(:,c-p:c+p)=0
    else if (potential_shape == 'young')then
        i = int(potential_position/xstep)                 
        j = potential_length
        p = potential_gap/2
        c = int(disc/2)
        a = potential_distance/2
        V=0
        V(i:i+j,:)=V0
        V(:,c-a-p:c-a+p)=0
        V(:,c+a-p:c+a+p)=0
    endif
  end subroutine pot
  

  ! RK4 Deriv subroutine
  subroutine deriv(t,psi,dx)
    use parametres
    implicit none
    real                        , intent(inout) :: t              
    complex , dimension(disc,disc)    :: w_x,w_y                  
    complex , dimension(disc,disc) , intent(inout) :: psi, dx     
    call deriv_2_esti(psi,w_x,1)                                 
    call deriv_2_esti(psi,w_y,2)
    dx=i0*(w_x+w_y-V*psi)
  end subroutine deriv

  
  ! Second spatial derivative estimation with finite difference method 
  subroutine deriv_2_esti(psi,w,a)
    use parametres
    implicit none
    complex(4) , dimension(disc,disc) , intent(in)      :: psi               
    complex(4) , dimension(disc,disc) , intent(inout)   :: w                  
    integer                                             :: i,a               
    w=0                                                                    
    if (a==1) then
        do i=2,disc-1,1                                                       
            w(i,:)=(psi(i+1,:)-2*psi(i,:)+psi(i-1,:))/(xstep**2)             
        enddo
    else if (a==2) then
        do i=2,disc-1,1                                                        
            w(:,i)=(psi(:,i+1)-2*psi(:,i)+psi(:,i-1))/(xstep**2)
        enddo
    endif
end subroutine deriv_2_esti

  ! File name generation function with integer input
  function file_name(i) result(fname)
    implicit none
    character(len=20) :: fname                            
    integer , intent(in) :: i                             ! File number
    write(fname,'(a,i0,a)')'data\file_',i,'.dat'          ! File name
  end function file_name
  
  ! Initialisation subroutine
  subroutine init(psi)
    use parametres
    use func
    implicit none
    complex(4) , dimension(disc)         :: x,y   
    complex(4) , dimension(disc,disc)    :: w_x,w_y
    complex(4) , dimension(disc,disc)    :: psi
    integer :: i,j           
    do j=1, disc,1
        y(j)=exp(-((j-1)*xstep-y0)**2/(4*sigma**2))*exp(i0*((j-1)*xstep-y0)*ky)  
        do i=1,disc,1  
            psi(i,j)=exp(-((i-1)*xstep-x0)**2/(4*sigma**2))*exp(i0*((i-1)*xstep-x0)*kx)*y(j)  
        enddo
    enddo
    call normalize(psi)
    call normal_check(psi)
    call deriv_2_esti(psi,w_x,1)
    call deriv_2_esti(psi,w_y,2)
    call write_in_file(1,psi)
  end subroutine init
  
  !RK4
  subroutine rk4(t,x,dt2,deriv)
    ! 4th order Runge-Kutta. See Numerical Recipes p. 701ff
    use parametres
    implicit none
    real              , intent(in)               :: t, dt2
    complex(4), dimension(disc,disc), intent(inout) :: x
    real                                         :: ddt
    complex, dimension(disc,disc)                   :: xp, k1, k2, k3, k4
    external :: deriv
    ddt = 0.5*dt2
    call deriv(t,x,k1)      ; xp = x + ddt*k1
    call deriv(+ddt,xp,k2) ;  xp = x + ddt*k2
    call deriv(t+ddt,xp,k3) ; xp = x +  dt*k3
    call deriv(t+dt,xp,k4)  ; x  = x +  dt*( k1 +2.0*k2 + 2.0*k3 + k4 )/6.0
  end subroutine rk4
