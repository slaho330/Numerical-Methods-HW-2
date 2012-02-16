program heat
!this code solves the 3D heat diffusion problem given set parameters.
!Options for solving the code include FTCS, Crank-Nicolson and ADI.
!FTCS and Crank-Nicolson both have constant and periodic b.c.'s.
!ADI only has constant b.c.'s

  implicit none
  real, parameter :: lx = 1
  real, parameter :: ly = 1
  real, parameter :: lz = 1
  integer, parameter :: nx = 10 !number of x divisions
  integer, parameter :: ny = 10 !number of y divisions
  integer, parameter :: nz = 10 !number of z divisions
  real, parameter :: alpha = 0.001 
  real, parameter :: dt = 0.005
  integer, parameter :: nsteps = 1
  real            :: dx, dy, dz, Cx, Cy, Cz
  real, dimension(nx+1,ny+1,nz+1)  :: T,S
  real            :: bc  !boundary condition (- temp = periodic)
  print*, "nsteps: ", nsteps

  dx = lx/nx  !distance between x grid points
  dy = ly/ny  !distance between y grid points
  dz = lz/nz  !distance between z grid points

  Cx = alpha * dt/dx**2
  Cy = alpha * dt/dy**2
  Cz = alpha * dt/dx**2

  !bc = -1  !periodic boundary conditions
  bc = 0   !constant boundary condition (set to desired (non negative) value)
  
  !Initializes T and the S term and then run desired method
  call initialize(dx, dy, dz, nx, ny, nz, T, S, dt, bc)
  call FTCS(nx, ny, nz, Cx, Cy, Cz, nsteps, T, S, dt, bc)

  call initialize(dx, dy, dz, nx, ny, nz, T, S, dt, bc)
  call CrankNicolson(nx, ny, nz, Cx, Cy, Cz, nsteps, T, S, dt, bc)
  
  call initialize(dx, dy, dz, nx, ny, nz, T, S, dt, bc)
  call ADI(nx, ny, nz, Cx, Cy, Cz, dt, T, S, bc, nsteps)
  
  print*, "done!"
end program heat

subroutine initialize(dx, dy, dz, nx, ny, nz, T, S, dt, bc)
  implicit none
  integer :: i, j, k, nx, ny, nz
  real    :: dx, dy, dz, dt, x, y, z
  real, dimension(nx+1, ny+1, nz+1) :: T,S
  real :: bc

  !initialize T
  do i = 1, nx+1
     do j = 1, ny+1
        do k = 1, nz+1
           x = dx * (i-1)
           y = dy * (j-1)
           z = dz * (k-1)
           !gaussian
           T(i,j,k) = exp(-(5*x - 2.5)**2) * exp(-(5*y - 2.5)**2) * exp(-(5*z - 2.5)**2)
           !initialize time independent source term matrix.
           S(i,j,k) = 0
        enddo
     enddo
  enddo
  
  !S(20,20,50) = 50

  !constant boundary conditions
  if (bc.ge.0) then
     do i = 1, ny+1
        do j = 1, nz+1
           T(1, i, j) = bc
           T(nx+1, i, j) = bc
        enddo
     enddo 
     do i = 1, nx+1
        do j = 1, nz+1
           T(i, 1, j) = bc
           T(i, ny+1, j) = bc
        enddo
     enddo
     do i = 1, nx+1
        do j = 1, ny+1
           T(i, j, 1) = bc
           T(i, j, nz+1) = bc
        enddo
     enddo
  endif
end subroutine initialize
   
subroutine FTCS(nx, ny, nz, Cx, Cy, Cz, nsteps, T, S, dt, bc)
  integer :: nsteps, n, i, j, k, nx, ny, nz, clock_rate, clock_max, t1, t2
  real    :: Cx, Cy, Cz, dt
  real, dimension(nx+1, ny+1, nz+1) :: T, Tnew, S
  real :: bc

  do n=1, nsteps
     print*, 'FTCS: ', nx
     call system_clock(t1,clock_rate,clock_max)

     !calcutlate the new T for all the internal elements
     do i=2, nx
        do j=2, ny
           do k=2, nz
              Tnew(i,j,k) = T(i,j,k) + Cx*(T(i+1,j,k) + T(i-1,j,k) - 2*T(i,j,k))&
                   + Cy*(T(i,j+1,k) + T(i,j-1,k) - 2*T(i,j,k))& 
                   + Cz*(T(i,j,k+1) + T(i,j,k-1) - 2*T(i,j,k)) + S(i,j,k)*dt
           enddo
        enddo
     enddo

     if (bc.lt.0) then !if periodic boundary conditions
        !deal with faces, no edges
        !x-faces
        do j = 2, ny
           do k = 2, nz
              Tnew(1,j,k)= T(1,j,k) + Cx*(T(2,j,k) + T(nx,j,k) - 2*T(1,j,k))&
                   + Cy*(T(1,j+1,k) + T(1,j-1,k) - 2*T(1,j,k))& 
                   + Cz*(T(1,j,k+1) + T(1,j,k-1) -2*T(1,j,k))
              Tnew(nx+1,j,k) = Tnew(1,j,k)
           enddo
        enddo
        
        !y-faces
        do i = 2, nx
           do k = 2, nz
              Tnew(i,1,k) = T(i,1,k) + Cx*(T(i+1,1,k) + T(i-1,1,k) - 2*T(i,1,k))&
                   + Cy*(T(i,2,k) + T(i,ny,k) - 2*T(i,1,k))& 
                   + Cz*(T(i,1,k+1) + T(i,1,k-1) -2*T(i,1,k))
              Tnew(i,ny+1,k) = Tnew(i,1,k)
           enddo
        enddo

        !z-faces
        do i = 2, nx
           do j = 2, ny
              Tnew(i,j,1) = T(i,j,1) + Cx*(T(i+1,j,1) + T(i-1,j,1) - 2*T(i,j,1))&
                   + Cy*(T(i,j+1,1) + T(i,j-1,1) - 2*T(i,j,1))& 
                   + Cz*(T(i,j,2) + T(i,j,nz) -2*T(i,j,1))
              Tnew(i,j,nz+1) = Tnew(i,j,1)
           enddo
        enddo

        !deal with edges, no corners
        !x-axis edges
        do i = 2, nx
           Tnew(i,1,1) = T(i,1,1) + Cx*(T(i+1,1,1) + T(i-1,1,1) - 2*T(i,1,1))&
                + Cy*(T(i,2,1) + T(i,ny,1) - 2*T(i,1,1))&
                + Cz*(T(i,1,2) + T(i,1,nz) - 2*T(i,1,1))
           Tnew(i,ny,1) = Tnew(i,1,1)
           Tnew(i,1,nz) = Tnew(i,1,1)
           Tnew(i,ny,nz) = Tnew(i,1,1)
        enddo
        
        !y-axis edges
        do j = 2, ny
           Tnew(1,j,1) = T(1,j,1) + Cx*(T(2,j,1) + T(nx,j,1) - 2*T(1,j,1))&
                + Cy*(T(1,j+1,1) + T(1,j-1,1) - 2*T(1,j,1))&
                + Cz*(T(1,j,2) + T(1,j,nz) - 2*T(1,j,1))
           Tnew(1,j,nz) = Tnew(1,j,1)
           Tnew(nx,j,1) = Tnew(1,j,1)
           Tnew(nx,j,nz) = Tnew(1,j,1)
        enddo

        !z-axis edges
        do k = 2, nz
           Tnew(1,1,k) = T(1,1,k) + Cx*(T(2,1,k) + T(nx,1,k) - 2*T(1,1,k))&
                + Cy*(T(1,2,k) + T(1,ny,k) - 2*T(1,1,k))&
                + Cz*(T(1,1,k+1) + T(1,1,k-1) - 2*T(1,1,k))
           Tnew(1,ny,k) = Tnew(1,1,k)
           Tnew(nx,1,k) = Tnew(1,1,k)
           Tnew(nx,ny,k) = Tnew(1,1,k)
        enddo

        !deal with corners
        Tnew(1,1,1) = T(1,1,1) + Cx*(T(2,1,1) + T(nx,1,1) - 2*T(1,1,1))&
             + Cy*(T(1,2,1) + T(1,ny,1) - 2*T(1,1,1))&
             + Cz*(T(1,1,2) + T(1,1,nz) - 2*T(1,1,1))
        Tnew(1,1,nz) = Tnew(1,1,1)
        Tnew(1,ny,1) = Tnew(1,1,1)
        Tnew(1,ny,nz) = Tnew(1,1,1)
        Tnew(nx,1,1) = Tnew(1,1,1)
        Tnew(nx,1,nz) = Tnew(1,1,1)
        Tnew(nx,ny,1) = Tnew(1,1,1)
        Tnew(nx,ny,nz) = Tnew(1,1,1)
     endif

     T=Tnew
     call system_clock(t2,clock_rate,clock_max)
     write(*,*) real(t2-t1)/real(clock_rate)

  enddo
 
end subroutine FTCS

subroutine CrankNicolson(nx, ny, nz, Cx, Cy, Cz, nsteps, T, S, dt, bc)
  integer :: nsteps, n, i, j, k, l, m, nx, ny, nz, t1, t2, clock_rate,clock_max, col, row
  real    :: Cx, Cy, Cz, dt
  real, dimension(nx+1, ny+1, nz+1) :: T,S
  real, dimension((nx+1)*(ny+1)*(nz+1)) :: Ti, Tf
  real, dimension((nx+1)*(ny+1)*(nz+1),(nx+1)*(ny+1)*(nz+1)) :: A
  real :: bc
  
  do n=1, nsteps
     print*,'Crank Nicholson: ', nx 
     call system_clock(t1,clock_rate,clock_max)
     
     !make A matrix, 7 diagonals
     do i = 1, (nx+1)*(ny+1)*(nz+1)
        do j=1, (nx+1)*(ny+1)*(nz+1)
           A(i,j) = 0
           if (i.eq.j) then 
              A(i,j) = 1+Cx+Cy+Cz
           endif
           if (abs(i-j).eq.1) then
              A(i,j) = -0.5*Cx
           endif
           if (abs(i-j).eq.(nx+1)) then
              A(i,j) = -0.5*Cy
           endif
           if (abs(i-j).eq.((nx+1)*(ny+1))) then
              A(i,j) = -0.5*Cz
           endif
        enddo
     enddo
     !make the "holes" in x and y diagonals
     do k=1, nz+1
        do j=1, ny+1
           do i=1, nx+1
              row = (i+(j-1)*(nx+1)+(k-1)*(nx+1)*(ny+1))
              if (i.eq.1 .and. row.ne.1) then
                 A(row, row-1) = 0
              endif
              if (i.eq.(nx+1) .and. row.ne.((nx+1)*(ny+1)*(nz+1))) then
                 A(row, row+1) = 0
              endif
              if (j.eq.1 .and. row.gt.(nx+1)) then
                 A(row, row-(nx+1)) = 0
              endif
              if (j.eq.(ny+1) .and. row.lt.((nx+1)*(ny+1)*(nz+1) - (nx+1))) then
                 A(row, row+(nx+1)) = 0
              endif
           enddo
        enddo
     enddo

     !periodic boundary conditions in the A matrix
     if (bc.lt.0) then
        do k = 1, nz+1
           do j = 1, ny+1
              do i = 1, nx+1
                 row = (i+(j-1)*(nx+1)+(k-1)*(nx+1)*(ny+1))
                 !x wrap around
                 if (i.eq.1) then
                    A(row,row+(nx-2)) = -0.5*Cx
                 endif
                 if (j.eq.1) then
                    A(row, row+(ny-2)*(nx+1)) = -0.5*Cy
                 endif
                 !y wrap around
                 if (k.eq.1) then
                    A(row, row+(nz-2)*(nx+1)*(ny+1)) = -0.5*Cz
                 endif
                 if (i.eq.(nx+1)) then
                    A(row, row-(nx-2)) = -0.5*Cx
                 endif
                 !z wrap around
                 if (j.eq.(ny+1)) then
                    A(row,row-(ny-2)*(nx+1)) = -0.5*Cy
                 endif
                 if (k.eq.(nz+1)) then
                    A(row,row-(nz-2)*(nx+1)*(ny+1)) = -0.5*Cz
                 endif
              enddo
           enddo
        enddo
     endif

     l=1
     do k=1, nz+1
        do j=1, ny+1
           do i=1, nx+1
              !make the T vector
              if (i.eq.1 .or. j.eq.1 .or. k.eq.1 .or. i.eq.(nx+1) .or. j.eq.(ny+1) .or. k.eq.(nz+1)) then
                 if (bc.ge.0) then
                    Ti(l) = bc !constant boundary conditions
                 else
                    !corners
                    if ((i.eq.1 .or. i.eq.(nx+1)) .and. (j.eq.1 .or. j.eq.(ny+1)) .and. (k.eq.1 .or. k.eq.(nz+1))) then  
                       Ti(l) = T(i,j,k)*(1-Cx-Cy-Cz) + 0.5*Cx*(T(nx,j,k) + T(2,j,k))&
                            + 0.5*Cy*(T(i,ny,k) + T(i,2,k)) + 0.5*Cz*(T(i,j,nz) + T(i,j,2))
                    !edges
                    !z edges
                    else if (((i.eq.1 .or. i.eq.(nx+1)) .and. (j.eq.1 .or. j.eq.(ny+1))) .and. k.ne.1 .and. k.ne.(nz+1)) then
                       Ti(l) = T(i,j,k) *(1-Cx-Cy-Cz) + 0.5*Cx*(T(nx,j,k) + T(2,j,k))&
                            + 0.5*Cy*(T(i,ny,k) + T(i,2,k)) + 0.5*Cz*(T(i,j,k-1) + T(i,j,k+1))
                    !y edges
                    else if (((i.eq.1 .or. i.eq.(nx+1)) .and. (k.eq.1 .or. k.eq.(nz+1))) .and. j.ne.1 .and. j.ne.(ny+1)) then
                       Ti(l) = T(i,j,k) *(1-Cx-Cy-Cz) + 0.5*Cx*(T(nx,j,k) + T(2,j,k))&
                            + 0.5*Cy*(T(i,j+1,k) + T(i,j-1,k)) + 0.5*Cz*(T(i,j,nz) + T(i,j,2))
                    !x edges
                    else if (((j.eq.1 .or. j.eq.(ny+1)) .and. (k.eq.1 .or. k.eq.(nz+1))) .and. i.ne.1 .and. i.ne.(nx+1)) then
                       Ti(l) = T(i,j,k) *(1-Cx-Cy-Cz) + 0.5*Cx*(T(i-1,j,k) + T(i+1,j,k))&
                            + 0.5*Cy*(T(i,ny,k) + T(i,2,k)) + 0.5*Cz*(T(i,j,nz) + T(i,j,2))
                    !faces
                    !constant x faces
                    else if ((i.eq.1 .or. i.eq.(nx+1)) .and. j.ne.1 .and. j.ne.(ny+1) .and. k.ne.1 .and. k.ne.(nz+1)) then
                       Ti(l) = T(i,j,k) *(1-Cx-Cy-Cz) + 0.5*Cx*(T(nx,j,k) + T(2,j,k))&
                            + 0.5*Cy*(T(i,j-1,k) + T(i,j+1,k)) + 0.5*Cz*(T(i,j,k-1) + T(i,j,k+1))
                    !constant y faces
                    else if ((j.eq.1 .or. j.eq.(ny+1)) .and. i.ne.1 .and. i.ne.(nx+1) .and. k.ne.1 .and. k.ne.(nz+1)) then
                       Ti(l) = T(i,j,k) *(1-Cx-Cy-Cz) + 0.5*Cx*(T(i-1,j,k) + T(i+1,j,k))&
                            + 0.5*Cy*(T(i,ny,k) + T(i,2,k)) + 0.5*Cz*(T(i,j,k-1) + T(i,j,k+1))
                    !constant z faces
                    else if ((k.eq.1 .or. k.eq.(nz+1)) .and. i.ne.1 .and. i.ne.(nx+1) .and. j.ne.1 .and. j.ne.(ny+1)) then
                       Ti(l) = T(i,j,k) *(1-Cx-Cy-Cz) + 0.5*Cx*(T(i+1,j,k) + T(i-1,j,k))&
                            + 0.5*Cy*(T(i,j+1,k) + T(i,j-1,k)) + 0.5*Cz*(T(i,j,nz) + T(i,j,2))
                    endif
                 endif
              else
                 Ti(l) = T(i,j,k)*(1-Cx-Cy-Cz) + 0.5*Cx*(T(i-1,j,k)+T(i+1,j,k))&
                      + 0.5*Cy*(T(i,j-1,k) + T(i,j+1,k)) + 0.5*Cz*(T(i,j,k-1)+T(i,j,k+1)) + S(i,j,k)*dt
              endif
              l=l+1
           enddo
        enddo
     enddo
     
     !solve A*Tf = Ti
     call ge_sub(A,Ti,Tf,nx,ny,nz)
     
     !turn back into cube
     l = 1
     do k=1,(nz+1)
        do j=1,(ny+1)
           do i=1,(nx+1)
              if (bc.ge.0) then
                 if (i.eq.1 .or. j.eq.1 .or. k.eq.1 .or. i.eq.(nx+1) .or. j.eq.(ny+1) .or. k.eq.(nz+1)) then
                    T(i,j,k) = 0 
                 else
                    T(i,j,k) = Tf(l)
                 endif
              else
                 T(i,j,k) = Tf(l)
              endif
              l=l+1
           enddo
        enddo
     enddo

     call system_clock(t2,clock_rate,clock_max)
     write(*,*) real(t2-t1)/real(clock_rate)
  
  enddo
  
end subroutine CrankNicolson

subroutine ADI(nx, ny, nz, Cx, Cy, Cz, dt, T, S, bc, nsteps)
  integer :: i, j, k, nx, ny, nz, n, m, nsteps, t1, t2, clock_rate, clock_max
  real :: Cx, Cy, Cz, dt
  real, dimension((nx+1),(ny+1),(nz+1)) :: T, S
  real, dimension(3) :: A
  real, dimension((ny+1)*(nz+1)) :: Tyzi,Tyzf
  real, dimension((nx+1)*(nz+1)) :: Txzi,Txzf
  real, dimension((nx+1)*(ny+1)) :: Txyi,Txyf
  real :: bc

  do n=1, nsteps
     print*, 'ADI: ', nx
     call system_clock(t1,clock_rate,clock_max)
     
     !x part
     !make A (3 element vector represents tridiagonal matrix)
     A(1) = -1/3.*0.5*Cx
     A(2) = 1 + 0.5*2./3.*Cx
     A(3) = -1/3.*0.5*Cx

     !put T into Ti vector for x direction
     do i=2, nx
        l=1
        do k=1, nz+1
           do j=1, ny+1
              !make the T vector
              !boundary
              if (i.eq.1 .or. j.eq.1 .or. k.eq.1 .or. i.eq.(nx+1) .or. j.eq.(ny+1) .or. k.eq.(nz+1)) then
                 Tyzi(l) = bc
              else
                 Tyzi(l) = T(i,j,k)*(1-0.25*(2./3.)*Cy-0.25*(2./3.)*Cz) + 0.25*(1./3.)*Cy*(T(i,j-1,k) + T(i,j+1,k))&
                      +0.25*(1./3.)*Cz*(T(i,j,k-1)+T(i,j,k+1)) + S(i,j,k)*dt/3
              endif
              l=l+1
           enddo
        enddo
        m=(ny+1)*(nz+1) !number of unknowns
        call tri_solve(A,Tyzi,Tyzf,m)
  
        l = 1
        !put calculated values in T
        do k=1,(nz+1)
           do j=1,(ny+1)
              !boundary
              if (i.eq.1 .or. j.eq.1 .or. k.eq.1 .or. i.eq.(nx+1) .or. j.eq.(ny+1) .or. k.eq.(nz+1)) then
                 T(i,j,k) = bc 
              else
                 T(i,j,k) = Tyzf(l)
              endif
              l=l+1
           enddo
        enddo
     enddo
     
     !y part
     !make A (tridiagonal)
     A(1) = -(1./3.)*0.5*Cy
     A(2) = 1 + 0.5*(2./3.)*Cy
     A(3) = -(1./3.)*0.5*Cy

     do j=2, ny
        l=1
        do k=1, nz+1
           do i=1, nx+1
              !make the T vector for the y direction at each y step
              !boundary
              if (i.eq.1 .or. j.eq.1 .or. k.eq.1 .or. i.eq.(nx+1) .or. j.eq.(ny+1) .or. k.eq.(nz+1)) then
                 Txzi(l) = bc
              else
                 Txzi(l) = T(i,j,k)*(1-0.25*(2./3.)*Cx-0.25*(2./3.)*Cz) + 0.25*(1./3.)*Cx*(T(i-1,j,k) + T(i+1,j,k))&
                      +0.25*(1./3.)*Cz*(T(i,j,k-1)+T(i,j,k+1)) + S(i,j,k)*dt/3
              endif
              l=l+1
           enddo
        enddo
        m=(nx+1)*(nz+1) !number of unknowns
        call tri_solve(A,Txzi,Txzf,m)
  
        l = 1
        !put calculated values in T
        do k=1,(nz+1)
           do i=1,(nx+1)
              !boundary
              if (i.eq.1 .or. j.eq.1 .or. k.eq.1 .or. i.eq.(nx+1) .or. j.eq.(ny+1) .or. k.eq.(nz+1)) then
                 T(i,j,k) = bc 
              else
                 T(i,j,k) = Txzf(l)
              endif
              l=l+1
           enddo
        enddo
     enddo
     
     !z part
     !make A (tridiagonal)
     A(1) = -(1./3.)*0.5*Cz
     A(2) = 1 + 0.5*(2./3.)*Cz
     A(3) = -(1./3.)*0.5*Cz

     !put T into Ti vector for each z step
     do k=2, nz
        l=1
        do j=1, ny+1
           do i=1, nx+1
              !make the T vector
              !boundary
              if (i.eq.1 .or. j.eq.1 .or. k.eq.1 .or. i.eq.(nx+1) .or. j.eq.(ny+1) .or. k.eq.(nz+1)) then
                 Txyi(l) = bc
              else
                 Txyi(l) = T(i,j,k)*(1-0.25*(2./3.)*Cy-0.25*(2./3.)*Cx) +0.25*(1./3.)*Cy*(T(i,j-1,k) + T(i,j+1,k))&
                      +0.25*(1./3.)*Cx*(T(i-1,j,k)+T(i+1,j,k)) + S(i,j,k)*dt/3.0
              endif
              l=l+1
           enddo
        enddo
        m=(ny+1)*(nx+1)
        call tri_solve(A,Txyi,Txyf,m)       
  
        l = 1
        !put the calculated values in T
        do j=1,(ny+1)
           do i=1,(nx+1)
              !boundary
              if (i.eq.1 .or. j.eq.1 .or. k.eq.1 .or. i.eq.(nx+1) .or. j.eq.(ny+1) .or. k.eq.(nz+1)) then
                 T(i,j,k) = bc
              else
                 T(i,j,k) = Txyf(l)
              endif
              l=l+1
           enddo
        enddo
     enddo
     
     call system_clock(t2,clock_rate,clock_max)
     write(*,*) real(t2-t1)/real(clock_rate)
  enddo

end subroutine ADI

subroutine tri_solve(A,Ti,Tf,n)
  implicit none
  integer :: nx, ny, nz, i, n
  real, dimension(3) :: A
  real, dimension(n) :: Tf, Ti, low, mid, high, midt, Tt 
  real :: scale

  !diagonal vectors
  do i=1, n
     low(i) = A(1)
     mid(i) = A(2)
     high(i) = A(3)
  enddo

  !remove lower diagonal
  do i = 2, n
     scale = low(i)/mid(i-1)
     mid(i) = mid(i) - scale*high(i-1)
     Ti(i) = Ti(i) - scale*Ti(i-1)
  enddo

  Tf(n) = Ti(n)/mid(n)

  !back sub for answers
  do i=n-1, 1, -1
     Tf(i) = (Ti(i) - high(i)*Tf(i+1))/mid(i)
  enddo
end subroutine tri_solve

subroutine ge_sub(A,Ti,Tf,nx,ny,nz)
  integer :: nx,ny,nz,m,i,j,k
  real :: sum
  real, dimension((nx+1)*(ny+1)*(nz+1)) :: Ti, Tf
  real, dimension((nx+1)*(ny+1)*(nz+1), (nx+1)*(ny+1)*(nz+1)) :: A
  
  m=(nx+1)*(ny+1)*(nz+1)
  !elimination
  do j=1,m-1
     do i = j+1,m
        if (A(i,j).ne.0) then 
           scale = A(i,j)/A(j,j) !scale factor for each row
           A(i,:) = A(i,:) - A(j,:)*scale !matrix values change 
           Ti(i) = Ti(i) - Ti(j)*scale !b vector changes
        endif
     enddo
  enddo
  !back substitution
  do i=1,m
     k=m-i+1
     !addition of lower elements
     sum = 0
     if (k.ne.m) then
        do j=k+1,m
           sum = sum + Tf(j) * A(k,j) 
        enddo
     endif
     Tf(k) = (Ti(k)-sum)/A(k,k) !fill in the x vector
  enddo
end subroutine ge_sub
