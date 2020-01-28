      ! export OMP_NUM_THREADS=3
      ! gfortran enzymekinetics2dopt.f -O3 -fopenmp -march=native -ffast-math
      PROGRAM ENZYMEKINETICS2D
!$    use omp_lib
      implicit none

      ! Parameters of the problem
      integer n                                     ! Number of spatial points
      integer m                                     ! Number of temporal points
      double precision eps1                         ! Diffusion coefficients
      double precision, allocatable::x(:)           ! x-grid
      double precision, allocatable::wold(:)        ! solution w=[u,v]
      double precision, allocatable::A1diag(:)      ! A1 for operator splitting
      double precision, allocatable::A1lower(:)
      double precision, allocatable::A1upper(:)
      double precision, allocatable::A2diag(:)      ! A2 for operator splitting
      double precision, allocatable::A2lower(:)
      double precision, allocatable::A2upper(:)
      double precision, allocatable::M1diag(:)
      double precision, allocatable::M1lower(:)
      double precision, allocatable::M1upper(:)
      double precision, allocatable::M2diag(:)
      double precision, allocatable::M2lower(:)
      double precision, allocatable::M2upper(:)
      double precision, allocatable::M3diag(:)
      double precision, allocatable::M3lower(:)
      double precision, allocatable::M3upper(:)
      double precision, allocatable::M11diag(:)
      double precision, allocatable::M11lower(:)
      double precision, allocatable::M11upper(:)
      double precision, allocatable::M22diag(:)
      double precision, allocatable::M22lower(:)
      double precision, allocatable::M22upper(:)
      double precision, allocatable::M33diag(:)
      double precision, allocatable::M33lower(:)
      double precision, allocatable::M33upper(:)
      double precision, allocatable::line(:)
      double precision, allocatable::Fold(:)
      double precision, allocatable::wstar(:)
      double precision, allocatable::cache1(:)
      double precision, allocatable::cache2(:)
      double precision, allocatable::Fstar(:)
      double precision, allocatable::wnext(:)
      double precision, allocatable::a1(:)
      double precision, allocatable::a2(:)
      double precision, allocatable::b1(:)
      double precision, allocatable::b2(:)
      double precision, allocatable::c1(:)
      double precision, allocatable::c2(:)
      double precision, allocatable::d1(:)
      double precision, allocatable::d2(:)
      double precision diff1
      double precision h                            ! spatial step size
      double precision dt                           ! temporal step size
      double precision pi
      double precision starttime,finaltime
      integer nn                                    ! n x n
      integer lowband,highband,dimprob
      integer i,j,ii

      n=399
      m=401
      lowband=1
      highband=n
      eps1=0.2d0
      h=1.0d0/(n+2-1.0d0) ! only interior nodes
      dt=1.0d0/(m-1.0d0)
      nn=n*n;
      dimprob=nn
      pi=4.0d0*atan(1.0d0)
      
      write(*,*) 'Spatial step size h=', h
      write(*,*) 'Temporal step size dt=', dt
      write(*,*) 'Dimension of the system nn=', nn

!$    starttime=omp_get_wtime();
      allocate(x(n))
      CALL LINSPACE(x,h,1.0d0-h,n)

      ! Initial condition for u and v
      allocate(wold(dimprob))
      j=1
      DO ii=1,n
         DO i=1,n
            wold(j)=1.0d0 !sin(PI*x(i))*sin(PI*x(ii)) ! v(x,y,0)=sin(pi*x)*sin(pi*y)
            j=j+1
         END DO
      END DO

      allocate(A1diag(dimprob))
      allocate(A1lower(dimprob))
      allocate(A1upper(dimprob))
      allocate(A2diag(dimprob))
      allocate(A2lower(dimprob))
      allocate(A2upper(dimprob))
      diff1=eps1*dt/(h*h)

      allocate(line(n))

      line=2.0d0
      A1diag=diff1*RESHAPE(SPREAD(line,1,n),(/nn/))
      A2diag=A1diag
      
      line=-1.0d0
      line(1)=0.0d0
      A1lower=diff1*RESHAPE(TRANSPOSE(SPREAD(line,1,n)),(/nn/))

      line=-1.0d0
      line(n)=0.0d0
      A1upper=diff1*RESHAPE(TRANSPOSE(SPREAD(line,1,n)),(/nn/))
      
      line=-1.0d0
      line(1)=0.0d0
      A2lower=diff1*RESHAPE(SPREAD(line,1,n),(/nn/))

      line=-1.0d0
      line(n)=0.0d0
      A2upper=diff1*RESHAPE(SPREAD(line,1,n),(/nn/))

      allocate(M1lower(dimprob))
      allocate(M1diag(dimprob))
      allocate(M1upper(dimprob))
      allocate(M2lower(dimprob))
      allocate(M2diag(dimprob))
      allocate(M2upper(dimprob))
      allocate(M3lower(dimprob))
      allocate(M3diag(dimprob))
      allocate(M3upper(dimprob))

      allocate(M11lower(dimprob))
      allocate(M11diag(dimprob))
      allocate(M11upper(dimprob))
      allocate(M22lower(dimprob))
      allocate(M22diag(dimprob))
      allocate(M22upper(dimprob))
      allocate(M33lower(dimprob))
      allocate(M33diag(dimprob))
      allocate(M33upper(dimprob))

      M1lower=A1lower
      M1diag=1.0d0+A1diag
      M1upper=A1upper
      M2lower=A1lower/3.0d0
      M2diag=1.0d0+A1diag/3.0d0
      M2upper=A1upper/3.0d0
      M3lower=A1lower/4.0d0
      M3diag=1.0d0+A1diag/4.0d0
      M3upper=A1upper/4.0d0

      M11lower=A2lower
      M11diag=1.0d0+A2diag
      M11upper=A2upper
      M22lower=A2lower/3.0d0
      M22diag=1.0d0+A2diag/3.0d0
      M22upper=A2upper/3.0d0
      M33lower=A2lower/4.0d0
      M33diag=1.0d0+A2diag/4.0d0
      M33upper=A2upper/4.0d0

      allocate(Fold(dimprob))
      allocate(wstar(dimprob))
      allocate(cache1(dimprob))
      allocate(cache2(dimprob))
      allocate(Fstar(dimprob))
      allocate(wnext(dimprob))
      allocate(a1(dimprob))
      allocate(a2(dimprob))
      allocate(b1(dimprob))
      allocate(b2(dimprob))
      allocate(c1(dimprob))
      allocate(c2(dimprob))
      allocate(d1(dimprob))
      allocate(d2(dimprob))

      DO i=2,m
         Fold=FUNC(dimprob,wold)
         
!$OMP PARALLEL DO DEFAULT(none) SHARED(wold,dt,Fold,wnext,dimprob)
         DO J=1,dimprob
            wnext(j)=wold(j)+dt*Fold(j)
         END DO
!$OMP END PARALLEL DO

         ! First block starts
!$OMP PARALLEL
!$OMP SECTIONS
!$OMP SECTION
         CALL FORBACK(dimprob,M11lower,M11diag,M11upper,wnext,
     *                wstar,highband)
!$OMP SECTION
         CALL FORBACK(dimprob,M2lower,M2diag,M2upper,wold,a1,lowband)
!$OMP SECTION
         CALL FORBACK(dimprob,M3lower,M3diag,M3upper,wold,b1,lowband)
!$OMP END SECTIONS
!$OMP END PARALLEL

!$OMP PARALLEL DO DEFAULT(none) SHARED(a1,b1,c1,dimprob)
         DO J=1,dimprob
            c1(j)=9.0d0*a1(j)-8.0d0*b1(j)
         END DO
!$OMP END PARALLEL DO

         ! second block starts
!$OMP PARALLEL
!$OMP SECTIONS
!$OMP SECTION
         CALL FORBACK(dimprob,M1lower,M1diag,M1upper,wstar,wstar,
     *                lowband)
!$OMP SECTION
         CALL FORBACK(dimprob,M2lower,M2diag,M2upper,Fold,a2,lowband)
!$OMP SECTION
         CALL FORBACK(dimprob,M3lower,M3diag,M3upper,Fold,b2,lowband)
!$OMP END SECTIONS
!$OMP END PARALLEL

!$OMP PARALLEL DO DEFAULT(none) SHARED(a2,b2,c2,dimprob)
         DO J=1,dimprob
            c2(j)=9.0d0*a2(j)-8.0d0*b2(j)
         END DO
!$OMP END PARALLEL DO

         Fstar=FUNC(dimprob,wstar)

!$OMP PARALLEL DO DEFAULT(none) SHARED(c1,c2,Fstar,dt,cache1,cache2)
!$OMP& SHARED(dimprob)         
         DO J=1,dimprob
            cache1(j)= 9.0d0*c1(j)+2.0d0*dt*c2(j)+dt*Fstar(j)
            cache2(j)=-8.0d0*c1(j)-1.5d0*dt*c2(j)-dt*Fstar(j)/2.0d0
         END DO
!$OMP END PARALLEL DO

         ! Third block starts
!$OMP PARALLEL
!$OMP SECTIONS
!$OMP SECTION
         CALL FORBACK(dimprob,M22lower,M22diag,M22upper,cache1,d1,
     *                highband)
!$OMP SECTION
         CALL FORBACK(dimprob,M33lower,M33diag,M33upper,cache2,d2,
     *                highband)
!$OMP END SECTIONS
!$OMP END PARALLEL

!$OMP PARALLEL DO DEFAULT(none) SHARED(d1,d2,wold,dimprob)
         DO J=1,dimprob
            wold(j)=d1(j)+d2(j)
         END DO
!$OMP END PARALLEL DO
      END DO
!$    finaltime=omp_get_wtime()
      write(*,*) 'Time of the loop',finaltime-starttime
      
      ! Write the solution to a file
      open(20,file='solutionenzyme2d.txt')
      DO I=1,dimprob
         write(20,*) wold(i)
      END DO
      close(20)

      deallocate(x)
      deallocate(wold)
      deallocate(A1diag)
      deallocate(A1lower)
      deallocate(A1upper)
      deallocate(A2diag)
      deallocate(A2lower)
      deallocate(A2upper)
      deallocate(M1lower)
      deallocate(M1diag)
      deallocate(M1upper)
      deallocate(M2lower)
      deallocate(M2diag)
      deallocate(M2upper)
      deallocate(M3lower)
      deallocate(M3diag)
      deallocate(M3upper)
      deallocate(M11lower)
      deallocate(M11diag)
      deallocate(M11upper)
      deallocate(M22lower)
      deallocate(M22diag)
      deallocate(M22upper)
      deallocate(M33lower)
      deallocate(M33diag)
      deallocate(M33upper)
      deallocate(line)
      deallocate(Fold)
      deallocate(wstar)
      deallocate(cache1)
      deallocate(cache2)
      deallocate(Fstar)
      deallocate(wnext)
      deallocate(a1)
      deallocate(a2)
      deallocate(b1)
      deallocate(b2)
      deallocate(c1)
      deallocate(c2)
      deallocate(d1)
      deallocate(d2)

      CONTAINS

      SUBROUTINE FORBACK(n,lower,diag,upper,b,x,band)
      implicit none
      integer i,n, band
      double precision, DIMENSION(n)::lower
      double precision, DIMENSION(n)::diag
      double precision, DIMENSION(n)::upper
      double precision, DIMENSION(n)::b
      double precision, DIMENSION(n)::x
      double precision, DIMENSION(n)::ci
      double precision, DIMENSION(n)::di

      DO I=1,band
         ci(i)=upper(i)/diag(i);
         di(i)=b(i)/diag(i);
      END DO

      DO I=band+1,n-1
         ci(i)=upper(i)/(diag(i)-ci(i-band)*lower(i));
      END DO

      DO I=band+1,n
         di(i)=(b(i)-di(i-band)*lower(i))/(diag(i)-ci(i-band)*lower(i));
      END DO

      DO I=1,band
         x(n-i+1)=di(n-i+1);
      END DO
      
      DO I=n-band,1,-1
         x(i)=di(i)-ci(i)*x(i+band);
      END DO

      END SUBROUTINE FORBACK

      FUNCTION FUNC(n,vec) RESULT(z)
      implicit none
      double precision vec(:)
      double precision z(n)
      double precision u
      integer n,i

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(u)
!$OMP& SHARED(n,z,vec)
      DO i=1,n
         u=vec(i)
         z(i)=-u/(1.0d0+u);
      END DO
!$OMP END PARALLEL DO

      END FUNCTION FUNC

      SUBROUTINE LINSPACE(vec,lower,upper,n)
      IMPLICIT NONE
      double precision, DIMENSION(n)::vec
      double precision lower, upper, h
      integer i, n

      h=1.0d0/(n+2-1.0d0)
      vec(1)=lower
      DO i=2,n-1
         vec(i)=vec(i-1)+h
      END DO
      vec(n)=upper
      RETURN
      END SUBROUTINE LINSPACE
      
      END PROGRAM ENZYMEKINETICS2D
