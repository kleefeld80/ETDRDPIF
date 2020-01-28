      ! export OMP_NUM_THREADS=3
      ! gfortran -c -O3 fftpack5.f90 
      ! gfortran ginzburg2dopt.f fftpack5.o -O3 -fopenmp -march=native -ffast-math
      PROGRAM GINZBURG2D
!$    use omp_lib
      implicit none

      ! Parameters of the problem
      integer n                                     ! Number of spatial points
      integer m                                     ! Number of temporal points
      double precision eps1                         ! Diffusion coefficients
      double precision, allocatable::x(:)           ! x-grid
      double complex, allocatable::wold(:)          ! solution w=[u,v]
      double complex, allocatable::wold2(:) 
      double complex, allocatable::M1col(:)
      double complex, allocatable::M2col(:)
      double complex, allocatable::M3col(:)
      double complex, allocatable::M11col(:)
      double complex, allocatable::M22col(:)
      double complex, allocatable::M33col(:)
      double complex, allocatable::Fold(:)
      double complex, allocatable::Fold2(:)
      double complex, allocatable::wstar(:)
      double complex, allocatable::wstar2(:)
      double complex, allocatable::cache1(:)
      double complex, allocatable::cache2(:)
      double complex, allocatable::Fstar(:)
      double complex, allocatable::wnext(:)
      double complex, allocatable::a1(:)
      double complex, allocatable::a2(:)
      double complex, allocatable::b1(:)
      double complex, allocatable::b2(:)
      double complex, allocatable::c1(:)
      double complex, allocatable::c2(:)
      double complex, allocatable::d1(:)
      double complex, allocatable::d2(:)
      double precision, allocatable::work(:)
      double precision, allocatable::wsave(:)
      double precision, allocatable::work2(:)
      double precision, allocatable::wsave2(:)
      double precision diff1
      double precision h                            ! spatial step size
      double precision dt                           ! temporal step size
      double precision starttime,finaltime
      integer nn                                    ! n x n
      integer lowband,highband,dimprob
      integer i,j,ii,ier,lensav,lenwrk,lensav2,lenwrk2

      n=400    !200
      m=2001   !2001    ! change number here
      lowband=1
      highband=n
      eps1=1.0d0
      h=200.0d0/n ! only interior nodes
      dt=100.0d0/(m-1.0d0)
      nn=n*n;
      dimprob=nn
      
      write(*,*) 'Spatial step size h=', h
      write(*,*) 'Temporal step size dt=', dt
      write(*,*) 'Dimension of the system nn=', nn

!$    starttime=omp_get_wtime();
      allocate(x(n))
      CALL LINSPACE(x,h,1.0d0-h,n)
      
      ! Initial condition for u and v
      allocate(wold(dimprob))
      allocate(wold2(dimprob))
      j=1
      DO ii=1,n
         DO i=1,n
            wold(j)=exp(-((x(i)-50.0d0)**2.0d0+
     *                  (x(ii)- 50.0d0)**2.0d0)/1000.0d0)
     *             -exp(-((x(i)-100.0d0)**2.0d0
     *                  +(x(ii)-100.0d0)**2.0d0)/1000.0d0)
     *             +exp(-((x(i)-100.0d0)**2.0d0
     *                  +(x(ii)- 50.0d0)**2.0d0)/1000.0d0)
            j=j+1
         END DO
      END DO
      
      diff1=eps1*dt/(h*h)

      allocate(Fold(dimprob))
      allocate(Fold2(dimprob))
      allocate(wstar(dimprob))
      allocate(wstar2(dimprob))
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

      lensav=2*n+int(log(real(n,kind = 8 )))+4
      lenwrk=2*n
      allocate(wsave(1:lensav))
      allocate(work(1:lenwrk))
      call zfft1i(n,wsave,lensav,ier)
      lensav2=2*dimprob+int(log(real(dimprob,kind = 8 )))+4
      lenwrk2=2*dimprob
      allocate(wsave2(1:lensav2))
      allocate(work2(1:lenwrk2))
      call zfft1i(dimprob,wsave2,lensav2,ier)
      
      allocate(M1col(n))
      M1col=0.0d0
      M1col(1)=1.0d0+diff1*2.0d0
      M1col(2)=-diff1
      M1col(n)=-diff1
      call zfft1f(n,1,M1col,n,wsave,lensav,work,lenwrk,ier)
      
      allocate(M2col(n))
      M2col=0.0d0
      M2col(1)=1.0d0+diff1*2.0d0/3.0d0
      M2col(2)=-diff1/3.0d0
      M2col(n)=-diff1/3.0d0
      call zfft1f(n,1,M2col,n,wsave,lensav,work,lenwrk,ier)
      
      allocate(M3col(n))
      M3col=0.0d0
      M3col(1)=1.0d0+diff1*2.0d0/4.0d0
      M3col(2)=-diff1/4.0d0
      M3col(n)=-diff1/4.0d0
      call zfft1f(n,1,M3col,n,wsave,lensav,work,lenwrk,ier)
      
      allocate(M11col(dimprob))
      M11col=0.0d0
      M11col(1)=1.0d0+2.0d0*diff1
      M11col(n+1)=-diff1
      M11col(n*(n-1)+1)=-diff1
      call zfft1f(dimprob,1,M11col,dimprob,wsave2,
     *     lensav2,work2,lenwrk2,ier)
      
      allocate(M22col(dimprob))
      M22col=0.0d0
      M22col(1)=1.0d0+2.0d0*diff1/3.0d0
      M22col(n+1)=-diff1/3.0d0
      M22col(n*(n-1)+1)=-diff1/3.0d0
      call zfft1f(dimprob,1,M22col,dimprob,wsave2,
     *     lensav2,work2,lenwrk2,ier)
      
      allocate(M33col(dimprob))
      M33col=0.0d0
      M33col(1)=1.0d0+diff1*2.0d0/4.0d0
      M33col(n+1)=-diff1/4.0d0
      M33col(n*(n-1)+1)=-diff1/4.0d0
      call zfft1f(dimprob,1,M33col,dimprob,wsave2,
     *     lensav2,work2,lenwrk2,ier)
      
      DO i=2,m
         Fold=FUNC(dimprob,wold)
!$OMP PARALLEL DO DEFAULT(none) SHARED(wold,wold2,Fold2,dt,Fold,wnext)
!$OMP& SHARED(dimprob) 
         DO J=1,dimprob
            Fold2(j)=Fold(j)
            wold2(j)=wold(j)
            wnext(j)=wold(j)+dt*Fold(j)
         END DO
!$OMP END PARALLEL DO
         
         ! First block starts
!$OMP PARALLEL DEFAULT(None)
!$OMP& FIRSTPRIVATE(dimprob,lensav2,lenwrk2,wsave2,work2)
!$OMP& FIRSTPRIVATE(n,lensav,lenwrk,wsave,work)
!$OMP& SHARED(wnext,wstar,wold,a1,M2col,wold2,M3col,b1,M11col)
!$OMP& PRIVATE(ier,j,II)         
!$OMP SECTIONS
!$OMP SECTION
         call zfft1f(dimprob,1,wnext,dimprob,wsave2,lensav2,
     *        work2,lenwrk2,ier)
         
         wstar=wnext/(dimprob*M11col)
         call zfft1b(dimprob,1,wstar,dimprob,wsave2,lensav2,
     *        work2,lenwrk2,ier)
!$OMP SECTION
         J=1
         DO II=1,n
            call zfft1f(n,1,wold(J:J+n-1),n,wsave,lensav,
     *           work,lenwrk,ier)
            a1(J:J+n-1)=wold(J:J+n-1)/(n*M2col)
            call zfft1b(n,1,a1(J:J+n-1),n,wsave,lensav,
     *           work,lenwrk,ier)
             J=J+n
         END DO
!$OMP SECTION
         J=1
         DO II=1,n
            call zfft1f(n,1,wold2(J:J+n-1),n,wsave,lensav,
     *           work,lenwrk,ier)
            b1(J:J+n-1)=wold2(J:J+n-1)/(n*M3col)
            call zfft1b(n,1,b1(J:J+n-1),n,wsave,lensav,
     *           work,lenwrk,ier)
            J=J+n
         END DO
!$OMP END SECTIONS
!$OMP END PARALLEL
 
!$OMP PARALLEL DO DEFAULT(none) SHARED(a1,b1,c1,dimprob)
         DO J=1,dimprob
            c1(j)=9.0d0*a1(j)-8.0d0*b1(j)
         END DO
!$OMP END PARALLEL DO

         ! second block starts
!$OMP PARALLEL DEFAULT(none)
!$OMP& FIRSTPRIVATE(n,lensav,lenwrk,wsave,work)
!$OMP& SHARED(wstar,wstar2,Fold,a2,M1col,Fold2,M2col,b2,M3col)
!$OMP& PRIVATE(ier,j,II)  
!$OMP SECTIONS
!$OMP SECTION
         J=1
         DO II=1,n
            call zfft1f(n,1,wstar(J:J+n-1),n,wsave,lensav,
     *           work,lenwrk,ier)
            wstar2(J:J+n-1)=wstar(J:J+n-1)/(n*M1col)
            call zfft1b(n,1,wstar2(J:J+n-1),n,wsave,lensav,
     *           work,lenwrk,ier)
             J=J+n
         END DO
!$OMP SECTION
         J=1
         DO II=1,n
            call zfft1f(n,1,Fold(J:J+n-1),n,wsave,lensav,
     *           work,lenwrk,ier)
            a2(J:J+n-1)=Fold(J:J+n-1)/(n*M2col)
            call zfft1b(n,1,a2(J:J+n-1),n,wsave,lensav,
     *           work,lenwrk,ier)
             J=J+n
         END DO
!$OMP SECTION
         J=1
         DO II=1,n
            call zfft1f(n,1,Fold2(J:J+n-1),n,wsave,lensav,
     *           work,lenwrk,ier)
            b2(J:J+n-1)=Fold2(J:J+n-1)/(n*M3col)
            call zfft1b(n,1,b2(J:J+n-1),n,wsave,lensav,
     *           work,lenwrk,ier)
            J=J+n
         END DO
!$OMP END SECTIONS
!$OMP END PARALLEL
         
!$OMP PARALLEL DO DEFAULT(none) SHARED(a2,b2,c2,dimprob)
         DO J=1,dimprob
            c2(j)=9.0d0*a2(j)-8.0d0*b2(j)
         END DO
!$OMP END PARALLEL DO

         Fstar=FUNC(dimprob,wstar2)

!$OMP PARALLEL DO DEFAULT(none) SHARED(c1,c2,Fstar,dt,cache1,cache2)
!$OMP& SHARED(dimprob)         
         DO J=1,dimprob
            cache1(j)= 9.0d0*c1(j)+2.0d0*dt*c2(j)+dt*Fstar(j)
            cache2(j)=-8.0d0*c1(j)-1.5d0*dt*c2(j)-dt*Fstar(j)/2.0d0
         END DO
!$OMP END PARALLEL DO

         ! Third block starts
!$OMP PARALLEL DEFAULT(none)
!$OMP& FIRSTPRIVATE(dimprob,lensav2,lenwrk2,wsave2,work2)
!$OMP& SHARED(cache1,cache2,d1,d2,M22col,M33col)
!$OMP& PRIVATE(ier,j,II) 
!$OMP SECTIONS
!$OMP SECTION
         call zfft1f(dimprob,1,cache1,dimprob,wsave2,lensav2,
     *        work2,lenwrk2,ier)
         
         d1=cache1/(dimprob*M22col)
         call zfft1b(dimprob,1,d1,dimprob,wsave2,lensav2,
     *        work2,lenwrk2,ier)
!$OMP SECTION
         call zfft1f(dimprob,1,cache2,dimprob,wsave2,lensav2,
     *        work2,lenwrk2,ier)
         
         d2=cache2/(dimprob*M33col)
         call zfft1b(dimprob,1,d2,dimprob,wsave2,lensav2,
     *        work2,lenwrk2,ier)
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
      open(20,file='solutionginzburg2d.txt')
      DO I=1,dimprob
         write(20,*) real(wold(i)), imag(wold(i))
      END DO
      close(20)

      deallocate(x)
      deallocate(wold)
      deallocate(wold2)
      deallocate(Fold)
      deallocate(Fold2)
      deallocate(wstar)
      deallocate(wstar2)
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
      deallocate(work)
      deallocate(wsave)
      deallocate(work2)
      deallocate(wsave2)
      deallocate(M1col)
      deallocate(M2col)
      deallocate(M3col)
      deallocate(M11col)
      deallocate(M22col)
      deallocate(M33col)
      CONTAINS

      FUNCTION FUNC(n,vec) RESULT(z)
      implicit none
      double complex vec(:)
      double complex z(n)
      double complex u
      double precision absu
      integer n,i

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(u)
!$OMP& SHARED(n,z,vec,absu)
      DO i=1,n
         u=vec(i)
         absu=abs(u)
         z(i)=u-dcmplx(1.0d0,1.3d0)*u*absu*absu;
      END DO
!$OMP END PARALLEL DO

      END FUNCTION FUNC

      SUBROUTINE LINSPACE(vec,lower,upper,n)
      IMPLICIT NONE
      double precision, DIMENSION(n)::vec
      double precision lower, upper, h
      integer i, n

      h=200.0d0/n
      vec(1)=0.0d0
      DO i=2,n
         vec(i)=vec(i-1)+h
      END DO
      
      RETURN
      END SUBROUTINE LINSPACE
      
      END PROGRAM GINZBURG2D
