      ! export OMP_NUM_THREADS=3
      ! gfortran -c -O3 fftpack5.f90 
      ! gfortran ginzburg3dopt.f fftpack5.o -O3 -fopenmp -march=native -ffast-math
      PROGRAM GINZBURG3D
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
      double complex, allocatable::M111col(:)
      double complex, allocatable::M222col(:)
      double complex, allocatable::M333col(:)
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
      double complex, allocatable::a3(:)
      double complex, allocatable::a4(:)
      double complex, allocatable::b1(:)
      double complex, allocatable::b2(:)
      double complex, allocatable::b3(:)
      double complex, allocatable::b4(:)
      double complex, allocatable::c1(:)
      double complex, allocatable::c1c(:)
      double complex, allocatable::c2(:)
      double complex, allocatable::c2c(:)
      double complex, allocatable::c3(:)
      double complex, allocatable::c4(:)
      double complex, allocatable::d1(:)
      double complex, allocatable::d2(:)
      double precision, allocatable::work(:)
      double precision, allocatable::wsave(:)
      double precision, allocatable::work2(:)
      double precision, allocatable::wsave2(:)
      double precision, allocatable::work3(:)
      double precision, allocatable::wsave3(:)
      double precision diff1
      double precision h                            ! spatial step size
      double precision dt                           ! temporal step size
      double precision starttime,finaltime
      integer nn                                    ! n x n
      integer lowband,highband,dimprob,nnn
      integer i,j,ii,iii,ier,lensav,lenwrk,lensav2,lenwrk2
      integer lensav3,lenwrk3

      n=50 !256
      m=2001   !301    ! change number here
      lowband=1
      highband=n
      eps1=1.0d0
      h=100.0d0/n ! only interior nodes
      dt=100.0d0/(m-1.0d0)
      nn=n*n;
      nnn=n*n*n
      dimprob=nnn
      
      write(*,*) 'Spatial step size h=', h
      write(*,*) 'Temporal step size dt=', dt
      write(*,*) 'Dimension of the system nnn=', nnn

!$    starttime=omp_get_wtime();
      allocate(x(n))
      CALL LINSPACE(x,h,1.0d0-h,n)

      ! Initial condition for u and v
      allocate(wold(dimprob))
      allocate(wold2(dimprob))
      j=1
      DO iii=1,n
         DO ii=1,n
            DO i=1,n
               wold(j)=exp(-((x(i)  -50.0d0)**2.0d0+
     *                       (x(ii) -50.0d0)**2.0d0+
     *                       (x(iii)-50.0d0)**2.0d0)/1000.0d0)
               j=j+1
            END DO
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
      allocate(a3(dimprob))
      allocate(a4(dimprob))
      allocate(b1(dimprob))
      allocate(b2(dimprob))
      allocate(b3(dimprob))
      allocate(b4(dimprob))
      allocate(c1(dimprob))
      allocate(c1c(dimprob))
      allocate(c2(dimprob))
      allocate(c2c(dimprob))
      allocate(c3(dimprob))
      allocate(c4(dimprob))
      allocate(d1(dimprob))
      allocate(d2(dimprob))

      lensav=2*n+int(log(real(n,kind = 8 )))+4
      lenwrk=2*n
      allocate(wsave(1:lensav))
      allocate(work(1:lenwrk))
      call zfft1i(n,wsave,lensav,ier)
      
      lensav2=2*n*n+int(log(real(n*n,kind = 8 )))+4
      lenwrk2=2*n*n
      allocate(wsave2(1:lensav2))
      allocate(work2(1:lenwrk2))
      call zfft1i(n*n,wsave2,lensav2,ier)

      lensav3=2*dimprob+int(log(real(dimprob,kind = 8 )))+4
      lenwrk3=2*dimprob
      allocate(wsave3(1:lensav3))
      allocate(work3(1:lenwrk3))
      call zfft1i(dimprob,wsave3,lensav3,ier)
      
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
      
      allocate(M11col(n*n))
      M11col=0.0d0
      M11col(1)=1.0d0+2.0d0*diff1
      M11col(n+1)=-diff1
      M11col(n*(n-1)+1)=-diff1
      call zfft1f(n*n,1,M11col,n*n,wsave2,
     *     lensav2,work2,lenwrk2,ier)
      
      allocate(M22col(n*n))
      M22col=0.0d0
      M22col(1)=1.0d0+2.0d0*diff1/3.0d0
      M22col(n+1)=-diff1/3.0d0
      M22col(n*(n-1)+1)=-diff1/3.0d0
      call zfft1f(n*n,1,M22col,n*n,wsave2,
     *     lensav2,work2,lenwrk2,ier)
      
      allocate(M33col(n*n))
      M33col=0.0d0
      M33col(1)=1.0d0+diff1*2.0d0/4.0d0
      M33col(n+1)=-diff1/4.0d0
      M33col(n*(n-1)+1)=-diff1/4.0d0
      call zfft1f(n*n,1,M33col,n*n,wsave2,
     *     lensav2,work2,lenwrk2,ier)

      allocate(M111col(dimprob))
      M111col=0.0d0
      M111col(1)=1.0d0+2.0d0*diff1
      M111col(n*n+1)=-diff1
      M111col(n*n*(n-1)+1)=-diff1
      call zfft1f(dimprob,1,M111col,dimprob,wsave3,
     *     lensav3,work3,lenwrk3,ier)
      
      allocate(M222col(dimprob))
      M222col=0.0d0
      M222col(1)=1.0d0+2.0d0*diff1/3.0d0
      M222col(n*n+1)=-diff1/3.0d0
      M222col(n*n*(n-1)+1)=-diff1/3.0d0
      call zfft1f(dimprob,1,M222col,dimprob,wsave3,
     *     lensav3,work3,lenwrk3,ier)
      
      allocate(M333col(dimprob))
      M333col=0.0d0
      M333col(1)=1.0d0+diff1*2.0d0/4.0d0
      M333col(n*n+1)=-diff1/4.0d0
      M333col(n*n*(n-1)+1)=-diff1/4.0d0
      call zfft1f(dimprob,1,M333col,dimprob,wsave3,
     *     lensav3,work3,lenwrk3,ier)

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
!$OMP& FIRSTPRIVATE(dimprob,lensav3,lenwrk3,wsave3,work3)
!$OMP& FIRSTPRIVATE(n,lensav,lenwrk,wsave,work)
!$OMP& SHARED(wnext,wstar,wold,a1,M2col,wold2,M3col,b1,M111col)
!$OMP& PRIVATE(ier,j,II)         
!$OMP SECTIONS
!$OMP SECTION
         call zfft1f(dimprob,1,wnext,dimprob,wsave3,lensav3,
     *        work3,lenwrk3,ier)
         
         wstar=wnext/(dimprob*M111col)
         call zfft1b(dimprob,1,wstar,dimprob,wsave3,lensav3,
     *        work3,lenwrk3,ier)
!$OMP SECTION
         J=1
         DO II=1,n*n
            call zfft1f(n,1,wold(J:J+n-1),n,wsave,lensav,
     *           work,lenwrk,ier)
            a1(J:J+n-1)=wold(J:J+n-1)/(n*M2col)
            call zfft1b(n,1,a1(J:J+n-1),n,wsave,lensav,
     *           work,lenwrk,ier)
             J=J+n
         END DO
!$OMP SECTION
         J=1
         DO II=1,n*n
            call zfft1f(n,1,wold2(J:J+n-1),n,wsave,lensav,
     *           work,lenwrk,ier)
            b1(J:J+n-1)=wold2(J:J+n-1)/(n*M3col)
            call zfft1b(n,1,b1(J:J+n-1),n,wsave,lensav,
     *           work,lenwrk,ier)
            J=J+n
         END DO
!$OMP END SECTIONS
!$OMP END PARALLEL

!$OMP PARALLEL DO DEFAULT(none) SHARED(a1,b1,c1,c1c,dimprob)
         DO J=1,dimprob
            c1(j)=9.0d0*a1(j)-8.0d0*b1(j)
            c1c(j)=c1(j)
         END DO
!$OMP END PARALLEL DO

         ! second block starts
!$OMP PARALLEL DEFAULT(none)
!$OMP& FIRSTPRIVATE(n,lensav,lenwrk,wsave,work)
!$OMP& FIRSTPRIVATE(lensav2,lenwrk2,wsave2,work2)
!$OMP& SHARED(wstar,wstar2,Fold,a2,M11col,Fold2,M2col,b2,M3col)
!$OMP& PRIVATE(ier,j,II)  
!$OMP SECTIONS
!$OMP SECTION
         J=1
         DO II=1,n
            call zfft1f(n*n,1,wstar(J:J+n*n-1),n*n,wsave2,lensav2,
     *           work2,lenwrk2,ier)
            wstar2(J:J+n*n-1)=wstar(J:J+n*n-1)/(n*n*M11col)
            call zfft1b(n*n,1,wstar2(J:J+n*n-1),n*n,wsave2,lensav2,
     *           work2,lenwrk2,ier)
             J=J+n*n
         END DO
!$OMP SECTION
         J=1
         DO II=1,n*n
            call zfft1f(n,1,Fold(J:J+n-1),n,wsave,lensav,
     *           work,lenwrk,ier)
            a2(J:J+n-1)=Fold(J:J+n-1)/(n*M2col)
            call zfft1b(n,1,a2(J:J+n-1),n,wsave,lensav,
     *           work,lenwrk,ier)
             J=J+n
         END DO
!$OMP SECTION
         J=1
         DO II=1,n*n
            call zfft1f(n,1,Fold2(J:J+n-1),n,wsave,lensav,
     *           work,lenwrk,ier)
            b2(J:J+n-1)=Fold2(J:J+n-1)/(n*M3col)
            call zfft1b(n,1,b2(J:J+n-1),n,wsave,lensav,
     *           work,lenwrk,ier)
            J=J+n
         END DO
!$OMP END SECTIONS
!$OMP END PARALLEL
         
!$OMP PARALLEL DO DEFAULT(none) SHARED(a2,b2,c2,c2c,dimprob)
         DO J=1,dimprob
            c2(j)=9.0d0*a2(j)-8.0d0*b2(j)
            c2c(j)=c2(j)
         END DO
!$OMP END PARALLEL DO

! third block starts
!$OMP PARALLEL DEFAULT(none)
!$OMP& FIRSTPRIVATE(n,lensav,lenwrk,wsave,work)
!$OMP& FIRSTPRIVATE(lensav2,lenwrk2,wsave2,work2)
!$OMP& SHARED(wstar,wstar2,c1,c1c,a3,M1col,Fold2,M22col,b3,M33col)
!$OMP& PRIVATE(ier,j,II)  
!$OMP SECTIONS
!$OMP SECTION
         J=1
         DO II=1,n*n
            call zfft1f(n,1,wstar2(J:J+n-1),n,wsave,lensav,
     *           work,lenwrk,ier)
            wstar(J:J+n-1)=wstar2(J:J+n-1)/(n*M1col)
            call zfft1b(n,1,wstar(J:J+n-1),n,wsave,lensav,
     *           work,lenwrk,ier)
             J=J+n
         END DO
!$OMP SECTION
         J=1
         DO II=1,n
            call zfft1f(n*n,1,c1(J:J+n*n-1),n*n,wsave2,lensav2,
     *           work2,lenwrk2,ier)
            a3(J:J+n*n-1)=c1(J:J+n*n-1)/(n*n*M22col)
            call zfft1b(n*n,1,a3(J:J+n*n-1),n*n,wsave2,lensav2,
     *           work2,lenwrk2,ier)
             J=J+n*n
         END DO
!$OMP SECTION
         J=1
         DO II=1,n
            call zfft1f(n*n,1,c1c(J:J+n*n-1),n*n,wsave2,lensav2,
     *           work2,lenwrk2,ier)
            b3(J:J+n*n-1)=c1c(J:J+n*n-1)/(n*n*M33col)
            call zfft1b(n*n,1,b3(J:J+n*n-1),n*n,wsave2,lensav2,
     *           work2,lenwrk2,ier)
             J=J+n*n
         END DO
!$OMP END SECTIONS
!$OMP END PARALLEL

!$OMP PARALLEL DO DEFAULT(none) SHARED(a3,b3,c3)
         DO J=1,dimprob
            c3(j)=9.0d0*a3(j)-8.0d0*b3(j)
         END DO
!$OMP END PARALLEL DO

! fourth block starts
!$OMP PARALLEL DEFAULT(none)
!$OMP& FIRSTPRIVATE(n,lensav2,lenwrk2,wsave2,work2)
!$OMP& SHARED(c2,c2c,a4,M22col,b4,M33col)
!$OMP& PRIVATE(ier,j,II)  
!$OMP SECTIONS
!$OMP SECTION
         J=1
         DO II=1,n
            call zfft1f(n*n,1,c2(J:J+n*n-1),n*n,wsave2,lensav2,
     *           work2,lenwrk2,ier)
            a4(J:J+n*n-1)=c2(J:J+n*n-1)/(n*n*M22col)
            call zfft1b(n*n,1,a4(J:J+n*n-1),n*n,wsave2,lensav2,
     *           work2,lenwrk2,ier)
             J=J+n*n
         END DO
!$OMP SECTION
         J=1
         DO II=1,n
            call zfft1f(n*n,1,c2c(J:J+n*n-1),n*n,wsave2,lensav2,
     *           work2,lenwrk2,ier)
            b4(J:J+n*n-1)=c2c(J:J+n*n-1)/(n*n*M33col)
            call zfft1b(n*n,1,b4(J:J+n*n-1),n*n,wsave2,lensav2,
     *           work2,lenwrk2,ier)
             J=J+n*n
         END DO
!$OMP END SECTIONS
!$OMP END PARALLEL

!$OMP PARALLEL DO DEFAULT(none) SHARED(a4,b4,c4)
         DO J=1,dimprob
            c4(j)=9.0d0*a4(j)-8.0d0*b4(j)
         END DO
!$OMP END PARALLEL DO
         
         Fstar=FUNC(dimprob,wstar)

!$OMP PARALLEL DO DEFAULT(none) SHARED(c3,c4,Fstar,dt,cache1,cache2)
!$OMP& SHARED(dimprob)         
         DO J=1,dimprob
            cache1(j)= 9.0d0*c3(j)+2.0d0*dt*c4(j)+dt*Fstar(j)
            cache2(j)=-8.0d0*c3(j)-1.5d0*dt*c4(j)-dt*Fstar(j)/2.0d0
         END DO
!$OMP END PARALLEL DO

         ! Fifth block starts
!$OMP PARALLEL DEFAULT(none)
!$OMP& FIRSTPRIVATE(dimprob,lensav3,lenwrk3,wsave3,work3)
!$OMP& SHARED(cache1,cache2,d1,d2,M222col,M333col)
!$OMP& PRIVATE(ier,j,II) 
!$OMP SECTIONS
!$OMP SECTION
         call zfft1f(dimprob,1,cache1,dimprob,wsave3,lensav3,
     *        work3,lenwrk3,ier)
         
         d1=cache1/(dimprob*M222col)
         call zfft1b(dimprob,1,d1,dimprob,wsave3,lensav3,
     *        work3,lenwrk3,ier)
!$OMP SECTION
         call zfft1f(dimprob,1,cache2,dimprob,wsave3,lensav3,
     *        work3,lenwrk3,ier)
         
         d2=cache2/(dimprob*M333col)
         call zfft1b(dimprob,1,d2,dimprob,wsave3,lensav3,
     *        work3,lenwrk3,ier)
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
      open(20,file='solutionginzburg3d.txt')
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
      deallocate(a3)
      deallocate(a4)
      deallocate(b1)
      deallocate(b2)
      deallocate(b3)
      deallocate(b4)
      deallocate(c1)
      deallocate(c1c)
      deallocate(c2)
      deallocate(c2c)
      deallocate(c3)
      deallocate(c4)
      deallocate(d1)
      deallocate(d2)
      deallocate(work)
      deallocate(wsave)
      deallocate(work2)
      deallocate(wsave2)
      deallocate(work3)
      deallocate(wsave3)
      deallocate(M1col)
      deallocate(M2col)
      deallocate(M3col)
      deallocate(M11col)
      deallocate(M22col)
      deallocate(M33col)
      deallocate(M111col)
      deallocate(M222col)
      deallocate(M333col)
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
      
      END PROGRAM GINZBURG3D
