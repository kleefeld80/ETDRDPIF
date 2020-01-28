% Create Table 2 third column
% First  row, use n= 99, m=101
% Second row, use n=199, m=201
% Third  row, use n=399, m=401
% Fourth row, use n=799, m=801
export OMP_NUM_THREADS=1
gfortran enzymekinetics2dopt.f -O3 -fopenmp -march=native -ffast-math -o enz.out
enz.out

% Create Table 2 fourth column
% First  row, use n= 99, m=101
% Second row, use n=199, m=201
% Third  row, use n=399, m=401
% Fourth row, use n=799, m=801
export OMP_NUM_THREADS=3
gfortran enzymekinetics2dopt.f -O3 -fopenmp -march=native -ffast-math -o enz.out
enz.out

% Create numbers for single precision
% First  number, use n= 99, m=101
% Second number, use n=199, m=201
% Third  number, use n=399, m=401
% Fourth number, use n=799, m=801
export OMP_NUM_THREADS=3
gfortran enzymekinetics2doptReal.f -O3 -fopenmp -march=native -ffast-math -o enz.out
enz.out

% Create Table 4 third column
% First  row, use n=101, m= 201
% Second row, use n=201, m= 401
% Third  row, use n=401, m= 801
% Fourth row, use n=801, m=1601
export OMP_NUM_THREADS=1
gfortran bruss2dopt.f -O3 -fopenmp -march=native -ffast-math -o bru.out
bru.out

% Create Table 4 fourth column
% First  row, use n=101, m= 201
% Second row, use n=201, m= 401
% Third  row, use n=401, m= 801
% Fourth row, use n=801, m=1601
export OMP_NUM_THREADS=3
gfortran bruss2dopt.f -O3 -fopenmp -march=native -ffast-math -o bru.out
bru.out

% 3D Brusselator parallelized version
export OMP_NUM_THREADS=1
gfortran bruss3dopt.f -O3 -fopenmp -march=native -ffast-math -o bru.out
bru.out
export OMP_NUM_THREADS=3
gfortran bruss3dopt.f -O3 -fopenmp -march=native -ffast-math -o bru.out
bru.out

% 2D Ginzburg (timings for Table 5)
% First  row, use n= 50, m= 2001
% Second row, use n=100, m= 2001
% Third  row, use n=200, m= 2001
% Fourth row, use n=400, m= 2001
gfortran -c -O3 fftpack5.f90
gfortran ginzburg2dopt.f fftpack5.o -O3 -ffast-math -march=native -fopenmp -o ginz.out
export OMP_NUM_THREADS=1
ginz.out
export OMP_NUM_THREADS=3
ginz.out

% 3D Ginzburg
gfortran -c -O3 fftpack5.f90
gfortran ginzburg3dopt.f fftpack5.o -O3 -ffast-math -march=native -fopenmp -o ginz.out
export OMP_NUM_THREADS=1
ginz.out
export OMP_NUM_THREADS=3
ginz.out

% 3D Ginzburg
gfortran -c -O3 fftpack5.f90 
gfortran ginzburg3dopt.f fftpack5.o -O3 -fopenmp -march=native -ffast-math -o ginz.out
export OMP_NUM_THREADS=1
ginz.out
export OMP_NUM_THREADS=3
ginz.out

% Timings for Table 6
% First  row, use n= 50, m= 2001
% Second row, use n=100, m= 2001
% Third  row, use n=200, m= 2001

gfortran -c -O3 fftpack5.f90 
gfortran ginzburg3dopt2.f fftpack5.o -O3 -fopenmp -march=native -ffast-math -o ginz.out
export OMP_NUM_THREADS=1
ginz.out
export OMP_NUM_THREADS=3
ginz.out

