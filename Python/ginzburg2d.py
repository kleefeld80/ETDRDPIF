# Ginzburg in 2d, single-threaded
# check pythons precision on your system
# sys.float_info.dig == 15  -> double precision

from math import exp
import time
from numpy.fft import fft, ifft
import numpy as np
import sys


def ginzburg2d(n, m):
    # for comparison, have a look at 'ginzburg2d.f'
    t_start = time.time()
    
    eps1 = 1  # diffusion coefficients
    h = 200 / n  # spatial step size # only interior nodes
    dt = 100 / (m - 1)  # temporal step size
    dimprob = n ** 2
    
    print("Spatial step size h =    ", h)
    print("Temporal step size dt =  ", dt)
    print("Dimension of the system n^2 = ", n ** 2)
    
    x = linspace(h, n)  # x-grid
    
    # initial condition for u and v
    wold = [
        exp(- ((x[b] - 50) ** 2 + (x[a] - 50) ** 2) / 1000)
        - exp(-((x[b] - 100) ** 2 + (x[a] - 100) ** 2) / 1000)
        + exp(-((x[b] - 100) ** 2 + (x[a] - 50) ** 2) / 1000)
        for a in range(n) for b in range(n)]
    wold = np.array(wold, dtype=complex)  # solution w=[u,v]
    
    diff1 = eps1 * dt / (h ** 2)
    
    m1col = [0] * n
    m1col[0] = 1 + 2 * diff1
    m1col[1] = -diff1
    m1col[n - 1] = -diff1
    m1col = fft(m1col)
    
    m2col = [0] * n
    m2col[0] = 1 + 2 / 3 * diff1
    m2col[1] = -diff1 / 3
    m2col[n - 1] = -diff1 / 3
    m2col = fft(m2col)
    
    m3col = [0] * n
    m3col[0] = 1 + 2 / 4 * diff1
    m3col[1] = -diff1 / 4
    m3col[n - 1] = -diff1 / 4
    m3col = fft(m3col)
    
    # dimprob = n**2
    m11col = [0] * dimprob
    m11col[0] = 1 + 2 * diff1
    m11col[n] = -diff1
    m11col[n * (n - 1)] = -diff1
    m11col = fft(m11col)
    
    m22col = [0] * dimprob
    m22col[0] = 1 + 2 / 3 * diff1
    m22col[n] = -diff1 / 3
    m22col[n * (n - 1)] = -diff1 / 3
    m22col = fft(m22col)
    
    m33col = [0] * dimprob
    m33col[0] = 1 + diff1 / 2
    m33col[n] = -diff1 / 4
    m33col[n * (n - 1)] = -diff1 / 4
    m33col = fft(m33col)
    
    t_prep = time.time()
    
    a1 = np.zeros(dimprob, dtype=complex)
    a2 = np.zeros(dimprob, dtype=complex)
    b1 = np.zeros(dimprob, dtype=complex)
    b2 = np.zeros(dimprob, dtype=complex)
    wstar2 = np.zeros(dimprob, dtype=complex)
    for i in range(2, m + 1):
        print("Iteration:           {}/{}".format(i - 1, m - 1))
        # Step 1
        fold = func(wold)

        # Step 2
        fold2 = np.copy(fold)
        wold2 = np.copy(wold)
        wnext = wold + dt * fold

        # Step 3
        wnext = fft(wnext) / dimprob
        wstar = wnext / m11col
        wstar = ifft(wstar) * dimprob
        
        # Step 4
        j = 0
        for _ in range(n):
            wold[j:j + n] = fft(wold[j:j + n]) / n
            a1[j:j + n] = wold[j:j + n] / m2col
            a1[j:j + n] = ifft(a1[j:j + n]) * n
            j += n
        
        j = 0
        for _ in range(n):
            wold2[j:j + n] = fft(wold2[j:j + n]) / n
            b1[j:j + n] = wold2[j:j + n] / m3col
            b1[j:j + n] = ifft(b1[j:j + n]) * n
            j += n
        
        # Step 5
        c1 = 9 * a1 - 8 * b1

        # Step 6
        j = 0
        for _ in range(n):
            wstar[j:j + n] = fft(wstar[j:j + n]) / n
            wstar2[j:j + n] = wstar[j:j + n] / m1col
            wstar2[j:j + n] = ifft(wstar2[j:j + n]) * n
            j += n

        j = 0
        for _ in range(n):
            fold[j:j + n] = fft(fold[j:j + n]) / n
            a2[j:j + n] = fold[j:j + n] / m2col
            a2[j:j + n] = ifft(a2[j:j + n]) * n
            j += n

        j = 0
        for _ in range(n):
            fold2[j:j + n] = fft(fold2[j:j + n]) / n
            b2[j:j + n] = fold2[j:j + n] / m3col
            b2[j:j + n] = ifft(b2[j:j + n]) * n
            j += n

        # Step 7
        c2 = 9 * a2 - 8 * b2

        # Step 8
        fstar = func(wstar2)

        # Step 9
        cache1 = 9 * c1 + 2 * dt * c2 + dt * fstar
        cache2 = -8 * c1 - 1.5 * dt * c2 - dt / 2 * fstar
        
        # Step 10
        cache1 = fft(cache1) / dimprob
        d1 = cache1 / m22col
        d1 = ifft(d1) * dimprob
        
        cache2 = fft(cache2) / dimprob
        d2 = cache2 / m33col
        d2 = ifft(d2) * dimprob
        
        # Step 11
        wold = d1 + d2
        # end for
    t_end = time.time()
    print("#" * 27)
    print('Total time:          {:5.3f}s'.format(t_end - t_start))
    print('Preparation time:    {:5.3f}s'.format(t_prep - t_start))
    print('Time of the loop:    {:5.3f}s'.format(t_end - t_prep))
    print("#" * 27)
    
    with open("pysolutionginzburg2d.txt", "w") as out:
        for val in wold:
            out.write(str(format(val.real, '.15e')) + "\t" + str(format(val.imag, '.15e')) + '\n')


def func(vec):
    return np.array([vec[i] - (1.0 + 1.3j) * vec[i] * (np.absolute(vec[i]) ** 2) for i in range(len(vec))], dtype=complex)


def linspace(h, n):
    return [h * i for i in range(n)]


if __name__ == "__main__":
    try:
        n = int(sys.argv[1])
        m = int(sys.argv[2])
        ginzburg2d(n, m)
    except:
        print("unknown parameters, try: 'python3 ginzburg2d.py <n> <m>' (without <>, n,m are integers)")
    # ginzburg2d(int(input("n: ")), int(input("m: ")))
