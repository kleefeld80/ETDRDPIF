# Ginzburg in 3d, single-threaded
# check pythons precision on your system
# sys.float_info.dig == 15  -> double precision

from math import exp
import time
from numpy.fft import fft, ifft
import numpy as np
import sys


def ginzburg3d(n, m):
    # for comparison, have a look at 'ginzburg3dopt2.f'
    t_start = time.time()
    
    eps1 = 1
    h = 200 / n
    dt = 100 / (m - 1)
    dimprob = n ** 3
    
    print("Spatial step size h =    ", h)
    print("Temporal step size dt =  ", dt)
    print("Dimension of the system n^3 = ", n ** 3)
    
    x = linspace(h, n)
    wold = [exp(-((x[i] - 50) ** 2 + (x[ii] - 50) ** 2 + (x[iii] - 50) ** 2) / 1000)
            - exp(-((x[i] - 100) ** 2 + (x[ii] - 100) ** 2 + (x[iii] - 100) ** 2) / 1000)
            for iii in range(n) for ii in range(n) for i in range(n)]
    wold = np.array(wold, dtype=complex)
    
    diff1 = eps1 * dt / (h ** 2)
    
    m1col = [0] * n
    m1col[0] = 1 + 2 * diff1
    m1col[1] = -diff1
    m1col[n - 1] = -diff1
    m1col = fft(m1col) / n
    
    m2col = [0] * n
    m2col[0] = 1 + 2 / 3 * diff1
    m2col[1] = -diff1 / 3
    m2col[n - 1] = -diff1 / 3
    m2col = fft(m2col)
    
    m3col = [0] * n
    m3col[0] = 1 + diff1 / 2
    m3col[1] = -diff1 / 4
    m3col[n - 1] = -diff1 / 4
    m3col = fft(m3col) / n
    
    m11col = [0] * n ** 2
    m11col[0] = 1 + 2 * diff1
    m11col[n] = -diff1
    m11col[n * (n - 1)] = -diff1
    m11col = fft(m11col) / n ** 2
    
    m22col = [0] * n ** 2
    m22col[0] = 1 + 2 / 3 * diff1
    m22col[n] = -diff1 / 3
    m22col[n * (n - 1)] = -diff1 / 3
    m22col = fft(m22col) / n ** 2
    
    m33col = [0] * n ** 2
    m33col[0] = 1 + diff1 / 2
    m33col[n] = -diff1 / 4
    m33col[n * (n - 1)] = -diff1 / 4
    m33col = fft(m33col) / n ** 2
    
    m111col = [0] * dimprob
    m111col[0] = 1 + 2 * diff1
    m111col[n ** 2] = -diff1
    m111col[n ** 2 * (n - 1)] = -diff1
    m111col = fft(m111col) / dimprob
    
    m222col = [0] * dimprob
    m222col[0] = 1 + 2 / 3 * diff1
    m222col[n ** 2] = -diff1 / 3
    m222col[n ** 2 * (n - 1)] = -diff1 / 3
    m222col = fft(m222col) / dimprob
    
    m333col = [0] * dimprob
    m333col[0] = 1 + diff1 / 2
    m333col[n ** 2] = -diff1 / 4
    m333col[n ** 2 * (n - 1)] = -diff1 / 4
    m333col = fft(m333col) / dimprob
    
    t_prep = time.time()
    
    a1 = np.zeros(dimprob, dtype=complex)
    a2 = np.zeros(dimprob, dtype=complex)
    a3 = np.zeros(dimprob, dtype=complex)
    a4 = np.zeros(dimprob, dtype=complex)
    b1 = np.zeros(dimprob, dtype=complex)
    b2 = np.zeros(dimprob, dtype=complex)
    b3 = np.zeros(dimprob, dtype=complex)
    b4 = np.zeros(dimprob, dtype=complex)
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
        wstar = wnext / m111col
        wstar = ifft(wstar)
        
        j = 0
        for _ in range(n ** 2):
            wold[j:j + n] = fft(wold[j:j + n]) / n
            a1[j:j + n] = wold[j:j + n] / m2col
            a1[j:j + n] = ifft(a1[j:j + n]) * n
            j += n
        
        j = 0
        for _ in range(n ** 2):
            wold2[j:j + n] = fft(wold2[j:j + n]) / n
            b1[j:j + n] = wold2[j:j + n] / m3col
            b1[j:j + n] = ifft(b1[j:j + n])
            j += n
        
        # Step 4
        c1 = 9 * a1 - 8 * b1
        c1c = np.copy(c1)
        
        # Step 5
        j = 0
        for _ in range(n):
            wstar[j:j + n ** 2] = fft(wstar[j:j + n ** 2]) / n ** 2
            wstar2[j:j + n ** 2] = wstar[j:j + n ** 2] / m11col
            wstar2[j:j + n ** 2] = ifft(wstar2[j:j + n ** 2])
            j += n ** 2
        
        j = 0
        for _ in range(n ** 2):
            fold[j:j + n] = fft(fold[j:j + n]) / n
            a2[j:j + n] = fold[j:j + n] / m2col
            a2[j:j + n] = ifft(a2[j:j + n]) * n
            j += n
        
        j = 0
        for _ in range(n ** 2):
            fold2[j:j + n] = fft(fold2[j:j + n]) / n
            b2[j:j + n] = fold2[j:j + n] / m3col
            b2[j:j + n] = ifft(b2[j:j + n])
            j += n
        
        # Step 6
        c2 = 9 * a2 - 8 * b2
        c2c = np.copy(c2)
        
        # Step 7
        j = 0
        for _ in range(n ** 2):
            wstar2[j:j + n] = fft(wstar2[j:j + n]) / n
            wstar[j:j + n] = wstar2[j:j + n] / m1col
            wstar[j:j + n] = ifft(wstar[j:j + n])
            j += n
        
        j = 0
        for _ in range(n):
            c1[j:j + n ** 2] = fft(c1[j:j + n ** 2]) / n ** 2
            a3[j:j + n ** 2] = c1[j:j + n ** 2] / m22col
            a3[j:j + n ** 2] = ifft(a3[j:j + n ** 2])
            j += n ** 2
        
        j = 0
        for _ in range(n):
            c1c[j:j + n ** 2] = fft(c1c[j:j + n ** 2]) / n ** 2
            b3[j:j + n ** 2] = c1c[j:j + n ** 2] / m33col
            b3[j:j + n ** 2] = ifft(b3[j:j + n ** 2])
            j += n ** 2
        
        # Step 8
        c3 = 9 * a3 - 8 * b3
        
        # Step 9
        j = 0
        for _ in range(n):
            c2[j:j + n ** 2] = fft(c2[j:j + n ** 2]) / n ** 2
            a4[j:j + n ** 2] = c2[j:j + n ** 2] / m22col
            a4[j:j + n ** 2] = ifft(a4[j:j + n ** 2])
            j += n ** 2
        
        j = 0
        for _ in range(n):
            c2c[j:j + n ** 2] = fft(c2c[j:j + n ** 2]) / n ** 2
            b4[j:j + n ** 2] = c2c[j:j + n ** 2] / m33col
            b4[j:j + n ** 2] = ifft(b4[j:j + n ** 2])
            j += n ** 2
        
        # Step 10
        c4 = 9 * a4 - 8 * b4
        
        # Step 11
        fstar = func(wstar)
        
        # Step 12
        cache1 = 9 * c3 + 2 * dt * c4 + dt * fstar
        cache2 = -8 * c3 - 1.5 * dt * c4 - dt / 2 * fstar
        
        # Step 13
        cache1 = fft(cache1) / dimprob
        d1 = cache1 / m222col
        d1 = ifft(d1)
        
        cache2 = fft(cache2) / dimprob
        d2 = cache2 / m333col
        d2 = ifft(d2)
        
        # Step 14
        wold = d1 + d2
        # end for
    t_end = time.time()
    print("#" * 27)
    print('Total time:          {:5.3f}s'.format(t_end - t_start))
    print('Preparation time:    {:5.3f}s'.format(t_prep - t_start))
    print('Time of the loop:    {:5.3f}s'.format(t_end - t_prep))
    print("#" * 27)
    
    with open("pysolutionginzburg3d.txt", "w") as out:
        for val in wold:
            out.write(str(format(val.real, '.15e')) + "\t" + str(format(val.imag, '.15e')) + '\n')


def func(vec):
    return np.array([vec[i] - (1 + 1.3j) * vec[i] * abs(vec[i]) ** 2 for i in range(len(vec))], dtype=complex)


def linspace(h, n):
    return [h * i for i in range(n)]


if __name__ == "__main__":
    try:
        n = int(sys.argv[1])
        m = int(sys.argv[2])
        ginzburg3d(n, m)
    except:
        print("unknown parameters, try: 'python3 ginzburg3d.py <n> <m>' (without <>, n,m are integers)")
    # ginzburg3d(int(input("n: ")), int(input("m: ")))
