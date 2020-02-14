# Brusselator in 2d, single-threaded
# check pythons precision on your system
# sys.float_info.dig == 15  -> double precision

from operations import *
import time
import sys


def bruss2d(n, m):
    # for comparison, have a look at 'bruss2dopt.f'
    t_start = time.time()
    
    lowband = 2
    highband = 2 * n
    eps1 = 2 * (10 ** (-3))
    eps2 = 2 * (10 ** (-3))
    h = 1 / (n - 1)
    dt = 2 / (m - 1)
    dimprob = 2 * n ** 2
    
    print("Spatial step size h =    ", h)
    print("Temporal step size dt =  ", dt)
    print("Dimension of the system n^2 = ", n ** 2)
    
    x = linspace(0, 1, n)
    
    wold = [0] * dimprob
    j = 0
    for i in range(0, n):
        for ii in range(0, n):
            wold[j] = 0.5 + x[i]
            wold[j + 1] = 1 + 5 * x[ii]
            j += 2
    
    diff1 = eps1 * dt / (h ** 2)
    diff2 = eps2 * dt / (h ** 2)
    
    a1lower = [0] * dimprob
    a1diag = [0] * dimprob
    a1upper = [0] * dimprob
    a2lower = [0] * dimprob
    a2diag = [0] * dimprob
    a2upper = [0] * dimprob
    
    line = [2 for _ in range(n)]
    ll = stretch(line)
    for i in range(0, dimprob, 2):
        a1diag[i] = ll[i // 2] * diff1
        a1diag[i + 1] = ll[i // 2] * diff2
    a2diag = a1diag[:]
    
    line = [-1 for _ in range(n)]
    line[0] = 0
    line[n - 1] = -2
    ll = duplicate(line)
    for i in range(0, dimprob, 2):
        a1lower[i] = ll[i // 2] * diff1
        a1lower[i + 1] = ll[i // 2] * diff2
    ll = stretch(line)
    for i in range(0, dimprob, 2):
        a2lower[i] = ll[i // 2] * diff1
        a2lower[i + 1] = ll[i // 2] * diff2
    
    line[0] = -2
    line[n - 1] = 0
    ll = duplicate(line)
    for i in range(0, dimprob, 2):
        a1upper[i] = ll[i // 2] * diff1
        a1upper[i + 1] = ll[i // 2] * diff2
    ll = stretch(line)
    for i in range(0, dimprob, 2):
        a2upper[i] = ll[i // 2] * diff1
        a2upper[i + 1] = ll[i // 2] * diff2
    
    m1lower = tuple(a1lower)
    m1diag = tuple(add(1, a1diag))
    m1upper = tuple(a1upper)
    m2lower = tuple(scalar(1 / 3, a1lower))
    m2diag = tuple(add(1, scalar(1 / 3, a1diag)))
    m2upper = tuple(scalar(1 / 3, a1upper))
    m3lower = tuple(scalar(1 / 4, a1lower))
    m3diag = tuple(add(1, scalar(1 / 4, a1diag)))
    m3upper = tuple(scalar(1 / 4, a1upper))
    
    m11lower = tuple(a2lower)
    m11diag = tuple(add(1, a2diag))
    m11upper = tuple(a2upper)
    m22lower = tuple(scalar(1 / 3, a2lower))
    m22diag = tuple(add(1, scalar(1 / 3, a2diag)))
    m22upper = tuple(scalar(1 / 3, a2upper))
    m33lower = tuple(scalar(1 / 4, a2lower))
    m33diag = tuple(add(1, scalar(1 / 4, a2diag)))
    m33upper = tuple(scalar(1 / 4, a2upper))
    
    t_prep = time.time()
    wnext = [0] * dimprob
    for i in range(2, m + 1):
        print("Iteration:           {}/{}".format(i - 1, m))
        # Step 1
        fold = func(wold)
        
        # Step 2
        for j in range(0, dimprob):
            wnext[j] = wold[j] + dt * fold[j]
            
        # Step 3
        wstar = forback(m11lower, m11diag, m11upper, wnext, highband)
        a1 = forback(m2lower, m2diag, m2upper, wold, lowband)
        b1 = forback(m3lower, m3diag, m3upper, wold, lowband)
        
        # Step 4
        c1 = addV(scalar(9, a1), scalar(-8, b1))
        
        # Step 5
        wstar = forback(m1lower, m1diag, m1upper, wstar, lowband)
        a2 = forback(m2lower, m2diag, m2upper, fold, lowband)
        b2 = forback(m3lower, m3diag, m3upper, fold, lowband)
        
        # Step 6
        c2 = addV(scalar(9, a2), scalar(-8, b2))
        
        # Step 7
        fstar = func(wstar)
        
        # Step 8
        f1 = addV(scalar(9, c1), scalar(2 * dt, c2), scalar(dt, fstar))  # cache1
        f2 = addV(scalar(-8, c1), scalar(-3 / 2 * dt, c2), scalar(-dt / 2, fstar))  # cache2
        
        # Step 9
        d1 = forback(m22lower, m22diag, m22upper, f1, highband)
        d2 = forback(m33lower, m33diag, m33upper, f2, highband)
        
        # Step 10
        wold = addV(d1, d2)
        # end for
    t_end = time.time()
    print("#" * 27)
    print('Total time:          {:5.3f}s'.format(t_end - t_start))
    print('Preparation time:    {:5.3f}s'.format(t_prep - t_start))
    print('Time of the loop:    {:5.3f}s'.format(t_end - t_prep))
    print("#" * 27)
    
    with open("pysolutionbruss2d.txt", "w") as out:
        for i in range(0, len(wold), 2):
            out.write(str(format(wold[i], '.15e')) + "\t\t" + str(format(wold[i + 1], '.15e')) + '\n')
    
    return t_end - t_start


def func(vec):
    n = len(vec)
    a = 1
    b = 3.4
    z = [0] * n
    for i in range(0, n - 1, 2):
        u = vec[i]
        v = vec[i + 1]
        cache1 = u ** 2 * v
        cache2 = b * u
        z[i] = a + cache1 - cache2 - u
        z[i + 1] = cache2 - cache1
    return z


def linspace(lower, upper, n):
    return [lower + (upper - lower) / (n - 1) * i for i in range(0, n)]


if __name__ == "__main__":
    try:
        n = int(sys.argv[1])
        m = int(sys.argv[2])
        bruss2d(n, m)
    except:
        print("unknown parameters, try: 'python3 bruss2d.py <n> <m>' (without <>, n,m are integers)")
    # bruss2d(int(input("n: ")), int(input("m: ")))
