# Brusselator in 3d, single-threaded
# check pythons precision on your system
# sys.float_info.dig == 15  -> double precision

from operations import *
from math import pi, sin
import time
import sys


def bruss3d(n, m):
    # for comparison, have a look at 'bruss3dopt.f'
    t_start = time.time()
    
    lowband = 2
    highband = 2 * n
    highestband = 2 * n ** 2
    eps1 = 0.02
    eps2 = 0.02
    h = 1 / (n - 1)
    dt = 5 / (m - 1)
    dimprob = 2 * n ** 3
    
    print("Spatial step size h =    ", h)
    print("Temporal step size dt =  ", dt)
    print("Dimension of the system n^3 = ", n ** 3)
    
    x = linspace(0, 1, n)
    
    pi2 = 2 * pi
    wold = [0] * dimprob
    
    i = 0
    for a in range(0, n):
        for b in range(0, n):
            for c in range(0, n):
                wold[i] = 1 + sin(pi2 * x[c]) * sin(pi2 * x[b]) * sin(pi2 * x[a])
                wold[i + 1] = 3
                i += 2
    
    diff1 = eps1 * dt / (h ** 2)
    diff2 = eps2 * dt / (h ** 2)
    
    a1lower = [0] * dimprob
    a2lower = [0] * dimprob
    a3lower = [0] * dimprob
    
    a1diag = [0] * dimprob
    # a2diag = [0] * dimprob    initialization not needed
    # a3diag = [0] * dimprob
    
    a1upper = [0] * dimprob
    a2upper = [0] * dimprob
    a3upper = [0] * dimprob
    
    line = [2 for _ in range(n)]
    longline = stretch(line, 3)
    j = 0
    for i in range(n ** 3):
        a1diag[j] = longline[i] * diff1
        a1diag[j + 1] = longline[i] * diff2
        j += 2
    a2diag = a1diag[:]
    a3diag = a1diag[:]
    
    line = [-1 for _ in range(n)]
    
    line[0] = 0
    line[n - 1] = -2
    longline = duplicate(line, 3)
    
    j = 0
    for i in range(0, n ** 3):
        a1lower[j] = longline[i] * diff1
        a1lower[j + 1] = longline[i] * diff2
        j += 2
    
    line[0] = -2
    line[n - 1] = 0
    longline = duplicate(line, 3)
    
    j = 0
    for i in range(0, n ** 3):
        a1upper[j] = longline[i] * diff1
        a1upper[j + 1] = longline[i] * diff2
        j += 2
    
    line[0] = 0
    line[n - 1] = -2
    longline = dupstretch(line)
    
    j = 0
    for i in range(0, n ** 3):
        a2lower[j] = longline[i] * diff1
        a2lower[j + 1] = longline[i] * diff2
        j += 2
    
    line[0] = -2
    line[n - 1] = 0
    longline = dupstretch(line)
    
    j = 0
    for i in range(0, n ** 3):
        a2upper[j] = longline[i] * diff1
        a2upper[j + 1] = longline[i] * diff2
        j += 2
    
    line[0] = 0
    line[n - 1] = -2
    longline = stretch(line, 3)
    
    j = 0
    for i in range(0, n ** 3):
        a3lower[j] = longline[i] * diff1
        a3lower[j + 1] = longline[i] * diff2
        j += 2
    
    line[0] = -2
    line[n - 1] = 0
    longline = stretch(line, 3)
    
    j = 0
    for i in range(0, n ** 3):
        a3upper[j] = longline[i] * diff1
        a3upper[j + 1] = longline[i] * diff2
        j += 2
    
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
    
    m111lower = tuple(a3lower)
    m111diag = tuple(add(1, a3diag))
    m111upper = tuple(a3upper)
    m222lower = tuple(scalar(1 / 3, a3lower))
    m222diag = tuple(add(1, scalar(1 / 3, a3diag)))
    m222upper = tuple(scalar(1 / 3, a3upper))
    m333lower = tuple(scalar(1 / 4, a3lower))
    m333diag = tuple(add(1, scalar(1 / 4, a3diag)))
    m333upper = tuple(scalar(1 / 4, a3upper))
    
    t_prep = time.time()
    
    for i in range(2, m + 1):
        print("Iteration:           {}/{}".format(i - 1, m))
        # Step 1
        fold = func(wold)
        
        # Step 2
        wnext = addV(wold, scalar(dt, fold))
        
        # Step 3
        wstar = forback(m111lower, m111diag, m111upper, wnext, highestband)
        a1 = forback(m2lower, m2diag, m2upper, wold, lowband)
        b1 = forback(m3lower, m3diag, m3upper, wold, lowband)
        
        # Step 4
        c1 = addV(scalar(9, a1), scalar(-8, b1))
        
        # Step 5
        wstar = forback(m11lower, m11diag, m11upper, wstar, highband)
        a2 = forback(m2lower, m2diag, m2upper, fold, lowband)
        b2 = forback(m3lower, m3diag, m3upper, fold, lowband)
        
        # Step 6
        c2 = addV(scalar(9, a2), scalar(-8, b2))
        
        # Step 7
        wstar = forback(m1lower, m1diag, m1upper, wstar, lowband)
        a3 = forback(m22lower, m22diag, m22upper, c1, highband)
        b3 = forback(m33lower, m33diag, m33upper, c1, highband)
        
        # Step 8
        c3 = addV(scalar(9, a3), scalar(-8, b3))
        
        # Step 9
        a4 = forback(m22lower, m22diag, m22upper, c2, highband)
        b4 = forback(m33lower, m33diag, m33upper, c2, highband)
        
        # Step 10
        c4 = addV(scalar(9, a4), scalar(-8, b4))
        
        # Step 11
        fstar = func(wstar)
        
        # Step 12
        f1 = addV(scalar(9, c3), scalar(2 * dt, c4), scalar(dt, fstar))
        f2 = addV(scalar(-8, c3), scalar(-1.5 * dt, c4), scalar(-dt / 2, fstar))
        
        # Step 13
        d1 = forback(m222lower, m222diag, m222upper, f1, highestband)
        d2 = forback(m333lower, m333diag, m333upper, f2, highestband)
        
        # Step 14
        wold = addV(d1, d2)
        # end for
    t_end = time.time()
    
    print("#" * 27)
    print('Total time:          {:5.3f}s'.format(t_end - t_start))
    print('Preparation time:    {:5.3f}s'.format(t_prep - t_start))
    print('Time of the loop:    {:5.3f}s'.format(t_end - t_prep))
    print("#" * 27)
    
    with open("pysolutionbruss3d.txt", "w") as out:
        for i in range(0, len(wold), 2):
            out.write(str(format(wold[i], '.15e')) + "\t\t" + str(format(wold[i + 1], '.15e')) + '\n')


def func(vec):
    a = 2
    b = 1
    z = [0] * len(vec)
    for i in range(0, len(vec) - 1, 2):
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
        bruss3d(n, m)
    except:
        print("unknown parameters, try: 'python3 bruss3d.py <n> <m>' (without <>, n,m are integers)")
    # bruss3d(int(input("n: ")), int(input("m: ")))
