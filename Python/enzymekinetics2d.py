# Enzymekinetics in 2d, single-threaded
# check pythons precision on your system
# sys.float_info.dig == 15  -> double precision

from operations import *
from multiprocessing import Pool
import math
import time
import sys


def enzymekinetics(n, m):
    # for comparison, have a look at 'enzymekinetics2dopt.f'
    t_start = time.time()
    
    lowband = 1
    highband = n
    eps1 = 0.2
    h = 1 / (n + 1)
    dt = 1 / (m - 1)  # only interior nodes
    pi = math.pi
    
    print("Spatial step size h =    ", h)
    print("Temporal step size dt =  ", dt)
    print("Dimension of the system n^2 = ", n ** 2)
    
    x = linspace(0, 1, n)
    # wold = [math.sin(pi * x[j]) * math.sin(pi * x[i]) for i in range(0, n) for j in range(0, n)]
    wold = [1] * n ** 2
    
    diff1 = eps1 * dt / (h ** 2)
    
    line = [2 for _ in range(0, n)]
    a1diag = scalar(diff1, stretch(line))
    a2diag = a1diag
    
    line = [-1 for _ in range(0, n)]
    line[0] = 0
    a1lower = scalar(diff1, duplicate(line))
    a2lower = scalar(diff1, stretch(line))
    
    line[0] = -1
    line[n - 1] = 0
    a1upper = scalar(diff1, duplicate(line))
    a2upper = scalar(diff1, stretch(line))
    
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
    # Algorithm
    for i in range(2, m + 1):
        print("Iteration:           {}/{}".format(i - 1, m - 1))
        # Step 1
        fold = func(wold)
        
        # Step 2
        wnext = addV(wold, scalar(dt, fold))
        
        # Step 3
        wstar = forback(m11lower, m11diag, m11upper, wnext, highband)
        a1 = forback(m2lower, m2diag, m2upper, wold, lowband)
        b1 = forback(m3lower, m3diag, m3upper, wold, lowband)
        
        # Step 4
        c1 = subV(scalar(9, a1), scalar(8, b1))
        
        # Step 5
        wstar = forback(m1lower, m1diag, m1upper, wstar, lowband)
        a2 = forback(m2lower, m2diag, m2upper, fold, lowband)
        b2 = forback(m3lower, m3diag, m3upper, fold, lowband)
        
        # Step 6
        c2 = subV(scalar(9, a2), scalar(8, b2))
        
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
    t_end = time.time()
    print("#" * 27)
    print('Total time:          {:5.3f}s'.format(t_end - t_start))
    print('Preparation time:    {:5.3f}s'.format(t_prep - t_start))
    print('Time of the loop:    {:5.3f}s'.format(t_end - t_prep))
    print("#" * 27)
    
    with open("pysolutionenzyme2d.txt.txt", "w") as out:
        for val in wold:
            out.write(str(format(val, '.15e')) + '\n')
    
    return t_end - t_start


def enzymekinetics_parallel(n, m):
    t_start = time.time()
    
    # Preparation
    lowband = 1
    highband = n
    eps1 = 0.2
    h = 1 / (n + 1)
    dt = 1 / (m - 1)  # only interior nodes
    pi = math.pi
    
    print("Spatial step size h =    ", h)
    print("Temporal step size dt =  ", dt)
    print("Dimension of the system nn = ", n ** 2)
    
    x = linspace(0, 1, n)
    wold = [math.sin(pi * x[j]) * math.sin(pi * x[i]) for i in range(0, n) for j in range(0, n)]
    
    diff1 = eps1 * dt / (h ** 2)  # float
    
    line = [2 for _ in range(0, n)]
    a1diag = scalar(diff1, stretch(line))
    a2diag = a1diag
    
    line = [-1 for _ in range(0, n)]
    line[0] = 0
    a1lower = scalar(diff1, duplicate(line))
    a2lower = scalar(diff1, stretch(line))
    
    line[0] = -1
    line[n - 1] = 0
    a1upper = scalar(diff1, duplicate(line))
    a2upper = scalar(diff1, stretch(line))
    
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
    
    p = Pool(processes=3)
    for i in range(2, m + 1):
        print("Iteration:           {}/{}".format(i - 1, m - 1))
        # Step 1
        fold = func(wold)
        
        # Step 2
        wnext = addV(wold, scalar(dt, fold))
        
        # Step 3
        (wstar, a1, b1) = p.starmap(forback, [(m11lower, m11diag, m11upper, wnext, highband),
                                              (m2lower, m2diag, m2upper, wold[:], lowband),
                                              (m3lower, m3diag, m3upper, wold[:], lowband)])
        
        # Step 4
        c1 = subV(scalar(9, a1), scalar(8, b1))
        
        # Step 5
        (wstar, a2, b2) = p.starmap(forback, [(m1lower, m1diag, m1upper, wstar, lowband),
                                              (m2lower, m2diag, m2upper, fold[:], lowband),
                                              (m3lower, m3diag, m3upper, fold[:], lowband)])
        
        # Step 6
        c2 = subV(scalar(9, a2), scalar(8, b2))
        
        # Step 7
        fstar = func(wstar)
        
        # Step 8
        (f1, f2) = p.starmap(cache, [(c1[:], 9, c2[:], 2 * dt, fstar[:], dt),
                                     (c1[:], -8, c2[:], -3 / 2 * dt, fstar[:], -dt / 2)])
        
        # Step 9
        (d1, d2) = p.starmap(forback,
                             [(m22lower, m22diag, m22upper, f1, highband),
                              (m33lower, m33diag, m33upper, f2, highband)])
        
        # Step 10
        wold = addV(d1, d2)
        # end for
    t_end = time.time()
    p.close()
    print("#" * 27)
    print('Total time:          {:5.3f}s'.format(t_end - t_start))
    print('Preparation time:    {:5.3f}s'.format(t_prep - t_start))
    print('Time of the loop:    {:5.3f}s'.format(t_end - t_prep))
    print("#" * 27)
    
    with open("pysolutionenzyme2d.txt", "w") as out:
        for val in wold:
            out.write(str(format(val, '.15e')) + '\n')
    
    return t_end - t_start


def func(vec):
    return [-vec[i] / (1 + vec[i]) for i in range(len(vec))]


def cache(arr1, coeff1, arr2, coeff2, arr3, coeff3):
    return addV(scalar(coeff1, arr1), scalar(coeff2, arr2), scalar(coeff3, arr3))


def linspace(lower, upper, n):
    return [lower + (upper - lower) / (n + 1) * i for i in range(1, n + 1)]


if __name__ == "__main__":
    try:
        n = int(sys.argv[1])
        m = int(sys.argv[2])
        enzymekinetics(n, m)
    except:
        print("unknown parameters, try: 'python3 enzymekinetics2d.py <n> <m>' (without <>, n,m are integers)")
    # enzymekinetics_parallel(int(input("n: ")), int(input("m: ")))
