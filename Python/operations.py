def forback(lower, diag, upper, b, band):
    assert (len(lower) == len(diag) and len(diag) == len(upper))
    n = len(lower)  # not global n !
    ci = [0] * n
    di = [0] * n
    x = [0] * n
    
    for i in range(0, band):
        ci[i] = upper[i] / diag[i]
        di[i] = b[i] / diag[i]
    for i in range(band, n - band + 1):
        ci[i] = upper[i] / (diag[i] - ci[i - band] * lower[i])
    for i in range(band, n):
        di[i] = (b[i] - di[i - band] * lower[i]) / (diag[i] - ci[i - band] * lower[i])
    for i in range(n - 1, n - 1 - band, -1):  # [n; n-band]
        x[i] = di[i]
    for i in range(n - 1 - band, -1, -1):  # [n-band-1; 0]
        x[i] = di[i] - ci[i] * x[i + band]
    return x


def stretch(vec, dim=2):
    """
    vec: one-dimensional list of length n
    return: one-dimensional list of length n^dim

    equal to 'RESHAPE(SPREAD(vec, 1, len(vec)), (/len(vec)**2/))' in Fortran
    example: stretch([1,2,3,4])
    return: [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4]
    """
    return [vec[i] for i in range(len(vec)) for _ in range(len(vec) ** (dim - 1))]


def duplicate(vec, dim=2):
    """
    vec: one-dimensional list of length n
    return: one-dimensional list of length n^dim

    equal to 'RESHAPE(TRANSPOSE(SPREAD(vec, 1, len(vec))), (/len(vec)**2/))' in Fortran
    example: duplicate([1,2,3,4])
    return: [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4]
    """
    return vec * len(vec) ** (dim - 1)


def dupstretch(vec):
    return [vec[i] for i in range(len(vec)) for _ in range(len(vec))] * len(vec)


def scalar(fac, vec):
    """
    fac: float
    vec: one-dimensional list of length n
    return: one-dimensional list of length n

    multiplicates vector with scalar
    """
    return [fac * vec[i] for i in range(len(vec))]


def add(val, vec):
    """
    val: float
    vec: one-dimensional list of length n
    return: one-dimensional list of length n
    """
    return [val + vec[i] for i in range(len(vec))]


def addV(vec1, vec2, vec3=None):
    """
    vec1: one-dimensional list of length n
    vec2: one-dimensional list of length n
    vec3: one-dimensional list of length n, default None
    """
    if vec3:
        assert (len(vec1) == len(vec2) and len(vec2) == len(vec3))
        return [vec1[i] + vec2[i] + vec3[i] for i in range(len(vec1))]
    assert (len(vec1) == len(vec2))
    return [vec1[i] + vec2[i] for i in range(len(vec1))]


def subV(vec1, vec2, vec3=None):
    """
    vec1: one-dimensional list of length n
    vec2: one-dimensional list of length n
    vec3: one-dimensional list of length n, default None
    """
    if vec3:
        assert (len(vec1) == len(vec2) and len(vec2) == len(vec3))
        return [vec1[i] - vec2[i] - vec3[i] for i in range(len(vec1))]
    assert (len(vec1) == len(vec2))
    return [vec1[i] - vec2[i] for i in range(len(vec1))]
