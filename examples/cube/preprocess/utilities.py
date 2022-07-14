def prime_factors(n):
    factor = 2
    while factor <= n:
        if n % factor == 0:
            n /= factor
            yield factor
        else:
            factor += 1


def even_repartition(n_cpu, dim):
    partition = [1] * dim

    for (i, prime) in enumerate(reversed(list(prime_factors(n_cpu)))):
        partition[i % len(partition)] *= prime

    return partition


def cubic_root(x):
    return int(round(x**(1/3))) if x >= 0 else -cubic_root(-x)


def is_exact_cube(x):
    return cubic_root(x)**3 == x or (cubic_root(x) + 1)**3 == x


def doubling_range(start=1):
    while True:
        yield start
        start *= 2
