"""Compare the ICK polygon kernel with closed-form exits from the pi/8 triangle."""
import hashlib
import itertools
import math
import subprocess
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path


def coefficients(rational, radical=0):
    return Fraction(rational), Fraction(radical)


def add(first, second):
    return first[0] + second[0], first[1] + second[1]


def subtract(first, second):
    return first[0] - second[0], first[1] - second[1]


def multiply(first, second):
    return (first[0] * second[0] + 2 * first[1] * second[1],
            first[0] * second[1] + first[1] * second[0])


def divide(first, second):
    norm = second[0] ** 2 - 2 * second[1] ** 2
    real, radical = multiply(first, (second[0], -second[1]))
    return real / norm, radical / norm


# Dyadic bounds distinguish every nonzero value in this finite corpus.
# A bound containing zero is a test failure, never a tie.
DENOMINATOR = 1 << 64
LOWER = Fraction(math.isqrt(2 * DENOMINATOR ** 2), DENOMINATOR)
UPPER = LOWER + Fraction(1, DENOMINATOR)


def sign(value):
    rational, radical = value
    if rational == radical == 0:
        return 0
    left = rational + radical * LOWER
    right = rational + radical * UPPER
    if left * right <= 0:
        raise ValueError(f"sign unresolved for {value}")
    return 1 if left > 0 else -1


def expected_exit(x, y, gap, dx, dy):
    """Take the earliest outward crossing of the three bounding half-planes."""
    candidates = []
    if dy < 0:
        candidates.append((0, divide(y, coefficients(-dy))))
    if dx > 0:
        candidates.append((1, divide(subtract(coefficients(1), x), coefficients(dx))))
    approach = subtract(coefficients(dy), multiply(coefficients(-1, 1), coefficients(dx)))
    if sign(approach) > 0:
        candidates.append((2, divide(gap, approach)))
    earliest = [(side, time) for side, time in candidates
                if all(sign(subtract(other, time)) >= 0 for _, other in candidates)]
    side, time = earliest[0]
    end_x = add(x, multiply(time, coefficients(dx)))
    end_y = add(y, multiply(time, coefficients(dy)))
    if len(earliest) == 2:
        pair = tuple(sorted(record[0] for record in earliest))
        vertex = {(0, 2): 0, (0, 1): 1, (1, 2): 2}[pair]
        return 2, vertex, (time, end_x, end_y, coefficients(0), coefficients(0))
    assert len(earliest) == 1
    reflections = ((coefficients(dx), coefficients(-dy)),
                   (coefficients(-dx), coefficients(dy)),
                   (coefficients(0, Fraction(dx + dy, 2)), coefficients(0, Fraction(dx - dy, 2))))
    return 1, side, (time, end_x, end_y, *reflections[side])


def decode(row):
    fields = list(map(int, row.split()))
    if len(fields) != 17:
        raise ValueError(f"wrong output field count: {row!r}")
    scalars = []
    for rational, radical, denominator in zip(fields[2::3], fields[3::3], fields[4::3]):
        if denominator <= 0 or math.gcd(rational, radical, denominator) != 1:
            raise ValueError(f"noncanonical coefficients: {row!r}")
        scalars.append((Fraction(rational, denominator), Fraction(radical, denominator)))
    return fields[0], fields[1], tuple(scalars)


def main(executable):
    inputs, expected = [], []
    for column, row in itertools.product(range(1, 16), range(1, 7)):
        x, y = coefficients(Fraction(column, 16)), coefficients(Fraction(row, 16))
        gap = subtract(multiply(x, coefficients(-1, 1)), y)
        if sign(gap) <= 0:
            continue
        for dx, dy in itertools.product(range(-4, 5), repeat=2):
            if dx == dy == 0:
                continue
            inputs.append(f"{column} {row} {dx} {dy}\n")
            expected.append(expected_exit(x, y, gap, dx, dy))
    result = subprocess.run([executable, "--trace-eighths"], input="".join(inputs).encode(),
                            stdout=subprocess.PIPE, check=True)
    rows = result.stdout.splitlines()
    if len(rows) != len(expected) or len(rows) != 3360:
        raise ValueError(f"expected {len(expected)} grid rows, got {len(rows)}")
    counts = Counter()
    for ray, wanted, line in zip(inputs, expected, rows):
        observed = decode(line)
        if observed != wanted:
            raise ValueError(f"{ray.strip()} produced {observed}; expected {wanted}")
        counts[observed[0]] += 1
    print(f"pi/8 oracle: {counts[1]} regular reflections and {counts[2]} vertex stops PASS")
    print("pi/8 trace SHA256", hashlib.sha256(result.stdout).hexdigest())


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
