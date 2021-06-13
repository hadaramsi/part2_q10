import math

import sp as sp

x = sp.symbols('x')
f = (x * math.exp(-x) - x + math.ln(x ** 2)) * (2 * x ** 2 + 2 * x ** 2 - 3 * x - 5)
start = 0
end = 1.5
