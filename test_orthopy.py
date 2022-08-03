import itertools
import orthopy
import sympy

x = sympy.Symbol("x")

evaluator = orthopy.c1.legendre.Eval(x, "classical")
for val in itertools.islice(evaluator, 5):
     print(sympy.expand(val))


