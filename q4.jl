using LinearAlgebra
using SymPy

@syms t::Real
D = Differential(t)

x1 = exp(2t) * (cos(2t) - sin(2t))
x2 = exp(2t) * (sin(2t) + cos(2t))
D(x2)

f = -sin(2t)
g = cos(2t)

subs(f + g, t, π / 2)

A = [5 -1 1; 1 3.0 0; -3 2 1]
A