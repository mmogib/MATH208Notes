### A Pluto.jl notebook ###
# v0.19.47

using Markdown
using InteractiveUtils

# ╔═╡ f2d4c2a5-f486-407b-b31b-d2efcc7476b3
begin
    using CommonMark
    using PlutoUI, PlutoExtras
    using Plots, PlotThemes, LaTeXStrings
    using Latexify
    using HypertextLiteral
    using Colors
    using LinearAlgebra, Random, Printf, SparseArrays
    # using Symbolics
    # using SymPy
    using QRCoders
    using PrettyTables
	# using Primes
    # using LinearSolve
    # using NonlinearSolve
    # using ForwardDiff
    # using Integrals
    # using OrdinaryDiffEq

end

# ╔═╡ 71bc54d5-d0ed-42d3-9bc1-48aa86e91d1d
TableOfContents(title="📚 MATH208: Differential Equations and Linear Algebra", indent=true, depth=4)

# ╔═╡ e414122f-b93a-4510-b8ae-026c303e0df9
begin
    struct LocalImage
        filename
    end

    function Base.show(io::IO, ::MIME"image/png", w::LocalImage)
        write(io, read(w.filename))
    end
end


# ╔═╡ cd269caf-ef81-43d7-a1a8-6668932b6363
# exportqrcode("https://www.mathmatize.com/")
LocalImage("./qrcode.png")

# ╔═╡ cb0eda5b-9dcb-421f-b6b6-ed4b51797df7
md"# 1.1 Differential Equations and Mathematical Models"

# ╔═╡ 20bd193d-c0f7-47e6-8b46-33e07338ad91
cm"""
__The study of differential equations has three principal goals:__
1. To discover the differential equation that describes a specified physical situation.
2. To find-either exactly or approximately-the appropriate solution of that equation.
3. To interpret the solution that is found.
"""

# ╔═╡ 555e1ddd-fcd2-4ff8-af06-f4caafc8dcff
md"## Mathematical Models"

# ╔═╡ d5117856-b9ca-4ed1-976e-9e20d613bfc1
md"## Examples and Terminology"

# ╔═╡ 55221819-6a40-4f85-93b0-76c5d05ed35e
cm"""
 - The __order__ of a differential equation is the __order of the highest derivative that appears in it__. For example
```math
y^{(4)}+x^2 y^{(3)}+x^5 y=\sin x \quad \text{is a fourth-order equation.}
```
- The most general form of an ``\boldsymbol{n}^t h``-order differential equation with independent variable ``x`` and unknown function or dependent variable ``y=y(x)`` is
$(texeq"F\left(x, y, y^{\prime}, y^{\prime \prime}, \ldots, y^{(n)}\right)=0\label{GDE}")
where ``F`` is a specific real-valued function of ``n+2`` variables.
- We say that the continuous function ``u=u(x)`` is a __solution__ of the differential equation in $(eqref("GDE")) on the interval ``I`` provided that the derivatives ``u^{\prime}, u^{\prime \prime}, \ldots, u^{(n)}`` exist on ``I`` and
$(texeq"F\left(x, u, u^{\prime}, u^{\prime \prime}, \ldots, u^{(n)}\right)=0")
for all ``x`` in ``I``. For the sake of brevity, we may say that ``u=u(x)`` satisfies the differential equation in $(eqref("GDE")) on ``I``.

"""

# ╔═╡ a9d46b08-d25a-49df-9f93-427e521055e0
md"# 1.2 Integrals as General and Particular Solutions"

# ╔═╡ 22715277-7264-44aa-a2be-0cf6bba6b24f
cm"""
__Second-order equations.__
```math
\frac{d^2y}{dx^2} = g(x),
```
"""

# ╔═╡ 07d56dc4-4a09-43eb-8ad6-8cb796b53c80
md"## Velocity and Acceleration"

# ╔═╡ 2d42c9ba-bc48-4836-aa56-bca41304a300
cm"""
- __The position function__ ``x=f(t)``.
- __The Velocity__ ``v(t)=f'(t)``, or ``v =\displaystyle\frac{dx}{dt}``.
- __The Accelaration__ ``a(t)=v'(t)=x''(t)``. or ``\displaystyle a = \frac{dv}{dt}=\frac{d^2x}{dt^2}``.
"""

# ╔═╡ b97a6f03-dbeb-40b9-be58-f9a9a50cac5a
md"## A Swimmer’s Problem"

# ╔═╡ a1d00dd2-59b0-427b-b2c0-25ec94e039a9
let
    dydx(x) = 3(1 - 4x^2)
    x = -0.5:0.01:0.5
    # y = x .+ dydx.(x)
    # 3x-4x.^3)
    y = zeros(length(x))
    p = plot(xlimits=(-0.55, 0.55), ylimits=(-10, 221),
        frame_style=:origin, legend=:topleft,
        xticks=-0.5:0.2:0.5
    )
    vline!(p, [-0.5, 0], c=:black, label=nothing)
    vline!(p, [0.5, 0], c=:black, label=nothing)
    annotate!(p, [(0.52, 10.5, text(L"x=\frac{1}{2}", 10))])
    annotate!(p, [(-0.52, 10.5, text(L"x=-\frac{1}{2}", 10))])
    anim = @animate for i in 2:length(x)
        y[i] = y[i-1] + dydx(x[i])
        plot(p, [(x[1:i], y[1:i]), ([x[i]], [y[i]])],
            seriestype=[:line :scatter],
            label=["y:: swimmer trajectory" "swimmer"],
            markershapes=[:none :star]
        )

    end
    gif(anim, "imgs/swimmer.gif", fps=10)

end

# ╔═╡ 8781e87a-146a-41b0-8981-243fec51bfa3
function vector_field(xs, ys, df; c=:black, args...)
    # xs = -2:0.3:2
    # ys = -2:0.2:2
    # df(x, y) = normalize([1, 1/x]) ./ 10

    xxs = repeat(xs', length(xs), length(ys))
    yys = reshape(repeat(ys, length(xs)), length(xs), length(ys))

    quiver(xxs, yys, quiver=df; framestyle=:origin, c=c, args...)

end

# ╔═╡ 261637c6-4da1-4b4a-9f01-3e2324c41b60
# let
# 	xs = -2:0.3:2
# 	ys = -2:0.2:2


# 	df(x, y) = normalize([1, 1/x]) ./ 10

# 	xxs = repeat(xs',length(xs),length(ys))
# 	yys = reshape(repeat(ys,length(xs)),length(xs),length(ys))

# quiver(xxs, yys, quiver=df, c=:black,framestyle=:origin)
# end
md"Vector Feilds"

# ╔═╡ 953e31c8-3e26-47cf-a675-2067802ba941
md"# 1.4 Separable Equations and Applications"

# ╔═╡ 75287b3b-3f90-4a0a-88d0-92130f84c0db
cm"""
The first-order differential equation ``\displaystyle \frac{dy}{dx}=f(x,y)`` is called __separable__ provided that ``f(x, y)`` can be written as the product of a function of ``x`` and a function of ``y`` :
```math
\frac{d y}{d x}=f(x, y)=g(x) k(y)=\frac{g(x)}{h(y)}
```
where ``h(y)=1 / k(y)``. 

In this case the variables ``x`` and ``y`` can be __separated__(isolated) on opposite sides of an equation-by writing informally the equation
```math
h(y) d y=g(x) d x
```
"""

# ╔═╡ 1790195b-439e-4472-96e8-6f87d8ff0601
md"## Implicit, General, and Singular Solution"

# ╔═╡ 085418be-7e83-4a6a-b3a8-fbd4d2ac2451
cm"""
- The equation ``K(x,y)=0`` is called an __implicit solution__ of a differential
equation if it is satisfied (on some interval) by some solution ``y=y(x)`` of the differential equation.
- A particular solution ``y=y(x)`` of ``K(x,y)`` may or may not satisfy a given initial condition.
- Not every possible algebraic solution ``y=y(x)`` of an implicit solution ``K(x,y)`` satisfies the same differential equation.
- Similarly, solutions of a given differential equation can be either gained or lost
when it is multiplied or divided by an algebraic factor.
- A solution of a differential equation that contains an “arbitrary constant” (like
the constant C) is commonly called __the general solution__ of the differential equation; any particular choice of a specific value for ``C`` yields a single particular solution of the equation.
- __A particular solution__ is one that is obtained by selecting a value for ``C``.
- A solution that cannot be obtained by selecting a value for ``C`` is called a __singular solution__.


"""

# ╔═╡ 7bacb885-372a-4ee9-ba9c-de505d332dfb
md"## Natural Growth and Decay"

# ╔═╡ 241b2797-3050-400b-a532-303a6feeb39e
cm"""
The differential equation
$(texeq"
\frac{d x}{d t}=k x \quad(k \text { a constant }) \label{}
")
serves as a mathematical model for a wide range of natural phenomena-any
"""

# ╔═╡ a5a9fda5-c077-4505-8fcf-2e23f9b3ad1f
md"###  POPULATION GROWTH"

# ╔═╡ 809d228a-65a3-4af5-b605-e030c41c4548
cm"""
```math
\frac{dP}{dt} = k P
```
where ``P(t)`` is the number of individuals in a
 population (of humans, or insects, or bacteria) having constant birth and death rates.
"""

# ╔═╡ 6b4ce7dc-39f3-4c4e-973f-b501916f0717
md"### COMPOUND INTEREST"

# ╔═╡ 0756191d-272f-4e62-8fed-1f731820325c
cm"""
```math
\frac{d A}{d t}=\lim _{\Delta t \rightarrow 0} \frac{\Delta A}{\Delta t}=r A
```
where  
- ``A(t)`` is the number of dollars in a savings account
 at time t (in years), and
- ``r`` is the interest compounded continuously at
 an annually.
"""

# ╔═╡ f832a1e4-17ff-43d8-97a4-c7e29bbe542f
md"### RADIOACTIVE DECAY"

# ╔═╡ a2cc59fc-01d3-45c6-80e7-3c8485632036
cm"""
```math
\frac{d N}{d t}=-k N
```
where 
- ``N(t)`` the number of atoms  of  certain radio active isotope at time ``t``.

The decay constant of a radioactive isotope is often specified in terms of another empirical constant, the __half-life__ of the isotope, because this parameter is more convenient. 

- The half-life ``\tau`` of a radioactive isotope is the time required for half of it to decay. 
- To find the relationship between ``k`` and ``\tau``, we set ``t=\tau`` and ``N=\frac{1}{2} N_0`` in the equation ``N(t)=N_0 e^{-k t}``, so that ``\frac{1}{2} N_0=N_0 e^{-k \tau}``. When we solve for ``\tau``, we find that
```math
\tau=\frac{\ln 2}{k}
```

- For example, the half-life of ``{ }^{14} \mathrm{C}`` is ``\tau \approx(\ln 2) /(0.0001216)``, approximately 5700
"""

# ╔═╡ 72cef7b7-f466-4872-bba0-6f46ac36de15
let
    k = 0.0001216
    t = -log(0.63) / k
end

# ╔═╡ 1cfb4072-1e58-4b7a-bcdc-672cd0183f75
md"###  DRUG ELIMINATION"

# ╔═╡ 07b46453-655c-4fc0-8590-3104ce5f6c99
cm"""
```math
\frac{d A}{d t}=-\lambda A
```
where 
- ``A(t)`` is the amount of a certain drug in the bloodstream, measured by the excess over the natural level of the drug, will decline at a rate proportional to the current excess amount. 
- ``\lambda`` is called the elimination constant of the drug.
"""

# ╔═╡ 27f89662-ef45-4469-bbc7-7a6b1949881c
let
    P₀ = 6
    Pprime0 = 0.000212 * (365.25)
    k = Pprime0 / P₀

    P_2051 = P₀ * exp(k * 51)

    P10 = floor(log(10) / k)
    1999 + Int(P10)
end

# ╔═╡ 124fc4b5-157b-46b1-8a77-06b119ee0a4d
md"##  Cooling and Heating"

# ╔═╡ c8d24a6a-b851-420f-9c4f-e42eeb934fc3
cm"""
__Newton's law of cooling__: the time rate of change of the temperature ``T(t)`` of a body immersed in a medium of constant temperature ``A`` is proportional to the difference ``A-T``. That is,
```math
\frac{d T}{d t}=k(A-T)
```
where ``k`` is a positive constant. This is an instance of the linear first-order differential equation with constant coefficients:
```math
\frac{d x}{d t}=a x+b
```

It includes the exponential equation as a special case ``(b=0)`` and is also easy to solve by separation of variables.
"""

# ╔═╡ c96bcf16-8467-4f8f-a55a-6769f73143da
let
    A = 375
    C = 375 - 50
    # (A-125)
    k = -log((A - 125) / 325) / 75
    T = 150
    t = -log((A - T) / C) / k
    # 105/60
    # (0.75*60)
end

# ╔═╡ e9d05d41-2d51-4f11-ad43-95571b99b466
md"# 1.5 Linear First-Order Equations"

# ╔═╡ e1c8a21a-8f19-4632-992f-014d4690baa6
cm"""
A __linear first-order equation__ is a differential equation of the form
```math
\frac{d y}{d x}+P(x) y=Q(x) \text {. }
```
"""

# ╔═╡ 16d641fc-590f-4ab2-9288-c0f03849caec
cm"""
multiplying by the __integrating factor__
```math
\rho(x) = e^{\int P(x)dx}
```

Solves the problem.

##### METHOD: SOLUTION OF LINEAR FIRST-ORDER EQUATIONS
1. Begin by calculating the integrating factor ``\rho(x)=e^{\int P(x) d x}``.
2. Then multiply both sides of the differential equation by ``\rho(x)``.
3. Next, recognize the left-hand side of the resulting equation as the derivative of a product:
```math
D_x[\rho(x) y(x)]=\rho(x) Q(x)
```
4. Finally, integrate this equation,
```math
\rho(x) y(x)=\int \rho(x) Q(x) d x+C
```
then solve for ``y`` to obtain the general solution of the original differential equation.

"""

# ╔═╡ ae81e505-f95a-4355-bf74-b308445d2830
let
    xs = -1.5:0.6:6
    ys = -4.5:0.6:2.5
    df(x, y) = normalize([1, y + (11 / 8) * exp(-x / 3)])
    p1 = vector_field(xs, ys, df)
    p1 = plot(p1, x -> (1 / 32) * (exp(x) - 33 * exp(-x / 3)), c=:blue,
        xlims=(min(xs...), max(xs...)),
        ylims=(min(ys...), max(ys...)),
        label=:none
    )
    plot(p1, [0], [-1], seriestype=:scatter, label=:none)
end

# ╔═╡ 5aa66e9f-5741-4306-aa36-8a28627bca0f
md"# 1.6 Substitution Methods and Exact Equations"

# ╔═╡ 419f6c6f-1aff-4b72-be83-3dc3ef01bc76
md"##  Homogeneous Equations"

# ╔═╡ cacc6877-5310-4f5a-a0c2-344cccd0ec8b
md"## Bernoulli Equations"

# ╔═╡ 564258db-41ca-41fd-90c8-13c8fd39fe93
cm"""
A first-order differential equation of the form
```math
\frac{d y}{d x}+P(x) y=Q(x) y^n
```
is called a __Bernoulli equation__. If either ``n=0`` or ``n=1``, then this equation is linear. Otherwise, we make the substitution
```math
v=y^{1-n}
```
transforming it the linear equation
```math
\frac{d v}{d x}+(1-n) P(x) v=(1-n) Q(x)
```
"""

# ╔═╡ 4ccc21f1-f253-4e7d-87e4-4c8f5e1db785
md"##  Exact Differential Equations"

# ╔═╡ 9ed8b8ed-1ac5-4e8b-8a74-17dff32f1ad5
let
    x = -2:0.1:2
    y = -2:0.1:2
    df(x, y) = 0.1 * normalize([1, (y^3 - 6 * x * y) / (4y + 3x^2 - 3x * y^2)])
    plt = vector_field(x, y, df; lw=0.01)
    contour!(x, y, (x, y) -> 3x^2 * y - x * y^3 + 2y^2, c=:red)
end

# ╔═╡ 80e20e84-51d2-4e6e-8eb3-b7a8dac22a8b
md"## Reducible Second-Order Equations"

# ╔═╡ 65b3d353-815c-4e5a-a887-24a82bdf005d
cm"""

A second-order differential equation involves the second derivative of the unknown function ``y(x)``, and thus has the general form
```math
F\left(x, y, y^{\prime}, y^{\prime \prime}\right)=0
```

If either the dependent variable ``y`` or the independent variable ``x`` is missing from a second-order equation, then it is easily reduced by a simple substitution to a firstorder equation that may be solvable by the methods of this chapter.
"""

# ╔═╡ 8d103106-40ac-49fa-86ae-1af3ecb8a3ec
let
    # 	@syms x::Real, y()
    # 	D = Differential(x)
    # 	Eq =2*x*exp(2*y(x))*D(y(x)) ~ 3x^4+exp(2*y(x))
    # 	dsolve(Eq)
end

# ╔═╡ c7cc5a14-f964-4cf9-82f0-f23904bbdace
let
    # @syms x(),y(), t::Real
    # ∂ = Differential(t)
    # Z =[x(t);y(t)]
    # S1 = ∂.(Z) .~ [4 1;6 -1]*Z
    # dsolve(S1)
    cm"Example ...."
end

# ╔═╡ 2678d172-8896-4681-b9af-45e3485f2312
md"# Chapter 3: Linear Systems and Matrices"

# ╔═╡ 27dbde5f-5e4d-4c67-bc67-ae0b9c40a774
md"## Gaussian Elimination"

# ╔═╡ a3567a11-b2c2-4545-b4a2-278b4b7e4b0e
cm"""
__We can use three elementary row operations only:__

1. Multiply one equation by a nonzero constant. (``\lambda R_i \to R_i``)
2. Interchange two equations. (`` R_i \leftrightarrow R_j``)
3. Add a constant multiple of (the terms of) one equation to (corresponding terms of) another equation. (``\lambda R_i +R_j\to R_j``)
"""

# ╔═╡ aa04cd10-62d5-49a9-8b60-7e5d928ec1ee
let
    # @syms a::Real, b::Real, c::Real
    # A = [2 -1 3;1 2 1;7 4 9]
    # B =[a;b;c]
    # Au=hcat(A,B)
    # Au[1,:],Au[2,:]=Au[2,:],Au[1,:]
    # Au[2,:]=-2Au[1,:]+Au[2,:]
    # Au[3,:]=-7Au[1,:]+Au[3,:]
    # Au[3,:]=-2Au[2,:]+Au[3,:]

    Au
end

# ╔═╡ 387953af-89c5-47ee-a1b1-af9b05baf949
md"## Inverse of a Matrix"

# ╔═╡ 747384a8-c06a-49b7-92f1-0cdd44043fb2
let
    A = [
        1 1 5
        1 4 13
        3 2 12
    ]
    B = inv(A)

end

# ╔═╡ a2a039ee-2097-4569-9d85-e96597ff47da
let
    A = [4 3 2; 5 6 3; 3 5 2]
    # inv(A)
    # Au = Rational.(hcat([4 3 2;5 6 3;3 5 2],I(3)))
    # Au[1,:]=-Au[3,:]+Au[1,:]
    # Au[2,:]=-5Au[1,:]+Au[2,:]
    # Au[3,:]=-3Au[1,:]+Au[3,:]
    # Au[2,:]=(1//16)Au[2,:]
    # Au[3,:]=-11Au[2,:]+Au[3,:]
    # Au[1,:]=2Au[2,:]+Au[1,:]
    # Au[3,:]=-16Au[3,:]
    # Au[1,:]=-(3//8)Au[3,:]+Au[1,:]
    # Au[2,:]=-(3//16)Au[3,:]+Au[2,:]
    # B=Au[:,4:end]
    # B*A
end

# ╔═╡ 25e77707-5852-48f8-9849-46e3f6841bb6
md"## Determinants"

# ╔═╡ f5bf2679-c765-4b80-bd27-021a82683611
md"### Cramer's Rule"

# ╔═╡ aff3151b-4b66-4d55-8e27-744a53f5a91e
md"## Row and Column Properties"

# ╔═╡ babc3fb7-0eaa-4af6-b7d1-cd0b0fd4c78b
let
    A1 = [1 2 3
        0 4 4
        0 6 7
    ]
    A2 = [1 2 3
        4 0 0
        0 6 7
    ]
    B = [1 2 3; 4 4 4; 0 6 7]
    det(B), det(A1), det(A2)
end

# ╔═╡ ea5ca7e4-2f30-45f2-8715-7e502a274a53
let
    A = [
        1 2 3 4
        0 5 6 7
        0 0 8 9
        2 4 6 9
    ]
end

# ╔═╡ e1199b37-1692-4d97-ba86-0b5239c316d1
md"## Inverses and the Adjoint Matrix"

# ╔═╡ ece946d9-b094-4e01-8dda-47740ec418b3
let
    A = [
        1 4 5
        4 2 5
        -3 3 -1
    ]
    # inv(Rational.(A))
end

# ╔═╡ 30795897-1e43-41ed-ba85-c3b89ecdef1f
md"# 4.1 The Vector Space ``\mathbb{R}^3``"

# ╔═╡ c233fecb-8525-493f-a088-011d5386ec1d
md"## The Vector Space ``\mathbb{R}^2``"

# ╔═╡ d5618ff1-3420-4c92-bcfa-238b4da78f17
cm"""
 ``\mathbb{R}^2`` is simply the set of all 3-dimensional vectors that have third
 component ``0``.
```math
\mathbb{R}^2 = \left\{(x,y,0)\in \mathbb{R}^3\;\;|\;\; x,y\in \mathbb{R}\right\}
```

__Two vectors ``\mathbf{u}`` and ``\mathbf{v}`` are collinear__-they lie on the same line through the origin and hence point either in the same direction  or in opposite directions-if and only if one is a scalar multiple of the other; that is, either
```math
\mathbf{u}=c \mathbf{v} \quad \text { or } \quad \mathbf{v}=c \mathbf{u}
```
for some scalar ``c``. If one of the relations holds for some scalar ``c``, then we say that the two vectors are __linearly dependent__. 
"""

# ╔═╡ 6a3e0f87-707c-4122-a5fd-dcef8fce8a7d
md"##  Linear Independence in ``\mathbb{R}^3``"

# ╔═╡ b7430962-79fb-427c-a23e-8cbb09ae8e46
md"##  Basis Vectors in ``\mathbb{R}^3``"

# ╔═╡ e0261d6d-2bc3-4eab-8eb2-ef799f0cece7
md"##  Subspaces of ``\mathbb{R}^3``"

# ╔═╡ 2af29504-36c2-4f20-a325-c32acf62d9bc
md"# 4.2 The Vector Space ``\mathbb{R}^n``"

# ╔═╡ 170f2255-078f-4c9e-86e9-0eb93b2ff87e
md"## Vector Spaces"

# ╔═╡ 6827e93a-43b1-4ee7-a650-573a3ce61676
md"## Subspaces"

# ╔═╡ 2f5bc349-791e-4ece-8933-f2baef3c59a1
md"# 4.3 Linear Combinations and Independence of Vectors"

# ╔═╡ ca3afed1-ded8-4792-8253-58a2b5bc0502
md"##  Linear Independence"

# ╔═╡ 4fc61c32-377d-475f-9587-fb884aaafcb6
let
    A = [1 2 3; 2 3 8; 2 4 8; 1 1 5]
    b = zeros(3)
    # A\b
    # solve()
end

# ╔═╡ 73ea96b9-be3c-4e26-8f54-914bd53d409b
let
    v1 = [2; 0; -3]
    v2 = [4, -5, -6]
    v3 = [-2, 1, 3]
    A = [v1 v2 v3]
    A[3, :] = (3 // 2)A[1, :] + A[3, :]
    A
    # det(A)
end

# ╔═╡ faf0e3cf-3973-4272-bda3-a3410c8e17d2
md"# 4.4 Bases and Dimension for Vector Spaces"

# ╔═╡ f911e13b-7830-449d-b8d2-ae6ecb388e5c
let
    A = [
        1 1 1 0
        -1 -1 -1 3
        -2 2 -3 -1
        -3 3 -2 2
    ]
    det(A)
end

# ╔═╡ dca3f82b-7cac-4849-88a2-d47805957566
md"##  Bases for Solution Spaces"

# ╔═╡ 03e906c2-29b0-487b-9223-6c7313c1d740
let
    A = [
        3 6 -1 -5 5
        2 4 -1 -3 2
        3//1 6 -2 -4 1]
    A[1, :] = A[1, :] - A[2, :]
    A[2, :] = -(A[2, 1] / A[1, 1]) * A[1, :] + A[2, :]
    A[3, :] = -(A[3, 1] / A[1, 1]) * A[1, :] + A[3, :]
    A[3, :] = -(A[3, 3] / A[2, 3]) * A[2, :] + A[3, :]

    A
end

# ╔═╡ ee16937c-617a-4f4a-a15e-c3853f0e8ca8
md"# 4.5 Row and Column Spaces"

# ╔═╡ d999ba91-2a73-4d8e-a99a-fd5187f24e6b
let
    A = [
        3 6 -1 -5 5
        2 4 -1 -3 2
        3 6 -2 -4 1]
    # A[1,:]=A[1,:]-A[2,:]
    # A[2,:]=-(A[2,1]/A[1,1])*A[1,:]+A[2,:]
    # A[3,:]=-(A[3,1]/A[1,1])*A[1,:]+A[3,:]
    # A[3,:]=-(A[3,3]/A[2,3])*A[2,:]+A[3,:]

    rank(A)
end

# ╔═╡ 28fb5089-e83b-47f5-93be-60f1f0ea1e30
md"# 5.1 Introduction: Second-Order Linear Equations"

# ╔═╡ 82a55302-3cf3-4c7e-861c-12f5f554f142
cm"""
__Homogeneous Second-Order Linear Equations__

Consider the general second-order linear equation
```math
A(x) y^{\prime \prime}+B(x) y^{\prime}+C(x) y=F(x)
```
where the coefficient functions ``A, B, C``, and ``F`` are continuous on the open interval ``I``. Here we assume in addition that ``A(x) \neq 0`` at each point of ``I``, so we can divide each term in Eq. (7) by ``A(x)`` and write it in the form

$(texeq"y^{\prime \prime}+p(x) y^{\prime}+q(x) y=f(x)")


We will discuss first the associated homogeneous equation

"""

# ╔═╡ c9ac1dc2-078e-47ae-b0cf-9f22eb294b39
texeq"y^{\prime \prime}+p(x) y^{\prime}+q(x) y=0 \label{hso}"

# ╔═╡ 98730a31-e9e9-45fc-ae12-5cf0832b5069
md"## LinearlyIndependentSolutions"

# ╔═╡ f9a9f5f5-5134-4f2f-8672-c2e665c546ce
md"## GeneralSolutions"

# ╔═╡ 06cbf42d-9811-40ae-bcd5-1bcc4ff60880
md"##  Linear Second-Order Equations with Constant Coefficients"

# ╔═╡ c60351f6-a57c-46db-a3ea-3deee1f9180b
cm"""
```math
a y^{\prime \prime}+b y^{\prime}+c y=0
```
with constant coefficients ``a, b``, and ``c``.
Because ``e^{r x}`` is never zero, we may conclude that ``y(x)=e^{r x}`` will satisfy the differential equation precisely when ``r`` is a root of the algebraic equation
$(texeq"
a r^2+b r+c=0 .
")

This quadratic equation is called the __characteristic equation__ of the homogeneous linear differential equation
$(texeq"
a y^{\prime \prime}+b y^{\prime}+c y=0 \label{eqhc}
")
"""

# ╔═╡ 381710b7-c487-4b65-9a07-88d773cdacf0
eqref("eqhc")

# ╔═╡ 90d2b2d2-c030-4afd-a529-cd940bc82688
md"# 5.2 General Solutions of Linear Equations"

# ╔═╡ 8695eb73-c264-4f4a-ab66-2ce879ff8c88
cm"""
Consider the __``n^{\text{th}}``-order linear__ differential equation
```math
P_0(x) y^{(n)}+P_1(x) y^{(n-1)}+\cdots+P_{n-1}(x) y^{\prime}+P_n(x) y=F(x).
```
which can be written as
$(texeq"\quad y^{(n)}+p_1(x) y^{(n-1)}+\cdots+p_{n-1}(x) y^{\prime}+p_n(x) y=f(x)\label{n_order_de}")

The homogeneous linear equation associated with (10) is

$(texeq"y^{(n)}+p_1(x) y^{(n-1)}+\cdots+p_{n-1}(x) y^{\prime}+p_n(x) y=0.")

"""

# ╔═╡ dc4b275f-9014-4014-bc59-03353d5de046
md"## Existence and Uniqueness of Solutions"

# ╔═╡ cc282ef0-67d4-4086-9916-49ee02ee8009
md"## Linearly Independent Solutions"

# ╔═╡ 0f267953-3b92-4550-af48-3bdc51364a26
md"## General Solutions"

# ╔═╡ 87c2749f-1390-4ce1-97f1-3a9609cba34e
md"## Nonhomogeneous Equations"

# ╔═╡ 42a14a32-71aa-4019-8eaf-e42f7f51e7a9
cm"""
We now consider the __nonhomogeneous ``n`` th-order linear__ differential equation

$(texeq"y^{(n)}+p_1(x) y^{(n-1)}+\cdots+p_{n-1}(x) y^{\prime}+p_n(x) y=f(x)")

with associated homogeneous equation

$(texeq"y^{(n)}+p_1(x) y^{(n-1)}+\cdots+p_{n-1}(x) y^{\prime}+p_n(x) y=0")

Suppose 
- ``y_p`` is a particular solution of the nonhomogeneous equation in (13)
- ``Y`` is any other solution of Eq. (13). 
- If ``y_c=Y-y_p``, then subsitution of ``y_c`` in the differential equation gives (using the linearity of differentiation)
```math
\begin{aligned}
& y_c^{(n)}+p_1 y_c^{(n-1)}+\cdots+p_{n-1} y_c^{\prime}+p_n y_c \\
&= {\left[\left(Y^{(n)}+p_1 Y^{(n-1)}+\cdots+p_{n-1} Y^{\prime}+p_n Y\right]\right.} \\
&-\left[\left(y_p^{(n)}+p_1 y_p^{(n-1)}+\cdots+p_{n-1} y_p^{\prime}+p_n y_p\right]\right. \\
&= f(x)-f(x)=0
\end{aligned}
```

Thus __``y_c=Y-y_p``__ is a solution of the associated homogeneous equation in (14). Then
```math
Y=y_c+y_p
```
and it follows from Theorem 4 that
```math
y_c=c_1 y_1+c_2 y_2+\cdots+c_n y_n
```
where ``y_1, y_2, \ldots, y_n`` are linearly independent solutions of the associated homogeneous equation. 

- We call ``y_c`` a __complementary function__ of the nonhomogeneous equation 
- a general solution of the nonhomogeneous equation in (13) is the sum of its complementary function ``y_c`` and a single particular solution ``y_p`` of Eq. (13).
"""

# ╔═╡ bcaf9ef2-6371-4de1-8254-b4ae57f77363
md"# 5.3 Homogeneous Equations with Constant Coefficients"

# ╔═╡ fb34e1b6-1cd2-4f0f-a8c2-988dd3d709d0
cm"""
We consider the equation of the form

$(texeq"
a_n y^{(n)}+a_{n-1} y^{(n-1)}+\cdots+a_2 y^{\prime \prime}+a_1 y^{\prime}+a_0 =0
\label{heqc}
")
where the coefficients ``a_0, a_1, a_2, \ldots, a_n`` are real constants with ``a_n \neq 0``.

and its corresponding __characteristic equation__ (__auxiliary equation__)
$(texeq"
a_n r^n+a_{n-1} r^{n-1}+\cdots+a_2 r^2+a_1 r+a_0=0
\label{heqc_c}
")
"""

# ╔═╡ fcd85a36-c264-43a4-8eb6-8bd025582ac2
cm"## Distinct Real Roots"

# ╔═╡ 92ebabd1-d3df-47e9-9322-5c415cf1dc7b
md"##  Polynomial Differential Operators"

# ╔═╡ f9d43fec-1397-4206-9136-adbedd8d7fd8
cm"""
We write Equation (15) as 
```math
Ly=0
```
where 
```math 
L=a_n \frac{d^n}{d x^n}+a_{n-1} \frac{d^{n-1}}{d x^{n-1}}+\cdots+a_2 \frac{d^2}{d x^2}+a_1 \frac{d}{d x}+a_0
```
called the __differential operator__ operating on ``y``.

We also use the notation ``D=\frac{d}{dx}``.

So for example
```math
Dy = y^{\prime} , \quad D^2 = y^{\prime\prime}, \cdots
```
Therefore,

```math 
L=a_n D+a_{n-1} D^{n-1}+\cdots+a_2 D^2+a_1 D+a_0
```
called __polynomial differential operator__ in ``D``.
"""

# ╔═╡ 3e22f633-ef1c-45b4-aa83-3e28c61a9365
md"## Repeated Real Roots"

# ╔═╡ a03566a1-b19b-4bfe-84e6-7850c401e299
md"## Complex-Valued Functions and Euler’s Formula"

# ╔═╡ eb6f270c-9e5b-41fa-bfd0-8244d8d811f1
cm"""
Recall __Euler Formula__
```math
e^{i\theta} = cos(\theta) + i\sin(\theta)
```
A __complex-valued function ``F``__ of the real variable ``x`` associates with each real number ``x`` (in its domain of definition) the complex number
```math
F(x)=f(x)+i g(x)
```

The real-valued functions ``f`` and ``g`` are called the __real__ and __imaginary__ parts, respectively, of ``F``. 

If they are differentiable, we define the derivative ``F^{\prime}`` of ``F`` by
```math
F^{\prime}(x)=f^{\prime}(x)+i g^{\prime}(x)
```

Thus we simply differentiate the real and imaginary parts of ``F`` separately.
"""

# ╔═╡ 96210d10-9d4f-4aea-9b76-a8e066f34e71
md"##  Complex Roots"

# ╔═╡ 44f2e591-8d5f-41a6-b0b2-66b3a7d83cf2
md"## Repeated Complex Roots

Same as repeated real Roots.
"

# ╔═╡ 08dec10f-1df7-4249-88a2-fc52be252e7e
cm"""
That is if ``a\pm ib`` has multiplicity ``k`` then
```math
\begin{gathered}
\left(A_1+A_2 x+\cdots+A_k x^{k-1}\right) e^{(a+b i) x}+\left(B_1+B_2 x+\cdots+B_k x^{k-1}\right) e^{(a-b i) x} \\
=\sum_{p=0}^{k-1} x^p e^{a x}\left(c_p \cos b x+d_p \sin b x\right)
\end{gathered}
```
"""

# ╔═╡ 9a964597-f45a-4e1e-a4a5-c6d0007b3e13
md"# 5.5 Nonhomogeneous Equations and Undetermined Coefficients"

# ╔═╡ c0461d53-1446-468f-9ae5-babc4b2bab32
cm"""
The __general nonhomogeneous ``n`` th-order linear equation with constant coefficients__ has the form

$(texeq"a_n y^{(n)}+a_{n-1} y^{(n-1)}+\cdots+a_1 y^{\prime}+a_0 y=f(x)")

A __general solution__ of this equation has the form

$(texeq"y=y_c+y_p")

where the __complementary function__ ``y_c(x)`` is a general solution of the __associated homogeneous equation__

$(texeq"a_n y^{(n)}+a_{n-1} y^{(n-1)}+\cdots+a_1 y^{\prime}+a_0 y=0")
"""

# ╔═╡ 0ee542fb-a1d2-46c8-9604-681485015268
md"## The method of undetermined coefficients "

# ╔═╡ 270a3760-25c6-4cb3-96c2-95a2b98506f9
md"## The Case of Duplication"

# ╔═╡ ae2c1de9-d67a-4c4f-a87d-90250e436fc8
cm"""

Now our problem is to find a particular solution of the equation 
```math
L y=f(x),
``` 
where ``f(x)`` is a __linear combination of products of the elementary functions listed in (above). 

Thus ``f(x)`` can be written as a sum of terms each of the form

$(texeq"P_m(x) e^{r x} \cos k x \quad \text { or } \quad P_m(x) e^{r x} \sin k x")

where ``P_m(x)`` is a polynomial in ``x`` of degree ``m``. 

"""

# ╔═╡ c0f7a0ed-0a89-4e5d-933d-f261c952a48c
cm"""
| ``\boldsymbol{f}(\boldsymbol{x})`` | ``\boldsymbol{y}_{\boldsymbol{p}}`` |
| :---: | :---: |
| ``P_m(x)=b_0+b_1 x+b_2 x^2+\cdots+b_m x^m`` | ``x^s\left(A_0+A_1 x+A_2 x^2+\cdots+A_m x^m\right)`` |
| ``a \cos k x+b \sin k x`` | ``x^s(A \cos k x+B \sin k x)`` |
| ``e^{r x}(a \cos k x+b \sin k x)`` | ``x^s e^{r x}(A \cos k x+B \sin k x)`` |
| ``P_m(x) e^{r x}`` | ``x^s\left(A_0+A_1 x+A_2 x^2+\cdots+A_m x^m\right) e^{r x}`` |
| ``P_m(x)(a \cos k x+b \sin k x)`` | ``x^s\left[\left(A_0+A_1 x+\cdots+A_m x^m\right) \cos k x\right.`` |
|  | ``\left.+\left(B_0+B_1 x+\cdots+B_m x^m\right) \sin k x\right]`` |
"""

# ╔═╡ b8db6d0c-7200-4c76-9fed-3a774ce8dae9
md"## Variation of Parameters"

# ╔═╡ 274e341c-66c4-434a-b002-079130d6ee81
md"# 6.1 Introduction to Eigenvalues"

# ╔═╡ 41799b6d-d000-40c9-8144-211645bf6a3e
md"##  The Characteristic Equation"

# ╔═╡ 50c06b0f-3a4e-4f91-ad85-5b814dffa8bb
md"## Eigenspaces"

# ╔═╡ 2e116819-87fb-4fef-892e-3fbd19f3862e
cm"""


Let ``\lambda`` be a fixed eigenvalue of the ``n \times n`` matrix ``\mathbf{A}``. Then the set of all eigenvectors associated with ``\mathbf{A}`` is the set of all nonzero solution vectors of the system
```math
(\mathbf{A}-\lambda \mathbf{I}) \mathbf{v}=\mathbf{0} .
```

The __solution space__ of this system is called the __eigenspace of ``\mathbf{A}`` associated with the eigenvalue ``\lambda``__. This subspace of ``\mathbf{R}^n`` consists of all eigenvectors associated with ``\lambda`` together with the zero vector. 
"""

# ╔═╡ 92297d90-f326-4f1e-80fc-5465621e8a6b
md"# 6.2 Diagonalization of Matrices"

# ╔═╡ 050d88a5-3b65-4940-b6b2-6750eaff6d66
let
	λ1 = -2
	v1 = [-1;1]
	λ2 = 3
	v2 = [-7;2]
	P = [v1 v2]
	Pinv= inv(P)
	D = Diagonal([-2;3])
	P*D*Pinv
end

# ╔═╡ b7e7238a-b2ec-4632-85b3-c069217b7ecc
md"## Similarity and Diagonalization"

# ╔═╡ 343d8a6b-c6bd-41ef-ab28-3872753864c6
let
	f(x)=x^3-7x^2+16x-12
	# map(x->(x,f(x)),factor(Vector,-12))
end

# ╔═╡ ef081dfa-b610-4c7a-a039-7258f4f6e80e
begin
    function add_space(n=1)
        repeat("&nbsp;", n)
    end
    function post_img(img::String, w=500)
        res = Resource(img, :width => w)
        cm"""
      <div class="img-container">

      $(res)

      </div>"""
    end
    function poolcode()
        cm"""
      <div class="img-container">

      $(Resource("https://www.dropbox.com/s/cat9ots4ausfzyc/qrcode_itempool.com_kfupm.png?raw=1",:width=>300))

      </div>"""
    end
    function define(t="")
        beginBlock("Definition", t)
    end
    function remark(t="")
        beginBlock("Remark", t)
    end
    function remarks(t="")
        beginBlock("Remarks", t)
    end
    function bbl(t)
        beginBlock(t, "")
    end
    function bbl(t, s)
        beginBlock(t, s)
    end
    ebl() = endBlock()
    function theorem(s)
        bth(s)
    end
    function bth(s)
        beginTheorem(s)
    end
    eth() = endTheorem()
    ex(n::Int; s::String="") = ex("Example $n", s)
    ex(t::Int, s::String) = example("Example $t", s)
    ex(t, s) = example(t, s)
    function beginBlock(title, subtitle)
        """<div style="box-sizing: border-box;">
       	<div style="display: flex;flex-direction: column;border: 6px solid rgba(200,200,200,0.5);box-sizing: border-box;">
       	<div style="display: flex;">
       	<div style="background-color: #FF9733;
       	    border-left: 10px solid #df7300;
       	    padding: 5px 10px;
       	    color: #fff!important;
       	    clear: left;
       	    margin-left: 0;font-size: 112%;
       	    line-height: 1.3;
       	    font-weight: 600;">$title</div>  <div style="olor: #000!important;
       	    margin: 0 0 20px 25px;
       	    float: none;
       	    clear: none;
       	    padding: 5px 0 0 0;
       	    margin: 0 0 0 20px;
       	    background-color: transparent;
       	    border: 0;
       	    overflow: hidden;
       	    min-width: 100px;font-weight: 600;
       	    line-height: 1.5;">$subtitle</div>
       	</div>
       	<p style="padding:5px;">
       """
    end
    function beginTheorem(subtitle)
        beginBlock("Theorem", subtitle)
    end
    function endBlock()
        """</p></div></div>"""
    end
    function endTheorem()
        endBlock()
    end
    ex() = example("Example", "")
    function example(lable, desc)
        """<div style="display:flex;">
       <div style="
       font-size: 112%;
           line-height: 1.3;
           font-weight: 600;
           color: #f9ce4e;
           float: left;
           background-color: #5c5c5c;
           border-left: 10px solid #474546;
           padding: 5px 10px;
           margin: 0 12px 20px 0;
           border-radius: 0;
       ">$lable:</div>
       <div style="flex-grow:3;
       line-height: 1.3;
           font-weight: 600;
           float: left;
           padding: 5px 10px;
           margin: 0 12px 20px 0;
           border-radius: 0;
       ">$desc</div>
       </div>"""
    end
    @htl("")
end

# ╔═╡ 8408e369-40eb-4f9b-a7d7-26cde3e34a74
begin
    text_book = post_img("https://www.pearson.com/store/medias/size-W650-desktop-A1030-00-23-46-A103000234623-A103000234623-Lrg.jpg?context=bWFzdGVyfGltYWdlc3w3NTU4NTd8aW1hZ2UvanBlZ3xzeXMtbWFzdGVyL2ltYWdlcy9oNjcvaDQyLzEzMDMzMzI2NTc1NjQ2L3NpemVfVzY1MF9kZXNrdG9wXy9BMTAzMC8wMC8yMy80Ni9BMTAzMDAwMjM0NjIzL0ExMDMwMDAyMzQ2MjNfTHJnLmpwZ3w2NWIzMWM3NzA5ZmMxOGUxYzA4MjcyMDQ3OTYzODZjNmRhYjY4NTEyMWE4ZmM2ODYyYTlkZTNkYTI5ZjgxY2M4", 200)
    md""" # Syllabus
    ## Syallbus
    See here [Term 241 - MATH208 - Syllabus](https://www.dropbox.com/scl/fi/5u3riz5gzbfuj9c5h8bn2/T241MATH208Syllabus.pdf?rlkey=pcjja6yvflfzazbxx1wkkyazo&raw=1)
    ## Textbook
    __Textbook: Edwards, C. H., Penney, D. E., and Calvis, D. T., Differential Equations and Linear Algebra, Fourth edition, Pearson, 2021__
    $text_book

    ## Office Hours
    I strongly encourage all students to make use of my office hours. These dedicated times are a valuable opportunity for you to ask questions, seek clarification on lecture material, discuss challenging problems, and get personalized feedback on your work. Engaging with me during office hours can greatly enhance your understanding of the course content and improve your performance. Whether you're struggling with a specific concept or simply want to delve deeper into the subject, I am here to support your learning journey. Don't hesitate to drop by; __your success is my priority__.

    | Day       | Time        |
    |-----------|-------------|
    | Sunday    | 11:00-11:50AM |
    | Monday | 12:00-12:50PM |
    Also you can ask for an online meeting through __TEAMS__.
    """
end

# ╔═╡ e6c9f44f-19dd-4922-bc05-8614f5441e80
cm"""
$(ex("Example 3","Rate of cooling"))
 Newton's law of cooling may be stated in this way: __The time rate of change (the rate of change with respect to time ``t`` ) of the temperature ``T(t)`` of a body is proportional to the difference between ``T`` and the temperature ``A`` of the surrounding medium (Fig. 1.1.1)__. That is,
```math
\frac{d T}{d t}=-k(T-A)
```
where 
- ``k`` is a positive constant. 
- Observe that if ``T > A``, then ``d T / d t < 0``, so the temperature is a decreasing function of ``t`` and the body is cooling. But if ``T < A``, then ``d T / d t >0``, so that ``T`` is increasing.

Thus the physical law is translated into a differential equation. If we are given the values of ``k`` and ``A``, we should be able to find an explicit formula for ``T(t)``, and then-with the aid of this formula-we can predict the future temperature of the body.
"""

# ╔═╡ a1cc7232-8dff-4fb4-908a-01c105e62797
cm"""
$(post_img("https://www.dropbox.com/scl/fi/l1ku5zzq9cartv6qck9e3/fig_1_1_1.png?rlkey=pleb3h5e9gyrjwq9zop53t3g5&dl=1",500))
"""

# ╔═╡ 0853d79c-5ae0-4690-9f09-9068505ba213
cm"""
$(ex("Example 5","Population growth")) 
__The time rate of change of a population ``P(t)`` with constant birth and death rates is, in many simple cases, proportional to the size of the population__. That is,
$(texeq"\frac{d P}{d t}=k P \label{pop_growth}")
where ``k`` is the constant of proportionality.

This equation has the solution

$(texeq"P(t)=Ce^{kt}\label{pop_growth_sol}")
"""

# ╔═╡ 2b5ca5c2-4232-48df-95ab-4a5765436995
cm"""
$(ex("Example 6","Population growth"))
Suppose that ``P(t)=C e^{k t}`` is the population of a colony of bacteria at time ``t``, that the population at time ``t=0`` (hours, h) was 1000 , and that the population doubled after 1 h . This additional information about ``P(t)`` yields the following equations:
```math
\begin{aligned}
& 1000=P(0)=C e^0=C \\
& 2000=P(1)=C e^k
\end{aligned}
```

It follows that ``C=1000`` and that ``e^k=2``, so ``k=\ln 2 \approx 0.693147``. With this value of ``k`` the differential equation in $(eqref("pop_growth")) is
```math
\frac{d P}{d t}=(\ln 2) P \approx(0.693147) P
```

Substitution of ``k=\ln 2`` and ``C=1000`` in Eq. $(eqref("pop_growth_sol")) yields the particular solution
```math
P(t)=1000 e^{(\ln 2) t}=1000\left(e^{\ln 2}\right)^t=1000 \cdot 2^t \quad\left(\text { because } e^{\ln 2}=2\right)
```
that satisfies the given conditions. We can use this particular solution to predict future populations of the bacteria colony. For instance, the predicted number of bacteria in the population after one and a half hours (when ``t=1.5`` ) is
```math
P(1.5)=1000 \cdot 2^{3 / 2} \approx 2828
```
"""

# ╔═╡ 3dab7c5b-363a-4b55-b8b0-342723d99ba5
cm"""
$(bbl("Remarks",""))

$(post_img("https://www.dropbox.com/scl/fi/9jkvevpmi4apwd5k73mgv/fig_1_1_4.png?rlkey=m2gu7plj5o7bjbudeqze2bzxq&dl=1",500))
$(ebl())
"""

# ╔═╡ a7abd027-69c5-43a9-a99e-164a9199a334
cm"""
$(ex(8))
Verify that the function ``y(x)=2 x^{1 / 2}-x^{1 / 2} \ln x`` satisfies the differential equation
```math
4 x^2 y^{\prime \prime}+y=0
```
for all ``x>0``.
"""

# ╔═╡ fbd68049-57ac-42fc-9c4c-5bd9f019dabb
cm"""
__The first-order equation__
\[
d y / d x=f(x, y)
\]
- __The general solution__
```math
y(x)=\int f(x) d x+C
```
- __A particular solution__ : a solution to the *initial-value* problem
```math
\frac{d y}{d x}=f(x), \quad y\left(x_0\right)=y_0
```

$(bbl("Remark",""))
> We will first find a general solution involving an arbitrary
constant ``C``. We can then attempt to obtain, by appropriate choice of ``C``, a particular
solution satisfying a given initial condition ``y(x_0) = y_0``.

$(ebl())
"""

# ╔═╡ 41344730-e8c5-4e2b-a98c-9995846c7244
cm"""
$(ex(1))
Solve the initial value problem 

```math
\frac{d y}{d x}=2 x+3, \quad y(1)=2
```
"""

# ╔═╡ ab0da87e-9b87-48eb-a235-44fdfd2f81f2
cm"""
$(post_img("https://www.dropbox.com/scl/fi/l6nrl55ye583z7vu2v55k/fig_1_2_5.png?rlkey=vpoky2h6hycj51kfemjjv0itg&dl=1",500))
"""

# ╔═╡ 6168170e-2fcb-4728-b8d8-ddc44992d3f9
cm"""
$(ex("Example 4","River crossing"))
Suppose that the river is ``1`` mile wide and that its midstream velocity is ``v_0=9 \mathrm{mi} / \mathrm{h}`` and the swimmer's velocity is ``v_S=3 \mathrm{mi} / \mathrm{h}``.
"""

# ╔═╡ 15531160-b5d7-4e55-848c-9b239d4f116c
cm"""
$(ex(1)) Solve the differential equation
```math
\frac{d y}{d x}=\frac{4-2 x}{3 y^2-5}
```
"""

# ╔═╡ 6d3c4d0e-a950-43bd-a2d7-46aec4417ab3
cm"""
$(ex(2)) Find all solutions of the differential equation
```math
\frac{d y}{d x}=6 x(y-1)^{2 / 3}
```

"""

# ╔═╡ a1af7fe8-cb18-4b36-a015-38df3d346b33
cm"""
$(ex("Example 3", "World population")) According to data listed at www. census.gov, the world's total population reached ``6`` billion persons in mid-``1999``, and was then increasing at the rate of about ``212`` thousand persons each day. Assuming that natural population growth at this rate continues, we want to answer these questions: 

<ol type="a">

<li> 

What is the annual growth rate ``k`` ? 
</li>

<li> What will be the world population at the middle of the 21 st century? </li>
<li> How long will it take the world population to increase tenfold-thereby reaching the 60 billion that some demographers believe to be the maximum for which the planet can provide adequate food supplies? </li>

</ol>
"""

# ╔═╡ ceb445e8-64ca-4d6f-86c2-98c96ad42fb2
cm"""
$(ex("Example 4","Radiometric dating"))
A specimen of charcoal found at Stonehenge turns out to contain ``63 \%`` as much ``{ }^{14} \mathrm{C}`` as a sample of present-day charcoal of equal mass. What is the age of the 
sample?
"""

# ╔═╡ 177d7ac6-2a24-4282-84a3-e42a1c03251f
cm"""
$(ex(5,s="Cooling"))
A ``4``-lb roast, initially at ``50^{\circ} \mathrm{F}``, is placed in a ``375^{\circ} \mathrm{F}`` oven at 5:00 P.M. After 75 minutes it is found that the temperature ``T(t)`` of the roast is ``125^{\circ} \mathrm{F}``. When will the roast be ``150^{\circ} \mathrm{F}`` (medium rare) ``?``
"""

# ╔═╡ cd7eb646-c23c-4ab8-9dc4-04c0db21d4bb
cm"""
$(ex(1)) 
Solve the initial value problem
```math
\frac{d y}{d x}-y=\frac{11}{8} e^{-x / 3}, \quad y(0)=-1
```
"""

# ╔═╡ 2e8d8477-a937-4cd4-b300-2cae6b9b2853
cm"""
$(ex(2)) Find a general solution of
```math
\left(x^2+1\right) \frac{d y}{d x}+3 x y=6 x
```
"""

# ╔═╡ e0c6b9b1-ab42-4e6f-a41e-440e68d68041
cm"""
$(bth("1 The Linear First-Order Equation"))
If the functions ``P(x)`` and ``Q(x)`` are continuous on the open interval ``I`` containing the point ``x_0``, then the initial value problem
```math
\frac{d y}{d x}+P(x) y=Q(x), \quad y\left(x_0\right)=y_0
```
has a unique solution ``y(x)`` on ``I``, given by the formula in 
```math
y(x)=e^{-\int P(x) d x}\left[\int\left(Q(x) e^{\int P(x) d x}\right) d x+C\right]
```
with an appropriate value of ``C``.
"""

# ╔═╡ 1d54e85b-f3e6-41f1-a818-ab483c918fee
cm"""
$(ex(3)) Solve the initial value problem
```math
x^2 \frac{d y}{d x}+x y=\sin x, \quad y(1)=y_0
```
"""

# ╔═╡ 4fbe0c02-a31f-4f37-b014-84e4bb914fff
cm"""
$(example("Example",""))
Solve
```math 
\left(x+y e^y\right) \frac{d y}{d x}=1
```
"""

# ╔═╡ c868de03-803a-49a9-b1f5-657f46ab8498
cm"""
$(ex(1)) Solve the differential equation
```math
\frac{d y}{d x}=(x+y+3)^2
```
"""

# ╔═╡ 90a2fcd1-1d38-4c12-8f7e-cb9f83f890f3
cm"""
$(define("Homogeneous DE")) A __homogeneous__ first-order differential equation is one that can be written in the form
```math
\frac{d y}{d x}=F\left(\frac{y}{x}\right)
```

If we make the substitutions
```math
v=\frac{y}{x}, \quad y=v x, \quad \frac{d y}{d x}=v+x \frac{d v}{d x}
```
"""

# ╔═╡ aca4f2f0-5c16-4e51-b83c-b76149f836a9
cm"""
$(ex(2))
Solve the differential equation
```math
2 x y \frac{d y}{d x}=4 x^2+3 y^2
```
"""

# ╔═╡ d743b750-2a29-4fb3-a729-07bec9637121
cm"""
$(ex(3)) Solve the initial value problem
```math
x \frac{d y}{d x}=y+\sqrt{x^2-y^2}, \quad y\left(x_0\right)=0
```
where ``x_0>0``.
"""

# ╔═╡ 7a405477-d020-4b31-b0be-036fd3f22322
cm"""
$(ex(5)) Solve the differential equation
```math
x \frac{d y}{d x}+6 y=3 x y^{4 / 3}
```
"""

# ╔═╡ 0e5e9fe4-0157-4497-b675-53071f4e6782
cm"""
$(ex(6))
Solve
```math
2 x e^{2 y} \frac{d y}{d x}=3 x^4+e^{2 y}
```

"""

# ╔═╡ a32e62d2-e65f-492a-a3f4-db8dec72db67
cm"""
- A general solution ``y(x)`` of a __first-order differential equation__ is often defined implicitly by an equation of the form
```math
F(x, y(x))=C\quad \text{where } C \text{ is a constant.}
```
- The original differentail equation is
```math
\frac{\partial F}{\partial x}+\frac{\partial F}{\partial y} \frac{d y}{d x}=0
```
$(add_space(10)) that is,
```math
M(x, y)+N(x, y) \frac{d y}{d x}=0
```
$(add_space(10))where ``M(x, y)=F_x(x, y)`` and ``N(x, y)=F_y(x, y)``
- We write 
```math
M(x, y) d x+N(x, y) d y=0
```
- This last form is called the __differential form__.
- If there exists a function F(x,y) such that 
```math
\frac{\partial F}{\partial x}=M \quad \text{and} \quad \frac{\partial F}{\partial y}=N
```
$(add_space(10))then the equation
```math
F(x, y)=C
```
implicitly defines a general solution. In this case, the equation
```math
M(x, y)+N(x, y) \frac{d y}{d x}=0
```
is called an __exact differential equation__ 
```math
(\text{the differential} dF=F_x d x+F_y d y \text{ of }  F(x, y) \text{ is exactly } M d x+N d y)
```

"""

# ╔═╡ 8792fe08-498e-4d18-aad2-de9bc3a7ade4
cm"""
$(bbl("Question 1"))
How can we determine whether the differential
 equation in 
```math
M(x, y) d x+N(x, y) d y=0
```
is exact?
$(ebl())

$(bbl("Question 2"))
 If it is exact, how can we find the function ``F`` such
 that 
```math
 F_x = M, \quad F_y = N?
```
$(ebl())
"""

# ╔═╡ c006d3aa-8cf1-4382-806a-4e9714f934f6
cm"""
$(bth("1 Criterion for Exactness"))
Suppose that the functions ``M(x, y)`` and ``N(x, y)`` are continuous and have continuous first-order partial derivatives in the open rectangle ``R`` : ``a < x < b, c < y < d``. Then the differential equation
```math
M(x, y) d x+N(x, y) d y=0
```
is exact in ``R`` if and only if
```math
\frac{\partial M}{\partial y}=\frac{\partial N}{\partial x}
```
at each point of ``R``. That is, there exists a function ``F(x, y)`` defined on ``R`` with ``\partial F / \partial x=M`` and ``\partial F / \partial y=N`` if and only if last equation holds on ``R``.
"""

# ╔═╡ 626fe783-d05a-49ca-8e3c-2ac76be27e34
cm"""
$(ex(8)) 
Show that the differential equation 
```math
\quad y^3 d x+3 x y^2 d y=0
```
is exact.
"""

# ╔═╡ ad29db84-3c84-4a6e-aeb0-5a51210ef05d
cm"""
$(bbl("Remarks",""))
- What happens if we divide by ``y^2`` both sides?
"""

# ╔═╡ 933c6345-fffc-4159-8b67-e1443b988f9f
cm"""
$(ex(9)) Solve the differential equation
```math
\left(6 x y-y^3\right) d x+\left(4 y+3 x^2-3 x y^2\right) d y=0
```
"""

# ╔═╡ a4256712-b454-4733-bcf5-c5764705b028
cm"""
$(ex(10))
Solve the equation ``x y^{\prime \prime}+2 y^{\prime}=6 x``.
"""

# ╔═╡ 87802187-2aeb-4e67-a2d7-1975cf588a17
cm"""
$(ex(11)) 
Solve the equation ``y y^{\prime \prime}=\left(y^{\prime}\right)^2``.
"""

# ╔═╡ 608e2182-ce7c-4626-bcda-bf4bd6f2f1c6
cm"""
## Using Julia
$(example("EXample",""))
Solve the differential equation
```math
2 x e^{2 y} \frac{d y}{d x}=3 x^4+e^{2 y}
```
using Julia
"""

# ╔═╡ 643d7441-1dbf-4b98-b76f-617babece0c0
cm"""
$(example("Example","")) Solve the linear system
```math
\begin{aligned}
x+2 y+z & =4 \\
3 x+8 y+7 z & =20 \\
2 x+7 y+9 z & =23
\end{aligned}
```
"""

# ╔═╡ 2aa7e38b-c48d-470c-943e-a6a73d053a70
cm"""
$(define("Echelon Matrix"))
The matrix ``\mathbf{E}`` is called an echelon matrix provided it has the following two properties:
1. Every row of ``\mathbf{E}`` that consists entirely of zeros (if any) lies beneath every row that contains a nonzero element.
2. In each row of ``\mathbf{E}`` that contains a nonzero element, the first nonzero element lies strictly to the right of the first (from left) nonzero element (called __leading entry__) in the preceding row (if there is a preceding row).
$(ebl())

__FOR EXAMPLE__
```math
\mathbf{E}=\left[\begin{array}{rrrrr}2 & -1 & 0 & 4 & 7 \\ 0 & 1 & 2 & 0 & -5 \\ 0 & 0 & 0 & 3 & 0 \\ 0 & 0 & 0 & 0 & 0\end{array}\right]
```
"""

# ╔═╡ 13748f1f-4092-4592-9a4f-8f17c38a875a
cm"""
$(define("Reduced Echelon Matrix"))
A reduced echelon matrix ``\mathbf{E}`` is an echelon matrix that has-in addition to Properties 1 and 2 -the following properties:

3. Each leading entry of ``\mathbf{E}`` is 1 .
4. Each leading entry of ``\mathbf{E}`` is the only nonzero element in its column.
"""

# ╔═╡ 2d3b73f6-8429-4551-96ac-b2db865189ed
cm"""
$(bbl("ALGORITHM Gauss-Jordan Elimination",""))
1. First transform ``\mathbf{A}`` into echelon form by Gaussian elimination.
2. Then divide each element of each nonzero row by its leading entry (to satisfy Property 3).
3. Finally, use each leading 1 to "clear out" any remaining nonzero elements in its column (to satisfy Property 4).
$(ebl())
> Every matrix is row equivalent to one and only one reduced echelon matrix.
"""

# ╔═╡ ab51c417-0482-4c68-ba07-38fd46e535cd
cm"""
$(ex())Find the reduced echelon form of the matrix
```math
\mathbf{A}=\left[\begin{array}{rrrr}
1 & 2 & 1 & 4 \\
3 & 8 & 7 & 20 \\
2 & 7 & 9 & 23
\end{array}\right]
```
"""

# ╔═╡ 3e062725-d0ac-4d5e-9baa-49b1b10de464
cm"""
$(ex())
Use Gauss-Jordan elimination to solve the linear system
```math
\begin{aligned}
x_1+x_2+x_3+x_4 & =12 \\
x_1+2 x_2+5 x_4 & =17 \\
3 x_1+2 x_2+4 x_3-x_4 & =31
\end{aligned}
```
"""

# ╔═╡ 34820013-08f9-4077-9403-206671bea915
cm"""
$(bth("The Three Possibilities")) 
A linear system of equations has either
- a unique solution, or
- no solution, or
- infinitely many solutions.
"""

# ╔═╡ 6f434bf9-93f6-4052-aed0-7ebe0c06b736
cm"""
$(bth(""))
Every homogeneous linear system with more variables than equations has infinitely many solutions.
"""

# ╔═╡ 22bb573e-8fa7-4f6c-b5f6-b96dee4c8507
cm"""
$(ex()) Determine the constants ``A`` and ``B`` so as to find a solution of the differential equation that satisfies the given initial conditions involving ``y(0)`` and ``y^{\prime}(0)``.
```math
\begin{aligned} & y^{\prime \prime}-25 y=0, y(x)=A e^{5 x}+B e^{-5 x} \\ & y(0)=10, y^{\prime}(0)=20\end{aligned}
```
"""

# ╔═╡ 6e605916-523b-466b-ad96-c206af3bf27c
cm"""
$(ex())
Under what condition on the constants ``a, b``, and ``c`` does the system
```math
\begin{array}{r}
2 x-y+3 z=a \\
x+2 y+z=b \\
7 x+4 y+9 z=c
\end{array}
```
have a unique solution? No solution? Infinitely many solutions?
"""

# ╔═╡ bf7f3424-92fd-4ece-8aeb-db8da6762056
cm"""
$(define("Invertible Matrix"))
The square matrix ``\mathbf{A}`` is called invertible if there exists a matrix ``\mathbf{B}`` such that
```math
\mathbf{A B}=\mathbf{B A}=\mathbf{I} .
```
"""

# ╔═╡ c35aca68-0b50-430f-b1e1-bc2ab187fa40
cm"""
$(ex())
Find ``\mathbf{A}^{-1}``. Then use ``\mathbf{A}^{-1}`` (as in Example 5) to solve the system ``\mathbf{A x}=\mathbf{b}``
```math
\mathbf{A}=\left[\begin{array}{ll}3 & 2 \\ 5 & 4\end{array}\right], \mathbf{b}=\left[\begin{array}{l}5 \\ 6\end{array}\right]
```
"""

# ╔═╡ eb2a3144-7dca-42d0-9fce-01c78728eac7
# ```math
# \mathbf{A}=\left[\begin{array}{lll}
# 4 & 3 & 2 \\
# 5 & 6 & 3 \\
# 3 & 5 & 2
# \end{array}\right]
# ```
cm"""
$(ex())Find the inverse of the ``3 \times 3`` matrix
```math
\mathbf{A}=\left[\begin{array}{rrr}1 & 1 & 5 \\ 1 & 4 & 13 \\ 3 & 2 & 12\end{array}\right]
```
"""

# ╔═╡ 6189764d-9da8-460b-b873-a040c413e721
cm"""
$(bth("Properties of Nonsingular Matrices"))
The following properties of an ``n \times n`` matrix ``\mathbf{A}`` are equivalent.
1. ``\mathbf{A}`` is invertible.
2. ``\mathbf{A}`` is row equivalent to the ``n \times n`` identity matrix ``\mathbf{I}``.
3. ``\mathbf{A x}=\mathbf{0}`` has only the trivial solution.
4. For every ``n``-vector ``\mathbf{b}``, the system ``\mathbf{A x}=\mathbf{b}`` has a unique solution.
5. For every ``n``-vector ``\mathbf{b}``, the system ``\mathbf{A x}=\mathbf{b}`` is consistent.
"""

# ╔═╡ 9d011957-4432-484f-99dc-da94395aeef1
cm"""
$(ex())
Find a matrix ``\mathbf{X}`` such that ``\mathbf{A X}=\mathbf{B}``.
```math
\mathbf{A}=\left[\begin{array}{ll}7 & 6 \\ 8 & 7\end{array}\right], \mathbf{B}=\left[\begin{array}{rrr}2 & 0 & 4 \\ 0 & 5 & -3\end{array}\right]
```
"""

# ╔═╡ fb9cd8c0-ace5-4123-a8a0-58ab32f849f1
cm"""
$(ex())
Apply Cramer's rule to solve the system
```math
\begin{aligned}
& 7 x+8 y=5 \\
& 6 x+9 y=4
\end{aligned}
```
"""

# ╔═╡ c63574f6-2f03-4f3d-b9eb-08ca321a0aad
cm"""
$(define("Minors and Cofactors"))
Let ``\mathbf{A}=\left[a_{i j}\right]`` be an ``n \times n`` matrix. The ``i j`` th minor of ``\mathbf{A}`` (also called the minor of ``a_{i j}`` ) is the determinant ``M_{i j}`` of the ``(n-1) \times(n-1)`` submatrix that remains after deleting the ``i`` th row and the ``j`` th column of ``\mathbf{A}``. The ``i j`` th cofactor ``A_{i j}`` of ``\mathbf{A}`` (or the cofactor of ``a_{i j}`` ) is defined to be
```math
A_{i j}=(-1)^{i+j} M_{i j}
```
$(ebl())
For the sign of the cofactors
```math
\left[\begin{array}{lll}+ & - & + \\ - & + & - \\ + & - & +\end{array}\right] \quad\text{and} \quad\left[\begin{array}{cccc}+ & - & + & - \\ - & + & - & + \\ + & - & + & - \\ - & + & - & +\end{array}\right]
```
So
```math
\begin{array}{llll}A_{11}=+M_{11}, & A_{12}=-M_{12}, & A_{13}=+M_{13}, & A_{14}=-M_{14} \\ A_{21}=-M_{21}, & A_{22}=+M_{22}, & A_{23}=-M_{23}, & A_{24}=+M_{24}\end{array}
```
"""

# ╔═╡ 737365c0-472d-4f0d-87c5-5d70a2fe027f
cm"""
$(define("Determinants"))
The determinant ``\operatorname{det} \mathbf{A}=\left|a_{i j}\right|`` of an ``n \times n`` matrix ``\mathbf{A}=\left[a_{i j}\right]`` is defined as
```math
\operatorname{det} \mathbf{A}=a_{11} A_{11}+a_{12} A_{12}+\cdots+a_{1 n} A_{1 n}
```

Thus we multiply each element of the first row of ``\mathbf{A}`` by its cofactor and then add these ``n`` products to get ``\operatorname{det} \mathbf{A}``.
"""

# ╔═╡ b583732f-ed42-4a4f-9dec-6369b5e85c41
cm"""
$(ex()) Evaluate the determinant of
```math
\mathbf{A}=\left[\begin{array}{rrrr}
2 & 0 & 0 & -3 \\
0 & -1 & 0 & 0 \\
7 & 4 & 3 & 5 \\
-6 & 2 & 2 & 4
\end{array}\right]
```
"""

# ╔═╡ 8dae3e8d-96e3-4fe0-8b9f-5074dff0f9e4
cm"""
$(bth("Cofactor Expansions of Determinants"))
The determinant of an ``n \times n`` matrix ``\mathbf{A}=\left[a_{i j}\right]`` can be obtained by expansion along any row or column. The cofactor expansion along the ``i`` th row is
```math
\operatorname{det} \mathbf{A}=a_{i 1} A_{i 1}+a_{i 2} A_{i 2}+\cdots+a_{i n} A_{i n} .
```

The cofactor expansion along the ``j`` th column is
```math
\operatorname{det} \mathbf{A}=a_{1 j} A_{1 j}+a_{2 j} A_{2 j}+\cdots+a_{n j} A_{n j}
```
"""

# ╔═╡ 1ba6a9f9-f0b3-4033-a20a-117b19d7c74b
cm"""
1. __Property 1__: If the ``n \times n`` matrix ``\mathbf{B}`` is obtained from ``\mathbf{A}`` by multiplying a single row (or a column) of ``\mathbf{A}`` by the constant ``k``, then ``\operatorname{det} \mathbf{B}=k \operatorname{det} \mathbf{A}``.
2. __Property 2__: If the ``n \times n`` matrix ``\mathbf{B}`` is obtained from ``\mathbf{A}`` by interchanging two rows (or two columns), then ``\operatorname{det} \mathbf{B}=-\operatorname{det} \mathbf{A}``.
3. __Property 3__: If two rows (or two columns) of the ``n \times n`` matrix ``\mathbf{A}`` are identical, then ``\operatorname{det} \mathbf{A}=0``.
4. __Property 4__: ``\quad`` Suppose that the ``n \times n`` matrices ``\mathbf{A}_1, \mathbf{A}_2``, and ``\mathbf{B}`` are identical except for their ``i`` th rows-that is, the other ``n-1`` rows of the three matrices are identicaland that the ``i`` th row of ``\mathbf{B}`` is the sum of the ``i`` th rows of ``\mathbf{A}_1`` and ``\mathbf{A}_2``. Then
```math
\operatorname{det} \mathbf{B}=\operatorname{det} \mathbf{A}_1+\operatorname{det} \mathbf{A}_2
```
$(add_space(10))This result also holds if columns are involved instead of rows.

5. __Property 5__: If the ``n \times n`` matrix ``\mathbf{B}`` is obtained by adding a constant multiple of one row (or column) of ``\mathbf{A}`` to another row (or column) of ``\mathbf{A}``, then ``\operatorname{det} \mathbf{B}=\operatorname{det} \mathbf{A}``.

6. __Property 6__: The determinant of a triangular matrix is equal to the product of its diagonal elements.
7. __Property 7__: If ``\mathbf{A}`` is a square matrix, then ``\operatorname{det}\left(\mathbf{A}^T\right)=\operatorname{det} \mathbf{A}``.

$(add_space(10))Note that:

$(add_space(20))``\text{(i)} \quad\left(\mathbf{A}^T\right)^T=\mathbf{A}``;

$(add_space(20))``\text{(ii)} \quad(\mathbf{A}+\mathbf{B})^T=\mathbf{A}^T+\mathbf{B}^T``;

$(add_space(20))``\text{(iii)} \quad(c \mathbf{A})^T=c \mathbf{A}^T``;

$(add_space(20))``\text{(iv)} \quad(\mathbf{A B})^T=\mathbf{B}^T \mathbf{A}^T``.

"""

# ╔═╡ ea830a4d-f11a-474e-9b42-dde6295e8d24
cm"""
$(ex()) Evaluate the determinant of
```math
\mathbf{A}=\left|\begin{array}{llll}1 & 2 & 3 & 4 \\ 0 & 5 & 6 & 7 \\ 0 & 0 & 8 & 9 \\ 2 & 4 & 6 & 9\end{array}\right|
```
"""

# ╔═╡ b0efeb12-e0fd-4a1e-80b9-24c230368c53
cm"""
$(bth("Determinants and Invertibility"))
The ``n \times n`` matrix ``\mathbf{A}`` is invertible if and only if ``\operatorname{det} \mathbf{A} \neq 0``.
"""

# ╔═╡ 42e01a69-848d-4d34-a5de-8c54ca3f3b83
cm"""
$(bth("The Inverse Matrix"))
The inverse of the invertible matrix ``\mathbf{A}`` is given by the formula
```math
\mathbf{A}^{-1}=\frac{\left[A_{i j}\right]^T}{|\mathbf{A}|}
```
where, as usual, ``A_{i j}`` denotes the ``i j`` th cofactor of ``\mathbf{A}``; that is, ``A_{i j}`` is the product of ``(-1)^{i+j}`` and the ``i j`` th minor determinant of ``\mathbf{A}``.
"""

# ╔═╡ e948494b-8a14-47b8-bfec-4390376e1bcd
cm"""
$(ex()) Apply the formula above to find the inverse of the matrix
```math
\mathbf{A}=\left[\begin{array}{rrr}
1 & 4 & 5 \\
4 & 2 & 5 \\
-3 & 3 & -1
\end{array}\right]
```
"""

# ╔═╡ 52d0d0d2-3079-4285-92de-c59961848921
cm"""
$(define("")) 
__Vector__ : A vector ``\mathbf{v}`` in 3 -space ``\mathbb{R}^3`` is simply an ordered triple ``(a, b, c)`` of real numbers. We write ``\mathbf{v}=(a, b, c)`` and call the numbers ``a, b``, and ``c`` the components (or coordinates) of the vector ``v``. We may also write 
```math
\mathbf{v} =\begin{bmatrix} a\\b\\c \end{bmatrix}.
```

__Addition of Vectors__: The sum ``\mathbf{u}+\mathbf{v}`` of the two vectors ``\mathbf{u}=\left(u_1, u_2, u_3\right)`` and ``\mathbf{v}=\left(v_1, v_2, v_3\right)`` is the vector
```math
\mathbf{u}+\mathbf{v}=\left(u_1+v_1, u_2+v_2, u_3+v_3\right)
```
that is obtained upon addition of respective components of ``\mathbf{u}`` and ``\mathbf{v}``.

__The geometric interpretation of vector addition__

$(post_img("https://www.dropbox.com/scl/fi/kvrnymjm7panzimoqe4ov/fig_4_1_2_and_3.png?rlkey=s6lpzojg7la8rybddv0t04gkf&raw=1",800))

__Multiplication of a Vector by a Scalar__: If ``\mathbf{v}=\left(v_1, v_2, v_3\right)`` is a vector and ``c`` is a real number, then the scalar multiple ``c \mathbf{v}`` is the vector
```math
c \mathbf{v}=\left(c v_1, c v_2, c v_3\right)
```
that is obtained upon multiplying each component of ``\mathbf{v}`` by ``c``.

__The length of a vector__: The length ``|\mathbf{v}|`` of the vector ``\mathbf{v}=(a, b, c)`` is defined to be the distance of the point ``P(a, b, c)`` from the origin,
```math
|\mathbf{v}|=\sqrt{a^2+b^2+c^2}
```
"""

# ╔═╡ aaae2aff-fc65-4b7f-9e28-b7577c5fdd51
cm"""
$(bth("3")) __``\mathbb{R}^3`` as a Vector Space__

If ``\mathbf{u}, \mathbf{v}``, and ``\mathbf{w}`` are vectors in ``\mathbf{R}^3``, and ``r`` and ``s`` are real numbers, then
1. ``\mathbf{u}+\mathbf{v}=\mathbf{v}+\mathbf{u}`` (commutativity)
2. ``\mathbf{u}+(\mathbf{v}+\mathbf{w})=(\mathbf{u}+\mathbf{v})+\mathbf{w}`` (associativity)
3. ``\mathbf{u}+\mathbf{0}=\mathbf{0}+\mathbf{u}=\mathbf{u}`` (zero element)
4. ``\mathbf{u}+(-\mathbf{u})=(-\mathbf{u})+\mathbf{u}=\mathbf{0}`` (additive inverse)
5. ``r(\mathbf{u}+\mathbf{v})=r \mathbf{u}+r \mathbf{v}`` (distributivity)
6. ``(r+s) \mathbf{u}=r \mathbf{u}+s \mathbf{u}``
7. ``r(s \mathbf{u})=(r s) \mathbf{u}``
8. ``1(\mathbf{u})=\mathbf{u}`` (multiplicative identity).
"""

# ╔═╡ ddda1623-3b11-46f7-9fda-fd6686dd4343
cm"""
$(bth("2 Two Linearly Dependent Vectors"))
The two vectors ``\mathbf{u}`` and ``\mathbf{v}`` are __linearly dependent__ if and only if there exist scalars ``a`` and ``b`` not both zero such that
```math
a \mathbf{u}+b \mathbf{v}=\mathbf{0} .
```
"""

# ╔═╡ c31e6c51-9f45-4581-a7b6-54bdf2e04ada
cm"""
$(bbl("Remark","Two Linearly Independent Vectors "))
The two vectors ``\mathbf{u}`` and ``\mathbf{v}`` are __linearly independent__ if and only if the relation
```math
a \mathbf{u}+b \mathbf{v}=\mathbf{0}
```
implies that ``a=b=0``.
"""

# ╔═╡ 856ff140-9c2b-4953-a033-849f1af34286
cm"""
$(ex(1))
Consider the vector 
```math
\mathbf{u}=(3,-2), \quad \mathbf{v}=(-6,4), \quad\text{and}\quad \mathbf{w}=(5,-7)
```
"""

# ╔═╡ 96a09993-9f1a-4908-9334-494b4e198455
cm"""
$(ex(2))
Express the vector ``\mathbf{w}=(11,4)`` as a linear combination of the vectors  ``\mathbf{u}=(3,-2)`` and ``\mathbf{v}= (-2,7)``

"""

# ╔═╡ dff6e998-431d-4493-bc34-6fa9b5de6f93
cm"""
$(define(""))
__Linearly Dependent Vectors in ``R^3``__

The three vectors ``\mathbf{u}, \mathbf{v}``, and ``\mathbf{w}`` in ``\mathbf{R}^3`` are said to be __linearly dependent__ provided that one of them is a linear combination of the other two-that is, either
```math
\begin{array}{ll}
\mathbf{w}=r \mathbf{u}+s \mathbf{v} & \text { or } \\
\mathbf{u}=r \mathbf{v}+s \mathbf{w} & \text { or } \\
\mathbf{v}=r \mathbf{u}+s \mathbf{w} &
\end{array}
```
for appropriate scalars ``r`` and ``s``.
"""

# ╔═╡ f61386e6-e32c-45f6-b1a9-dec901195e96
cm"""
$(bth("3 Three Linearly Dependent Vectors"))
The three vectors ``\mathbf{u}, \mathbf{v}``, and ``\mathbf{w}`` in ``\mathbf{R}^3`` are __linearly dependent__ if and only if there exist scalars ``a, b``, and ``c`` not all zero such that
```math
a \mathbf{u}+b \mathbf{v}+c \mathbf{w}=\mathbf{0}
```
"""

# ╔═╡ 20916835-0034-40b5-ab98-1f24fb78247c
cm"""
$(bbl("Remark","Linear Independence"))
The vectors ``\mathbf{u}, \mathbf{v}``, and ``\mathbf{w}`` are __linearly independent__ if and only if the relation
```math
a \mathbf{u}+b \mathbf{v}+c \mathbf{w}=\mathbf{0}
```
implies that ``a=b=c=0``.
"""

# ╔═╡ f868d290-9a4e-4a6f-933e-ee642deabeac
cm"""
$(bth("4 Three Linearly Independent Vectors"))
The vectors ``\mathbf{u}=\left(u_1, u_2, u_3\right), \mathbf{v}=\left(v_1, v_2, v_3\right)``, and ``\mathbf{w}=\left(w_1, w_2, w_3\right)`` are linearly independent if and only if
```math
\left|\begin{array}{lll}
u_1 & v_1 & w_1 \\
u_2 & v_2 & w_2 \\
u_3 & v_3 & w_3
\end{array}\right| \neq 0
```
"""

# ╔═╡ 20e79ac6-f0f4-40e6-89e0-2eaa44e0ff2b
cm"""
$(ex(3))
Check linear independence for

```math
\mathbf{u}=(1,2,-3), \quad \mathbf{v}=(3,1,-2), \quad \text{and}\quad  \mathbf{w}=(5,-5,6)
```
"""

# ╔═╡ f6653849-9370-4f84-a07a-ad924983fda6
cm"""
__The basic unit vectors__
```math
\mathbf{i}=(1,0,0), \quad \mathbf{j}=(0,1,0), \quad \text { and } \quad \mathbf{k}=(0,0,1)
```

The expression
```math
\mathbf{v}=a \mathbf{i}+b \mathbf{j}+c \mathbf{k}=(a, b, c)
```
shows both that
- the three vectors ``\mathbf{i}, \mathbf{j}``, and ``\mathbf{k}`` are linearly independent (because ``\mathbf{v}=\mathbf{0}`` immediately implies ``a=b=c=0`` ), and that
- any vector in ``\mathbf{R}^3`` can be expressed as a linear combination of ``\mathbf{i}, \mathbf{j}``, and ``\mathbf{k}``.

- A __basis__ for ``\mathbf{R}^3`` is a triple ``\mathbf{u}, \mathbf{v}, \mathbf{w}`` of vectors such that every vector ``\mathbf{t}`` in ``\mathbf{R}^3`` can be expressed as a linear combination
```math
\mathbf{t}=a \mathbf{u}+b \mathbf{v}+c \mathbf{w}
```
$(add_space(10))of them. 

That is, given any vector ``\mathbf{t}`` in ``\mathbf{R}^3``, there exist scalars ``a, b, c`` such that Eq. 
```math
\mathbf{t}=a \mathbf{u}+b \mathbf{v}+c \mathbf{w}
```
holds
"""

# ╔═╡ 4faab984-e4a5-4620-9c7e-61abd525a1fa
cm"""
$(bth("5"))
__Basis for ``\mathrm{R}^3``__

If the vectors ``\mathbf{u}, \mathbf{v}``, and ``\mathbf{w}`` in ``\mathbf{R}^3`` are linearly independent, then they constitute a basis for ``\mathbf{R}^3``.
"""

# ╔═╡ cf4fc0c1-c9d8-4597-b1d7-a5d65b59421d
cm"""
$(ex(4))
Express the vector ``\mathbf{t}=(4,20,23)`` as a combination of the linearly independent vectors ``\mathbf{u}=(1,3,2), \mathbf{v}=(2,8,7)``, and ``\mathbf{w}=(1,7,9)``
"""

# ╔═╡ d76bbff8-1730-41cb-93d4-ddb35f66b19b
cm"""
$(bbl("Subspaces",""))
A nonempty subset ``V`` of ``\mathbf{R}^3`` is a __subspace__ of ``\mathbf{R}^3`` if and only if it satisfies the following two conditions:

1. If ``\mathbf{u}`` and ``\mathbf{v}`` are vectors in ``V``, then ``\mathbf{u}+\mathbf{v}`` is also in ``V``  __(closure under addition)__.
2. If ``\mathbf{u}`` is a vector in ``V`` and ``c`` is a scalar, then ``c \mathbf{u}`` is in ``V`` __(closure under multiplication by scalars)__.
"""

# ╔═╡ d711bb6f-f627-4568-97f5-5d7690088663
cm"""
$(bbl("Remarks",""))

- The subspaces ``\{\boldsymbol{0}\}`` and ``\mathbf{R}^3`` are sometimes called the __trivial subspaces__ of ``\mathbf{R}^3`` (because the verification that they are subspaces is quite trivial). 
- All subspaces other than ``\{0\}`` and ``\mathbf{R}^3`` itself are called __proper subspaces__ of ``\mathbf{R}^3``.
- The proper subspaces of ``\mathbf{R}^3`` are what we customarily call __lines__ and __planes__ through the origin. 
"""

# ╔═╡ 263e7015-3130-4597-b440-a996e52967b9
cm"""
$(ex(5))
Let ``V`` be the set of all vectors ``(x, y)`` in ``\mathbf{R}^2`` such that ``y=x``. That is
```math
V = \left\{(x,y)\in \mathbb{R}^2 \;|\; y=x\right\}
````

Show that the set ``V`` is a subspace of ``\mathbf{R}^2``.
"""

# ╔═╡ 64e0c10c-a4e8-4435-86fc-60bff5e476d4
cm"""
$(ex(6))
Let ``V`` be the set of all vectors ``(x, y)`` in ``\mathbf{R}^2`` such that ``x+y=1``. That is
```math
V = \left\{(x,y)\in \mathbb{R}^2 \;|\; x+y=1\right\}
````

Show that the set ``V`` is a NOT subspace of ``\mathbf{R}^2``.
"""

# ╔═╡ d9a91130-0cf8-414e-b8a8-dd352090f7a6
cm"""
$(define(""))
__``n``-Space ``R^n``__

The ``\boldsymbol{n}``-dimensional space ``\mathbf{R}^{\boldsymbol{n}}`` is the set of all ``n``-tuples ``\left(x_1, x_2, x_3, \ldots, x_n\right)`` of real numbers.
"""

# ╔═╡ f1464582-faf2-4902-96a6-bb1cf17c1102
cm"""
$(define("Vector Space"))
Let ``V`` be a set of elements called vectors, in which the operations of addition of vectors and multiplication of vectors by scalars are defined. That is, given vectors ``\mathbf{u}`` and ``\mathbf{v}`` in ``V`` and a scalar ``c``, the vectors ``\mathbf{u}+\mathbf{v}`` and ``c \mathbf{u}`` are also in ``V`` (so that ``V`` is closed under vector addition and multiplication by scalars). Then, with these operations, ``V`` is called a vector space provided that—given any vectors ``\mathbf{u}, \mathbf{v}``, and ``\mathbf{w}`` in ``V`` and any scalars ``a`` and ``b``-the following properties hold true:
1. ``\mathbf{u}+\mathbf{v}=\mathbf{v}+\mathbf{u}``
(commutativity)
2. ``\mathbf{u}+(\mathbf{v}+\mathbf{w})=(\mathbf{u}+\mathbf{v})+\mathbf{w}``
(associativity)
3. ``\mathbf{u}+\mathbf{0}=\mathbf{0}+\mathbf{u}=\mathbf{u}``
(zero element)
4. ``\mathbf{u}+(-\mathbf{u})=(-\mathbf{u})+\mathbf{u}=\mathbf{0}``
(additive inverse)
5. ``a(\mathbf{u}+\mathbf{v})=a \mathbf{u}+a \mathbf{v}``
(distributivity)
6. ``(a+b) \mathbf{u}=a \mathbf{u}+b \mathbf{u}``
7. ``a(b \mathbf{u})=(a b) \mathbf{u}``
8. ``(1) \mathbf{u}=\mathbf{u}``
"""

# ╔═╡ f3d58858-f9ab-42bb-a822-14d03916984e
cm"""
$(ex(1)) 
Let ``\mathscr{F}`` be the set of all real-valued functions defined on the real number line ``\mathbb{R}`` is a vector space with function addition and scalar multiplication.
"""

# ╔═╡ 4b2763bd-cf78-4913-8d5d-ac10d08b1031
cm"""
$(define("Subspace"))
Let ``W`` be a nonempty subset of the vector space ``V``. Then ``W`` is a subspace of ``V`` provided that ``W`` itself is a vector space with the operations of addition and multiplication by scalars as defined in ``V``.
"""

# ╔═╡ abefd094-b57a-46c7-a314-eb836e58c6a0
cm"""
$(bth("1 Conditions for a Subspace"))
The nonempty subset ``W`` of the vector space ``V`` is a subspace of ``V`` if and only if it satisfies the following two conditions:
- (i) If ``\mathbf{u}`` and ``\mathbf{v}`` are vectors in ``W``, then ``\mathbf{u}+\mathbf{v}`` is also in ``W``.
- (ii) If ``\mathbf{u}`` is in ``W`` and ``c`` is a scalar, then the vector ``c \mathbf{u}`` is also in ``W``.
"""

# ╔═╡ f8cd8873-f4fb-4300-8aa2-4aee9990cf5e
cm"""
$(ex(2))
Let ``W`` be the subset of ``\mathbf{R}^n`` consisting of all those vectors ``\left(x_1, x_2, \ldots, x_n\right)`` whose coordinates satisfy the single homogeneous linear equation
```math
a_1 x_1+a_2 x_2+\cdots+a_n x_n=0
```
where the given coefficients ``a_1, a_2, \ldots, a_n`` are not all zero. Show that ``W`` is a subspace of ``\mathbb{R}^n``.
"""

# ╔═╡ e17f0845-12a9-4da8-b81d-bf3e61d95260
cm"""
$(ex(4))
Let ``W`` be the set of all those vectors ``\left(x_1, x_2, x_3, x_4\right)`` in ``\mathbb{R}^4`` such that ``x_1 x_4=0``. 

Is ``W`` a subspace of ``\mathbb{R}^4``?
"""

# ╔═╡ b13e61f1-3d38-4f87-838f-5c448ce1193a
cm"""
$(bth("2 Solution Subspaces"))
If ``\mathbf{A}`` is a (constant) ``m \times n`` matrix, then the solution set of the homogeneous linear system
```math
\mathbf{A x}=\mathbf{0}
```
is a subspace of ``\mathbf{R}^n``. This subspace is called __solution space__ of the system.
"""

# ╔═╡ add8e82e-02cb-4262-b8e4-ac7a8c317257
cm"""
$(ex(5))
In Example 4 of Section 3.4 we considered the homogeneous system
```math
\begin{aligned}
x_1+3 x_2-15 x_3+7 x_4 & =0 \\
x_1+4 x_2-19 x_3+10 x_4 & =0 \\
2 x_1+5 x_2-26 x_3+11 x_4 & =0
\end{aligned}
```
The reduced echelon form of the coefficient matrix of this system is
```math
\left[\begin{array}{rrrr}
1 & 0 & -3 & -2 \\
0 & 1 & -4 & 3 \\
0 & 0 & 0 & 0
\end{array}\right]
```

What is the solution space of the system?
"""

# ╔═╡ 56aee5a8-caf4-4cba-9e1e-7e319c21047b
cm"""
$(define("Linear Combination"))
The vector ``\mathbf{w}`` is called a __linear combination__ of the vectors ``\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k`` provided that there exist scalars ``c_1, c_2, \ldots, c_k`` such that
```math
\mathbf{w}=c_1 \mathbf{v}_1+c_2 \mathbf{v}_2+\cdots+c_k \mathbf{v}_k
```
"""

# ╔═╡ f110354f-ab42-481d-9057-fd467fa9a664
cm"""
$(define("Spanning Set"))
Suppose that ``\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k`` are vectors in a vector space ``V``. Then we say that the vectors ``\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k`` __span__ the vector space ``V`` provided that every vector in ``V`` is a linear combination of these ``k`` vectors. 

We may also say that the set ``S=\left\{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k\right\}`` of vectors is a __spanning set__ for ``V``.
"""

# ╔═╡ 4e16d3d9-668d-4565-95f2-f8c20e73839c
cm"""
$(bth("1 The Span of a Set of Vectors"))
Let ``\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k`` be vectors in the vector space ``V``. Then the set ``W`` of all linear combinations of ``\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k`` is a subspace of ``V``.
$(ebl())

- We sometimes write
```math
W=\operatorname{span}(S)=\operatorname{span}\left\{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k\right\}
```

- ``\mathbb{R}^3=\operatorname{span}\{\mathbf{i}, \mathbf{j}, \mathbf{k}\}``. 
"""

# ╔═╡ 4a6e9768-71f4-4a2b-9746-d27c64d31314
cm"""
$(define("Linear Independence"))
The vectors ``\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k`` in a vector space ``V`` are said to be __linearly independent__ provided that the equation
```math
c_1 \mathbf{v}_1+c_2 \mathbf{v}_2+\cdots+c_k \mathbf{v}_k=\mathbf{0}
```
has only the trivial solution ``c_1=c_2=\cdots=c_k=0``. That is, the only linear combination of ``\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k`` that represents the zero vector ``\mathbf{0}`` is the trivial combination 
```math
0 \mathbf{v}_1+0 \mathbf{v}_2+\cdots+0 \mathbf{v}_k.
```
"""

# ╔═╡ d14cb002-bfcd-49d4-8318-62ca9a2e6521
cm"""
$(ex(4)) 
__The standard unit vectors__

```math
\begin{aligned}
\mathbf{e}_1 & =(1,0,0, \ldots, 0) \\
\mathbf{e}_2 & =(0,1,0, \ldots, 0) \\
& \vdots \\
\mathbf{e}_n & =(0,0,0, \ldots, 1)
\end{aligned}
```
in ``\mathbb{R}^n`` are __linearly independent__. 
"""

# ╔═╡ 9221b160-0335-42e3-830f-35be3df0a1ed
cm"""
$(bbl("Remark","Linear Dependent"))
A set of vectors is called __linearly dependent__ provided it is not linearly independent. Hence the vectors ``\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k`` are linearly dependent if and only if there exist scalars ``c_1, c_2, \ldots, c_k`` not all zero such that
```math
c_1 \mathbf{v}_1+c_2 \mathbf{v}_2+\cdots+c_k \mathbf{v}_k=\mathbf{0}
```

In short, a (finite) set of vectors is linearly dependent provided that some nontrivial linear combination of them equals the zero vector.
"""

# ╔═╡ 23159722-3535-468c-8501-987e3786dcdb
cm"""
$(ex(5))
Determine whether the vectors ``\mathbf{v}_1=(1,2,2,1), \mathbf{v}_2=(2,3,4,1)``, and ``\mathbf{v}_3=(3,8,7,5)`` in ``\mathbb{R}^4`` are linearly independent.
"""

# ╔═╡ 26755d8b-8b1f-43df-8117-a38840b9349b
cm"""
$(remarks())
- Observe that linear independence of the vectors ``\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k`` actually is a property of the set ``S=\left\{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k\right\}`` whose elements are these vectors. Occasionally the phraseology "the set ``S=\left\{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k\right\}`` is linearly independent" is more convenient. 

- Any subset of a linearly independent set ``S=`` ``\left\{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k\right\}`` is a linearly independent set of vectors.

- the coefficients in a linear combination of the linearly independent vectors ``\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k`` are unique. If both
"""

# ╔═╡ 57369806-889e-42cf-9da3-7325ff1176f5
cm"""
$(ex(6))
Let ``\mathbf{v}_1=(2,1,3), \mathbf{v}_2=(5,-2,4), \mathbf{v}_3=(3,8,-6)``, and ``\mathbf{v}_4=(2,7,-4)``. Show that these vectors are linealy depedent.

Then the equation ``c_1 \mathbf{v}_1+c_2 \mathbf{v}_2+c_3 \mathbf{v}_3+c_4 \mathbf{v}_4=\mathbf{0}`` is equivalent to the linear system
```math
\begin{aligned}
2 c_1+5 c_2+3 c_3+2 c_4 & =0 \\
c_1-2 c_2+8 c_3+7 c_4 & =0 \\
3 c_1+4 c_2-6 c_3-4 c_4 & =0
\end{aligned}
```
of three equations in four unknowns. 
"""

# ╔═╡ 87f79eab-0b30-4bc8-8a47-fe027d91e476
cm"""
$(remark())
Any set of more than ``n`` vectors in ``\mathbb{R}^n`` is linearly dependent.
"""

# ╔═╡ 7a6a7ba0-cd23-4d36-b6d7-6af7514d9bb6
cm"""
$(bth("2"))
__Independence of ``\boldsymbol{n}`` Vectors in ``\mathbb{R}^n``__

The ``n`` vectors ``\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_n`` in ``\mathbf{R}^n`` are linearly independent if and only if the ``n \times n`` matrix
```math
\mathbf{A}=\left[\begin{array}{llll}
\mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_n
\end{array}\right]
```
having them as its column vectors has nonzero determinant.
"""

# ╔═╡ 46874b7b-0d64-4a7f-bbbc-4744aef58e55
cm"""
$(bth("3"))
__Independence of Fewer Than ``\boldsymbol{n}`` Vectors in ``\mathbb{R}^n``__

Consider ``k`` vectors ``\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k`` in ``\mathbf{R}^n``, with ``k<n``. Let
```math
\mathbf{A}=\left[\begin{array}{llll}
\mathbf{v}_1 & \mathbf{v}_2 & \cdots & \mathbf{v}_k
\end{array}\right]
```
be the ``n \times k`` matrix having them as its column vectors. Then the vectors ``\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k`` are linearly independent if and only if some ``k \times k`` submatrix of A has nonzero determinant.
"""

# ╔═╡ 6cd1ff75-1846-45d4-a24e-dc5a8d4b5176
cm"""
$(ex())
If theyare linearly independent, showthis; otherwise find a
 non trivial linear combination of them that is equal to the zero
 vector.
```math
\mathbf{v}_1=(2,0,-3), \mathbf{v}_2=(4,-5,-6), \mathbf{v}_3=(-2,1,3)
```
"""

# ╔═╡ e0a788cb-4029-4673-8220-4edf08026ca9
cm"""
$(ex())
 The vectors ``\{v_1,v_2\}`` are known to be linearly
 independent. Apply the definition of linear independence to
 show that the vectors fuig are also linearly independent.
```math
\mathbf{u}_1=\mathbf{v}_1+\mathbf{v}_2, \mathbf{u}_2=2 \mathbf{v}_1+3 \mathbf{v}_2
```

"""

# ╔═╡ de7140da-31e6-45a6-ba94-83f98cfd9dbf
cm"""
$(define("Basis"))
A finite set ``S`` of vectors in a vector space ``V`` is called a basis for ``V`` provided that
- (a) the vectors in ``S`` are linearly independent, and
- (b) the vectors in ``S`` span ``V``.
"""

# ╔═╡ 14f1a87d-56e1-4606-bfed-8b0fa558ab93
cm"""
$(ex(1))
The standard basis for ``\mathbb{R}^n`` consists of the unit vectors
```math
\mathbf{e}_1=(1,0,0, \ldots, 0), \mathbf{e}_2=(0,1,0, \ldots, 0), \ldots, \mathbf{e}_n=(0,0,0, \ldots, 1)
```
"""

# ╔═╡ 215b534d-7de6-4e6e-a20c-7ca936bf3af9
cm"""
$(bbl("Remark",""))

Any set of ``n`` linearly independent vectors in ``\mathbb{R}^n`` is a __basis__ for ``\mathbb{R}^n``.


$(ebl())
"""

# ╔═╡ 53b4b372-cc7f-48d7-930e-8bf007820dae
cm"""
$(ex(3))
Let ``\mathbf{v}_1=(1,-1,-2,-3), \mathbf{v}_2=(1,-1,2,3), \mathbf{v}_3=(1,-1,-3,-2)``, and ``\mathbf{v}_4=(0,3,-1,2)``. Is the set ``\{v_1,v_2,v_3,v_4\}`` a basis for ``\mathbb{R}^4``?
"""

# ╔═╡ 91007be4-b072-44bc-a484-5ace7c30a174
cm"""
$(bth("1 Bases as Maximal Linearly Independent Sets"))
Let ``S=\left\{\mathbf{v}_{\mathbf{1}}, \mathbf{v}_2, \ldots, \mathbf{v}_n\right\}`` be a basis for the vector space ``V``. Then any set of more than ``n`` vectors in ``V`` is linearly dependent.
"""

# ╔═╡ 83e6ef8d-5c04-4805-a9fe-5ef096b1c5a6
cm"""
$(bth("The Dimension of a Vector Space"))
Any two bases for a vector space consist of the same number of vectors.
"""

# ╔═╡ 390905b4-84be-430a-be26-da43007d8de2
cm"""
$(define("Dimension of a Vector Space"))
A nonzero vector space ``V`` is called __finite dimensional__ provided that there exists a basis for ``V`` consisting of a __finite number of vectors__ from ``V``. In this case the number ``n`` of vectors in each basis for ``V`` is called the __dimension of ``V``, denoted by ``n=\operatorname{dim} V``__. 
"""

# ╔═╡ 36b5cd1f-7872-491c-95c9-e0503b7d5ca6
cm"""
$(bbl("Remarks"))
- Note that the zero vector space ``\{\boldsymbol{0}\}`` has no basis because it contains no linearly independent set of vectors. (Sometimes it is convenient to adopt the convention that the null set is a basis for ``\{\boldsymbol{0}\}``.) Here we define ``\operatorname{dim}\{\boldsymbol{0}\}=0``. 
- A nonzero vector space that has no finite basis is called __infinite dimensional__. 
"""

# ╔═╡ a35044ca-ff7d-411c-ad19-b1067dc96b0d
cm"""
$(ex(4))
Let ``\mathcal{P}`` be the set of all polynomials of the form 
```math
p(x)=a_0+a_1 x+a_2 x^2+\cdots+a_n x^n
```
where the largest exponent ``n \geq 0`` that appears is the degree of the polynomial ``p(x)``, and the coefficients ``a_0, a_1, a_2, \ldots, a_n`` are real numbers. 

``\mathcal{P}`` is an infinite-dimensional space.
"""

# ╔═╡ 77fbd41a-f3ac-477e-b4e0-01cb0049362f
cm"""
$(bth("3 Independent Sets, Spanning Sets, and Bases"))
Let ``V`` be an ``n``-dimensional vector space and let ``S`` be a subset of ``V``. Then
- (a) If ``S`` is linearly independent and consists of ``n`` vectors, then ``S`` is a basis for ``V``;
- (b) If ``S`` spans ``V`` and consists of ``n`` vectors, then ``S`` is a basis for ``V``;
- (c) If ``S`` is linearly independent, then ``S`` is contained in a basis for ``V``;
- (d) If ``S`` spans ``V``, then ``S`` contains a basis for ``V``.
"""

# ╔═╡ 50e51cd6-6337-47ee-b405-3b0ec86a9ad3
cm"""
$(bbl("ALGORITHM A Basis for the Solution Space"))
To find a basis for the solution space ``W`` of the homogeneous linear system ``\mathbf{A x}=`` ``\mathbf{0}``, carry out the following steps.
1. Reduce the coefficient matrix ``\mathbf{A}`` to echelon form.
2. Identify the ``r`` leading variables and the ``k=n-r`` free variables. If ``k=0``, then ``W=\{0\}``.
3. Set the free variables equal to parameters ``t_1, t_2, \ldots, t_k``, and then solve by back substitution for the leading variables in terms of these parameters.
4. Let ``\mathbf{v}_j`` be the solution vector obtained by setting ``t_j`` equal to 1 and the other parameters equal to zero. Then ``\left\{\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k\right\}`` is a basis for ``W``.
"""

# ╔═╡ 0a7fda1b-7041-481d-8935-e25e91040e3a
cm"""
$(ex(5))
Find a basis for the solution space of the homogeneous linear system
```math
\begin{aligned}
& 3 x_1+6 x_2-x_3-5 x_4+5 x_5=0 \\
& 2 x_1+4 x_2-x_3-3 x_4+2 x_5=0 \\
& 3 x_1+6 x_2-2 x_3-4 x_4+x_5=0
\end{aligned}
```
"""

# ╔═╡ cd3dcd43-5aa6-4f0d-9e26-2a7e5d5db706
cm"""
$(define("Row Sapce of a Matrix"))
Let ``A\in \mathbb{R}^{m\times n}`` be an ``m\times n`` matrix. The __row space__ of ``A`` (denoted by ``\textrm{Row}(A)``) is the subspace of ``\mathbb{R}^n`` spanned by the ``m`` rows of ``A``. 
The dimension of ``\textrm{Row}(A)`` is called the __row rank__ of ``A``.
$(ebl())

$(define("Column Sapce of a Matrix"))
Let ``A\in \mathbb{R}^{m\times n}`` be an ``m\times n`` matrix. The __column space__ of ``A`` (denoted by ``\textrm{Col}(A)``) is the subspace of ``\mathbb{R}^m`` spanned by the ``n`` columns of ``A``. 
The dimension of ``\textrm{Col}(A)`` is called the __column rank__ of ``A``.
$(ebl())

"""

# ╔═╡ 131a8093-a933-4ffe-b5b9-396cbe6b03e2
cm"""
$(bth("Rank of a Matrix"))
Let ``A\in \mathbb{R}^{m\times n}`` be an ``m\times n`` matrix. The __column rank__ of ``A`` is equal to the __row rank__ of ``A``. So we define the __rank of ``A``__ as
```math
\textrm{rank}(A) = \textrm{dim}(\textrm{Col}(A)) =\textrm{dim}(\textrm{Row}(A))
```
"""

# ╔═╡ 87933565-4591-492c-a819-6136c1b1b547
cm"""
$(ex())
Find ``\textrm{rank}(A)`` of 
```math
A = \begin{bmatrix}
		3& 6 &-1& -5& 5 \\
		2& 4& -1& -3& 2 \\
		3& 6& -2& -4& 1\\

	\end{bmatrix}
```
	
"""

# ╔═╡ 55fb1d51-1582-4465-ae90-d21a4322281b
cm"""
$(bth("1 Principle of Superposition for Homogeneous Equations"))
Let ``y_1`` and ``y_2`` be two solutions of the homogeneous linear equation in (7) on the interval ``I``. If ``c_1`` and ``c_2`` are constants, then the linear combination
```math
y=c_1 y_1+c_2 y_2
```
is also a solution of Eq. (7) on ``I``.
"""

# ╔═╡ 03dcf152-ab37-478f-8130-b53a8047c7b0
cm"""

$(bth("2 Existence and Uniqueness for Linear Equations"))
Suppose that the functions ``p, q``, and ``f`` are continuous on the open interval ``I`` containing the point ``a``. Then, given any two numbers ``b_0`` and ``b_1``, the equation
```math
y^{\prime \prime}+p(x) y^{\prime}+q(x) y=f(x)
```
has a unique (that is, one and only one) solution on the entire interval ``I`` that satisfies the initial conditions
```math
y(a)=b_0, \quad y^{\prime}(a)=b_1
```
"""

# ╔═╡ 21aea514-7f39-41b8-9e51-4139e5d1825d
cm"""
$(ex(2)) Verify that the functions
```math
y_1(x)=e^x \quad \text { and } \quad y_2(x)=x e^x
```
are solutions of the differential equation
```math
y^{\prime \prime}-2 y^{\prime}+y=0
```
"""

# ╔═╡ b35479c0-8873-42d5-9296-8e22323f7e96
cm"""
$(define("Linear Independence of Two Functions"))
Two functions defined on an open interval ``I`` are said to be linearly independent on ``I`` provided that neither is a constant multiple of the other.
"""

# ╔═╡ d3371b50-9c33-4df0-91b2-1c3bd13b7093
cm"""
$(bth("3 Wronskians of Solutions"))
Suppose that ``y_1`` and ``y_2`` are two solutions of the homogeneous second-order linear equation (Eq. (7))
```math
y^{\prime \prime}+p(x) y^{\prime}+q(x) y=0
```
on an open interval ``I`` on which ``p`` and ``q`` are continuous.
- (a) If ``y_1`` and ``y_2`` are linearly dependent, then ``W\left(y_1, y_2\right) \equiv 0`` on ``I``.
- (b) If ``y_1`` and ``y_2`` are linearly independent, then ``W\left(y_1, y_2\right) \neq 0`` at each point of ``I``.

Where ``W(y_1,y_2)`` is called the __Wronskian__ of ``y_1`` and ``y_2`` and is defined as the determinant
```math
W=\left|\begin{array}{cc}
y_1 & y_2 \\
y_1^{\prime} & y_2^{\prime}
\end{array}\right|=y_1 y_2^{\prime}-y_1^{\prime} y_2
```
"""

# ╔═╡ 48739d88-721f-471c-afa6-7cbb2e92ce3d
cm"""
$(bth("4 General Solutions of Homogeneous Equations"))
Let ``y_1`` and ``y_2`` be two linearly independent solutions of the homogeneous equation (Eq. (7))
```math
y^{\prime \prime}+p(x) y^{\prime}+q(x) y=0
```
with ``p`` and ``q`` continuous on the open interval ``I``. If ``Y`` is any solution whatsoever of Eq. (7) on ``I``, then there exist numbers ``c_1`` and ``c_2`` such that
```math
Y(x)=c_1 y_1(x)+c_2 y_2(x)
```
for all ``x`` in ``I``.
"""

# ╔═╡ a7f45d96-911d-470e-a370-7037820a927b
cm"""
$(bth("5 Distinct Real Roots"))
If the roots ``r_1`` and ``r_2`` of the characteristic equation in (8) are real and distinct, then
```math
y(x)=c_1 e^{r_1 x}+c_2 e^{r_2 x}
```
is a general solution of Eq. (9). Thus the solution space of the equation ``a y^{\prime \prime}+`` ``b y^{\prime}+c y=0`` has basis ``\left\{e^{r_1 x}, e^{r_2 x}\right\}``.
"""

# ╔═╡ dbbdb46b-5cb6-4ee2-831e-251f8b39f9b8
cm"""
$(ex(5)) Find the general solution of ``\quad 2 y^{\prime \prime}-7 y^{\prime}+3 y=0``
"""

# ╔═╡ 0354703b-ee11-4d5c-8707-6b8559ca3fac
cm"""
$(bth("6 Repeated Roots"))
If the characteristic equation in (8) has equal (necessarily real) roots ``r_1=r_2``, then
```math
y(x)=\left(c_1+c_2 x\right) e^{r_1 x}
```
is a general solution of Eq. (9). In this case the solution space of the equation ``a y^{\prime \prime}+b y^{\prime}+c y=0`` has basis ``\left\{e^{r_1 x}, x e^{r_1 x}\right\}``.
"""

# ╔═╡ e15959e9-60c8-4084-9adb-012ca2ab1c2e
cm"""
$(ex(7))
To solve the initial value problem
```math
\begin{aligned}
& y^{\prime \prime}+2 y^{\prime}+y=0 \\
& y(0)=5, \quad y^{\prime}(0)=-3
\end{aligned}
```
"""

# ╔═╡ adea2ccd-4182-4317-b949-cd76dc7cddcf
cm"""
$(bth(" 1 Principle of Superposition for Homogeneous Equations"))
Let ``y_1, y_2, \ldots, y_n`` be ``n`` solutions of the homogeneous linear equation in (11) on the interval ``I``. If ``c_1, c_2, \ldots, c_n`` are constants, then the linear combination
```math
y=c_1 y_1+c_2 y_2+\cdots+c_n y_n
```
is also a solution of Eq. (11) on ``I``.
"""

# ╔═╡ 1cd59ab0-6e97-469b-9723-466105b77ef7
cm"""
$(bth(" 2 Existence and Uniqueness for Linear Equations"))
Suppose that the functions ``p_1, p_2, \ldots, p_n``, and ``f`` are continuous on the open interval ``I`` containing the point ``a``. Then, given ``n`` numbers ``b_0, b_1, \ldots, b_{n-1}``, the ``n`` th-order linear equation (Eq. (10))
```math
y^{(n)}+p_1(x) y^{(n-1)}+\cdots+p_{n-1}(x) y^{\prime}+p_n(x) y=f(x)
```
has a unique (that is, one and only one) solution on the entire interval ``I`` that satisfies the ``n`` initial conditions
```math
y(a)=b_0, \quad y^{\prime}(a)=b_1, \quad \ldots, \quad y^{(n-1)}(a)=b_{n-1} .
```
"""

# ╔═╡ 3288d394-9746-450d-aab1-8dc3689b6c4f
cm"""
$(define("Linear Dependence of Functions"))
The ``n`` functions ``f_1, f_2, \ldots, f_n`` are said to be linearly dependent on the interval ``I`` provided that there exist constants ``c_1, c_2, \ldots, c_n`` not all zero such that
```math
c_1 f_1+c_2 f_2+\cdots+c_n f_n=0
```
on ``I``; that is,
```math
c_1 f_1(x)+c_2 f_2(x)+\cdots+c_n f_n(x)=0
```
for all ``x`` in ``I``.
"""

# ╔═╡ 576d55ef-48c6-4785-b431-3b36324650fd
cm"""
$(bbl("Reamrks"))
To show that the functions ``f_1, f_2, \ldots, f_n`` are __linearly independent__ on the interval ``I``, it suffices to show that their Wronskian is nonzero at just one point of ``I``. i.e.
```math
W=\left|\begin{array}{cccc}f_1 & f_2 & \cdots & f_n \\ f_1^{\prime} & f_2^{\prime} & \cdots & f_n^{\prime} \\ \vdots & \vdots & & \vdots \\ f_1^{(n-1)} & f_2^{(n-1)} & \cdots & f_n^{(n-1)}\end{array}\right| \not = 0
```
"""

# ╔═╡ ae990960-82a6-41f6-93e5-62a399b22288
cm"""
$(bth("3 Wronskians of Solutions"))
Suppose that ``y_1, y_2, \ldots, y_n`` are ``n`` solutions of the homogeneous ``n`` th-order linear equation
```math
y^{(n)}+p_1(x) y^{(n-1)}+\cdots+p_{n-1}(x) y^{\prime}+p_n(x) y=0
```
on an open interval ``I``, where each ``p_i`` is continuous. Let
```math
W=W\left(y_1, y_2, \ldots, y_n\right)
```
- (a) If ``y_1, y_2, \ldots, y_n`` are linearly dependent, then ``W \equiv 0`` on ``I``.
- (b) If ``y_1, y_2, \ldots, y_n`` are linearly independent, then ``W \neq 0`` at each point of ``I``.

Thus there are just two possibilities: Either ``W=0`` everywhere on ``I``, or ``W \neq 0`` everywhere on ``I``.
"""

# ╔═╡ f3a35115-28e5-4d79-a3d5-a7e4370ad976
cm"""
$(ex())
Show that the functions ``y_1(x)=e^{-3 x}, y_2(x)=\cos 2 x``, and ``y_3(x)=\sin 2 x``  are linearly independent.
"""

# ╔═╡ 6e905945-6a6a-4bc0-8693-0206e1d920c0
cm"""
$(bth("4 General Solutions of Homogeneous Equations"))
Let ``y_1, y_2, \ldots, y_n`` be ``n`` linearly independent solutions of the homogeneous equation

$(texeq"y^{(n)}+p_1(x) y^{(n-1)}+\cdots+p_{n-1}(x) y^{\prime}+p_n(x) y=0")

on an open interval ``I`` where the ``p_i`` are continuous. If ``Y`` is any solution whatsoever of Eq. (12), then there exist numbers ``c_1, c_2, \ldots, c_n`` such that
```math
Y(x)=c_1 y_1(x)+c_2 y_2(x)+\cdots+c_n y_n(x)
```
for all ``x`` in ``I``.
"""

# ╔═╡ 2c586df8-7dc5-4692-ab09-ec0d42e22904
cm"""
$(bth(" 5 Solutions of Nonhomogeneous Equations"))
Let ``y_p`` be a particular solution of the nonhomogeneous equation in (13) on an open interval ``I`` where the functions ``p_i`` and ``f`` are continuous. Let ``y_1, y_2, \ldots, y_n`` be linearly independent solutions of the associated homogeneous equation in (14). If ``Y`` is any solution whatsoever of Eq. (13) on ``I``, then there exist numbers ``c_1, c_2, \ldots``, ``c_n`` such that
```math
Y(x)=c_1 y_1(x)+c_2 y_2(x)+\cdots+c_n y_n(x)+y_p(x)
```
for all ``x`` in ``I``.
"""

# ╔═╡ 91ccaa9f-65c3-4527-8ff2-48e02381a8d6
cm"""
$(ex())
Solve ``y^{\prime \prime}-4 y=12 x``.

"""

# ╔═╡ f7081366-8968-4b03-96ad-cb13a7da4cb2
cm"""
$(bth(" 1 Distinct Real Roots"))
If the roots ``r_1, r_2, \ldots, r_n`` of the characteristic equation in (3) are real and distinct, then
```math
y(x)=c_1 e^{r_1 x}+c_2 e^{r_2 x}+\cdots+c_n e^{r_n x}
```
is a general solution of Eq. (1). Thus the ``n`` linearly independent functions ``\left\{e^{r_1 x}\right.``, ``\left.e^{r_2 x}, \ldots, e^{r_n x}\right\}`` constitute a basis for the ``n``-dimensional solution space of Eq. (15).
"""

# ╔═╡ 18389ff8-adf5-4d98-98ed-a7589c09c908
cm"""
$(ex(1)) Solve the initial value problem
```math
\begin{aligned}
& y^{(3)}+3 y^{\prime \prime}-10 y^{\prime}=0 \\
& y(0)=7, \quad y^{\prime}(0)=0, \quad y^{\prime \prime}(0)=70
\end{aligned}
```
"""

# ╔═╡ e9a04d99-5857-4bd5-bd7f-1681af50cd2e
cm"""
$(bth("2 Repeated Roots"))
If the characteristic equation in (16) has a repeated root ``r`` of multiplicity ``k``, then the part of a general solution of the differential equation in (15) corresponding to ``r`` is of the form
```math
\left(c_1+c_2 x+c_3 x^2+\cdots+c_k x^{k-1}\right) e^{r x}
```
"""

# ╔═╡ a3e7d794-182d-41df-a914-1b425a4998ff
cm"""
$(ex(2))
Find a general solution of the fifth-order differential equation
```math
9 y^{(5)}-6 y^{(4)}+y^{(3)}=0 \text {. }
```
"""

# ╔═╡ 369dcfa6-ebfd-4a9a-85fc-0a6ecfc059a6
cm"""
$(bth("3 Complex Roots"))
If the characteristic equation in (16) has an unrepeated pair of complex conjugate roots ``a \pm b i`` (with ``b \neq 0`` ), then the corresponding part of a general solution of Eq. (15) has the form
```math
e^{a x}\left(c_1 \cos b x+c_2 \sin b x\right)
```

Thus the linearly independent solutions ``e^{a x} \cos b x`` and ``e^{a x} \sin b x`` (corresponding to the complex conjugate characteristic roots ``a \pm b i`` ) generate a 2-dimensional subspace of the solution space of the differential equation.
"""

# ╔═╡ e343be8d-2697-4b5e-9787-659f9ddc53bb
cm"""
$(ex()) Find the particular solution of
```math
y^{\prime \prime}-4 y^{\prime}+5 y=0
```
for which ``y(0)=1`` and ``y^{\prime}(0)=5``.
"""

# ╔═╡ 4d647c9b-50da-4998-b2f2-cebce39ff874
cm"""
$(ex())
Find a general solution of ``y^{(4)}+4 y=0``.
"""

# ╔═╡ abf14a2e-e06e-47cc-a000-f78eae4ec377
cm"""
$(ex())Find a general solution of ``\left(D^2+6 D+13\right)^2 y=0``.
"""

# ╔═╡ b0e322b7-1ec6-4f22-b5ed-204cd0d25f2c
cm"""
$(ex())
The roots of the characteristic equation of a certain differential equation are 

```math
3,-5,0,0,0,0, \quad -5,2 \pm 3 i, \quad \text{and}\quad 2 \pm 3 i.
```

Write a general solution of this homogeneous differential equation.
"""

# ╔═╡ 0afe6221-5946-498b-9c1a-6b34d731f403
cm"""
$(ex(1))
Find a particular solution of ``y^{\prime \prime}+3 y^{\prime}+4 y=3 x+2``.
"""

# ╔═╡ ef47af20-62c9-4d29-bb94-c076e3d0c87e
cm"""
$(ex(2))
Find a particular solution of ``y^{\prime \prime}-4 y=2 e^{3 x}``.
"""

# ╔═╡ a3e3bcd6-63f1-4e1f-9f33-45d6567ac47b
cm"""
$(ex(3))
Find a particular solution of ``3 y^{\prime \prime}+y^{\prime}-2 y=2 \cos x``.
"""

# ╔═╡ 8fc8c2f1-efd8-43e0-b0a3-3f107c78add9
cm"""
$(ex(4)) 
Find a particular solution of ``y^{\prime \prime}-4 y=2 e^{2 x}``.
"""

# ╔═╡ 15ed1ab6-c1e4-4be8-bc71-7e2b37c678b4
cm"""
$(bbl("Remark","The General Approach"))
The method of undetermined coefficients applies whenever the function ``f(x)`` in Eq. (17) is a linear combination of (finite) products of functions of the following three types:
1. A polynomial in ``x``;
2. An exponential function ``e^{r x}``;
3. ``\cos k x`` or ``\sin k x``.

Any such function, for example,
```math
f(x)=\left(3-4 x^2\right) e^{5 x}-4 x^3 \cos 10 x
```
"""

# ╔═╡ a53c0017-5bd9-4730-b272-6f5b64978c30
cm"""
$(bbl("RULE 1","Method of Undetermined Coefficients"))

Suppose that 

__no term appearing either in ``f(x)`` or in any of its derivatives satisfies the associated homogeneous equation ``L y=0``__.

Then take as a trial solution for ``y_p`` a linear combination of all linearly independent such terms and their derivatives. Then determine the coefficients by substitution of this trial solution into the nonhomogeneous equation ``L y=f(x)``.
"""

# ╔═╡ 3f7b2e7d-ea10-487c-a19c-a5769848cbd3
cm"""
$(ex(5))
Find a particular solution of ``y^{\prime \prime}+4 y=3 x^3``.
"""

# ╔═╡ b925ecfd-7745-4843-9245-a1191e1505b6
cm"""
$(ex(6))
Solve the initial value problem
```math
\begin{aligned}
& y^{\prime \prime}-3 y^{\prime}+2 y=3 e^{-x}-10 \cos 3 x \\
& y(0)=1, \quad y^{\prime}(0)=2
\end{aligned}
```
"""

# ╔═╡ 5889801d-3c1a-45e3-89f0-f6d1e52fbd2b
cm"""
$(ex(7))
Find the general form of a particular solution of
```math
y^{(3)}+9 y^{\prime}=x \sin x+x^2 e^{2 x}
```
"""

# ╔═╡ 81dd2f88-0e8c-4b67-beb7-f234984c22f3
cm"""
$(bbl("RULE 2","Method of Undetermined Coefficients"))

If the function ``f(x)`` is of either form in (20), take as the trial solution
```math
\begin{aligned}
y_p(x)= & x^s\left[\left(A_0+A_1 x+A_2 x^2+\cdots+A_m x^m\right) e^{r x} \cos k x\right. \\
& \left.+\left(B_0+B_1 x+B_2 x^2+\cdots+B_m x^m\right) e^{r x} \sin k x\right]
\end{aligned}
```
where ``s`` is the smallest nonnegative integer such that no term in ``y_p`` duplicates a term in the complementary function ``y_c``. Then determine the coefficients in Eq. (15) by substituting ``y_p`` into the nonhomogeneous equation.

"""

# ╔═╡ 5e05ad97-9f9a-45d7-b211-3056ddee07f3
cm"""
$(ex(8))
Find a particular solution of
```math
y^{(3)}+y^{\prime \prime}=3 e^x+4 x^2 \text {. }
```
"""

# ╔═╡ 5405ebd1-dbd9-42d1-a542-1d8c9f4b9aaa
cm"""
$(ex(9))
Determine the appropriate form for a particular solution of
```math
y^{\prime \prime}+6 y^{\prime}+13 y=e^{-3 x} \cos 2 x
```
"""

# ╔═╡ fa4b457e-6833-49b2-a4ee-db18905061d7
cm"""
$(ex(10))
Determine the appropriate form for a particular solution of the fifth-order equation
```math
(D-2)^3\left(D^2+9\right) y=x^2 e^{2 x}+x \sin 3 x
```
"""

# ╔═╡ 0fec7bee-9c5a-4288-b8e0-568cbe0b0466
cm"""
$(bth("1 Variation of Parameters"))
If the nonhomogeneous equation ``y^{\prime \prime}+P(x) y^{\prime}+Q(x) y=f(x)`` has complementary function ``y_c(x)=c_1 y_1(x)+c_2 y_2(x)``, then a particular solution is given by
```math
y_p(x)=-y_1(x) \int \frac{y_2(x) f(x)}{W(x)} d x+y_2(x) \int \frac{y_1(x) f(x)}{W(x)} d x,
```
where ``W=W\left(y_1, y_2\right)`` is the Wronskian of the two independent solutions ``y_1`` and ``y_2`` of the associated homogeneous equation.
"""

# ╔═╡ 35574f43-874e-4c68-9830-85bbf452296d
cm"""
$(ex(11))
Find a particular solution of the equation ``y^{\prime \prime}+y=\tan x``.
"""

# ╔═╡ bc149a7f-b9ff-4724-814f-e7e18e3dea86
cm"""
$(define("Eigenvalues and Eigenvectors"))
The number ``\lambda`` is said to be an __eigenvalue__ of the ``n \times n`` matrix ``\mathbf{A}`` provided there exists a __nonzero vector__ ``\mathbf{v}`` such that
```math
\mathbf{A v}=\lambda \mathbf{v},
```
in which case the vector ``\mathbf{v}`` is called an __eigenvector__ of the matrix ``\mathbf{A}``. We also say that the eigenvector ``\mathbf{v}`` is associated with the eigenvalue ``\lambda``, or that the eigenvalue ``\lambda`` corresponds to the eigenvector ``\mathbf{v}``.
"""

# ╔═╡ 82353e45-9faf-4818-8c04-c2c41cc69866
cm"""
$(bbl("Remarks",""))
The system 
```math
\mathbf{A v}=\lambda \mathbf{v},
```
has a nontrivial solution ``\mathbf{v} \neq \mathbf{0}`` if and only if the determinant
```math
\operatorname{det}(\mathbf{A}-\lambda \mathbf{I})=|\mathbf{A}-\lambda \mathbf{I}|
```
of its coefficient matrix is zero. The equation 
```math 
|\mathbf{A}-\lambda \mathbf{I}|=0
``` 
is called the __characteristic equation__ of the square matrix ``\mathbf{A}``, and we have proved that there exists an eigenvector ``\mathbf{v}`` associated with ``\lambda`` if and only if ``\lambda`` satisfies this equation.
"""

# ╔═╡ 8ce55646-dfa7-47df-8164-d90fc824d477
cm"""
$(bth("1 The Characteristic Equation"))
The number ``\lambda`` is an eigenvalue of the ``n \times n`` matrix ``\mathbf{A}`` if and only if ``\lambda`` satisfies the characteristic equation
```math
|\mathbf{A}-\lambda \mathbf{I}|=0 .
```
"""

# ╔═╡ b87b52a1-dfce-4891-abb4-8583cdd23744
cm"""
$(ex())
Find the eigenvalues and associated eigenvectors of the matrix
```math
\mathbf{A}=\left[\begin{array}{rr}
5 & 7 \\
-2 & -4
\end{array}\right] .
```
"""

# ╔═╡ 2718bf8e-7b0c-42b3-913e-a18f16d22d0b
cm"""
$(ex())
Find the eigenvalues and associated eigenvectors of the matrix
```math
\mathbf{A}=\left[\begin{array}{rr}0 & 8 \\ -2 & 0\end{array}\right]
```
"""

# ╔═╡ 9a981f7a-9c28-4af0-ad2a-ff6f0c3ec2ea
cm"""
$(ex()) 
Find the eigenvalues and associated eigenvectors of the ``I_2``.
"""

# ╔═╡ 16d4239c-a7c0-42ce-8640-8afeffecd47d
cm"""
$(ex())
Find the eigenvalues and associated eigenvectors of the matrix
```math
\mathbf{A}=\left[\begin{array}{ll}2 & 3 \\ 0 & 2\end{array}\right]
```
"""

# ╔═╡ a4ac83bb-c85a-4a6f-86b7-7594043e9294
cm"""
$(bbl("Remark",""))
These examples illustrate the four possibilities for a ``2 \times 2`` matrix ``\mathbf{A}``. It can have either
- two distinct real eigenvalues, each corresponding to a single eigenvector;
- a single real eigenvalue corresponding to a single eigenvector;
- a single real eigenvalue corresponding to two linearly independent eigenvectors; or
- two complex conjugate eigenvalues corresponding to complex conjugate eigenvectors.

$(bbl("Remark",""))
Substitution of ``\lambda=0`` in the characteristic equation ``|\mathbf{A}-\lambda \mathbf{I}|=0`` yields ``|\mathbf{A}|=0``. Therefore, ``\lambda=0`` is an eigenvalue of the matrix ``\mathbf{A}`` if and only if ``\mathbf{A}`` is singular: ``|\mathbf{A}|=0``.
"""

# ╔═╡ 31061f53-97f2-456c-b5b5-8e4f3229e0f4
cm"""
$(ex())Find bases for the eigenspaces of the matrix
```math
\mathbf{A}=\left[\begin{array}{rrr}
4 & -2 & 1 \\
2 & 0 & 1 \\
2 & -2 & 3
\end{array}\right]
```
"""

# ╔═╡ a78f4ee9-0ddc-47b8-8252-4c6a1fe0073d
cm"""
Recall:

$(ex())
Find the eigenvalues and associated eigenvectors of the matrix
```math
\mathbf{A}=\left[\begin{array}{rr}
5 & 7 \\
-2 & -4
\end{array}\right] .
```
"""

# ╔═╡ 8d5ea72f-66be-4a22-811f-1405a6b724c8
cm"""
$(define("Similar Matrices"))
The ``n \times n`` matrices ``\mathbf{A}`` and ``\mathbf{B}`` are called similar provided that there exists an invertible matrix ``\mathbf{P}`` such that
```math
\mathbf{B}=\mathbf{P}^{-1} \mathbf{A P}
```
"""

# ╔═╡ 122b69e6-2e6d-4338-81bb-c7c050d0cc8a
cm"""
$(bbl("Remarks",""))

An ``n \times n`` matrix ``\mathbf{A}`` is called __diagonalizable__ if it is similar to a diagonal matrix ``\mathbf{D}``; that is, there exist a diagonal matrix ``\mathbf{D}`` and an invertible matrix ``\mathbf{P}`` such that ``\mathbf{A}=\mathbf{P D P}^{-1}``, and so
```math
\mathbf{P}^{-1} \mathbf{A P}=\mathbf{D} .
```

The process of finding the __diagonalizing__ matrix ``\mathbf{P}`` and the diagonal matrix ``\mathbf{D}`` is called __diagonalization__ of the matrix A. 

In Example 1 we showed that the matrix ``A`` is __diagonalizable__.
"""

# ╔═╡ 9ededb67-b1a0-4ee9-ac2f-f2bfbdb9e772
cm"""
$(bth("1 Criterion for Diagonalizability"))
The ``n \times n`` matrix ``\mathbf{A}`` is diagonalizable if and only if it has ``n`` linearly independent eigenvectors.
"""

# ╔═╡ 6ec677b8-aae1-4a49-861f-ee85e94f872d
cm"""
$(ex())
The matrix 
```math
\mathbf{A}=\left[\begin{array}{ll}2 & 3 \\ 0 & 2\end{array}\right]
```
is NOT diagonalizable. 
"""

# ╔═╡ d7f27e6b-b367-4375-a4ee-1ffaee828fba
cm"""
$(ex())Find bases for the eigenspaces of the matrix
```math
\mathbf{A}=\left[\begin{array}{rrr}
4 & -2 & 1 \\
2 & 0 & 1 \\
2 & -2 & 3
\end{array}\right]
```
Is ``A`` __diagonalizable__?
"""

# ╔═╡ 4569a153-4e04-4a6b-b379-adf9f3ded21a
cm"""
$(bth("2 Eigenvectors Associated with Distinct Eigenvalues"))
Suppose that the eigenvectors ``\mathbf{v}_1, \mathbf{v}_2, \ldots, \mathbf{v}_k`` are associated with the distinct eigenvalues ``\lambda_1, \lambda_2, \ldots, \lambda_k`` of the matrix ``\mathbf{A}``. Then these ``k`` eigenvectors are linearly independent.
"""

# ╔═╡ 1cbf11da-3e9f-4bfb-b298-8212bdedd95b
cm"""
$(bth("3"))
If the ``n \times n`` matrix ``\mathbf{A}`` has ``n`` distinct eigenvalues, then it is diagonalizable.
"""

# ╔═╡ 6c44073c-b617-4d7f-83f6-458d8830328d
cm"""
$(bth("4 Complete Independence of Eigenvectors"))
Let ``\lambda_1, \lambda_2, \ldots, \lambda_k`` be the distinct eigenvalues of the ``n \times n`` matrix ``\mathbf{A}``. For each ``i=1,2, \ldots, k``, let ``S_i`` be a basis for the eigenspace associated with ``\lambda_i``. Then the union ``S`` of the bases ``S_1, S_2, \ldots, S_k`` is a linearly independent set of eigenvectors of ``\mathbf{A}``.
"""

# ╔═╡ c447fe70-1e3d-4c90-aebd-7b6ed096bf70
cm"""
$(ex())
Is the matrix
```math
\mathbf{A}=\left[\begin{array}{rrr}3 & 0 & 0 \\ -4 & 6 & 2 \\ 16 & -15 & -5\end{array}\right]
```
diagonalizable?
"""

# ╔═╡ da9230a6-088d-4735-b206-9514c12dd223
initialize_eqref()

# ╔═╡ 107407c8-5da0-4833-9965-75a82d84a0fb
@htl("""
<style>
@import url("https://mmogib.github.io/math102/custom.css");

ul {
  list-style: none;
}

ul li:before {
  content: '💡 ';
}

.p40 {
	padding-left: 40px;
}
</style>

</style>
""")

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
CommonMark = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlotThemes = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoExtras = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
QRCoders = "f42e9828-16f3-11ed-2883-9126170b272d"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[compat]
Colors = "~0.12.11"
CommonMark = "~0.8.15"
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.3.1"
Latexify = "~0.16.5"
PlotThemes = "~3.2.0"
Plots = "~1.40.8"
PlutoExtras = "~0.7.13"
PlutoUI = "~0.7.60"
PrettyTables = "~2.4.0"
QRCoders = "~1.4.5"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.1"
manifest_format = "2.0"
project_hash = "9d4264e62aea4a6dbff5896fd89dd0d573ffedaa"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

    [deps.AbstractFFTs.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "009060c9a6168704143100f36ab08f06c2af4642"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.2+1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "bce6804e5e6044c6daab27bb533d1295e4a2e759"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.6"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b5278586822443594ff615963b0c09755771b3e0"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.26.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.CommonMark]]
deps = ["Crayons", "PrecompileTools"]
git-tree-sha1 = "3faae67b8899797592335832fccf4b3c80bb04fa"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.15"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "ea32b83ca4fefa1768dc84e504cc0a94fb1ab8d1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.2"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"
weakdeps = ["IntervalSets", "LinearAlgebra", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.Extents]]
git-tree-sha1 = "81023caa0021a41712685887db1fc03db26f41f5"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.4"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "53ebe7511fa11d33bec688a9178fac4e49eeee00"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "82d8afa92ecf4b52d78d869f038ebfb881267322"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.3"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "b5c7fe9cea653443736d264b85466bad8c574f4a"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.9.9"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "532f9126ad901533af1d4f5c198867227a7bb077"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+1"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "629693584cef594c3f6f99e76e7a7ad17e60e8d5"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.7"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a8863b69c2a0859f2c2c87ebdc4c6712e88bdf0d"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.7+0"

[[deps.GeoFormatTypes]]
git-tree-sha1 = "59107c179a586f0fe667024c5eb7033e81333271"
uuid = "68eda718-8dee-11e9-39e7-89f7f65f511f"
version = "0.4.2"

[[deps.GeoInterface]]
deps = ["Extents", "GeoFormatTypes"]
git-tree-sha1 = "2f6fce56cdb8373637a6614e14a5768a88450de2"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.7"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "b62f2b2d76cee0d61a2ef2b3118cd2a3215d3134"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.11"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "43ba3d3c82c18d88471cfd2924931658838c9d8f"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.0+4"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "674ff0db93fffcd11a3573986e550d66cd4fd71f"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.5+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "d1d712be3164d61d1fb98e7ce9bcbc6cc06b45ed"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.8"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "401e4f3f30f43af2c8478fc008da50096ea5240f"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.3.1+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "2e4520d67b0cef90865b3ef727594d2a58e0e1f8"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.11"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "437abb322a41d527c197fa800455f79d414f0a3c"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.8"

[[deps.ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils", "Libdl", "Pkg", "Random"]
git-tree-sha1 = "5bc1cb62e0c5f1005868358db0692c994c3a13c6"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.2.1"

[[deps.ImageMagick_jll]]
deps = ["Artifacts", "Ghostscript_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d65554bad8b16d9562050c67e7223abf91eaba2f"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "6.9.13+0"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "355e2b974f2e3212a75dfb60519de21361ad3cb7"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.9"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0936ba688c6d201805a83da835b55c61a180db52"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.11+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "d1b1b796e47d94588b3757fe84fbf65a5ec4a80d"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "39d64b09147620f5ffbf6b2d3255be3c901bec63"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.8"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "f389674c99bfcde17dc57454011aa44d5a260a40"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "fa6d0bcff8583bac20f1ffa708c3913ca605c611"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.5"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "25ee0be4d43d0269027024d75a24c24d6c6e590c"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.4+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "854a9c268c43b77b0a27f22d7fab8d33cdb3a731"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+1"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll"]
git-tree-sha1 = "fa7fd067dca76cadd880f1ca937b4f387975a9f5"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.16.0+0"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.MarchingCubes]]
deps = ["PrecompileTools", "StaticArrays"]
git-tree-sha1 = "27d162f37cc29de047b527dab11a826dd3a650ad"
uuid = "299715c1-40a9-479a-aaf9-4a633d36f717"
version = "0.1.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
git-tree-sha1 = "1a27764e945a152f7ca7efa04de513d473e9542e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.1"

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

    [deps.OffsetArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "8292dd5c8a38257111ada2174000a33745b06d4e"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.2.4+0"

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "libpng_jll"]
git-tree-sha1 = "f4cb457ffac5f5cf695699f82c537073958a6a6c"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.5.2+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7493f61f55a6cce7325f197443aa80d32554ba10"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.15+1"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "67186a2bc9a90f9f85ff3cc8277868961fb57cbd"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.3"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e127b609fb9ecba6f201ba7ab753d5a605d53801"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.1+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "6e55c6841ce3411ccb3457ee52fc48cb698d6fb0"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.2.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "45470145863035bb124ca51b320ed35d071cc6c2"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.8"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoExtras]]
deps = ["AbstractPlutoDingetjes", "DocStringExtensions", "HypertextLiteral", "InteractiveUtils", "Markdown", "PlutoUI", "REPL", "Random"]
git-tree-sha1 = "681f89bdd5c1da76b31a524af798efb5eb332ee9"
uuid = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
version = "0.7.13"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "8f6bc219586aef8baf0ff9a5fe16ee9c70cb65e4"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.2"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QRCoders]]
deps = ["FileIO", "ImageCore", "ImageIO", "ImageMagick", "StatsBase", "UnicodePlots"]
git-tree-sha1 = "b3e5fcc7a7ade2d43f0ffd178c299b7a264c268a"
uuid = "f42e9828-16f3-11ed-2883-9126170b272d"
version = "1.4.5"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["PrecompileTools"]
git-tree-sha1 = "98ca7c29edd6fc79cd74c61accb7010a4e7aee33"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.6.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eeafab08ae20c62c44c8399ccb9354a04b80db50"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.7"

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

    [deps.StaticArrays.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a6b1675a536c5ad1a60e5a5153e1fee12eb146e3"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.0"

[[deps.StructArrays]]
deps = ["ConstructionBase", "DataAPI", "Tables"]
git-tree-sha1 = "f4dc295e983502292c4c3f951dbb4e985e35b3be"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.18"

    [deps.StructArrays.extensions]
    StructArraysAdaptExt = "Adapt"
    StructArraysGPUArraysCoreExt = "GPUArraysCore"
    StructArraysSparseArraysExt = "SparseArrays"
    StructArraysStaticArraysExt = "StaticArrays"

    [deps.StructArrays.weakdeps]
    Adapt = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "657f0a3fdc8ff4a1802b984872468ae1649aebb3"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.10.1"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.UnicodePlots]]
deps = ["ColorTypes", "Contour", "Crayons", "Dates", "FileIO", "FreeTypeAbstraction", "LazyModules", "LinearAlgebra", "MarchingCubes", "NaNMath", "Printf", "SparseArrays", "StaticArrays", "StatsBase", "Unitful"]
git-tree-sha1 = "ae67ab0505b9453655f7d5ea65183a1cd1b3cfa0"
uuid = "b8865327-cd53-5732-bb35-84acbb429228"
version = "2.12.4"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d95fe458f26209c66a187b1114df96fd70839efd"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.21.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "1165b0443d0eca63ac1e32b8c0eb69ed2f4f8127"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "555d1076590a6cc2fdee2ef1469451f872d8b41b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "936081b536ae4aa65415d869287d43ef3cb576b2"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.53.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "b70c870239dc3d7bc094eb2d6be9b73d27bef280"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.44+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "7dfa0fd9c783d3d0cc43ea1af53d69ba45c447df"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─71bc54d5-d0ed-42d3-9bc1-48aa86e91d1d
# ╟─e414122f-b93a-4510-b8ae-026c303e0df9
# ╟─8408e369-40eb-4f9b-a7d7-26cde3e34a74
# ╟─cd269caf-ef81-43d7-a1a8-6668932b6363
# ╟─cb0eda5b-9dcb-421f-b6b6-ed4b51797df7
# ╟─20bd193d-c0f7-47e6-8b46-33e07338ad91
# ╟─e6c9f44f-19dd-4922-bc05-8614f5441e80
# ╟─a1cc7232-8dff-4fb4-908a-01c105e62797
# ╟─0853d79c-5ae0-4690-9f09-9068505ba213
# ╟─2b5ca5c2-4232-48df-95ab-4a5765436995
# ╟─555e1ddd-fcd2-4ff8-af06-f4caafc8dcff
# ╟─3dab7c5b-363a-4b55-b8b0-342723d99ba5
# ╟─d5117856-b9ca-4ed1-976e-9e20d613bfc1
# ╟─55221819-6a40-4f85-93b0-76c5d05ed35e
# ╟─a7abd027-69c5-43a9-a99e-164a9199a334
# ╟─a9d46b08-d25a-49df-9f93-427e521055e0
# ╟─fbd68049-57ac-42fc-9c4c-5bd9f019dabb
# ╟─41344730-e8c5-4e2b-a98c-9995846c7244
# ╟─22715277-7264-44aa-a2be-0cf6bba6b24f
# ╟─07d56dc4-4a09-43eb-8ad6-8cb796b53c80
# ╟─2d42c9ba-bc48-4836-aa56-bca41304a300
# ╟─b97a6f03-dbeb-40b9-be58-f9a9a50cac5a
# ╟─ab0da87e-9b87-48eb-a235-44fdfd2f81f2
# ╟─6168170e-2fcb-4728-b8d8-ddc44992d3f9
# ╠═a1d00dd2-59b0-427b-b2c0-25ec94e039a9
# ╠═8781e87a-146a-41b0-8981-243fec51bfa3
# ╠═261637c6-4da1-4b4a-9f01-3e2324c41b60
# ╟─953e31c8-3e26-47cf-a675-2067802ba941
# ╟─75287b3b-3f90-4a0a-88d0-92130f84c0db
# ╟─15531160-b5d7-4e55-848c-9b239d4f116c
# ╟─1790195b-439e-4472-96e8-6f87d8ff0601
# ╟─085418be-7e83-4a6a-b3a8-fbd4d2ac2451
# ╟─6d3c4d0e-a950-43bd-a2d7-46aec4417ab3
# ╟─7bacb885-372a-4ee9-ba9c-de505d332dfb
# ╟─241b2797-3050-400b-a532-303a6feeb39e
# ╟─a5a9fda5-c077-4505-8fcf-2e23f9b3ad1f
# ╟─809d228a-65a3-4af5-b605-e030c41c4548
# ╟─a1af7fe8-cb18-4b36-a015-38df3d346b33
# ╟─6b4ce7dc-39f3-4c4e-973f-b501916f0717
# ╟─0756191d-272f-4e62-8fed-1f731820325c
# ╟─f832a1e4-17ff-43d8-97a4-c7e29bbe542f
# ╟─a2cc59fc-01d3-45c6-80e7-3c8485632036
# ╟─ceb445e8-64ca-4d6f-86c2-98c96ad42fb2
# ╟─72cef7b7-f466-4872-bba0-6f46ac36de15
# ╟─1cfb4072-1e58-4b7a-bcdc-672cd0183f75
# ╟─07b46453-655c-4fc0-8590-3104ce5f6c99
# ╠═27f89662-ef45-4469-bbc7-7a6b1949881c
# ╟─124fc4b5-157b-46b1-8a77-06b119ee0a4d
# ╟─c8d24a6a-b851-420f-9c4f-e42eeb934fc3
# ╟─177d7ac6-2a24-4282-84a3-e42a1c03251f
# ╠═c96bcf16-8467-4f8f-a55a-6769f73143da
# ╟─e9d05d41-2d51-4f11-ad43-95571b99b466
# ╟─e1c8a21a-8f19-4632-992f-014d4690baa6
# ╟─16d641fc-590f-4ab2-9288-c0f03849caec
# ╟─cd7eb646-c23c-4ab8-9dc4-04c0db21d4bb
# ╠═ae81e505-f95a-4355-bf74-b308445d2830
# ╟─2e8d8477-a937-4cd4-b300-2cae6b9b2853
# ╟─e0c6b9b1-ab42-4e6f-a41e-440e68d68041
# ╟─1d54e85b-f3e6-41f1-a818-ab483c918fee
# ╟─4fbe0c02-a31f-4f37-b014-84e4bb914fff
# ╟─5aa66e9f-5741-4306-aa36-8a28627bca0f
# ╟─c868de03-803a-49a9-b1f5-657f46ab8498
# ╟─419f6c6f-1aff-4b72-be83-3dc3ef01bc76
# ╟─90a2fcd1-1d38-4c12-8f7e-cb9f83f890f3
# ╟─aca4f2f0-5c16-4e51-b83c-b76149f836a9
# ╟─d743b750-2a29-4fb3-a729-07bec9637121
# ╟─cacc6877-5310-4f5a-a0c2-344cccd0ec8b
# ╟─564258db-41ca-41fd-90c8-13c8fd39fe93
# ╟─7a405477-d020-4b31-b0be-036fd3f22322
# ╟─0e5e9fe4-0157-4497-b675-53071f4e6782
# ╟─4ccc21f1-f253-4e7d-87e4-4c8f5e1db785
# ╟─a32e62d2-e65f-492a-a3f4-db8dec72db67
# ╟─8792fe08-498e-4d18-aad2-de9bc3a7ade4
# ╟─c006d3aa-8cf1-4382-806a-4e9714f934f6
# ╟─626fe783-d05a-49ca-8e3c-2ac76be27e34
# ╟─ad29db84-3c84-4a6e-aeb0-5a51210ef05d
# ╟─933c6345-fffc-4159-8b67-e1443b988f9f
# ╠═9ed8b8ed-1ac5-4e8b-8a74-17dff32f1ad5
# ╟─80e20e84-51d2-4e6e-8eb3-b7a8dac22a8b
# ╟─65b3d353-815c-4e5a-a887-24a82bdf005d
# ╟─a4256712-b454-4733-bcf5-c5764705b028
# ╟─87802187-2aeb-4e67-a2d7-1975cf588a17
# ╟─608e2182-ce7c-4626-bcda-bf4bd6f2f1c6
# ╠═8d103106-40ac-49fa-86ae-1af3ecb8a3ec
# ╟─c7cc5a14-f964-4cf9-82f0-f23904bbdace
# ╟─2678d172-8896-4681-b9af-45e3485f2312
# ╟─27dbde5f-5e4d-4c67-bc67-ae0b9c40a774
# ╟─a3567a11-b2c2-4545-b4a2-278b4b7e4b0e
# ╟─643d7441-1dbf-4b98-b76f-617babece0c0
# ╟─2aa7e38b-c48d-470c-943e-a6a73d053a70
# ╟─13748f1f-4092-4592-9a4f-8f17c38a875a
# ╟─2d3b73f6-8429-4551-96ac-b2db865189ed
# ╟─ab51c417-0482-4c68-ba07-38fd46e535cd
# ╟─3e062725-d0ac-4d5e-9baa-49b1b10de464
# ╟─34820013-08f9-4077-9403-206671bea915
# ╟─6f434bf9-93f6-4052-aed0-7ebe0c06b736
# ╟─22bb573e-8fa7-4f6c-b5f6-b96dee4c8507
# ╟─6e605916-523b-466b-ad96-c206af3bf27c
# ╠═aa04cd10-62d5-49a9-8b60-7e5d928ec1ee
# ╟─387953af-89c5-47ee-a1b1-af9b05baf949
# ╟─bf7f3424-92fd-4ece-8aeb-db8da6762056
# ╟─c35aca68-0b50-430f-b1e1-bc2ab187fa40
# ╟─eb2a3144-7dca-42d0-9fce-01c78728eac7
# ╠═747384a8-c06a-49b7-92f1-0cdd44043fb2
# ╠═a2a039ee-2097-4569-9d85-e96597ff47da
# ╟─6189764d-9da8-460b-b873-a040c413e721
# ╟─9d011957-4432-484f-99dc-da94395aeef1
# ╟─25e77707-5852-48f8-9849-46e3f6841bb6
# ╟─f5bf2679-c765-4b80-bd27-021a82683611
# ╟─fb9cd8c0-ace5-4123-a8a0-58ab32f849f1
# ╟─c63574f6-2f03-4f3d-b9eb-08ca321a0aad
# ╟─737365c0-472d-4f0d-87c5-5d70a2fe027f
# ╟─b583732f-ed42-4a4f-9dec-6369b5e85c41
# ╟─8dae3e8d-96e3-4fe0-8b9f-5074dff0f9e4
# ╟─aff3151b-4b66-4d55-8e27-744a53f5a91e
# ╟─1ba6a9f9-f0b3-4033-a20a-117b19d7c74b
# ╠═babc3fb7-0eaa-4af6-b7d1-cd0b0fd4c78b
# ╟─ea830a4d-f11a-474e-9b42-dde6295e8d24
# ╠═ea5ca7e4-2f30-45f2-8715-7e502a274a53
# ╟─b0efeb12-e0fd-4a1e-80b9-24c230368c53
# ╟─e1199b37-1692-4d97-ba86-0b5239c316d1
# ╟─42e01a69-848d-4d34-a5de-8c54ca3f3b83
# ╟─e948494b-8a14-47b8-bfec-4390376e1bcd
# ╠═ece946d9-b094-4e01-8dda-47740ec418b3
# ╟─30795897-1e43-41ed-ba85-c3b89ecdef1f
# ╟─52d0d0d2-3079-4285-92de-c59961848921
# ╟─aaae2aff-fc65-4b7f-9e28-b7577c5fdd51
# ╟─c233fecb-8525-493f-a088-011d5386ec1d
# ╟─d5618ff1-3420-4c92-bcfa-238b4da78f17
# ╟─ddda1623-3b11-46f7-9fda-fd6686dd4343
# ╟─c31e6c51-9f45-4581-a7b6-54bdf2e04ada
# ╟─856ff140-9c2b-4953-a033-849f1af34286
# ╟─96a09993-9f1a-4908-9334-494b4e198455
# ╟─6a3e0f87-707c-4122-a5fd-dcef8fce8a7d
# ╟─dff6e998-431d-4493-bc34-6fa9b5de6f93
# ╟─f61386e6-e32c-45f6-b1a9-dec901195e96
# ╟─20916835-0034-40b5-ab98-1f24fb78247c
# ╟─f868d290-9a4e-4a6f-933e-ee642deabeac
# ╟─20e79ac6-f0f4-40e6-89e0-2eaa44e0ff2b
# ╟─b7430962-79fb-427c-a23e-8cbb09ae8e46
# ╟─f6653849-9370-4f84-a07a-ad924983fda6
# ╟─4faab984-e4a5-4620-9c7e-61abd525a1fa
# ╟─cf4fc0c1-c9d8-4597-b1d7-a5d65b59421d
# ╟─e0261d6d-2bc3-4eab-8eb2-ef799f0cece7
# ╟─d76bbff8-1730-41cb-93d4-ddb35f66b19b
# ╟─d711bb6f-f627-4568-97f5-5d7690088663
# ╟─263e7015-3130-4597-b440-a996e52967b9
# ╟─64e0c10c-a4e8-4435-86fc-60bff5e476d4
# ╟─2af29504-36c2-4f20-a325-c32acf62d9bc
# ╟─d9a91130-0cf8-414e-b8a8-dd352090f7a6
# ╟─170f2255-078f-4c9e-86e9-0eb93b2ff87e
# ╟─f1464582-faf2-4902-96a6-bb1cf17c1102
# ╟─f3d58858-f9ab-42bb-a822-14d03916984e
# ╟─6827e93a-43b1-4ee7-a650-573a3ce61676
# ╟─4b2763bd-cf78-4913-8d5d-ac10d08b1031
# ╟─abefd094-b57a-46c7-a314-eb836e58c6a0
# ╟─f8cd8873-f4fb-4300-8aa2-4aee9990cf5e
# ╟─e17f0845-12a9-4da8-b81d-bf3e61d95260
# ╟─b13e61f1-3d38-4f87-838f-5c448ce1193a
# ╟─add8e82e-02cb-4262-b8e4-ac7a8c317257
# ╟─2f5bc349-791e-4ece-8933-f2baef3c59a1
# ╟─56aee5a8-caf4-4cba-9e1e-7e319c21047b
# ╟─f110354f-ab42-481d-9057-fd467fa9a664
# ╟─4e16d3d9-668d-4565-95f2-f8c20e73839c
# ╟─ca3afed1-ded8-4792-8253-58a2b5bc0502
# ╟─4a6e9768-71f4-4a2b-9746-d27c64d31314
# ╟─d14cb002-bfcd-49d4-8318-62ca9a2e6521
# ╟─9221b160-0335-42e3-830f-35be3df0a1ed
# ╟─23159722-3535-468c-8501-987e3786dcdb
# ╠═4fc61c32-377d-475f-9587-fb884aaafcb6
# ╟─26755d8b-8b1f-43df-8117-a38840b9349b
# ╟─57369806-889e-42cf-9da3-7325ff1176f5
# ╟─87f79eab-0b30-4bc8-8a47-fe027d91e476
# ╟─7a6a7ba0-cd23-4d36-b6d7-6af7514d9bb6
# ╟─46874b7b-0d64-4a7f-bbbc-4744aef58e55
# ╟─6cd1ff75-1846-45d4-a24e-dc5a8d4b5176
# ╠═73ea96b9-be3c-4e26-8f54-914bd53d409b
# ╟─e0a788cb-4029-4673-8220-4edf08026ca9
# ╟─faf0e3cf-3973-4272-bda3-a3410c8e17d2
# ╟─de7140da-31e6-45a6-ba94-83f98cfd9dbf
# ╟─14f1a87d-56e1-4606-bfed-8b0fa558ab93
# ╟─215b534d-7de6-4e6e-a20c-7ca936bf3af9
# ╟─53b4b372-cc7f-48d7-930e-8bf007820dae
# ╠═f911e13b-7830-449d-b8d2-ae6ecb388e5c
# ╟─91007be4-b072-44bc-a484-5ace7c30a174
# ╟─83e6ef8d-5c04-4805-a9fe-5ef096b1c5a6
# ╟─390905b4-84be-430a-be26-da43007d8de2
# ╟─36b5cd1f-7872-491c-95c9-e0503b7d5ca6
# ╟─a35044ca-ff7d-411c-ad19-b1067dc96b0d
# ╟─77fbd41a-f3ac-477e-b4e0-01cb0049362f
# ╟─dca3f82b-7cac-4849-88a2-d47805957566
# ╟─50e51cd6-6337-47ee-b405-3b0ec86a9ad3
# ╟─0a7fda1b-7041-481d-8935-e25e91040e3a
# ╟─03e906c2-29b0-487b-9223-6c7313c1d740
# ╟─ee16937c-617a-4f4a-a15e-c3853f0e8ca8
# ╟─cd3dcd43-5aa6-4f0d-9e26-2a7e5d5db706
# ╟─131a8093-a933-4ffe-b5b9-396cbe6b03e2
# ╟─87933565-4591-492c-a819-6136c1b1b547
# ╠═d999ba91-2a73-4d8e-a99a-fd5187f24e6b
# ╟─28fb5089-e83b-47f5-93be-60f1f0ea1e30
# ╟─82a55302-3cf3-4c7e-861c-12f5f554f142
# ╟─c9ac1dc2-078e-47ae-b0cf-9f22eb294b39
# ╟─55fb1d51-1582-4465-ae90-d21a4322281b
# ╟─03dcf152-ab37-478f-8130-b53a8047c7b0
# ╟─21aea514-7f39-41b8-9e51-4139e5d1825d
# ╟─98730a31-e9e9-45fc-ae12-5cf0832b5069
# ╟─b35479c0-8873-42d5-9296-8e22323f7e96
# ╟─f9a9f5f5-5134-4f2f-8672-c2e665c546ce
# ╟─d3371b50-9c33-4df0-91b2-1c3bd13b7093
# ╟─48739d88-721f-471c-afa6-7cbb2e92ce3d
# ╟─06cbf42d-9811-40ae-bcd5-1bcc4ff60880
# ╟─c60351f6-a57c-46db-a3ea-3deee1f9180b
# ╠═381710b7-c487-4b65-9a07-88d773cdacf0
# ╟─a7f45d96-911d-470e-a370-7037820a927b
# ╟─dbbdb46b-5cb6-4ee2-831e-251f8b39f9b8
# ╟─0354703b-ee11-4d5c-8707-6b8559ca3fac
# ╟─e15959e9-60c8-4084-9adb-012ca2ab1c2e
# ╟─90d2b2d2-c030-4afd-a529-cd940bc82688
# ╟─8695eb73-c264-4f4a-ab66-2ce879ff8c88
# ╟─adea2ccd-4182-4317-b949-cd76dc7cddcf
# ╟─dc4b275f-9014-4014-bc59-03353d5de046
# ╟─1cd59ab0-6e97-469b-9723-466105b77ef7
# ╟─cc282ef0-67d4-4086-9916-49ee02ee8009
# ╟─3288d394-9746-450d-aab1-8dc3689b6c4f
# ╟─576d55ef-48c6-4785-b431-3b36324650fd
# ╟─ae990960-82a6-41f6-93e5-62a399b22288
# ╟─f3a35115-28e5-4d79-a3d5-a7e4370ad976
# ╟─0f267953-3b92-4550-af48-3bdc51364a26
# ╟─6e905945-6a6a-4bc0-8693-0206e1d920c0
# ╟─87c2749f-1390-4ce1-97f1-3a9609cba34e
# ╟─42a14a32-71aa-4019-8eaf-e42f7f51e7a9
# ╟─2c586df8-7dc5-4692-ab09-ec0d42e22904
# ╟─91ccaa9f-65c3-4527-8ff2-48e02381a8d6
# ╟─bcaf9ef2-6371-4de1-8254-b4ae57f77363
# ╟─fb34e1b6-1cd2-4f0f-a8c2-988dd3d709d0
# ╟─fcd85a36-c264-43a4-8eb6-8bd025582ac2
# ╟─f7081366-8968-4b03-96ad-cb13a7da4cb2
# ╟─18389ff8-adf5-4d98-98ed-a7589c09c908
# ╟─92ebabd1-d3df-47e9-9322-5c415cf1dc7b
# ╟─f9d43fec-1397-4206-9136-adbedd8d7fd8
# ╟─3e22f633-ef1c-45b4-aa83-3e28c61a9365
# ╟─e9a04d99-5857-4bd5-bd7f-1681af50cd2e
# ╟─a3e7d794-182d-41df-a914-1b425a4998ff
# ╟─a03566a1-b19b-4bfe-84e6-7850c401e299
# ╟─eb6f270c-9e5b-41fa-bfd0-8244d8d811f1
# ╟─96210d10-9d4f-4aea-9b76-a8e066f34e71
# ╟─369dcfa6-ebfd-4a9a-85fc-0a6ecfc059a6
# ╟─e343be8d-2697-4b5e-9787-659f9ddc53bb
# ╟─4d647c9b-50da-4998-b2f2-cebce39ff874
# ╟─44f2e591-8d5f-41a6-b0b2-66b3a7d83cf2
# ╟─08dec10f-1df7-4249-88a2-fc52be252e7e
# ╟─abf14a2e-e06e-47cc-a000-f78eae4ec377
# ╟─b0e322b7-1ec6-4f22-b5ed-204cd0d25f2c
# ╟─9a964597-f45a-4e1e-a4a5-c6d0007b3e13
# ╟─c0461d53-1446-468f-9ae5-babc4b2bab32
# ╟─0ee542fb-a1d2-46c8-9604-681485015268
# ╟─0afe6221-5946-498b-9c1a-6b34d731f403
# ╟─ef47af20-62c9-4d29-bb94-c076e3d0c87e
# ╟─a3e3bcd6-63f1-4e1f-9f33-45d6567ac47b
# ╟─8fc8c2f1-efd8-43e0-b0a3-3f107c78add9
# ╟─15ed1ab6-c1e4-4be8-bc71-7e2b37c678b4
# ╟─a53c0017-5bd9-4730-b272-6f5b64978c30
# ╟─3f7b2e7d-ea10-487c-a19c-a5769848cbd3
# ╟─b925ecfd-7745-4843-9245-a1191e1505b6
# ╟─5889801d-3c1a-45e3-89f0-f6d1e52fbd2b
# ╟─270a3760-25c6-4cb3-96c2-95a2b98506f9
# ╟─ae2c1de9-d67a-4c4f-a87d-90250e436fc8
# ╟─81dd2f88-0e8c-4b67-beb7-f234984c22f3
# ╟─c0f7a0ed-0a89-4e5d-933d-f261c952a48c
# ╟─5e05ad97-9f9a-45d7-b211-3056ddee07f3
# ╟─5405ebd1-dbd9-42d1-a542-1d8c9f4b9aaa
# ╟─fa4b457e-6833-49b2-a4ee-db18905061d7
# ╟─b8db6d0c-7200-4c76-9fed-3a774ce8dae9
# ╟─0fec7bee-9c5a-4288-b8e0-568cbe0b0466
# ╟─35574f43-874e-4c68-9830-85bbf452296d
# ╟─274e341c-66c4-434a-b002-079130d6ee81
# ╟─bc149a7f-b9ff-4724-814f-e7e18e3dea86
# ╟─41799b6d-d000-40c9-8144-211645bf6a3e
# ╟─82353e45-9faf-4818-8c04-c2c41cc69866
# ╟─8ce55646-dfa7-47df-8164-d90fc824d477
# ╟─b87b52a1-dfce-4891-abb4-8583cdd23744
# ╟─2718bf8e-7b0c-42b3-913e-a18f16d22d0b
# ╟─9a981f7a-9c28-4af0-ad2a-ff6f0c3ec2ea
# ╟─16d4239c-a7c0-42ce-8640-8afeffecd47d
# ╟─a4ac83bb-c85a-4a6f-86b7-7594043e9294
# ╟─50c06b0f-3a4e-4f91-ad85-5b814dffa8bb
# ╟─2e116819-87fb-4fef-892e-3fbd19f3862e
# ╟─31061f53-97f2-456c-b5b5-8e4f3229e0f4
# ╟─92297d90-f326-4f1e-80fc-5465621e8a6b
# ╟─a78f4ee9-0ddc-47b8-8252-4c6a1fe0073d
# ╠═050d88a5-3b65-4940-b6b2-6750eaff6d66
# ╟─b7e7238a-b2ec-4632-85b3-c069217b7ecc
# ╟─8d5ea72f-66be-4a22-811f-1405a6b724c8
# ╟─122b69e6-2e6d-4338-81bb-c7c050d0cc8a
# ╟─9ededb67-b1a0-4ee9-ac2f-f2bfbdb9e772
# ╟─6ec677b8-aae1-4a49-861f-ee85e94f872d
# ╟─d7f27e6b-b367-4375-a4ee-1ffaee828fba
# ╠═343d8a6b-c6bd-41ef-ab28-3872753864c6
# ╟─4569a153-4e04-4a6b-b379-adf9f3ded21a
# ╟─1cbf11da-3e9f-4bfb-b298-8212bdedd95b
# ╟─6c44073c-b617-4d7f-83f6-458d8830328d
# ╟─c447fe70-1e3d-4c90-aebd-7b6ed096bf70
# ╠═f2d4c2a5-f486-407b-b31b-d2efcc7476b3
# ╟─ef081dfa-b610-4c7a-a039-7258f4f6e80e
# ╟─da9230a6-088d-4735-b206-9514c12dd223
# ╟─107407c8-5da0-4833-9965-75a82d84a0fb
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
