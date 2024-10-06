### A Pluto.jl notebook ###
# v0.19.46

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
	using SymPy
    using QRCoders
    using PrettyTables
	using LinearSolve
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
	dydx(x) = 3(1-4x^2)
	x = -0.5:0.01:0.5
	# y = x .+ dydx.(x)
	# 3x-4x.^3)
	y = zeros(length(x));
	p = plot(xlimits=(-0.55,0.55),ylimits=(-10,221), 
		frame_style=:origin, legend=:topleft,
		xticks=-0.5:0.2:0.5
	)
	vline!(p,[-0.5,0],c=:black,label=nothing)
	vline!(p,[0.5,0],c=:black,label=nothing)
	annotate!(p,[(0.52,10.5,text(L"x=\frac{1}{2}",10))])
	annotate!(p,[(-0.52,10.5,text(L"x=-\frac{1}{2}",10))])
	anim = @animate for i in 2:length(x)
		y[i] = y[i-1] + dydx(x[i])
		plot(p,[(x[1:i],y[1:i]),([x[i]],[y[i]])], 
			seriestype=[:line  :scatter],
			label=["y:: swimmer trajectory" "swimmer"],
			markershapes = [:none :star]
		)
		
	end
	gif(anim,"imgs/swimmer.gif",fps=10)
	
end

# ╔═╡ 8781e87a-146a-41b0-8981-243fec51bfa3
function vector_field(xs,ys,df; c=:black, args...)
	# xs = -2:0.3:2
	# ys = -2:0.2:2
	# df(x, y) = normalize([1, 1/x]) ./ 10

	xxs = repeat(xs',length(xs),length(ys))
	yys = reshape(repeat(ys,length(xs)),length(xs),length(ys))

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
	t = -log(0.63)/k
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
	P₀=6
	Pprime0 = 0.000212 * (365.25)
	k = Pprime0/P₀

	P_2051 = P₀*exp(k*51)

	P10=floor(log(10)/k)
	1999+Int(P10)
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
	C = 375-50
	# (A-125)
	k = -log((A-125)/325)/75
	T = 150
	t = -log((A-T)/C)/k
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
	df(x,y) = normalize([1,y+(11/8)*exp(-x/3)])
	p1 = vector_field(xs,ys,df)
	p1 = plot(p1,x->(1/32)*(exp(x)-33*exp(-x/3)),c=:blue,
		xlims=(min(xs...),max(xs...)),
		ylims=(min(ys...),max(ys...)),
		label=:none
	)
	plot(p1,[0],[-1], seriestype=:scatter,label=:none)
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
	df(x,y) = 0.1*normalize([1,(y^3-6*x*y)/(4y+3x^2-3x*y^2)])
	plt = vector_field(x,y,df; lw=0.01)
	contour!(x,y,(x,y)->3x^2*y-x*y^3+2y^2, c=:red)
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
	@syms x::Real, y()
	D = Differential(x)
	Eq =2*x*exp(2*y(x))*D(y(x)) ~ 3x^4+exp(2*y(x))
	dsolve(Eq)
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
	@syms a::Real, b::Real, c::Real
	A = [2 -1 3;1 2 1;7 4 9]
	B =[a;b;c]
	Au=hcat(A,B)
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
		1  1  5 
		1  4  13 
		3  2  12
	]
	B=inv(A)

end

# ╔═╡ a2a039ee-2097-4569-9d85-e96597ff47da
let
	A = [4 3 2;5 6 3;3 5 2]
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
	A1 =[1 2 3
		 0 4 4
		0 6 7
	]
	A2 =[1 2 3
		 4 0 0
		0 6 7
	]
	B = [1 2 3;4 4 4;0 6 7]
	det(B),det(A1),det(A2)
end

# ╔═╡ ea5ca7e4-2f30-45f2-8715-7e502a274a53
let
	A =[
		1  2  3  4 
		0  5  6  7 
		0  0  8  9 
		2  4  6  9
	]
end

# ╔═╡ e1199b37-1692-4d97-ba86-0b5239c316d1
md"## Inverses and the Adjoint Matrix"

# ╔═╡ ece946d9-b094-4e01-8dda-47740ec418b3
let
	A = [
		1  4  5 
		4  2  5
		-3  3  -1
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
	A=[1 2 3;2 3 8;2 4 8;1 1 5]
	b= zeros(3)
	# A\b
	solve()
end

# ╔═╡ 73ea96b9-be3c-4e26-8f54-914bd53d409b
let
	v1=[2;0;-3]
	v2=[4,-5,-6]
	v3=[-2,1,3]
	A=[v1 v2 v3]
	A[3,:]=(3//2)A[1,:]+A[3,:]
	A
end

# ╔═╡ faf0e3cf-3973-4272-bda3-a3410c8e17d2
md"# 4.4 Bases and Dimension for Vector Spaces"

# ╔═╡ ef081dfa-b610-4c7a-a039-7258f4f6e80e
begin
	function add_space(n=1)
		repeat("&nbsp;",n)
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
	ex() = example("Example","")
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
# ╠═f2d4c2a5-f486-407b-b31b-d2efcc7476b3
# ╟─ef081dfa-b610-4c7a-a039-7258f4f6e80e
# ╟─da9230a6-088d-4735-b206-9514c12dd223
# ╟─107407c8-5da0-4833-9965-75a82d84a0fb
