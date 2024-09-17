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
	Au[1,:],Au[2,:]=Au[2,:],Au[1,:]
	Au[2,:]=-2Au[1,:]+Au[2,:]
	Au[3,:]=-7Au[1,:]+Au[3,:]
	Au[3,:]=-2Au[2,:]+Au[3,:]

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
	inv(A)
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

# ╔═╡ e1199b37-1692-4d97-ba86-0b5239c316d1
md"## Inverses and the Adjoint Matrix"

# ╔═╡ ece946d9-b094-4e01-8dda-47740ec418b3
let
	A = [
		1  4  5 
		4  2  5
		-3  3  -1
	]
	inv(Rational.(A))
end

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
SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[compat]
Colors = "~0.12.11"
CommonMark = "~0.8.12"
HypertextLiteral = "~0.9.5"
LaTeXStrings = "~1.3.1"
Latexify = "~0.16.5"
PlotThemes = "~3.2.0"
Plots = "~1.40.5"
PlutoExtras = "~0.7.12"
PlutoUI = "~0.7.59"
PrettyTables = "~2.3.2"
QRCoders = "~1.4.5"
SymPy = "~2.2.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.5"
manifest_format = "2.0"
project_hash = "f4703e0dff485bac3fb1464187004fe2ccfdbd03"

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
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

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
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

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

[[deps.CommonEq]]
git-tree-sha1 = "6b0f0354b8eb954cdba708fb262ef00ee7274468"
uuid = "3709ef60-1bee-4518-9f2f-acd86f176c50"
version = "0.2.1"

[[deps.CommonMark]]
deps = ["Crayons", "JSON", "PrecompileTools", "URIs"]
git-tree-sha1 = "532c4185d3c9037c0237546d817858b23cf9e071"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.12"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

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

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "b19db3927f0db4151cb86d073689f2428e524576"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.10.2"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d8a9c0b6ac2d9081bf76324b39c78ca3ce4f0c98"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.6"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
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
git-tree-sha1 = "94997910aca72897524d2237c41eb852153b0f65"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.3"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

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
git-tree-sha1 = "3f74912a156096bd8fdbef211eff66ab446e7297"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+0"

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

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "801aef8228f7f04972e596b09d4dba481807c913"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.4"

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
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

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
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

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
git-tree-sha1 = "a53ebe394b71470c7f97c2e7e170d51df21b17af"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.7"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

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
git-tree-sha1 = "c84a835e1a09b289ffcd2271bf2a337bbdda6637"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.3+0"

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
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

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
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

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

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

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

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

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
version = "0.3.23+4"

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
git-tree-sha1 = "a028ee3cb5641cccc4c24e90c36b0a4f7707bdf5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

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
git-tree-sha1 = "cb5a2ab6763464ae0f19c86c56c63d4a2b0f5bda"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.52.2+0"

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
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

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
git-tree-sha1 = "082f0c4b70c202c37784ce4bfbc33c9f437685bf"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.5"

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

[[deps.PlutoDevMacros]]
deps = ["AbstractPlutoDingetjes", "DocStringExtensions", "HypertextLiteral", "InteractiveUtils", "MacroTools", "Markdown", "Pkg", "Random", "TOML"]
git-tree-sha1 = "c3839362a712e6d9c2845d179edafe74371cb77b"
uuid = "a0499f29-c39b-4c5c-807c-88074221b949"
version = "0.7.4"

[[deps.PlutoExtras]]
deps = ["AbstractPlutoDingetjes", "HypertextLiteral", "InteractiveUtils", "Markdown", "PlutoDevMacros", "PlutoUI", "REPL"]
git-tree-sha1 = "93d8c75734da9192d0639406fe6fb446be0fba4f"
uuid = "ed5d0301-4775-4676-b788-cf71e66ff8ed"
version = "0.7.12"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

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
git-tree-sha1 = "66b20dd35966a748321d3b2537c4584cf40387c7"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.3.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "8f6bc219586aef8baf0ff9a5fe16ee9c70cb65e4"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.2"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "9816a3826b0ebf49ab4926e2b18842ad8b5c8f04"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.96.4"

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
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

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
git-tree-sha1 = "2803cab51702db743f3fda07dd1745aadfbf43bd"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.5.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

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

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

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
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

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
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

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

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.SymPy]]
deps = ["CommonEq", "CommonSolve", "LinearAlgebra", "PyCall", "SpecialFunctions", "SymPyCore"]
git-tree-sha1 = "d35b297be048dfac05bcff29e55d6106808e3c5a"
uuid = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
version = "2.2.0"

[[deps.SymPyCore]]
deps = ["CommonEq", "CommonSolve", "Latexify", "LinearAlgebra", "Markdown", "RecipesBase", "SpecialFunctions"]
git-tree-sha1 = "bef92ec4c31804bdc9c44cb00eaf0348eac383fb"
uuid = "458b697b-88f0-4a86-b56b-78b75cfb3531"
version = "0.2.5"

    [deps.SymPyCore.extensions]
    SymPyCoreTermInterfaceExt = "TermInterface"

    [deps.SymPyCore.weakdeps]
    TermInterface = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"

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

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "SIMD", "UUIDs"]
git-tree-sha1 = "bc7fd5c91041f44636b2c134041f7e5263ce58ae"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.10.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "e84b3a11b9bece70d14cce63406bbc79ed3464d2"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.2"

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

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

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

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

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
git-tree-sha1 = "d9717ce3518dc68a99e6b96300813760d887a01d"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.1+0"

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
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

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
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

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
version = "1.52.0+1"

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
# ╟─747384a8-c06a-49b7-92f1-0cdd44043fb2
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
# ╟─ea830a4d-f11a-474e-9b42-dde6295e8d24
# ╟─b0efeb12-e0fd-4a1e-80b9-24c230368c53
# ╟─e1199b37-1692-4d97-ba86-0b5239c316d1
# ╟─42e01a69-848d-4d34-a5de-8c54ca3f3b83
# ╟─e948494b-8a14-47b8-bfec-4390376e1bcd
# ╟─ece946d9-b094-4e01-8dda-47740ec418b3
# ╠═f2d4c2a5-f486-407b-b31b-d2efcc7476b3
# ╟─ef081dfa-b610-4c7a-a039-7258f4f6e80e
# ╟─da9230a6-088d-4735-b206-9514c12dd223
# ╟─107407c8-5da0-4833-9965-75a82d84a0fb
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
