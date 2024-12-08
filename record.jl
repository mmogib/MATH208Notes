L(xs)=(x,i)->prod((x-xs[k])/(xs[i]-xs[k]) for k in 1:length(xs) if k!=i)
xs=[0, t, 1]
@syms t::Real
xs=[0, t, 1]
xs=[0, t, 2]
f(x)=sqrt(x-2x^2)
L(xs)=(x,i)->prod((x-xs[k])/(xs[i]-xs[k]) for k in 1:length(xs) if k!=i)
P2 = P(xs,f)
P(xs,f)=x->sum(Li(x,i)*f(xs[i]) for i in 1:length(xs))
P2 = P(xs,f)
Eq = f(1)-P2(1) ~ -0.5
f(x)=sqrt(2x-x^2)
P2 = P(xs,f)
Eq = f(1)-P2(1) ~ -0.5
Li=L(xs)
Eq = f(1)-P2(1) ~ -0.5
print(ans)
Eq = f(1)-P2(1) ~ -0.5
Eq = f(1)-P2(1) ~ -1
print(ans)
function get_REPL_as_history()           history = Base.active_repl.interface.modes[1].hist           history.history[history.start_idx+1:end] end
function get_REPL_as_history()
history = Base.active_repl.interface.modes[1].hist
history.history[history.start_idx+1:end]
end
function save_REPL_history(file)
isfile(file) && error("file already exists")
open(file,"w") do io
txt = join(get_REPL_as_history(),"\n")
write(io,txt)
end
end
get_REPL_as_history()
save_REPL_history("./record.jl")