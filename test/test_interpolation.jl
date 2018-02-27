using Plots
# test cubic Interpolation
include("monotone_cubic_interpolation.jl")

xs = linspace(0., 1., 11)
func(x) = (x-0.5).^2
ys = func(xs)

xeval = collect(linspace(0.,1.,101))

res = interpolateCubicHermite(xeval, xs, ys,"FritschCarlson", .5)

truevals = func(xeval)
dist = res - truevals
plot(xeval, truevals)
plot!(xeval, res)
