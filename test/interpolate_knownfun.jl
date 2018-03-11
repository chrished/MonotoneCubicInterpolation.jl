"""
:linear                yuck
:finitedifference      classic cubic interpolation, no tension parameter
:cardinal              cubic cardinal splines, uses tension parameter which must be between [0,1]
:fritschcarlson        monotonic - tangents are first initialized, then adjusted if they are not monotonic
:fritschbutland        monotonic - faster algorithm (only requires one pass) but somewhat higher apparent "tension"
"""
methods = [:linear, :finitedifference, :cardinal, :fritschcarlson, :fritschbutland]

tension = 0.

xs = linspace(0., 1., 11)
func(x) = (x-0.5)
ys = func(xs)

xeval = collect(linspace(0.,1.,101))
truevals = func(xeval)

for method in methods
    res = interpolateCubicHermite(xeval, xs, ys, method, tension)
    dist = maximum(abs.(res - truevals))
    @test isapprox(dist,0., atol = eps())
end

func2(x) = (x-0.5).^2
ys = func(xs)
xeval = collect(linspace(0.,1.,101))
truevals = func(xeval)

for method in methods[2:end]
    res = interpolateCubicHermite(xeval, xs, ys, method, tension)
    dist = maximum(abs.(res - truevals))
    @test isapprox(dist,0., atol = eps())
end


func3(x) = (x-0.5).^3
ys = func(xs)
xeval = collect(linspace(0.,1.,101))
truevals = func(xeval)

for method in methods[2:end]
    res = interpolateCubicHermite(xeval, xs, ys, method, tension)
    dist = maximum(abs.(res - truevals))
    @test isapprox(dist,0., atol = eps())
end
