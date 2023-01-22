using Revise
using Plots
include(joinpath(@__DIR__,"../interpolation.jl"))

xs = range(0, 2π, length=10)
y = zeros(28)
xq = 1:3:30
y[1:3:30] = sin.(xs)
sticks(1:length(y), y)
ys = interp(Linear(), y, xq, exterp=true)

sticks(0:length(y)-1, y, marker=(:circle, 8))
plot!(0:length(ys)-1, ys)

ys = interp(Square(), y, xq, exterp=true)
sticks(range(0,1,length=length(y)), y, marker=(:circle, 8))
plot!(range(0,1,length=length(y)), ys)