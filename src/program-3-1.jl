# Program 3.1 1DbareBones.c: Bare-bones one-dimensional simulation with a hard source.
# John B. Schneider
# M.K AbdElRahman

using Plots
using LaTeXStrings

const Real = Float32

const η₀ = Real(377.0)
const inv_η₀ = inv(η₀)


NCells = 200
NTimeSteps = 500



ez = zeros(Real, NCells)
hy = zeros(Real, NCells)

animation = @animate for t = 1:NTimeSteps
    # update magnetic field
    @inbounds for x = 1:NCells-1 
        hy[x] = hy[x] + (ez[x+1] - ez[x])
        # hy[NCells] remains zero --> PMC
    end 
    # update electeric field
    @inbounds for x = 2:NCells
        # ez[1] remains zero --> PEC
        ez[x] = ez[x] + (hy[x] - hy[x-1])
    end
    # hardwire a source node
    ez[1] = exp(-(t - 30)^2 / 100)

    plot(
        ez,
        label = L"E_z",
        xlabel = L"x",
        color = :red,
        linewidth = 3,
        yrange = [-1.5, 1.5],
        framestyle = :box,
        widen = false,
        title = "t = $(t)",
    )
    plot!(hy, label = L"H_y", color = :blue, linewidth = 3)
end every 5

gif(animation, "animations/hard_source.gif", fps = 15)