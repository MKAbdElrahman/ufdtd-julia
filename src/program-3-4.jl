# Program 3.4 1Dadditive.c: One-dimensional FDTD program with an additive source.
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
    # use additive source at node 50
    ez[50] += exp(-(t - 30)^2 / 100)

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

gif(animation, "animations/Bare-bones 1D FDTD simulation with an additive source.gif", fps = 15)