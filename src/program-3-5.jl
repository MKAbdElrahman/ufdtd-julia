# Program 3.5 : One-dimensional FDTD program with an TF/SF source and ABC.
# John B. Schneider
# M.K AbdElRahman
using Plots
using LaTeXStrings

const Real = Float32

const η₀ = Real(377.0)
const inv_η₀ = inv(η₀)


NCells = 200
NTimeSteps = 500
TFSFSourceNode = 50
PulseOffset = 30
PulseWidth = 100


ez = zeros(Real, NCells)
hy = zeros(Real, NCells)

animation = @animate for t = 1:NTimeSteps
    # update magnetic field

    hy[NCells] = hy[NCells-1] # right ABC
    @inbounds for x = 1:NCells-1
        hy[x] = hy[x] + (ez[x+1] - ez[x])
    end
    # correction for Hy adjacent to TFSF boundary 
    hy[TFSFSourceNode] -= exp(-(t - PulseOffset)^2 / PulseWidth)
    
    # update electeric field

    ez[1] = ez[2] # left ABC
    @inbounds for x = 2:NCells
        ez[x] = ez[x] + (hy[x] - hy[x-1])
    end

    # correction for Ez adjacent to TFSF boundary
    ez[TFSFSourceNode+1] += exp(-(t + 1.0 - PulseOffset)^2 / PulseWidth)

    plot(ez,
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

gif(
    animation,
    "animations/TFSF_source_and_ABC.gif",
    fps = 15)
