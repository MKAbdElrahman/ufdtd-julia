# Program 3.6 : One-dimensional FDTD program with an TF/SF source, ABC and nonhomogenous media .
# John B. Schneider
# M.K AbdElRahman

# New Changes:
# 1. The Domain begins and ends with an electeric field, thus the magnetic field array is
# one element smaller than the electric field array.
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
hy = zeros(Real, NCells - 1)

ϵᵣ = ones(Real, NCells)
μᵣ = ones(Real, NCells - 1)

# build device
for x = 1:NCells
    if(x < 100)
        ϵᵣ[x] = 1.0
    else
        ϵᵣ[x] = 9.0
    end
end

ϵᵣ_inv = inv.(ϵᵣ)
μᵣ_inv = inv.(μᵣ)

plot(ϵᵣ,
    label = L"εᵣ",
    xlabel = L"x",
    color = :black,
    linewidth = 3,
    yrange = [0, 10],
    framestyle = :box,
    widen = false,
)


function apply_ABC!(field)
    field[1] = field[2] # left ABC
    field[end] = field[end-1] # right ABC
end

function update_hy!(hy,ez)
    @inbounds for x in eachindex(hy)
        hy[x] = hy[x] + (ez[x+1] - ez[x]) * μᵣ_inv[x]
    end
end
function update_ez!(ez,hy)
    @inbounds for x = 2:length(ez)-1
        ez[x] = ez[x] + (hy[x] - hy[x-1]) * ϵᵣ_inv[x]
    end
end

animation = @animate for t = 1:NTimeSteps

    update_hy!(hy,ez)
    hy[TFSFSourceNode] -= exp(-(t - PulseOffset)^2 / PulseWidth)

    apply_ABC!(ez)

    update_ez!(ez,hy)
    ez[TFSFSourceNode+1] += exp(-(t + 1.0 - PulseOffset)^2 / PulseWidth)


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
gif(animation, "animations/TFSF_source_and_ABC_NonHomgenous.gif", fps = 15)
