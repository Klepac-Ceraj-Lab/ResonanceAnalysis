using Resonance

kneadread = Resonance.load_knead()

DataFrames.transform!(kneadread, 
    AsTable(Cols(r"raw")) => ByRow(sum) => :raw,
    AsTable(Cols(r"trimmed")) => ByRow(sum) => :trimmed,
    AsTable(Cols(r"decontaminated")) => ByRow(sum) => :decontaminated,
    AsTable(Cols(r"final")) => ByRow(sum) => :final
)

DataFrames.transform!(kneadread, 
    AsTable([:raw, :trimmed]) => ByRow(x-> ismissing(x[2]) ? 0 : (x[1] - x[2]) / x[1] * 100) => :perc_trimmed,
    AsTable([:raw, :decontaminated]) => ByRow(x-> ismissing(x[2]) ? 0 : (x[1] - x[2]) / x[1] * 100) => :perc_decontaminated
)

kneadlong = stack(select(kneadread, [:Sample, :raw, :trimmed, :final]), [:raw, :trimmed, :final])

##

using CairoMakie
using AlgebraOfGraphics

##

fig = Figure(; resolution=(1100, 600))
lims = extrema(skipmissing(kneadread.raw))

ax1 = Axis(fig[1,1]; xlabel="Raw reads", ylabel="Final reads", yscale=log10, xscale=log10)
ax2 = Axis(fig[1,3]; xlabel="Raw reads", ylabel="Final reads", yscale=log10, xscale=log10)

sc1 = scatter!(ax1, kneadread.raw, kneadread.final; color = kneadread.perc_trimmed)
lines!(ax1, Point2f.(zip(lims, lims)); color = :gray)

sc2 = scatter!(ax2, kneadread.raw, kneadread.final; color = kneadread.perc_decontaminated, colormap=:plasma)
lines!(ax2, Point2f.(zip(lims, lims)); color = :gray)

Colorbar(fig[1,2], sc1; label = "Percent trimmed")
Colorbar(fig[1,4], sc2; label = "Percent human")

save("figures/kneaddata_counts.png", fig)
fig

#-

