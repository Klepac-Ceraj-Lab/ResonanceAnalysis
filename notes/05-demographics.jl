using Resonance
omni, tps = startup(; dfs = [:omni, :tps])

using CairoMakie
using AlgebraOfGraphics

hist(collect(skipmissing(tps.ageMonths)))

##

labels = [
    "White",
    "Black",
    "Asian",
    "Mixed",
    "Other",
    "Unknown"
]

counts = [
    "White"   =>   1142,
    "Black"   =>   290,
    "Asian"   =>   1,
    "Mixed"   =>   343,
    "Other"   =>   30,
    "Unknown" =>   594
]

colors = [
    :navajowhite3,
    :coral,
    :darkorange3,
    :darkseagreen3,
    :cyan4,
    :dimgray
]

fig, ax, p = pie(counts, color = colors, axis=(;aspect=1))
hidedecorations!(ax)
hidespines!(ax)

Legend(fig[2,1], [MarkerElement(color=c, marker=:circle) for c in colors],
                labels, "Race",
        orientation=:horizontal,
        tellheight=true, tellwidth=false, nbanks=3
)

# save("figures/demo_race.png", fig)
fig
##


ed = skipmissing(unique(tps, [:subject]).ed) |> collect

fig, ax, p = pie(map(levels(ed)) do lev
        count(==(lev), ed)
    end;
    color = Makie.resample_cmap(:viridis, 6), axis=(;aspect=1)
)
hidedecorations!(ax)
hidespines!(ax)

Legend(fig[2,1],
    [MarkerElement(color=c, marker=:circle) for c in Makie.resample_cmap(:viridis, 6)],
    levels(ed),
    "Maternal Education";
    orientation=:horizontal,
    tellheight=true, tellwidth=false, nbanks=3
)

save("figures/demo_education.png", fig)
fig


##


rce = skipmissing(unique(tps, [:subject]).rce) |> collect

fig, ax, p = pie(map(levels(rce)) do lev
        count(==(lev), skipmissing(rce))
    end;
    color = [:navajowhite3,:coral,:darkorange3,:darkseagreen3,:cyan4,:dimgray],
    axis  =(; aspect=1)
)

hidedecorations!(ax)
hidespines!(ax)

Legend(fig[2,1],
    [MarkerElement(color=c, marker=:circle) for c in  [:navajowhite3,:coral,:darkorange3,:darkseagreen3,:cyan4,:dimgray]],
    levels(rce),
    "Race";
    orientation=:horizontal,
    tellheight=true, tellwidth=false, nbanks=3
)

save("figures/demo_race.png", fig)
fig


##

unique!(omni, :sample)

dyads = @chain omni begin
    groupby(:subject)
    combine(:Mother_Child => (mc-> length(unique(mc)) == 2) => :hasdyad)
    @rsubset!(:hasdyad)
end

toplot = @chain omni begin
    @rsubset(:subject in dyads.subject)
    groupby(:subject)
    combine(:Mother_Child=> (mc-> count(==("C"), mc))=> :nchild, :Mother_Child=> (mc -> count(==("M"), mc)) => :nmother)
end

fig, ax1, h1 = hist(toplot.nchild, axis=(; xlabel="number of child samples", ylabel="number of dyads"))
ax2 = Axis(fig[1,2], xlabel="number of mother samples")
h2 = hist!(ax2, toplot.nmother)

linkyaxes!(ax1, ax2)

save("figures/dyads.png", fig)
fig