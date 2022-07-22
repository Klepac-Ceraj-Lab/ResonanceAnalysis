using Resonance
omni, etoh, tps, complete_brain, metabolites, species = startup()

#-


unique!(tps, ["subject", "timepoint"])

set!(species, leftjoin(select(omni, ["subject", "timepoint", "sample"]), tps, on=[:subject, :timepoint]))
set!(metabolites, leftjoin(select(etoh, ["subject", "timepoint", "sample"]), tps, on=[:subject, :timepoint]))

#-

data_overlaps = DataFrame(subject = union(omni.subject, tps.subject))
data_overlaps.stool = let
    subjects = get(species, :subject)
    map(data_overlaps.subject) do s
        count(==(s), subjects)
    end
end

data_overlaps.stool_kids = let
    subjects = get(species, :subject)
    kids = get(species, :Mother_Child) .== "C"
    map(data_overlaps.subject) do s
        count(==(s), subjects[kids])
    end
end

data_overlaps.stool_moms = let
    subjects = get(species, :subject)
    moms = get(species, :Mother_Child) .== "M"
    map(data_overlaps.subject) do s
        count(==(s), subjects[moms])
    end
end

data_overlaps.brains = let
    gdf = groupby(tps, "subject")
    map(data_overlaps.subject) do s
        count(!ismissing, gdf[(;subject=s)]."Left_Thalamus")
    end
end

data_overlaps.ages_months = let
    gdf = groupby(tps, "subject")
    map(data_overlaps.subject) do s
        if all(ismissing, gdf[(;subject=s)]."ageMonths")
            ()
        else
            Tuple(skipmissing(gdf[(;subject=s)]."ageMonths"))
        end
    end
end

data_overlaps.ages_stool = let
    gdf = groupby(tps, "subject")
    stps = collect(zip(get(species, :subject), get(species, :timepoint)))
    map(data_overlaps.subject) do subject
        g = gdf[(; subject)]
        count(row-> (row.subject, row.timepoint) in stps, eachrow(g))
    end
end

CSV.write(datafiles("sample_summaries.csv"), data_overlaps)

#-

using CairoMakie
using AlgebraOfGraphics

hist(collect(skipmissing(tps.ageMonths)))

# ## Old numbers
#
# These are numbers sent from a presentation given by Rhode Island groupby.
# Deprecating in favor of pulling from metadata table.

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

fig, ax, p = pie(map(p-> p[2], counts), color = colors, axis=(;aspect=1))
hidedecorations!(ax)
hidespines!(ax)

Legend(fig[2,1], [MarkerElement(color=c, marker=:circle) for c in colors],
                labels, "Race",
        orientation=:horizontal,
        tellheight=true, tellwidth=false, nbanks=3
)

## save("figures/demo_race.png", fig)
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