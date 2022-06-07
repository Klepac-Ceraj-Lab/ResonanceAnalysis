using Resonance
omni, etoh, tps, complete_brain, metabolites, species = startup()

## 

using CairoMakie
using AlgebraOfGraphics
using Microbiome.Distances
using Microbiome.MultivariateStats

@chain DataFrame(subject=get(species, :subject), timepoint = get(species, :timepoint), age = get(species, :ageMonths)) begin
    groupby(:subject)
    combine(:age => (a-> all(ismissing, a) ? 0 : maximum(skipmissing(a))) => :maxage)
    subset(:maxage => ByRow(a-> a > 24))
end

met_pcoa = pcoa(metabolites)
spec_pcoa = pcoa(species)

kidsidx = get(species, :Mother_Child) .== "C"
kidsspecies = species[:, kidsidx]
kidsspecies = kidsspecies[vec(prevalence(kidsspecies) .> 0), :]
ginisimpson!(kidsspecies)
shannon!(kidsspecies)

kids_pcoa = pcoa(kidsspecies)

using Microbiome.Distances
using MultivariateStats

brain_dm = pairwise(Euclidean(), Matrix(tps[complete_brain, Resonance.brainmeta]), dims=1)
brain_pcoa = fit(MDS, brain_dm, distances=true)

## 
function getcolors(vals, clip=(0.1, 30); highclip=Makie.to_color(:white), lowclip=Makie.to_color(:black), colormap=:plasma)
    kidscolor = Makie.to_colormap(colormap)

    map(vals) do val
        if ismissing(val)
            return Makie.to_color(:gray)
        elseif val < clip[1]
            return lowclip
        elseif val > clip[2]
            return highclip
        else
            Makie.interpolated_getindex(kidscolor, val, clip)
        end
    end
end

##

fig = Figure(resolution=(800, 800))
ax1 = Axis(fig[1,1], xlabel="MDS1 ($(round(varexplained(met_pcoa)[1] * 100, digits=2))%)",
                     ylabel="MDS2 ($(round(varexplained(met_pcoa)[2] * 100, digits=2))%)")

scatter!(ax1, Resonance.loadings(met_pcoa)[:,1], Resonance.loadings(met_pcoa)[:,2],
        color=[ismissing(x) ? :grey : x == "M" ? :dodgerblue : :orange for x in get(metabolites, :Mother_Child)])

leg = Legend(fig[1,2], [MarkerElement(color=:orange, marker=:circle),
                        MarkerElement(color=:dodgerblue, marker=:circle)],
                        ["Kids", "Moms"])

ax2 = Axis(fig[2,1],
            xlabel="MDS1 ($(round(varexplained(met_pcoa)[1] * 100, digits=2))%)",
            ylabel="MDS2 ($(round(varexplained(met_pcoa)[2] * 100, digits=2))%)"
)

scatter!(ax2, Resonance.loadings(met_pcoa)[:,1], Resonance.loadings(met_pcoa)[:,2],
         color = getcolors(get(metabolites, :ageMonths), (0,25)))
cleg = Colorbar(fig[2, 2], limits=(0.1, 25), colormap=:plasma, highclip=:white, lowclip=:black, label="Age in Months")
save("figures/metabolites_pcoa.png", fig)

fig

##

fig = Figure(resolution=(800, 800))
ax1 = Axis(fig[1,1], xlabel="MDS1 ($(round(varexplained(spec_pcoa)[1] * 100, digits=2))%)",
                     ylabel="MDS2 ($(round(varexplained(spec_pcoa)[2] * 100, digits=2))%)")
scatter!(ax1, Resonance.loadings(spec_pcoa)[:,1], Resonance.loadings(spec_pcoa)[:,2],
        color=[ismissing(x) ? :grey : x == "M" ? :dodgerblue : :orange for x in get(species, :Mother_Child)],
        strokewidth=0.5
)

leg = Legend(fig[1,2], [MarkerElement(color=:orange, marker=:circle, strokewidth=0.5),
                        MarkerElement(color=:dodgerblue, marker=:circle, strokewidth=0.5)],
                        ["Kids", "Moms"])

ax2 = Axis(fig[2,1],
            xlabel="MDS1 ($(round(varexplained(spec_pcoa)[1] * 100, digits=2))%)",
            ylabel="MDS2 ($(round(varexplained(spec_pcoa)[2] * 100, digits=2))%)"
)

scatter!(ax2, Resonance.loadings(spec_pcoa)[:,1], Resonance.loadings(spec_pcoa)[:,2],
         color = getcolors(get(species, :ageMonths) ./ 12, (0,10)),
         strokewidth=0.5
)
cleg = Colorbar(fig[2, 2], limits=(0.1, 10), colormap=:plasma, highclip=:white, lowclip=:black, label="Age in Years")
save("figures/species_pcoa.png", fig)

fig

##

fig = Figure(resolution=(600, 600))
ax1 = Axis(fig[1,1], xlabel="Age (years)",
                     ylabel="MDS1 ($(round(varexplained(kids_pcoa)[1] * 100, digits=2))%)"
)

sc = scatter!(ax1, get(kidsspecies, :ageMonths) ./12, Resonance.loadings(kids_pcoa)[:,1],
        color = get(kidsspecies, :shannon),
        strokewidth=0.5
)

lb = Label(fig[1,2], "α diveristy (Shannon)", rotation=π/2, tellheight=false, padding=(-10,-10,0,0))
cleg = Colorbar(fig[1, 3], sc)

save("figures/species_pcoa_kids_age.png", fig)

fig
##

fig = Figure(resolution=(800, 800))
ax1 = Axis(fig[1,1],
            xlabel="PCA1 ($(round(varexplained(brain_pcoa)[1] * 100, digits=2))%)",
            ylabel="PCA2 ($(round(varexplained(brain_pcoa)[2] * 100, digits=2))%)"
)

scatter!(ax1, Resonance.loadings(brain_pcoa)[:,1], Resonance.loadings(brain_pcoa)[:,2],
         color = getcolors(tps[complete_brain, :ageMonths] ./ 12, (0,10)),
         strokewidth=0.5
)
cleg = Colorbar(fig[1, 2], limits=(0.1, 10), colormap=:plasma, highclip=:white, lowclip=:black, label="Age in Years")

ax2 = Axis(fig[2,1],
        xlabel="PCA1 ($(round(varexplained(brain_pcoa)[1] * 100, digits=2))%)",
        ylabel="Age in Years"
)

scatter!(ax2, Resonance.loadings(brain_pcoa)[:,1], tps[complete_brain, :ageMonths] ./ 12,
        color = getcolors(tps[complete_brain, :cogScore], (65,135), colormap=:viridis),
        strokewidth=0.5
)

cleg = Colorbar(fig[2, 2], limits=(65, 135), colormap=:viridis, highclip=:white, lowclip=:black, label="Cognitive function score")

save("figures/brain_pcoa.png", fig)

fig

##

srt_metab, srt_spec = stp_overlap(collect(zip(get(metabolites, :subject), get(metabolites, :timepoint))),
                                  collect(zip(get(species, :subject), get(species, :timepoint)))
)

##

fig = Figure(resolution=(800, 800))
ax1 = Axis(fig[1,1], 
            xlabel="Species MDS1 ($(round(varexplained(spec_pcoa)[1] * 100, digits=2))%)",
            ylabel="Metabolites MDS1 ($(round(varexplained(met_pcoa)[1] * 100, digits=2))%)"
)

scatter!(ax1, Resonance.loadings(spec_pcoa)[srt_spec,1], Resonance.loadings(met_pcoa)[srt_metab,1],
        color=[ismissing(x) ? :grey : x == "M" ? :dodgerblue : :orange for x in get(metabolites, :Mother_Child)[srt_metab]],
        strokewidth=0.5
)

leg = Legend(fig[1,2], [MarkerElement(color=:orange, marker=:circle),
                        MarkerElement(color=:dodgerblue, marker=:circle)],
                        ["Kids", "Moms"])

ax2 = Axis(fig[2,1],
            xlabel="Species MDS1 ($(round(varexplained(spec_pcoa)[1] * 100, digits=2))%)",
            ylabel="Metabolites MDS1 ($(round(varexplained(met_pcoa)[1] * 100, digits=2))%)"
)

scatter!(ax2, Resonance.loadings(spec_pcoa)[srt_spec,1], Resonance.loadings(met_pcoa)[srt_metab,1],
         color = getcolors(get(metabolites, :ageMonths)[srt_metab], (0,25)),
         strokewidth=0.5
)
cleg = Colorbar(fig[2, 2], limits=(0.1, 25), colormap=:plasma, highclip=:white, lowclip=:black, label="Age in Months")

fig

##

srt_brain, srt_spec = stp_overlap(collect(zip(tps[complete_brain, :subject], tps[complete_brain, :timepoint])),
                                  collect(zip(get(species, :subject), get(species, :timepoint)))
)

##

fig = Figure(resolution=(800, 800))
ax1 = Axis(fig[1,1], 
            xlabel="Species MDS1 ($(round(varexplained(spec_pcoa)[1] * 100, digits=2))%)",
            ylabel="Brain PCA1 ($(round(varexplained(brain_pcoa)[1] * 100, digits=2))%)"
)

scatter!(ax1, Resonance.loadings(spec_pcoa)[srt_spec,1], Resonance.loadings(brain_pcoa)[srt_brain,1],
        color=[ismissing(x) ? :grey : x == "M" ? :dodgerblue : :orange for x in get(species, :Mother_Child)[srt_spec]],
        strokewidth=0.5
)

leg = Legend(fig[1,2], [MarkerElement(color=:orange, marker=:circle),
                        MarkerElement(color=:dodgerblue, marker=:circle)],
                        ["Kids", "Moms"])

ax2 = Axis(fig[2,1],
            xlabel="Species MDS1 ($(round(varexplained(spec_pcoa)[1] * 100, digits=2))%)",
            ylabel="Brain MDS1 ($(round(varexplained(brain_pcoa)[1] * 100, digits=2))%)"
)

scatter!(ax2, Resonance.loadings(spec_pcoa)[srt_spec,1], Resonance.loadings(brain_pcoa)[srt_brain,1],
         color = getcolors(get(species, :ageMonths)[srt_spec], (0,120)),
         strokewidth=0.5
)
cleg = Colorbar(fig[2, 2], limits=(0.1, 120), colormap=:plasma, highclip=:white, lowclip=:black, label="Age in Months")

fig

##

gaba = findfirst(f-> commonname(f) === "gamma-Aminobutyric acid", features(metabolites))
glutamate = findfirst(f-> commonname(f) === "Glutamic acid", features(metabolites))

momsidx = string.(get(metabolites, :Mother_Child)) .=== "M"
kidsidx = string.(get(metabolites, :Mother_Child)) .=== "C"

##

fig = Figure(resolution=(800,800))
ax1 = Axis(fig[1:4,1:4], xlabel="GABA (log abundance)", ylabel="Glutamate (log abundance)")
ax2 = Axis(fig[0,1:4], height=200, xlabels=false, xticklabelsvisible=false)
ax3 = Axis(fig[2:5, 5], width=200, ylabels=false, yticklabelsvisible=false)


scmom = scatter!(ax1, log.(vec(abundances(metabolites[gaba, momsidx]))), log.(vec(abundances(metabolites[glutamate,momsidx]))))
hist!(ax2, log.(vec(metabolites[gaba, momsidx] |> abundances)))
hist!(ax3, log.(vec(metabolites[glutamate, momsidx] |> abundances)), direction=:x)

sckid = scatter!(ax1, log.(vec(metabolites[gaba,kidsidx] |> abundances)), log.(vec(metabolites[glutamate, kidsidx] |> abundances)))
hist!(ax2, log.(vec(metabolites[gaba, kidsidx] |> abundances)))
hist!(ax3, log.(vec(metabolites[glutamate, kidsidx] |> abundances)), direction=:x)


leg = Legend(fig[1,5], [scmom, sckid], ["Moms", "Kids"], tellwidth = false, tellheight = false)

save("figures/gaba-glutamate.png", fig)
fig


##

fig = Figure(resolution=(1200,1200))
ax1 = Axis(fig[1:2,1:2], xlabel="GABA (log)", ylabel="Glutamate (log)")
ax2 = Axis(fig[0,1:2], height=200)
ax3 = Axis(fig[2:3, 3], width=200)


scmom = scatter!(ax1, log.(vec(abundances(metabolites[gaba,momsidx]))), 
                      log.(vec(abundances(metabolites[glutamate, momsidx]))), color=:gray)
histmom = hist!(ax2, log.(vec(abundances(metabolites[gaba,momsidx]))), color=:gray)
hist!(ax3, log.(vec(abundances(metabolites[glutamate, momsidx]))), direction=:x, color=:gray)


sckid = scatter!(ax1, log.(vec(abundances(metabolites[gaba,kidsidx]))), 
                      log.(vec(abundances(metabolites[glutamate, kidsidx]))),
                      color = getcolors(get(metabolites[:, kidsidx], :ageMonths), (0,25)),
                      strokewidth=1)
histkid = hist!(ax2, log.(vec(abundances(metabolites[gaba,kidsidx]))))
hist!(ax3, log.(vec(abundances(metabolites[glutamate,kidsidx]))), direction=:x)



cleg = Colorbar(fig[2:3, 4], limits=(0.1, 25), colormap=:plasma, highclip=:white, lowclip=:black, label="Age in Months")
leg = Legend(fig[1,3], [histmom, histkid], ["Moms", "Kids"], tellwidth = false, tellheight = false)

save("figures/gaba-glutamate_age.png", fig)
fig

##

allmeta = leftjoin(tps, select(omni, [:subject, :timepoint, :sample]), on=[:subject, :timepoint])
unirefs = Resonance.read_gfs_arrow()
set!(unirefs, allmeta)
##

neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(unirefs)))

overlap = intersect(Set(zip(get(unirefs, :subject), get(unirefs, :timepoint))), Set(zip(get(metabolites, :subject), get(metabolites, :timepoint))))

metaboverlap = metabolites[:, map(s-> (s.subject, s.timepoint) in overlap, samples(metabolites))]

##

using GLM
using MixedModels
using AlgebraOfGraphics



genemetab = DataFrame(
    gabasynth = log.(1 .+ map(sum, eachcol(abundances(unirefs[neuroactive["GABA synthesis"], [esmap[s] for s in samplenames(metaboverlap)]])))),
    gabadegr  = log.(1 .+ map(sum, eachcol(abundances(unirefs[neuroactive["GABA degradation"], [esmap[s] for s in samplenames(metaboverlap)]])))),
    gabagut   = log.(vec(abundances(metaboverlap[gaba, :]))),
    glutsynth = log.(1 .+ map(sum, eachcol(abundances(unirefs[neuroactive["Glutamate synthesis"], [esmap[s] for s in samplenames(metaboverlap)]])))),
    glutdegr  = log.(1 .+ map(sum, eachcol(abundances(unirefs[neuroactive["Glutamate degradation"], [esmap[s] for s in samplenames(metaboverlap)]])))),
    glutgut   = log.(vec(abundances(metaboverlap[glutamate, :]))),
    mc        = get(metaboverlap, :Mother_Child))

gabasynthlm = lm(@formula(gabagut ~ gabasynth + mc), genemetab)
gabadegrlm = lm(@formula(gabagut ~ gabadegr + mc), genemetab)
glutsynthlm = lm(@formula(glutgut ~ glutsynth + mc), genemetab)
glutdegrlm = lm(@formula(glutgut ~ glutdegr + mc), genemetab)

##

pred = DataFrame(gabasynth = repeat(range(extrema(genemetab.gabasynth)..., length=50), outer=2),
                 gabadegr  = repeat(range(extrema(genemetab.gabadegr)..., length=50), outer=2),
                 gabagut   = zeros(100),
                 glutsynth = repeat(range(extrema(genemetab.glutsynth)..., length=50), outer=2),
                 glutdegr  = repeat(range(extrema(genemetab.glutdegr)..., length=50), outer=2),
                 glutgut   = zeros(100),
                 mc        = repeat(["M", "C"], inner=50))

pred.gabasynth_pred = predict(gabasynthlm, pred)
pred.gabadegr_pred  = predict(gabadegrlm, pred)
pred.glutsynth_pred = predict(glutsynthlm, pred)
pred.glutdegr_pred  = predict(glutdegrlm, pred)

predgrp = groupby(pred, :mc)

##

fig = Figure(resolution=(800,800))
ax1 = Axis(fig[1,1], xlabel="GABA Synthesis log(RPKM)", ylabel="GABA (log abundance)")
ax2 = Axis(fig[2,1], xlabel="GABA Degradation log(RPKM)", ylabel="GABA (log abundance)")

scatter!(ax1, log.(1+genemetab.gabasynth), log.(1+genemetab.gabagut),
              color=[ismissing(x) ? :gray : x == "M" ? :dodgerblue : :orange for x in genemetab.mc])
              
for grp in groupby(pred, :mc)
    c = first(grp.mc) == "M" ? :dodgerblue : :orange
    lines!(ax1, grp.gabasynth, grp.gabasynth_pred, color=c)
end
            

scatter!(ax2, log.(1+genemetab.gabadegr), log.(1+genemetab.gabagut),
              color=[ismissing(x) ? :gray : x == "M" ? :dodgerblue : :orange for x in genemetab.mc])


for grp in groupby(pred, :mc)
    c = first(grp.mc) == "M" ? :dodgerblue : :orange
    lines!(ax2, grp.gabadegr, grp.gabadegr_pred, color=c)
end  

##

m1 = MarkerElement(color=:dodgerblue, marker=:circle)
l1 = LineElement(color=:dodgerblue)
m2 = MarkerElement(color=:orange, marker=:circle)
l2 = LineElement(color=:orange)
Legend(fig[1:2, 2], [[m1, l1], [m2, l2]], ["moms", "kids"])

save("figures/gaba_genes_metabolites.png", fig)
fig

##

##

fig = Figure(resolution=(850, 1100))
ax1 = Axis(fig[1,1], xlabel="Glutamate Synthesis log(RPKM)", ylabel="Glutamate (log abundance)")
ax2 = Axis(fig[2,1], xlabel="Glutamate Degradation log(RPKM)", ylabel="Glutamate (log abundance)")


scatter!(ax1, log.(1+genemetab.glutsynth), log.(1+genemetab.glutgut),
              color=[ismissing(x) ? :gray : x == "M" ? :dodgerblue : :orange for x in genemetab.mc])
              
for grp in groupby(pred, :mc)
    c = first(grp.mc) == "M" ? :dodgerblue : :orange
    lines!(ax1, grp.glutsynth, grp.glutsynth_pred, color=c)
end
            

scatter!(ax2, log.(1+genemetab.glutdegr), log.(1+genemetab.glutgut),
              color=[ismissing(x) ? :gray : x == "M" ? :dodgerblue : :orange for x in genemetab.mc])


for grp in groupby(pred, :mc)
    c = first(grp.mc) == "M" ? :dodgerblue : :orange
    lines!(ax2, grp.glutdegr, grp.glutdegr_pred, color=c)
end  

##

m1 = MarkerElement(color=:dodgerblue, marker=:circle)
l1 = LineElement(color=:dodgerblue)
m2 = MarkerElement(color=:orange, marker=:circle)
l2 = LineElement(color=:orange)
Legend(fig[4, 1], [[m1, l1], [m2, l2]], ["moms", "kids"], orientation=:horizontal, tellheight=true, tellwidth=false)

save("figures/glutamate_genes_metabolites.png", fig)
fig

##

scatter(log.(1 .+ genemetab.glutdegr), log.(1 .+ genemetab.glutsynth), genemetab.glutgut,
    axis=(; xlabel="degradataion", ylabel="synthesis", zlabel="concentration"))

## 

using PERMANOVA

