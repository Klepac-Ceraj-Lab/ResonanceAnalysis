using Resonance
using AlgebraOfGraphics

omni = CSV.read("data/wrangled/omnisamples.csv", DataFrame)
etoh = CSV.read("data/wrangled/etohsamples.csv", DataFrame)

##

tps  = CSV.read("data/wrangled/timepoints.csv", DataFrame)
tps."Left-Thalamus" = map(eachrow(tps)) do row
    (t, p) = (row."Left-Thalamus", row."Left-Thalamus-Proper")
    all(ismissing, (t,p)) && return missing
    return max(coalesce(t, 0), coalesce(p, 0))
end

tps."Right-Thalamus" = map(eachrow(tps)) do row
    (t, p) = (row."Right-Thalamus", row."Right-Thalamus-Proper")
    all(ismissing, (t,p)) && return missing
    return max(coalesce(t, 0), coalesce(p, 0))
end

mainmeta = [
    "ageMonths",
    "age0to3mo",
    "age3to6mo",
    "age6to12mo",
    "age12moplus",
    "mother_HHS_Education",
    "simple_race",
    "cogScore",
    "has_segmentation"
]

brainmeta = ["CortexVol",
            #  "CorticalWhiteMatterVol",
             "SubCortGrayVol",
             "TotalGrayVol",
             "BrainSegVol-to-eTIV",
            #  "CerebralWhiteMatterVol",
             "EstimatedTotalIntraCranialVol",
            #  "lhCorticalWhiteMatterVol",
            #  "lhCerebralWhiteMatterVol",
             "lhCortexVol",
             "Left-Thalamus",
             "Left-Lateral-Ventricle",
             "Left-Cerebellum-White-Matter",
             "Left-Cerebellum-Cortex",
             "Left-Caudate",
             "Left-Putamen",
             "Left-Pallidum",
             "Left-Hippocampus",
             "Left-Amygdala",
             "Left-Accumbens-area",
             "Left-VentralDC",
             "Left-choroid-plexus",
            #  "rhCorticalWhiteMatterVol",
            #  "rhCerebralWhiteMatterVol",
             "rhCortexVol",
             "Right-Thalamus",
             "Right-Lateral-Ventricle",
             "Right-Cerebellum-White-Matter",
             "Right-Cerebellum-Cortex",
             "Right-Caudate",
             "Right-Putamen",
             "Right-Pallidum",
             "Right-Hippocampus",
             "Right-Amygdala",
             "Right-Accumbens-area",
             "Right-VentralDC",
             "Right-choroid-plexus",
             "Brain-Stem",
             "CSF"
]


select!(tps, ["subject", "timepoint", mainmeta..., brainmeta...])
rename!(tps, Dict(k=> replace(k, "-"=>"_") for k in brainmeta))
foreach(i-> (brainmeta[i] = replace(brainmeta[i], "-"=>"_")), eachindex(brainmeta))

for m in brainmeta
    c = count(x-> !ismissing(x) && x != 0, tps[!, m])
    @info "Not missing or 0 `$m`: $c"
end

##

metabolites = CSV.read("data/wrangled/metabolites.csv", DataFrame)
ms = [Resonance.Metabolite(row[:uid], row[:Metabolite], row[:MZ], row[:RT]) for row in eachrow(metabolites)]
metabolites = CommunityProfile(Matrix(metabolites[!, 9:end]), ms, MicrobiomeSample.(names(metabolites)[9:end]))
set!(metabolites, leftjoin(etoh, tps, on=[:subject, :timepoint], makeunique=true))

species = CSV.read("data/wrangled/species.csv", DataFrame)
species = CommunityProfile(Matrix(species[!, 2:end]), Taxon.(species[!, 1]), MicrobiomeSample.(names(species)[2:end]))
set!(species, leftjoin(omni, tps, on=[:subject, :timepoint], makeunique=true))


## 
using Microbiome.Distances
using Microbiome.MultivariateStats

met_pcoa = pcoa(metabolites)
spec_pcoa = pcoa(species)

brain_dm = pairwise(Euclidean(), Matrix(tps[completecases(tps[:, brainmeta]), brainmeta]), dims=1)
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

fig = Figure(resolution=(800, 800))
ax1 = Axis(fig[1,1],
            xlabel="PCA1 ($(round(varexplained(brain_pcoa)[1] * 100, digits=2))%)",
            ylabel="PCA2 ($(round(varexplained(brain_pcoa)[2] * 100, digits=2))%)"
)

scatter!(ax1, Resonance.loadings(brain_pcoa)[:,1], Resonance.loadings(brain_pcoa)[:,2],
         color = getcolors(tps[completecases(tps[:, brainmeta]), :ageMonths] ./ 12, (0,10)),
         strokewidth=0.5
)
cleg = Colorbar(fig[1, 2], limits=(0.1, 10), colormap=:plasma, highclip=:white, lowclip=:black, label="Age in Years")

ax2 = Axis(fig[2,1],
        xlabel="PCA1 ($(round(varexplained(brain_pcoa)[1] * 100, digits=2))%)",
        ylabel="Age in Years"
)

scatter!(ax2, Resonance.loadings(brain_pcoa)[:,1], tps[completecases(tps[:, brainmeta]), :ageMonths] ./ 12,
        color = getcolors(tps[completecases(tps[:, brainmeta]), :cogScore], (65,135), colormap=:viridis),
        strokewidth=0.5
)

cleg = Colorbar(fig[2, 2], limits=(65, 135), colormap=:viridis, highclip=:white, lowclip=:black, label="Cognitive function score")

save("figures/brain_pcoa.png", fig)

fig

##

overlap = zip(get(metabolites, :subject), get(metabolites, :timepoint)) ∩
          zip(get(species, :subject), get(species, :timepoint))

sort!(overlap, lt=(x,y)-> x[1] == y[1] ? x[2] < y[2] : x < y)

keepetoh = [x ∈ overlap for x in zip(get(metabolites, :subject), get(metabolites, :timepoint))]
keepomni = [!any(ismissing, x) && x ∈ overlap for x in zip(get(species, :subject), get(species, :timepoint))]
sum(keepomni)
##

gaba = findfirst(f-> commonname(f) === "gamma-Aminobutyric acid", features(metabolites))
glutamate = findfirst(f-> commonname(f) === "Glutamic acid", features(metabolites))

momsidx = string.(get(metabolites, :Mother_Child)) .=== "M"
kidsidx = string.(get(metabolites, :Mother_Child)) .=== "C"

##

fig = Figure(resolution=(1200,1200))
ax1 = Axis(fig[1:2,1:2], xlabel="GABA (log)", ylabel="Glutamate (log)")
ax2 = Axis(fig[0,1:2], height=200)
ax3 = Axis(fig[2:3, 3], width=200)


scmom = scatter!(ax1, log.(vec(abundances(metabolites[gaba, momsidx]))), log.(vec(abundances(metabolites[glutamate,momsidx]))))
hist!(ax2, log.(vec(metabolites[gaba, momsidx] |> abundances)))
hist!(ax3, log.(vec(metabolites[glutamate, momsidx] |> abundances)), direction=:x)

sckid = scatter!(ax1, log.(vec(metabolites[gaba,kidsidx] |> abundances)), log.(vec(metabolites[glutamate, kidsidx] |> abundances)))
hist!(ax2, log.(vec(metabolites[gaba, kidsidx] |> abundances)))
hist!(ax3, log.(vec(metabolites[glutamate, kidsidx] |> abundances)), direction=:x)


leg = Legend(fig[1,3], [scmom, sckid], ["Moms", "Kids"], tellwidth = false, tellheight = false)

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


unirefs = Resonance.load_genefamilies()
set!(unirefs, allmeta)

neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(unirefs)))

metagrp = groupby(allmeta, [:subject, :timepoint])
esmap = Dict()
for grp in metagrp
    ss = unique(skipmissing(grp.sample))
    nrow(grp) > 1 || continue
    any(s-> contains(s, "FE"), ss) || continue
    for s in ss
        contains(s, "FE") || continue
        fgs = filter(s2-> contains(s2, "FG"), ss)
        isempty(fgs) && @info ss
        esmap[s] = isempty(fgs) ? missing : first(fgs)
    end
end

metagrp[(; subject=774, timepoint=2)]

for s in samples(metabolites)
    @show s.subject, s.timepoint
    break
end

metaboverlap = metabolites[:, findall(s-> haskey(esmap, s) && esmap[s] ∈ samplenames(unirefs), samplenames(metabolites))]

##

using GLM
using MixedModels
using AlgebraOfGraphics



genemetab = DataFrame(
    gabasynth = map(sum, eachcol(abundances(unirefs[neuroactive["GABA synthesis"], [esmap[s] for s in samplenames(metaboverlap)]]))),
    gabadegr  = map(sum, eachcol(abundances(unirefs[neuroactive["GABA degradation"], [esmap[s] for s in samplenames(metaboverlap)]]))),
    gabagut   = log.(vec(abundances(metaboverlap[gaba, :]))),
    glutsynth = map(sum, eachcol(abundances(unirefs[neuroactive["Glutamate synthesis"], [esmap[s] for s in samplenames(metaboverlap)]]))),
    glutdegr  = map(sum, eachcol(abundances(unirefs[neuroactive["Glutamate degradation"], [esmap[s] for s in samplenames(metaboverlap)]]))),
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
ax1 = Axis(fig[1,1], xlabel="GABA Synthesis (RPKM)", ylabel="GABA (log abundance)")
ax2 = Axis(fig[2,1], xlabel="GABA Degradation (RPKM)", ylabel="GABA (log abundance)")

scatter!(ax1, genemetab.gabasynth, genemetab.gabagut,
              color=[ismissing(x) ? :gray : x == "M" ? :dodgerblue : :orange for x in genemetab.mc])
              
for grp in groupby(pred, :mc)
    c = first(grp.mc) == "M" ? :dodgerblue : :orange
    lines!(ax1, grp.gabasynth, grp.gabasynth_pred, color=c)
end
            

scatter!(ax2, genemetab.gabadegr, genemetab.gabagut,
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
ax1 = Axis(fig[1,1], xlabel="Glutamate Synthesis (RPKM)", ylabel="Glutamate (log abundance)")
ax2 = Axis(fig[2,1], xlabel="Glutamate Degradation (RPKM)", ylabel="Glutamate (log abundance)")
ax3 = Axis(fig[3,1], xlabel="Glutamate Degradation (RPKM)", ylabel="Glutamate Glutamate Synthesis (RPKM)")


scatter!(ax1, genemetab.glutsynth, genemetab.glutgut,
              color=[ismissing(x) ? :gray : x == "M" ? :dodgerblue : :orange for x in genemetab.mc])
              
for grp in groupby(pred, :mc)
    c = first(grp.mc) == "M" ? :dodgerblue : :orange
    lines!(ax1, grp.glutsynth, grp.glutsynth_pred, color=c)
end
            

scatter!(ax2, genemetab.glutdegr, genemetab.glutgut,
              color=[ismissing(x) ? :gray : x == "M" ? :dodgerblue : :orange for x in genemetab.mc])


for grp in groupby(pred, :mc)
    c = first(grp.mc) == "M" ? :dodgerblue : :orange
    lines!(ax2, grp.glutdegr, grp.glutdegr_pred, color=c)
end  

scatter!(ax3, log.(1 .+ genemetab.glutdegr), log.(1 .+ genemetab.glutsynth),
              color=[ismissing(x) ? :gray : x == "M" ? :dodgerblue : :orange for x in genemetab.mc])
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

