using Resonance
using ProgressLogging
using CategoricalArrays
using CairoMakie
using AlgebraOfGraphics
using Microbiome.Distances
using Microbiome.MultivariateStats
using GLM
using MixedModels
using MultipleTesting
using AlgebraOfGraphics

omni, etoh, tps, metabolites, species = startup()
genes = Resonance.read_gfs_arrow()
ecs = Resonance.read_gfs_arrow(kind="ecs_rename")
ecs = filter(!hastaxon, ecs)
kos = Resonance.read_gfs_arrow(kind="kos_rename")
kos = filter(!hastaxon, kos)

unique!(tps, ["subject", "timepoint"])

set!(species, leftjoin(select(omni, ["subject", "timepoint", "sample"]), tps, on=[:subject, :timepoint]))
set!(genes, leftjoin(select(omni, ["subject", "timepoint", "sample"]), tps, on=[:subject, :timepoint]))
set!(ecs, leftjoin(select(omni, ["subject", "timepoint", "sample"]), tps, on=[:subject, :timepoint]))
set!(kos, leftjoin(select(omni, ["subject", "timepoint", "sample"]), tps, on=[:subject, :timepoint]))
set!(metabolites, leftjoin(select(etoh, ["subject", "timepoint", "sample"]), tps, on=[:subject, :timepoint]))

## 

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

complete_brain = map(!ismissing, tps."Left-Thalamus")

brain_dm = pairwise(Euclidean(), Matrix(tps[complete_brain, replace.(Resonance.brainmeta, "-"=>"_")]), dims=1)
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


(genesoverlap, metaboverlap) = comm_overlap(genes, metabolites)
neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(genesoverlap)))

gabaol = findfirst(f-> commonname(f) === "gamma-Aminobutyric acid", features(metaboverlap))
glutamateol = findfirst(f-> commonname(f) === "Glutamic acid", features(metaboverlap))


##

genemetab = DataFrame(
    gabasynth = log.(1 .+ map(sum, eachcol(abundances(genesoverlap[neuroactive["GABA synthesis"], :])))),
    gabadegr  = log.(1 .+ map(sum, eachcol(abundances(genesoverlap[neuroactive["GABA degradation"], :])))),
    gabagut   = log.(vec(abundances(metaboverlap[gabaol, :]))),
    glutsynth = log.(1 .+ map(sum, eachcol(abundances(genesoverlap[neuroactive["Glutamate synthesis"], :])))),
    glutdegr  = log.(1 .+ map(sum, eachcol(abundances(genesoverlap[neuroactive["Glutamate degradation"], :])))),
    glutgut   = log.(vec(abundances(metaboverlap[glutamateol, :]))),
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

scatter!(ax1, log.(1 .+ genemetab.gabasynth), log.(1 .+ genemetab.gabagut),
              color=[ismissing(x) ? :gray : x == "M" ? :dodgerblue : :orange for x in genemetab.mc])
              
for grp in groupby(pred, :mc)
    c = first(grp.mc) == "M" ? :dodgerblue : :orange
    lines!(ax1, log.(1 .+ grp.gabasynth), log.(1 .+ grp.gabasynth_pred), color=c)
end
            

scatter!(ax2, log.(1 .+ genemetab.gabadegr), log.(1 .+ genemetab.gabagut),
              color=[ismissing(x) ? :gray : x == "M" ? :dodgerblue : :orange for x in genemetab.mc])


for grp in groupby(pred, :mc)
    c = first(grp.mc) == "M" ? :dodgerblue : :orange
    lines!(ax2, log.(1 .+ grp.gabadegr), log.(1 .+ grp.gabadegr_pred), color=c)
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


scatter!(ax1, log.(1 .+ genemetab.glutsynth), log.(1 .+ genemetab.glutgut),
              color=[ismissing(x) ? :gray : x == "M" ? :dodgerblue : :orange for x in genemetab.mc])
              
for grp in groupby(pred, :mc)
    c = first(grp.mc) == "M" ? :dodgerblue : :orange
    lines!(ax1, log.(1 .+ grp.glutsynth), log.(1 .+ grp.glutsynth_pred), color=c)
end
            

scatter!(ax2, log.(1 .+ genemetab.glutdegr), log.(1 .+ genemetab.glutgut),
              color=[ismissing(x) ? :gray : x == "M" ? :dodgerblue : :orange for x in genemetab.mc])


for grp in groupby(pred, :mc)
    c = first(grp.mc) == "M" ? :dodgerblue : :orange
    lines!(ax2, log.(1 .+ grp.glutdegr), log.(1 .+ grp.glutdegr_pred), color=c)
end  

##

m1 = MarkerElement(color=:dodgerblue, marker=:circle)
l1 = LineElement(color=:dodgerblue)
m2 = MarkerElement(color=:orange, marker=:circle)
l2 = LineElement(color=:orange)
Legend(fig[3, 1], [[m1, l1], [m2, l2]], ["moms", "kids"], orientation=:horizontal, tellheight=true, tellwidth=false)

save("figures/glutamate_genes_metabolites.png", fig)
fig

##

# ## Linear models
#
# For each feature type (`species`, `ECs`, `KOs`),
# we will run the following linear mixed-effects models:
# 
# - All samples: feature ~ (1|subject) + ageMonths + maternal_ed + race + sex + read_depth
# - just kids < 6mo: feature ~ (1|subject) + cogScore + ageMonths + cesarean + breastfeeding + maternal_ed + race + sex + read_depth
# - just kids > 1 year : feature ~ (1|subject) + cogScore + ageMonths + maternal_ed + race + sex + read_depth
# - For kids that have a stool sample < 6 months and a future cogscore: feature_u6mo ~ cogScore_later + age_stool + maternal_ed + race + sex + read_depth
# 

# ### Read count
#
# All linear models will include read-depth as a covariate

reads = CSV.read(datafiles("read_counts.csv"), DataFrame)

for col in names(reads, Not("Sample"))
    if eltype(reads[!, col]) <: Union{Missing, <:AbstractString}
        reads[!, col] = [ismissing(x) ? missing : parse(Float64, x) for x in reads[!, col]]
    end
end

select!(reads,
    "Sample"=> ByRow(s-> replace(s, r"_S\d+_kneaddata"=> "")) => "sample",
    AsTable(r"final") => ByRow(r-> sum(values(r))) => "reads"
)


set!(species, reads)
set!(kos, reads)
set!(ecs, reads)

# ### Species lms

speclm_in = DataFrame(
    sample      = samplenames(species),
    subject     = get(species, :subject),
    read_depth  = get(species, :reads),
    ageMonths   = get(species, :ageMonths),
    race        = get(species, :rce),
    maternal_ed = levelcode.(get(species, :ed)),
    cogScore    = get(species, :cogScore)
)

speclm_in = subset(speclm_in, :ageMonths   => ByRow(!ismissing),
                              :read_depth  => ByRow(!ismissing),
                              :maternal_ed => ByRow(!ismissing),
                              :cogScore    => ByRow(!ismissing)
)

disallowmissing!(speclm_in, [:ageMonths, :read_depth, :maternal_ed, :cogScore])

let comm = species[:, [s in speclm_in.sample for s in samplenames(species)]]
    for sp in features(comm)[vec(prevalence(comm)) .> 0.05]
        speclm_in[!, name(sp)] = vec(abundances(comm[sp, :]))
    end
end

speclms = DataFrame()

@progress for sp in names(speclm_in, Not(["sample", "subject", "read_depth", "ageMonths", "race", "maternal_ed"]))
    ssp = Symbol(sp)
    f = @eval @formula($ssp ~ (1|subject) + read_depth + ageMonths + race + maternal_ed)
    
    modl = lm(f, speclm_in)

    df = DataFrame(coeftable(modl))
    df.species .= sp
    df.model .= "Complete"
    rename!(df, "Pr(>|t|)"=>"pvalue", "Coef."=> "coef")
    append!(speclms, df)
end

let indf = subset(speclm_in, :ageMonths => ByRow(<(6)))
    @progress for sp in names(speclm_in, Not(["sample", "subject", "read_depth", "ageMonths", "race", "maternal_ed", "cogScore"]))
        count(>(0), indf[:, sp]) / size(indf, 1) > 0.05 || continue

        ssp = Symbol(sp)
        
        f = @eval @formula($ssp ~ (1|subject) + cogScore + read_depth + ageMonths + race + maternal_ed)
        
        modl = lm(f, speclm_in)

        df = DataFrame(coeftable(modl))
        df.species .= sp
        df.model .= "Under 6mo"
        rename!(df, "Pr(>|t|)"=>"pvalue", "Coef."=> "coef")
        append!(speclms, df)
    end
end

let indf = subset(speclm_in, :ageMonths => ByRow(>(12)))
    @progress for sp in names(speclm_in, Not(["sample", "subject", "read_depth", "ageMonths", "race", "maternal_ed", "cogScore"]))
        count(>(0), indf[:, sp]) / size(indf, 1) > 0.05 || continue

        ssp = Symbol(sp)
        
        f = @eval @formula($ssp ~ (1|subject) + cogScore + read_depth + ageMonths + race + maternal_ed)
        
        modl = lm(f, speclm_in)

        df = DataFrame(coeftable(modl))
        df.species .= sp
        df.model .= "Over 12mo"
        rename!(df, "Pr(>|t|)"=>"pvalue", "Coef."=> "coef")
        append!(speclms, df)
    end
end

subset!(speclms, :Name=> ByRow(n-> any(m-> contains(n, m), ("cogScore", "maternal_ed", "race"))), :pvalue=> ByRow(!isnan))
speclms.qvalue = adjust(speclms.pvalue, BenjaminiHochberg())
sort!(speclms, :qvalue)

isdir(datafiles("lms")) || mkpath(datafiles("lms"))
CSV.write(datafiles("lms", "species_lms.csv"), speclms)

scatter(speclm_in.maternal_ed, speclm_in.cogScore)

# ### KOs lms

kolm_in = DataFrame(
    sample      = samplenames(kos),
    subject     = categorical(get(kos, :subject)),
    read_depth  = get(kos, :reads),
    ageMonths   = get(kos, :ageMonths),
    race        = get(kos, :rce),
    maternal_ed = levelcode.(get(kos, :ed)),
    cogScore    = get(kos, :cogScore)
)

kolm_in = subset(kolm_in, :ageMonths   => ByRow(!ismissing),
                          :read_depth  => ByRow(!ismissing),
                          :maternal_ed => ByRow(!ismissing),
                          :cogScore => ByRow(!ismissing)
)

disallowmissing!(kolm_in, [:ageMonths, :read_depth, :maternal_ed, :cogScore])

let comm = kos[:, [s in kolm_in.sample for s in samplenames(kos)]]
    for sp in features(comm)[vec(prevalence(comm)) .> 0.1]
        kolm_in[!, name(sp)] = vec(abundances(comm[sp, :]))
    end
end

#-

kos_cols = names(kolm_in, Not(["sample", "subject", "read_depth", "ageMonths", "race", "maternal_ed"]))
kos_sp_dfs = [DataFrame() for _ in kos_cols]

Threads.@threads for i in eachindex(kos_cols)
    sp = kos_cols[i]

    ssp = Symbol(sp)
    terms = sum(term(t) for t in (1, :read_depth, :ageMonths, :race, :maternal_ed))
    resp = term(ssp)
    sub = (term(1) | term(:subject))

    f = resp ~ sub + terms
    
    modl = fit(MixedModel, f, kolm_in)

    df = DataFrame(coeftable(modl))
    df.kos .= sp
    df.model .= "Complete"
    rename!(df, "Pr(>|z|)"=>"pvalue", "Coef."=> "coef")
    subset!(df, :z=> ByRow(!isnan))
    kos_sp_dfs[i] = df
end

complete = reduce(vcat, kos_sp_dfs)

#-

kos_cols = names(kolm_in, Not(["sample", "subject", "read_depth", "ageMonths", "race", "maternal_ed", "cogScore"]))
kos_sp_dfs = [DataFrame() for _ in kos_cols]

indf = subset(kolm_in, :ageMonths => ByRow(<(6)))
indf = unique(indf, :subject)


Threads.@threads for i in eachindex(kos_cols)
    sp = kos_cols[i]
    count(x-> x > 0, indf[!, sp]) / size(indf, 1) > 0.05 || continue

    ssp = Symbol(sp)
    terms = sum(term(t) for t in (1, :cogScore, :read_depth, :ageMonths,  :maternal_ed))
    resp = term(ssp)

    f = resp ~ terms
    
    try
        modl = fit(MixedModel, f, indf)
    catch e
        continue
    end

    df = DataFrame(coeftable(modl))
    df.kos .= sp
    df.model .= "Complete"
    rename!(df, "Pr(>|z|)"=>"pvalue", "Coef."=> "coef")
    subset!(df, :z=> ByRow(!isnan))
    kos_sp_dfs[i] = df
end

u6mo = reduce(vcat, kos_sp_dfs)

#-

kos_cols = names(kolm_in, Not(["sample", "subject", "read_depth", "ageMonths", "race", "maternal_ed", "cogScore"]))
kos_sp_dfs = [DataFrame() for _ in kos_cols]

indf = subset(kolm_in, :ageMonths => ByRow(>(12)))
indf = unique(indf, :subject)

Threads.@threads for i in eachindex(kos_cols)
    sp = kos_cols[i]
    count(x-> x > 0, indf[!, sp]) / size(indf, 1) > 0.05 || continue

    ssp = Symbol(sp)
    terms = sum(term(t) for t in (1, :cogScore, :read_depth, :ageMonths,  :maternal_ed))
    resp = term(ssp)

    f = resp ~ terms
    
    try
        modl = fit(LinearModel, f, indf)
    catch e
        continue
    end

    df = DataFrame(coeftable(modl))
    df.kos .= sp
    df.model .= "Complete"
    rename!(df, "Pr(>|t|)"=>"pvalue", "Coef."=> "coef")
    subset!(df, :t=> ByRow(!isnan))
    kos_sp_dfs[i] = df
end

o12mo = reduce(vcat, kos_sp_dfs)

#-

kolms = vcat(complete, u6mo, o12mo)

subset!(kolms, :Name=> ByRow(n-> any(m-> contains(n, m), ("cogScore", "maternal_ed", "race"))), :pvalue=> ByRow(!isnan))
kolms.qvalue = adjust(kolms.pvalue, BenjaminiHochberg())
sort!(kolms, :qvalue)

isdir(datafiles("lms")) || mkpath(datafiles("lms"))
CSV.write(datafiles("lms", "kos_lms.csv"), kolms)

scatter(kolm_in.maternal_ed, kolm_in.cogScore)