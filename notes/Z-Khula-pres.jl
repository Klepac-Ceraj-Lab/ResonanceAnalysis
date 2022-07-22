using Resonance
omni, etoh, tps, metabolites, species = startup()
genes = Resonance.read_gfs_arrow()

specu1 = species[:, [!ismissing(a) && 2 < a < 10 for a in get(species, :ageMonths)]]
shannon!(specu1)

set!(genes, leftjoin(tps, select(omni, [:sample, :subject, :timepoint]), on = [:subject, :timepoint]))
genesu1 = genes[:, [!ismissing(a) && 2 < a < 10 for a in get(genes, :ageMonths)]]
shannon!(genesu1)

set!(metabolites, leftjoin(tps, select(etoh, [:sample, :subject, :timepoint]), on = [:subject, :timepoint]))
metabolitesu1 = metabolites[:, [!ismissing(a) && 2 < a < 10 for a in get(metabolites, :ageMonths)]]

pcspec = pcoa(specu1)

pcgenes = pcoa(genesu1)
#-

using CairoMakie
using AlgebraOfGraphics

#-

fig = Figure(resolution=(1000, 1000))
ax1 = Axis(fig[1,1]; ylabel=mdsaxis(pcspec, 2))
ax2 = Axis(fig[2,1]; xlabel=mdsaxis(pcspec, 1), ylabel=mdsaxis(pcspec, 2))
ax3 = Axis(fig[1,2]; ylabel=mdsaxis(pcgenes, 2))
ax4 = Axis(fig[2,2]; xlabel=mdsaxis(pcgenes, 1), ylabel=mdsaxis(pcgenes, 2))

sc1 = scatter!(ax1, loadings(pcspec, 1), loadings(pcspec,2);
                    xlabel=mdsaxis(pcspec, 1),
                    ylabel=mdsaxis(pcspec, 2), 
                    color=get(specu1, :ageMonths)
)
sc2 = scatter!(ax2, loadings(pcspec, 1), loadings(pcspec,2);
                    color=get(specu1, :shannon),
                    colormap=:plasma
)

sc3 = scatter!(ax3, loadings(pcgenes, 1), loadings(pcgenes,2);
                    xlabel=mdsaxis(pcgenes, 1),
                    ylabel=mdsaxis(pcgenes, 2), 
                    color=get(genesu1, :ageMonths)
)
sc4 = scatter!(ax4, loadings(pcgenes, 1), loadings(pcgenes,2);
                    color=get(genesu1, :shannon),
                    colormap=:plasma
)
Colorbar(fig[1,3], sc1; label="Age (months)")
Colorbar(fig[2,3], sc2; label="Shannon diversity")

Label(fig[0, 1], "Taxonomic Ordinations"; textsize=20, tellwidth=false)
Label(fig[1, 2], "Functional Ordinations"; textsize=20, tellwidth=false)

fig

#-

using PERMANOVA

preds = DataFrame(sample = samplenames(specu1),
                  ageMonths = get(specu1, :ageMonths),
                  cogScore=get(specu1, :cogScore),
                  race=get(specu1, :simple_race)
)

hascog = findall(!ismissing, preds.cogScore)
hasrace = findall(!ismissing, preds.race)

p_spec_age = permanova(select(preds, [:ageMonths]), abundances(specu1)', BrayCurtis, @formula(1~ageMonths))
p_spec_cog = permanova(select(preds, [:cogScore])[hascog, :], abundances(specu1)[:, hascog]', BrayCurtis, @formula(1~cogScore))
p_spec_race = permanova(select(preds, [:race])[hasrace, :], abundances(specu1)[:, hasrace]', BrayCurtis, @formula(1~race))

#- 

preds = DataFrame(sample = samplenames(genesu1),
                  ageMonths = get(genesu1, :ageMonths),
                  cogScore=get(genesu1, :cogScore),
                  race=get(genesu1, :simple_race)
)

hascog = findall(!ismissing, preds.cogScore)
hasrace = findall(!ismissing, preds.race)

p_genes_age = permanova(select(preds, [:ageMonths]), abundances(genesu1)', BrayCurtis, @formula(1~ageMonths))
p_genes_cog = permanova(select(preds, [:cogScore])[hascog, :], abundances(genesu1)[:, hascog]', BrayCurtis, @formula(1~cogScore))
p_genes_race = permanova(select(preds, [:race])[hasrace, :], abundances(genesu1)[:, hasrace]', BrayCurtis, @formula(1~race))

#- 

preds = DataFrame(sample = samplenames(metabolitesu1),
                  ageMonths = get(metabolitesu1, :ageMonths),
                  cogScore=get(metabolitesu1, :cogScore),
                  race=get(metabolitesu1, :simple_race)
)

hascog = findall(!ismissing, preds.cogScore)
hasrace = findall(!ismissing, preds.race)

p_metabolites_age = permanova(select(preds, [:ageMonths]), abundances(metabolitesu1)', BrayCurtis, @formula(1~ageMonths))
p_metabolites_cog = permanova(select(preds, [:cogScore])[hascog, :], abundances(metabolitesu1)[:, hascog]', BrayCurtis, @formula(1~cogScore))
p_metabolites_race = permanova(select(preds, [:race])[hasrace, :], abundances(metabolitesu1)[:, hasrace]', BrayCurtis, @formula(1~race))

varexpl(p::PERMANOVA.PSummary) = p.results[1, 3] * 100
pvalue(p::PERMANOVA.PSummary) = p.results[1, 5]

vmat = varexpl.([
    p_spec_age  p_genes_age  p_metabolites_age;
    p_spec_cog  p_genes_cog  p_metabolites_cog;
    p_spec_race p_genes_race p_metabolites_race;
])


pmat = pvalue.([
    p_spec_age  p_genes_age  p_metabolites_age;
    p_spec_cog  p_genes_cog  p_metabolites_cog;
    p_spec_race p_genes_race p_metabolites_race;
])

heatmap(vmat; colormap=:blues, colorrange=(0, 10))

#-

using GLM
using MixedModels
using AlgebraOfGraphics
using MultipleTesting

prev = prevalence(specu1, 0.01)

specdf = DataFrame(age=get(specu1, :ageMonths), cogScore=get(specu1, :cogScore))
for sp in features(specu1)[findall(p-> p > 0.1, vec(prev))]
    specdf[:, name(sp)] = vec(abundances(specu1[sp, :]))
end

specdf


lmres = DataFrame()

for sp in Symbol.(name.(features(specu1)[findall(p-> p > 0.1, vec(prev))]))
    f = @eval @formula($sp ~ cogScore+age)
    modl = lm(f, specdf)
    df = DataFrame(coeftable(modl))
    df.species .= string(sp)
    rename!(df, "Pr(>|t|)"=>"pvalue", "Coef."=> "coef")
    append!(lmres, df)
end

subset!(lmres, "Name"=> ByRow(!=("(Intercept)")))
lmres.qvalue = adjust(lmres.pvalue, BenjaminiHochberg())

sort!(lmres, :qvalue)

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