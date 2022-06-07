using Resonance
using PERMANOVA

omni, etoh, tps, complete_brain, metabolites, species = startup()
genes = Resonance.read_gfs_arrow()



#-

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

using CairoMakie
using AlgebraOfGraphics