using Resonance
using CairoMakie


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

#-

specu6mo = species[:, map(a-> !ismissing(a) && a < 6, get(species, :ageMonths))]
specu6mo = specu6mo[vec(prevalence(specu6mo)) .> 0, :]

genesu6mo = genes[:, map(a-> !ismissing(a) && a < 6, get(genes, :ageMonths))]
genesu6mo = genesu6mo[vec(prevalence(genesu6mo)) .> 0, :]

ecsu6mo = ecs[:, map(a-> !ismissing(a) && a < 6, get(ecs, :ageMonths))]
ecsu6mo = ecsu6mo[vec(prevalence(ecsu6mo)) .> 0, :]

kosu6mo = kos[:, map(a-> !ismissing(a) && a < 6, get(kos, :ageMonths))]
kosu6mo = kosu6mo[vec(prevalence(kosu6mo)) .> 0, :]

metabu6mo = metabolites[:, map(a-> !ismissing(a) && a < 6, get(metabolites, :ageMonths))]
metabu6mo = metabu6mo[vec(prevalence(metabu6mo)) .> 0, :]


#-

fig, ax, hm = plot_permanovas([specu6mo, genesu6mo, ecsu6mo, kosu6mo, metabu6mo], [:ageMonths, :cogScore, :rce, :ed];
    commlabels = ["species", "genes", "ecs", "kos", "metabolites"],
    mdlabels = ["age", "cogScore", "race", "education"]
)

Label(fig[0,1], "Kids under 6 mo"; tellwidth=false)
save("figures/kids_under6mo_permanovas.png", fig)

fig

#-

speco1y = species[:, map(a-> !ismissing(a) && a > 12, get(species, :ageMonths))]
speco1y = speco1y[:, unique(DataFrame(metadata(speco1y)), [:subject, :timepoint]).sample]
speco1y = speco1y[vec(prevalence(speco1y)) .> 0, :]

geneso1y = genes[:, map(a-> !ismissing(a) && a > 12, get(genes, :ageMonths))]
geneso1y = geneso1y[:, unique(DataFrame(metadata(geneso1y)), [:subject, :timepoint]).sample]
geneso1y = geneso1y[vec(prevalence(geneso1y)) .> 0, :]

ecso1y = ecs[:, map(a-> !ismissing(a) && a > 12, get(ecs, :ageMonths))]
ecso1y = ecso1y[:, unique(DataFrame(metadata(ecso1y)), [:subject, :timepoint]).sample]
ecso1y = ecso1y[vec(prevalence(ecso1y)) .> 0, :]

koso1y = kos[:, map(a-> !ismissing(a) && a > 12, get(kos, :ageMonths))]
koso1y = koso1y[:, unique(DataFrame(metadata(koso1y)), [:subject, :timepoint]).sample]
koso1y = koso1y[vec(prevalence(koso1y)) .> 0, :]


#-

fig, ax, hm = plot_permanovas([speco1y, geneso1y, ecso1y, koso1y], [:ageMonths, :cogScore, :rce, :ed];
    commlabels = ["species", "genes", "ecs", "kos"],
    mdlabels = ["age", "cogScore", "race", "education"]
)

Label(fig[0,1], "Kids over 1 year"; tellwidth=false)
save("figures/kids_over1y_permanovas.png", fig)

fig

#-

specallkids = species[:, map(a-> !ismissing(a), get(species, :ageMonths))]
specallkids = specallkids[:, unique(DataFrame(metadata(specallkids)), ["subject", "timepoint"]).sample]
specallkids = specallkids[vec(prevalence(specallkids)) .> 0, :]

genesallkids = genes[:, map(a-> !ismissing(a), get(genes, :ageMonths))]
genesallkids = genesallkids[:, unique(DataFrame(metadata(genesallkids)), ["subject", "timepoint"]).sample]
genesallkids = genesallkids[vec(prevalence(genesallkids)) .> 0, :]

ecsallkids = ecs[:, map(a-> !ismissing(a), get(ecs, :ageMonths))]
ecsallkids = ecsallkids[:, unique(DataFrame(metadata(ecsallkids)), ["subject", "timepoint"]).sample]
ecsallkids = ecsallkids[vec(prevalence(ecsallkids)) .> 0, :]

kosallkids = kos[:, map(a-> !ismissing(a), get(kos, :ageMonths))]
kosallkids = kosallkids[:, unique(DataFrame(metadata(kosallkids)), ["subject", "timepoint"]).sample]
kosallkids = kosallkids[vec(prevalence(kosallkids)) .> 0, :]

metabsallkids = metabolites[:, map(a-> !ismissing(a), get(metabolites, :ageMonths))]
metabsallkids = metabsallkids[:, unique(DataFrame(metadata(metabsallkids)), ["subject", "timepoint"]).sample]
metabsallkids = metabsallkids[vec(prevalence(metabsallkids)) .> 0, :]

#-

fig, ax, hm = plot_permanovas([specallkids, genesallkids, ecsallkids, kosallkids], [:ageMonths, :cogScore, :rce, :ed];
    commlabels = ["species", "genes", "ecs", "kos"],
    mdlabels = ["age", "cogScore", "race", "education"]
)

Label(fig[0,1], "All kids"; tellwidth=false)
save("figures/kids_all_permanovas.png", fig)

fig

#-

## Mantel tests

using Combinatorics

manteldf = DataFrame()

comms = [specallkids, genesallkids, ecsallkids, kosallkids, metabsallkids]
labels = ["species", "genes", "ecs", "kos", "metabolites"]

for ((c1, l1), (c2, l2)) in combinations(collect(zip(comms, labels)), 2)
    
    c1, c2 = comm_overlap(c1, c2)
    c1 = c1[featurenames(c1) .!= "UNGROUPED", :]
    c2 = c2[featurenames(c2) .!= "UNGROUPED", :]

    dm1 = braycurtis(c1)
    dm2 = braycurtis(c2)

    m, p = mantel(dm1, dm2)
    push!(manteldf, (; stat = m, pvalue = p, thing1 = l1, thing2 = l2))
      
end

manteldf

#-

n = length(comms)
vmat = zeros(n, n)
pmat = ones(n, n)

manteldf.idx = [CartesianIndex(i, j) for (i,j) in combinations(1:n, 2)]
for row in eachrow(manteldf)
    vmat[row.idx] = row.stat
    pmat[row.idx] = row.pvalue
end

#-

fig = Figure()
ax = Axis(fig[1,1], title="Mantel tests")
hm = heatmap!(ax, vmat[1:n-1, 2:n]'; colormap=:viridis, colorrange = (0.01, 1), lowclip=:lightgray)

for ci in manteldf.idx
    c = vmat[ci] < 0.5 ? :lightgray : :black
    text!(string(round(vmat[ci], digits=4)); position=(ci[2]-1,ci[1]), align=(:center, :center), color=c)
    p = pmat[ci]
    stars = p < 0.001 ? "***" : p < 0.01 ? "**" : p < 0.05 ? "*" : ""
    text!(stars; position=(ci[2]-1,ci[1]), align=(:center, :bottom), color=c)
end

ax.xticks = (1:n-1, isempty(labels) ? ["comm$i" for i in 1:n-1] : labels[2:end])
ax.yticks = (1:n-1, isempty(labels) ? ["comm$i" for i in 1:n] : labels[1:end-1])

save("figures/mantel.png", fig)
fig