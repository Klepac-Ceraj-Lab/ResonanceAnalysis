include("startup.jl")

using CairoMakie
using GLM
using MixedModels
using AlgebraOfGraphics
using Statistics
using MultipleTesting
using CategoricalArrays


lq, uq = quantile(skipmissing(get(species, :cogScore)), [0.25, 0.75])
for smp in samplenames(species)
    cs = get(species, smp, :cogScore)
    qrt = ismissing(cs) ? missing : 
          cs <= lq      ? "lower" :
          cs >= uq      ? "upper" : 
                          "middle"
                  
    set!(species, smp, :cogQuartile, qrt)
end

specdf = metadata(species) |> DataFrame

##

prevspecies = prevalence(species) |> vec .> 0.1
totest = map(x-> !ismissing(x) && x >=18, get(species, :ageMonths)) .&
         map(x-> !ismissing(x), get(species, :cogScore))

@info count(totest)

testcomm = species[prevspecies, totest]
lm_results = DataFrame()

##

for spc in features(testcomm)
    @info name(spc)
    ab = vec(abundances(testcomm[spc, :]))
    over0 = ab .> 0

    df = DataFrame(
        age = get(testcomm, :ageMonths)[over0],
        cogScore = get(testcomm, :cogScore)[over0],
        cogQuartile = categorical(get(testcomm, :cogQuartile)[over0], levels=["lower", "middle", "upper"]),
        bug = log.(collect(ab[over0]))
    )

    mod = lm(@formula(bug ~ age + cogScore), df)
    ct = DataFrame(coeftable(mod))
    ct.species .= name(spc)
    ct.testvar .= "cogScore"
    append!(lm_results, ct)
    
    mod = lm(@formula(bug ~ age + levelcode(cogQuartile)), df)
    ct = DataFrame(coeftable(mod))
    ct.species .= name(spc)
    ct.testvar .= "cogQuartile"
    append!(lm_results, ct)

end    

select!(lm_results, Cols(:species, :testvar, :))
rename!(lm_results, "Pr(>|t|)"=>"pvalue")

##

adjust = MultipleTesting.adjust

@chain lm_results begin
    groupby(:testvar)
    transform!(:pvalue => (col-> adjust(collect(col), BenjaminiHochberg())) => :qvalue)
end


sort!(lm_results, :qvalue)

CSV.write("data/lm_results/species_cogscore.csv", lm_results)
##

cl = ["lower"=> :seagreen, "middle"=> :lightgray, "upper"=> :dodgerblue]

dt = data(specdf[completecases(specdf, [:ageMonths, :cogScore]), :]) * mapping(:ageMonths => "Age (months)", :cogScore=> "Cognitive function score", color = :cogQuartile)

out = draw(dt * visual(Scatter, strokewidth=0.5), palettes=(;color=cl))

save("figures/cogScore_quartiles.png", out)

##


bugs = [
    "s__Bifidobacterium_catenulatum",
    "s__Coprococcus_eutactus", # butyrate producer
    "s__Streptococcus_thermophilus",
    "s__Eubacterium_eligens",
    "s__Ruminococcus_gnavus",
    "s__Veillonella_parvula"
]

bugdf = DataFrame()
for bug in bugs
    ab = vec(abundances(testcomm[bug, :]))
    over0 = ab .> 0
    df = DataFrame(
        age = get(testcomm, :ageMonths)[over0],
        cogScore = get(testcomm, :cogScore)[over0],
        cogQuartile = categorical(get(testcomm, :cogQuartile)[over0], levels=["lower", "middle", "upper"]),
        bug = log.(collect(ab[over0])),
        regbug = collect(ab[over0])
    )
    
    df.bugname .= replace(bug, r"s__([A-Z])(\w+)_(\w+)"=> s"\1. \3")
    append!(bugdf, df)
end

##

fig = Figure(resolution=(800,1200))
draw!(fig[1,1], 
      data(bugdf) *
      mapping(:cogQuartile=>"Quartile",
              :bug=> "Log(relative abundance)";
              color=:cogQuartile, row=:bugname
              ) * 
      (visual(BoxPlot) + visual(Scatter, strokewidth=0.5)),
      palettes= (; color=cl),
      axis = (; titlevisible=false)
)


draw!(fig[1,2], 
      data(bugdf) *
      mapping(:cogScore=>"Cognitive Function Score",
              :bug=> "Log(relative abundance)";
              color=:cogQuartile, row=:bugname
              );
      palettes= (; color=cl),
      axis=(; ylabel="", yticklabelsvisible=false)
)
                   
save("figures/lms.png", fig)
fig

##

##

fig = Figure(resolution=(900,600))
# draw!(fig[1,1], 
#       data(subset(bugdf, :bugname => ByRow(b-> b in replace.(bugs[1:3], r"s__([A-Z])(\w+)_(\w+)"=> s"\1. \3")))) *
#       mapping(:cogQuartile=>"Quartile",
#               :bug=> "Log(relative abundance)";
#               color=:cogQuartile, row=:bugname
#               ) * 
#       (visual(BoxPlot) + visual(Scatter, strokewidth=0.5)),
#       palettes= (; color=cl),
#       axis = (; titlevisible=false, xticklabelsize=14)
# )


draw!(fig[1,1], 
      data(subset(bugdf, :bugname => ByRow(b-> b in replace.(bugs[1:3], r"s__([A-Z])(\w+)_(\w+)"=> s"\1. \3")))) *
      (mapping(:cogScore=>"Cognitive Function Score",
              :bug=> "Log(relative abundance)";
              color=:cogQuartile, row=:bugname
              ) +
       mapping(:cogScore=>"Cognitive Function Score",
              :bug=> "Log(relative abundance)";
              row=:bugname
              ) * linear());
      palettes= (; color=cl),
      axis=(; xticklabelsize=14)
)
                   
# draw!(fig[1,3], 
#       data(subset(bugdf, :bugname => ByRow(b-> b in replace.(bugs[4:6], r"s__([A-Z])(\w+)_(\w+)"=> s"\1. \3")))) *
#       mapping(:cogQuartile=>"Quartile",
#               :bug=> "Log(relative abundance)";
#               color=:cogQuartile, row=:bugname
#               ) * 
#       (visual(BoxPlot) + visual(Scatter, strokewidth=0.5)),
#       palettes= (; color=cl),
#       axis = (; titlevisible=false, ylabelvisible = false, yticklabelsvisible=false, xticklabelsize=14)
# )


draw!(fig[1,2], 
      data(subset(bugdf, :bugname => ByRow(b-> b in replace.(bugs[4:6], r"s__([A-Z])(\w+)_(\w+)"=> s"\1. \3")))) *
      (mapping(:cogScore=>"Cognitive Function Score",
              :bug=> "Log(relative abundance)";
              color=:cogQuartile, row=:bugname
              ) +
       mapping(:cogScore=>"Cognitive Function Score",
              :bug=> "Log(relative abundance)";
              row=:bugname
              ) * linear());
      palettes= (; color=cl),
)
                   
save("figures/lms-split-2.png", fig)
fig

##

################
# Longitudinal #
################

sort!(tps, [:subject, :timepoint])
grps = groupby(tps, :subject)
fecals = Set(zip(get(species, :subject), get(species, :timepoint)))

##

fig = Figure()
ax = Axis(fig[1,1], xlabel="timepoint", ylabel="cogScore")

for grp in grps
    any(!ismissing, grp.cogScore) || continue
    keep = subset(grp, :cogScore => ByRow(!ismissing))
    scatter!(ax, keep.timepoint, keep.cogScore)
    lines!(ax, keep.timepoint, keep.cogScore)
end 

fig

##

diffs = Float64[]

for key in keys(grps)
    grp = grps[key]
    scores = collect(skipmissing(grp.cogScore))
    n = length(scores)
    n > 1 || continue
    d = collect(Iterators.map(x-> abs(x[1] - x[2]), Iterators.product(scores, scores)))
    for i in 1:n, j in 1:n
        i < j || continue
        push!(diffs, d[i] - d[j])
    end
end

##

hist(abs.(diffs), axis=(; xlabel="abs(score difference)", ylabel="count"))

##

tpsfecal = leftjoin(tps, DataFrame(subject=get(species, :subject), timepoint=get(species, :timepoint), sample=samplenames(species)), on=[:subject, :timepoint])
sort!(tpsfecal, [:subject, :timepoint])
grps = groupby(tpsfecal, :subject)

##

youngfecal = Dict()
oldfecal = Dict()

for grp in grps
    subject = first(grp.subject)
    issorted(grp.timepoint) || error()

    firstyoung = findfirst(((a, f),) -> !ismissing(f) && !ismissing(a) && a < 12, collect(zip(grp.ageMonths, grp.sample)))
    firstold   = findfirst(((a, f),) -> !ismissing(f) && !ismissing(a) && a > 18, collect(zip(grp.ageMonths, grp.sample)))

    if !isnothing(firstyoung)
        laterscores = map(((i, s),) -> i > firstyoung && !ismissing(s), enumerate(grp.cogScore))
        if sum(laterscores) > 0
            lscr = grp[laterscores, :cogScore]
            ages = grp[laterscores, :ageMonths]
            youngfecal[subject] = (;fecal = grp[firstyoung, :sample], cogscores = collect(skipmissing(lscr[.!ismissing.(lscr)])), ages = collect(skipmissing(ages[.!ismissing.(lscr)])))
        end
    end

    if !isnothing(firstold)
        laterscores = map(((i, s),) -> i > firstold && !ismissing(s), enumerate(grp.cogScore))
        if sum(laterscores) > 0
            lscr = grp[laterscores, :cogScore]
            ages = grp[laterscores, :ageMonths]
            oldfecal[subject] = (;fecal = grp[firstold, :sample], cogscores = collect(skipmissing(lscr[.!ismissing.(lscr)])), ages = collect(skipmissing(ages[.!ismissing.(lscr)])))
        end
    end
end

##

yfdf = DataFrame(sample = [youngfecal[k][:fecal] for k in keys(youngfecal)])
yfdf.age = [get(species, s, :ageMonths) for s in yfdf.sample]
yfdf.hhs = [get(species, s, :mother_HHS_Education) for s in yfdf.sample]
yfdf.firstscore = [first(youngfecal[k][:cogscores]) for k in keys(youngfecal)]
yfdf.meanscore = [mean(youngfecal[k][:cogscores]) for k in keys(youngfecal)]
yfdf.allscores = [youngfecal[k][:cogscores] for k in keys(youngfecal)]
yfdf.allages = [youngfecal[k][:ages] for k in keys(youngfecal)]


##

youngsamples = map(a-> !ismissing(a) && a <= 12, get(species, :ageMonths))
youngprev = vec((prevalence(species[:, youngsamples]) .> 0.1))

youngcomm = species[youngprev, youngsamples]

for sp in features(youngcomm)
    yfdf[!, Microbiome.name(sp)] = collect(vec(abundances(species[sp, yfdf.sample])))
end

CSV.write("data/lm_results/under12m_cogScores_taxa.csv", yfdf)

##

lm_results = DataFrame()

for spc in featurenames(youngcomm)
    @info spc
    df = yfdf[:, ["sample", "age", "firstscore", "meanscore", spc]]
    subset!(df, spc => ByRow(!=(0)))
    size(df, 1) > 5 || continue

    df.bug = log.(df[!, spc])

    mod = lm(@formula(bug ~ age + firstscore), df)
    ct = DataFrame(coeftable(mod))
    ct.species .= spc
    ct.testvar .= "firstscore"
    append!(lm_results, ct)
    
    mod = lm(@formula(bug ~ age + meanscore), df)
    ct = DataFrame(coeftable(mod))
    ct.species .= spc
    ct.testvar .= "meanscore"
    append!(lm_results, ct)

end    

select!(lm_results, Cols(:species, :testvar, :))
rename!(lm_results, "Pr(>|t|)"=>"pvalue")

##

@chain lm_results begin
    groupby(:testvar)
    transform!(:pvalue => (col-> adjust(collect(col), BenjaminiHochberg())) => :qvalue)
    @rsubset!(:Name in ("firstscore", "meanscore"))
end
lm_results

sort!(lm_results, :qvalue)
CSV.write("data/lm_results/young_species_predict_cogscore.csv", lm_results)

##

##

ofdf = DataFrame(sample = [oldfecal[k][:fecal] for k in keys(oldfecal)])
ofdf.age = [get(species, s, :ageMonths) for s in ofdf.sample]
ofdf.hhs = [get(species, s, :mother_HHS_Education) for s in ofdf.sample]
ofdf.firstscore = [first(oldfecal[k][:cogscores]) for k in keys(oldfecal)]
ofdf.meanscore = [mean(oldfecal[k][:cogscores]) for k in keys(oldfecal)]
ofdf.allscores = [oldfecal[k][:cogscores] for k in keys(oldfecal)]
ofdf.allages = [oldfecal[k][:ages] for k in keys(oldfecal)]


##

oldsamples = map(a-> !ismissing(a) && a > 12, get(species, :ageMonths))
oldprev = vec((prevalence(species[:, oldsamples]) .> 0.1))

oldcomm = species[oldprev, oldsamples]

for sp in features(oldcomm)
    ofdf[!, Microbiome.name(sp)] = collect(vec(abundances(species[sp, ofdf.sample])))
end

CSV.write("data/lm_results/over12m_cogScores_taxa.csv", ofdf)

##

lm_results = DataFrame()

for spc in featurenames(oldcomm)
    @info spc
    df = ofdf[:, ["sample", "age", "firstscore", "meanscore", spc]]
    subset!(df, spc => ByRow(!=(0)))
    size(df, 1) > 5 || continue

    df.bug = log.(df[!, spc])

    mod = lm(@formula(bug ~ age + firstscore), df)
    ct = DataFrame(coeftable(mod))
    ct.species .= spc
    ct.testvar .= "firstscore"
    append!(lm_results, ct)
    
    mod = lm(@formula(bug ~ age + meanscore), df)
    ct = DataFrame(coeftable(mod))
    ct.species .= spc
    ct.testvar .= "meanscore"
    append!(lm_results, ct)

end    

select!(lm_results, Cols(:species, :testvar, :))
rename!(lm_results, "Pr(>|t|)"=>"pvalue")

##

@chain lm_results begin
    groupby(:testvar)
    transform!(:pvalue => (col-> adjust(collect(col), BenjaminiHochberg())) => :qvalue)
    @rsubset!(:Name in ("firstscore", "meanscore"))
end
lm_results

sort!(lm_results, :qvalue)
CSV.write("data/lm_results/old_species_predict_cogscore.csv", lm_results)

##

using MLJ

yfdf.quartile = categorical([x < lq ? "lower" : x > uq ? "upper" : "middle" for x in yfdf.firstscore], levels=["lower", "middle", "upper"])

y, X = unpack(yfdf, ==(:quartile), (x-> !in(x, (:sample, :meanscore, :firstscore, :age))))

Tree = @load DecisionTreeClassifier pkg=DecisionTree
tree = Tree()

fprs = Float64[]
tprs = Float64[]

for _ in 1:100
    res = MLJ.evaluate(tree, X, y,
        resampling=Holdout(; fraction_train=0.9, shuffle=true),
                measures=[Accuracy(), MulticlassTruePositiveRate(), MulticlassFalsePositiveRate()],
                verbosity=2
    )

    (ac, tpr, fpr) = res.measurement
    push!(fprs, fpr)
    push!(tprs, tpr)
end

srt = sortperm(fprs)
scatter(fprs, tprs, axis=(; limits=(0,1,0,1)))


##

tree3 = Tree(max_depth=3)
tree5 = Tree(max_depth=5)
tree8 = Tree(max_depth=8)
tree10 = Tree(max_depth=10)

mach3 = machine(tree3, X, y)
mach5 = machine(tree5, X, y)
mach8 = machine(tree8, X, y)
mach10 = machine(tree10, X, y)
fit!.([mach3,mach5,mach8,mach10])

yhat3 = MLJ.predict(mach3)
yhat5 = MLJ.predict(mach5)
yhat8 = MLJ.predict(mach8)
yhat10 = MLJ.predict(mach10)


fprs3, tprs3, ts3 = roc_curve(yhat3, y)
fprs5, tprs5, ts5 = roc_curve(yhat5, y)
fprs8, tprs8, ts8 = roc_curve(yhat8, y)
fprs10, tprs10, ts10 = roc_curve(yhat10, y)

##
fig = Figure()
ax = Axis(fig[1,1], ylabel="true positive rate (sensitivity)", xlabel="false positive rate (specificity)")
noinf = lines!(ax, [0, 1], [0, 1], linestyle=:dash, color=:black)



sc3 = scatter!(ax, fprs3, tprs3)
ln3 = lines!(ax, fprs3, tprs3; color=sc3.attributes[:color][])
sc5 = scatter!(ax, fprs5, tprs5)
ln5 = lines!(ax, fprs5, tprs5; color=sc5.attributes[:color][])
sc8 = scatter!(ax, fprs8, tprs8)
ln8 = lines!(ax, fprs8, tprs8; color=sc8.attributes[:color][])
sc10 = scatter!(ax, fprs10, tprs10)
ln10 = lines!(ax, fprs10, tprs10; color=sc10.attributes[:color][])

Legend(fig[1,2], [[sc3, ln3], [sc5, ln5], [sc8, ln8], [sc10, ln10], noinf],
                 ["Depth = 3", "Depth = 5", "Depth = 8", "Depth = 10", "no information"])

fig

