using Resonance
using ProgressLogging
using CategoricalArrays
using CairoMakie
using AlgebraOfGraphics
using Microbiome.Distances
using Microbiome.MultivariateStats
using HypothesisTests
using MultipleTesting
using AlgebraOfGraphics
using Statistics
using DataFrames.InvertedIndices
using ThreadsX

omni, etoh, tps, complete_brain, metabolites, species = startup()
genes = Resonance.read_gfs_arrow()
unique!(tps, ["subject", "timepoint"])

set!(species, leftjoin(select(omni, ["subject", "timepoint", "sample"]), tps, on=[:subject, :timepoint]))
set!(genes, leftjoin(select(omni, ["subject", "timepoint", "sample"]), tps, on=[:subject, :timepoint]))
set!(metabolites, leftjoin(select(etoh, ["subject", "timepoint", "sample"]), tps, on=[:subject, :timepoint]))

neuroactive = Resonance.getneuroactive(map(f-> replace(f, "UniRef90_"=>""), featurenames(genes)))

#-

gsgenes = genes[:, .!ismissing.(get(genes, :cogScore))]

cscor = cor(abundances(gsgenes), get(gsgenes, :cogScore), dims=2)

fsdf = DataFrame(
    ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        isempty(ixs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cscor[ixs])
        isempty(cs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cscor[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)

        return (; geneset = gs, U = mwu.U, median = mwu.median, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end
)

subset!(fsdf, :pvalue=> ByRow(!isnan))
fsdf.qvalue = adjust(fsdf.pvalue, BenjaminiHochberg())
sort!(fsdf, :qvalue)
CSV.write("data/fsea_all.csv")

#- 

gs = "Menaquinone synthesis"
ixs = neuroactive[gs]
cs = filter(!isnan, cscor[ixs])
acs = filter(!isnan, cscor[Not(ixs)])

Resonance.plot_fsea(acs, cs)

#-

nodupes_samples = unique(DataFrame(metadata(gsgenes)), :subject).sample

gsnodupes = gsgenes[:, map(s-> name(s) ∈ nodupes_samples, samples(gsgenes))]

cscor = cor(abundances(gsnodupes), get(gsnodupes, :cogScore), dims=2)

fsdf = DataFrame(
    ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        isempty(ixs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cscor[ixs])
        isempty(cs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cscor[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)

        return (; geneset = gs, U = mwu.U, median = mwu.median, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end
)

subset!(fsdf, :pvalue=> ByRow(!isnan))
fsdf.qvalue = adjust(fsdf.pvalue, BenjaminiHochberg())
sort!(fsdf, :qvalue)
CSV.write("data/fsea_nodupe.csv")

#-

u6_samples = unique(subset(DataFrame(metadata(gsgenes)), :ageMonths => ByRow(<(6))), :subject).sample

gsu6 = gsgenes[:, map(s-> name(s) ∈ u6_samples, samples(gsgenes))]

cscor = cor(abundances(gsu6), get(gsu6, :cogScore), dims=2)

fsdf = DataFrame(
    ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        isempty(ixs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cscor[ixs])
        isempty(cs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cscor[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)

        return (; geneset = gs, U = mwu.U, median = mwu.median, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end
)

subset!(fsdf, :pvalue=> ByRow(!isnan))
fsdf.qvalue = adjust(fsdf.pvalue, BenjaminiHochberg())
sort!(fsdf, :qvalue)
CSV.write("data/fsea_u6.csv")


#-

o12_samples = unique(subset(DataFrame(metadata(gsgenes)), :ageMonths => ByRow(>(12))), :subject).sample

gso12 = gsgenes[:, map(s-> name(s) ∈ o12_samples, samples(gsgenes))]

cscor = cor(abundances(gso12), get(gso12, :cogScore), dims=2)

fsdf = DataFrame(
    ThreadsX.map(collect(keys(neuroactive))) do gs
        ixs = neuroactive[gs]
        isempty(ixs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        cs = filter(!isnan, cscor[ixs])
        isempty(cs) && return (; geneset = gs, U = NaN, median = NaN, mu = NaN, sigma = NaN, pvalue = NaN)

        acs = filter(!isnan, cscor[Not(ixs)])
        mwu = MannWhitneyUTest(cs, acs)

        return (; geneset = gs, U = mwu.U, median = mwu.median, mu = mwu.mu, sigma = mwu.sigma, pvalue=pvalue(mwu))
    end
)

subset!(fsdf, :pvalue=> ByRow(!isnan))
fsdf.qvalue = adjust(fsdf.pvalue, BenjaminiHochberg())
sort!(fsdf, :qvalue)
CSV.write("data/fsea_o12.csv")