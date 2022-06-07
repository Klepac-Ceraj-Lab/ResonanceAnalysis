_letter2number(l) = findfirst(==(only(l)), "abcdefghijklmnopqrstuvwxyz")


function split_tp(id)
    m = match(r"(^\d+)(\w)?$", id)
    isnothing(m) && throw(ArgumentError("`$id` does not match expected pattern: `^\\d+\\w?\$`"))

    sub = parse(Int, m[1])
    tp = isnothing(m[2]) ? 1 : _letter2number(m[2])
    return (subject=sub, timepoint=tp)
end

function brain_ingest(file; label=nothing)
    df = CSV.read(file, DataFrame)
    df = hcat(df, DataFrame(split_tp.(df[!, 1])))

    !isnothing(label) && (df[:, "has_$label"] = fill(true, nrow(df)))
    return df
end

"""
Given vector of timepoints and a boolean vector
indicating if matching timepoint has a stool sample,
returns a boolean vector indicating if each timepoint
has an earlier timepoint that has a stool sample.
"""
function findprevstool(timepoints, has_stool; rev=false)
    length(timepoints) == 1 && return [false]
    srt = sortperm(timepoints; rev)
    fs = findfirst(has_stool[srt])
    isnothing(fs) && return fill(false, length(timepoints))
    return (timepoints .> fs) .& has_stool
end

countmap(v) = Dict(k=> count(x -> !ismissing(x) && x == k, v) for k in unique(v))

## Breastfeeding

_codebreastfeeding(::Missing) = missing
_codebreastfeeding(bf::Bool) = bf
_codebreastfeeding(bf::AbstractString) = lowercase(bf) == "yes" ? true : lowercase(bf) == "no" ? false : error("invalid Breastfeeding code: $bf")

function _codebreastfeeding(bf::Real)
    0 <= bf <= 100 || error("invalid breastfeeding percent")
    if bf == 100
        return "breastmilk"
    elseif bf == 0
        return "formula"
    else
        return "mixed"
    end
end

function _processbfpercent(bfp)
    n = length(bfp)
    f = findfirst(!ismissing, bfp)
    isnothing(f) && return fill(missing, n)
    
    codes = _codebreastfeeding.(bfp)
    code = codes[f]
    codes[1:f] .= code
    if code == "mixed"
        return fill("mixed", n)
    elseif f == n
        return fill(code, n)
    else
        f2 = findfirst(x-> !ismissing(x) && x != code, codes[f+1:end])
        if isnothing(f2)
            return fill(code, n)
        else
            f2 += f
            codes[f:f2-1] .= code
            codes[f2:end] .= "mixed"
            return String.(codes)
        end
    end
end

function _extractbf(col)
    vals = unique(skipmissing(col))
    isempty(vals) && return missing
    return only(vals)
end

function _codebreastfeeding(dfg)
    dfg = DataFrame(dfg)
    n = nrow(dfg)
    eff = _extractbf(dfg."BreastfeedingDone::exclusiveFormulaFed")
    ebf = _extractbf(dfg."BreastfeedingDone::exclusivelyBottlefedBreastmilk")
    en = _extractbf(dfg."BreastfeedingDone::exclusivelyNursed")
    
    @assert issorted(dfg.timepoint)

    bfp = _processbfpercent(dfg."BreastfeedingStill::breastFedPercent")

    if all(ismissing, [eff, ebf, en])
        bfe = fill(missing, n)
    elseif !any(x-> !ismissing(x) && x, [eff, ebf, en])
        @warn "$(only(unique(dfg.subject))) has odd breastfeeding information"
        bfe = fill("mixed", n)
    elseif !ismissing(eff) && eff
        if any(x->!ismissing(x) && x, (ebf, en))
            @warn "$(only(unique(dfg.subject))) has contradictory breastfeeding information"
            bfe = fill("mixed", n)
        else
            bfe = fill("formula", n)
        end
    else
        bfe = fill("breastmilk", n)
    end
    return map(x-> coalesce(x...), zip(bfp, bfe))
end

function codebreastfeeding!(df::DataFrame)
    bfcols = ["BreastfeedingDone::exclusiveFormulaFed",
          "BreastfeedingDone::exclusivelyBottlefedBreastmilk",
          "BreastfeedingDone::exclusivelyNursed",
          "BreastfeedingStill::breastFedPercent"
    ]

    for col in bfcols[1:3]
        df[!, col] = _codebreastfeeding.(df[!, col])
    end

    sort!(df, [:subject, :timepoint])
    gdf = groupby(df, :subject)
    DataFrames.transform!(gdf, AsTable(["subject", "timepoint", bfcols...]) => _codebreastfeeding => :bfcalculated)
end

function _keep_unique(ds1, ds2, srt)
    recip = Set(ds2)
    keep = Bool[]
    used = Set(eltype(ds1)[])
    for d in ds1[srt]
        if d ∈ used
            push!(keep, false)
        else
            @info d maxlog = 10
            push!(used, d)
            @info d ∈ recip maxlog = 10
            push!(keep, d ∈ recip)
        end
    end
    @info count(keep)
    return keep
end

"""
Given 2 vectors of (:subject, :timepoint) tuples,
find the overlap, and put them in the same order.
"""
function stp_overlap(ds1, ds2; lt = (x,y)-> x[1] == y[1] ? x[2] < y[2] : x[1] < y[1]) 
    srt1 = sortperm(ds1; lt)
    srt2 = sortperm(ds2; lt)
    
    keep1 = _keep_unique(ds1, ds2, srt1)
    keep2 = _keep_unique(ds2, ds1, srt2)

    return srt1[keep1], srt2[keep2]
end