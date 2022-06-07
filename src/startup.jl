const mainmeta = [
        "ageMonths",
        "age0to3mo",
        "age3to6mo",
        "age6to12mo",
        "age12moplus",
        "mother_HHS_Education",
        "race",
        "childGender",
        "cogScore",
        "has_segmentation",
        "ECHOTPCoded",
        "assessmentDate",
        "ed",
        "rce"
    ]

    
const brainmeta = [
            #  "CortexVol",
            #  "CorticalWhiteMatterVol",
            #  "SubCortGrayVol",
            #  "TotalGrayVol",
            #  "BrainSegVol-to-eTIV",
            #  "CerebralWhiteMatterVol",
            #  "lhCorticalWhiteMatterVol",
            #  "lhCerebralWhiteMatterVol",
            #  "lhCortexVol",
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
            #  "rhCortexVol",
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

"""
Generates some data objects useful in many scripts.
Requires files the following files,
generated in notebooks 1-3:

- "data/wrangled/omnisamples.csv"
- "data/wrangled/etohsamples.csv"
- "data/wrangled/timepoints.csv"
- "data/wrangled/metabolites.csv"
- "data/wrangled/species.csv"
"""
function startup(; dfs=[:omni, :etoh, :tps, :complete_brain, :metabolites, :species])
    omni = CSV.read(joinpath(@__DIR__, "..", "data/wrangled/omnisamples.csv"), DataFrame)
    etoh = CSV.read(joinpath(@__DIR__, "..", "data/wrangled/etohsamples.csv"), DataFrame)
    
    tps, complete_brain = _gentps()
    
    metabolites = _genmetabolites(etoh, tps)
    
   species = _genspecies(omni, tps)

   return (; omni, etoh, tps, complete_brain, metabolites, species)[dfs]
end

function _gentps()
    tps  = CSV.read(joinpath(@__DIR__, "..", "data/wrangled/timepoints.csv"), DataFrame)

    DataFrames.transform!(groupby(tps, :subject), :mother_HHS_Education => (r->coalesce(r...)) => :hhs)
    tps.ed = categorical(tps.hhs; levels=[-8 , 2:7...], ordered=true)
    tps.ed = recode(tps.ed,
        -8 => missing,
        2 => "Junior high school",
        3 => "Some high school",
        4 => "High school grad",
        5 => "Some college",
        6 => "College grad",
        7 => "Grad/professional school")

    DataFrames.transform!(groupby(tps, :subject), :race => (r->coalesce(r...)) => :race)
    tps.rce = categorical(tps.race; ordered=true)
    tps.rce = recode(tps.rce, 
        "American Indian or Alaska Native"=> "Other",
        "Some other race"                 => "Other",
        "Asian Indian"                    => "Asian",
        "Black or African American"       => "Black",
        missing                           => "Unknown"
    )
    levels!(tps.rce, ["White","Black","Asian","Mixed","Other","Unknown"])

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
    brainmeta = Resonance.brainmeta
    for m in brainmeta
        tps[!, m] .= tps[!, m] ./ map(x-> ismissing(x) || x == 0 ? missing : x, tps."EstimatedTotalIntraCranialVol")
    end


    
    select!(tps, ["subject", "timepoint", mainmeta..., brainmeta...])
    rename!(tps, Dict(k=> replace(k, "-"=>"_") for k in brainmeta))
    brainmeta = map(i-> replace(brainmeta[i], "-"=>"_"), eachindex(brainmeta))
    
    complete_brain = completecases(tps[:, brainmeta])
    return tps, complete_brain
end

function _genmetabolites(etoh_samples, tps)
    metabolites = CSV.read(joinpath(@__DIR__, "..", "data/wrangled/metabolites.csv"), DataFrame)
    ms = [Resonance.Metabolite(row[:uid], row[:Metabolite], row[:MZ], row[:RT]) for row in eachrow(metabolites)]
    metabolites = CommunityProfile(Matrix(metabolites[!, 9:end]), ms, MicrobiomeSample.(names(metabolites)[9:end]))
    set!(metabolites, leftjoin(etoh_samples, tps[!, ["subject",  "timepoint", mainmeta...]], on=[:subject, :timepoint], makeunique=true))
    
    return metabolites
end

function _genspecies(omni_samples, tps)
    species = CSV.read(joinpath(@__DIR__, "..", "data/wrangled/species.csv"), DataFrame)
    species = CommunityProfile(Matrix(species[!, 2:end]), taxon.(species[!, 1]), MicrobiomeSample.(names(species)[2:end]))
    set!(species, leftjoin(omni_samples, tps[!, ["subject", "timepoint", mainmeta...]], on=[:subject, :timepoint], makeunique=true))
    species = species[:, map(!ismissing, get(species, :subject)) .& map(!ismissing, get(species, :timepoint))]

    return species
end
