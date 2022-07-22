using Resonance, CategoricalArrays

#include(joinpath(pathof(Resonance), "..", "..", "notes", "01-metadata_wrangle.jl"))
#include(joinpath(pathof(Resonance), "..", "..", "notes", "02-species.jl"))
#include(joinpath(pathof(Resonance), "..", "..", "notes", "03-metabolites.jl"))
#include(joinpath(pathof(Resonance), "..", "..", "notes", "04-functions.jl"))
include(joinpath(pathof(Resonance), "..", "..", "src", "startup.jl"))

omni, etoh, tps, complete_brain, metabolites, species = startup()

for col in names(tps)
    tps[!, col] = collect(tps[!, col])
end

tps_specific = stack(tps, ["ageMonths", "cogScore", "ECHOTPCoded", brainmeta_underscore...], [:subject, :timepoint])
subj_specific = stack(tps, [:mother_HHS_Education, :rce], [:subject]) # changed simple_race to rce
grp = groupby(subj_specific, [:subject, :variable])
subj_specific = DataFrames.combine(grp, :value=> (v-> coalesce(v...))=> :value)

##

tps_specific = vcat(tps_specific, DataFrame(subject   = get(species, :subject), 
                                timepoint = get(species, :timepoint),
                                variable  = fill("sample", nsamples(species)),
                                value     = samplenames(species))
)

for s in samples(species)
    subject   = s.subject
    timepoint = s.timepoint
    for f in features(species)
        variable = name(f)
        value = species[f, s]
        push!(tps_specific, (; subject, timepoint, variable, value))
    end
end

CSV.write(datafiles("wrangled", "tidy_timepoints_with_brain.csv"), tps_specific)
CSV.write(datafiles("wrangled", "tidy_subjects.csv"), subj_specific)


##