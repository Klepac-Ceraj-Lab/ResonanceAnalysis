# #Metadata wrangling

# ## Get sample-specific metadata from airtable database

using Resonance

samplemeta = airtable_metadata() # having set ENV["AIRTABLE_KEY"]
# Now to add important metadata to samples

stools = MicrobiomeSample[]

let fields = [
    "timepoint",
    "subject",
    "sid_old",
    "CovidCollectionNumber",
    "Mgx_batch",
    "Metabolomics_batch",
    "MaternalID"]

    for row in eachrow(samplemeta)
        s = MicrobiomeSample(row.sample)
        set!(s, NamedTuple(Symbol(k) => v for (k, v) in pairs(row[Not(:sample)]) if !ismissing(v)))
        push!(stools, s)
    end
end

stool_all = DataFrame()
for stool in stools
    push!(stool_all, NamedTuple(pairs(metadata(stool))), cols=:union)
end


stool_all = @chain stool_all begin
    groupby([:subject, :timepoint])
    transform(:has_metabolomics = any(!ismissing, :Metabolomics_batch))
end

##

# ## Metadata stored in filemaker pro database
#
# Read in the different tables as `DataFrame`s,
# then normalize certain columns.
# First, samples

fmp_sample = DataFrame(XLSX.readtable("data/Sample_Centric_10252021.xlsx", "Sheet1", infer_eltypes=true)...)
rename!(fmp_sample, Dict(:SampleID=> :sample, :studyID=>:subject, :collectionNum=> :timepoint))
@subset!(fmp_sample, @byrow !startswith(:sample, "FE"))

# Then, subject-specific data

fmp_subject = DataFrame(XLSX.readtable("data/Subject_Centric_10252021.xlsx", "Sheet1", infer_eltypes=true)...)
rename!(fmp_subject, Dict(:studyID=>:subject))

# Then, timepoint-specific data

fmp_timepoint = DataFrame(XLSX.readtable("data/Timepoint_Centric_10252021.xlsx", "Sheet1", infer_eltypes=true)...)
rename!(fmp_timepoint, Dict(:studyID=>:subject))

# and COVID-specific samples

fmp_covid = DataFrame(XLSX.readtable("data/COVID_Fecal_10252021.xlsx", "Sheet1", infer_eltypes=true)...)
rename!(fmp_covid, Dict(:studyID=>:subject))

# ## Getting data joined together

fmp_timed = outerjoin(fmp_timepoint,
                      @subset(fmp_sample, @byrow !ismissing(:timepoint)), # filter out covid samples
                      on=[:subject, :timepoint])

fmp_alltp = leftjoin(fmp_timed, fmp_subject, on=[:subject])
unique!(fmp_alltp, [:subject, :timepoint])

# Add info about brain data (just if it's there)

brain = let 
    fcleft = brain_ingest("data/freesurfer_curvature_leftHemi_oct2021.csv"; label="curvature_leftHemi")
    fcright = brain_ingest("data/freesurfer_curvature_rightHemi_oct2021.csv"; label="curvature_rightHemi")
    ftleft = brain_ingest("data/freesurfer_thickness_leftHemi_oct2021.csv"; label="thickness_leftHemi")
    ftright = brain_ingest("data/freesurfer_thickness_rightHemi_oct2021.csv"; label="thickness_rightHemi")
    seg = brain_ingest("/home/kevin/Repos/Resonance/data/segmentationVolumeMeasurements_oct2021.csv"; label="segmentation")

    outerjoin(fcleft, fcright, ftleft, ftright, seg, on=[:subject, :timepoint], makeunique=true)
end

## Validation

for row in eachrow(brain)
    all(h-> !ismissing(h) && h, row[r"has_"]) || @info row[r"has_"]
end

fmp_alltp = leftjoin(fmp_alltp, brain, on=[:subject, :timepoint])
fmp_alltp = leftjoin(fmp_alltp, stool_all, on=[:subject, :timepoint], makeunique=true)

sort!(fmp_alltp, [:subject, :timepoint])

# ## Getting values for presence of data types
#
# We want to plot set intersections for having various kinds of data.
# In some cases, additional wrangling is necessary

fmp_alltp.simple_race = map(fmp_alltp.simple_race) do r
    ismissing(r) && return missing
    r == "Unknown" && return missing
    r == "Decline to Answer" && return missing
    r ∈ ("Mixed", "Mixed Race") && return "Mixed"
    return r
end

unique!(fmp_alltp, [:subject, :timepoint])

fmp_alltp.has_race = .!ismissing.(fmp_alltp."simple_race")

##

subj = groupby(fmp_alltp, :subject)
transform!(subj, nrow => :n_samples; ungroup=false)

transform!(subj, AsTable(r"subject|breast"i) => (s -> begin
    if any(!ismissing, s.breastFedPercent)
        return fill(true, length(s[1]))
    else
        return fill(false, length(s[1]))
    end
    fill(0, length(s[1]))
end) => :has_bfperc_subj; ungroup=false)

transform!(subj, AsTable(r"subject|breast"i) => (s -> begin
    if any(!ismissing, s.breastFedPercent)
        return fill(true, length(s[1]))
    else
        return fill(false, length(s[1]))
    end
    fill(0, length(s[1]))
end) => :has_bfperc_subj; ungroup=false)

fmp_alltp = transform(subj, AsTable([:timepoint, :sample]) => (s -> begin
    has_stool = .!ismissing.(s.sample)
    has_prevstool = fill(false, length(s[1]))
    for i in 1:length(s[1])
        i == 1 && continue
        any(!ismissing, s.sample[1:i-1]) && (has_prevstool[i] = true)
    end
    (; has_stool, has_prevstool)
end)=> [:has_stool, :has_prevstool])

fmp_alltp.has_everbreast = .!ismissing.(fmp_alltp."everBreastFed")
fmp_alltp.has_bfperc = .!ismissing.(fmp_alltp."breastFedPercent")
fmp_alltp.has_everbreast = .!ismissing.(fmp_alltp."Metabolomics_batch")


# ## Dealing with ages

fmp_alltp.ageMonths = map(eachrow(fmp_alltp)) do row
    if ismissing(row.scanAgeMonths)
        (row.assessmentAgeDays ./ 365 .* 12) .+ row.assessmentAgeMonths
    else
        (row.scanAgeDays ./ 365 .* 12) .+ row.scanAgeMonths
    end
end

fmp_alltp.age0to3mo   = map(a-> !ismissing(a) && a < 3,       fmp_alltp.ageMonths)
fmp_alltp.age3to6mo   = map(a-> !ismissing(a) && 3 <= a < 6,  fmp_alltp.ageMonths)
fmp_alltp.age6to12mo  = map(a-> !ismissing(a) && 6 <= a < 12, fmp_alltp.ageMonths)
fmp_alltp.age12moplus = map(a-> !ismissing(a) && 12 <= a,     fmp_alltp.ageMonths)

CSV.write("data/wrangled.csv", fmp_alltp)

