# #Metadata wrangling

# ## Get sample-specific metadata from airtable database

using Resonance

samplemeta = airtable_metadata() # having set ENV["AIRTABLE_KEY"]
# Now to add important metadata to samples

@chain samplemeta begin
    groupby([:subject, :timepoint])
    @transform!(:has_metabolomics = any(!ismissing, :Metabolomics_batch))
    sort!([:subject, :timepoint])
end


fmp_samples = DataFrame(XLSX.readtable("data/resonance_fmp/Fecal_All_033022.xlsx", "Sheet1", infer_eltypes=true)...)
rename!(fmp_samples, Dict(:studyID=>:subject, :collectionNum=> :timepoint))
@rsubset! fmp_samples begin
    :subject in samplemeta.subject 
    !ismissing(:timepoint)
end

samplemeta = leftjoin(samplemeta, select(fmp_samples, [:subject, :timepoint, :collectionDate]), on=[:subject, :timepoint])
unique!(samplemeta)

##

# ## Metadata stored in filemaker pro database
#
# Read in the different tables as `DataFrame`s,
# then normalize certain columns.
# First, subject-specific data

fmp_subject = DataFrame(XLSX.readtable("data/resonance_fmp/Subject_Centric_060222.xlsx", "Sheet1", infer_eltypes=true)...)
rename!(fmp_subject, Dict(:studyID=>:subject))
@rsubset! fmp_subject :subject in samplemeta.subject

# Then, timepoint-specific data

fmp_timepoint = DataFrame(XLSX.readtable("data/resonance_fmp/Timepoint_Centric_033122.xlsx", "Sheet1", infer_eltypes=true)...)
rename!(fmp_timepoint, Dict(:studyID=>:subject))

# # and COVID-specific samples

fmp_covid = DataFrame(XLSX.readtable("data/resonance_fmp/COVID_Fecal_040122.xlsx", "Sheet1", infer_eltypes=true)...)
rename!(fmp_covid, Dict(:studyID=>:subject))
fmp_covid = leftjoin(fmp_covid, fmp_subject, on=:subject)

# ## Getting data joined together
@rsubset!(fmp_timepoint, !ismissing(:subject))
fmp_alltp = leftjoin(fmp_timepoint, fmp_subject, on=[:subject])


codebreastfeeding!(fmp_alltp)

count(!ismissing, fmp_alltp.bfcalculated)

# Add info about brain data (just if it's there)

brain = let 
    fcleft = brain_ingest("data/brain/freesurfer_curvature_leftHemi_oct2021.csv"; label="curvature_leftHemi")
    fcright = brain_ingest("data/brain/freesurfer_curvature_rightHemi_oct2021.csv"; label="curvature_rightHemi")
    ftleft = brain_ingest("data/brain/freesurfer_thickness_leftHemi_oct2021.csv"; label="thickness_leftHemi")
    ftright = brain_ingest("data/brain/freesurfer_thickness_rightHemi_oct2021.csv"; label="thickness_rightHemi")
    seg = brain_ingest("data/brain/segmentationVolumeMeasurements_oct2021.csv"; label="segmentation")

    outerjoin(fcleft, fcright, ftleft, ftright, seg, on=[:subject, :timepoint], makeunique=true)
end
brain.has_segmentation .= true

## Validation

for row in eachrow(brain)
    all(h-> !ismissing(h) && h, row[r"has_"]) || @info row[r"has_"]
end

fmp_alltp = leftjoin(fmp_alltp, brain, on=[:subject, :timepoint])
fmp_alltp.has_segmentation = [ismissing(x) ? false : x for x in fmp_alltp.has_segmentation]

# transform!(fmp_alltp, :sid_old => ByRow(id-> ismissing(id) ? id : replace(id, r"_(\d+)F_"=>s"_\1E_")) => :sid_old_etoh)
# etoh_map = Dict((old=>new for (old, new) in zip(samplemeta.sid_old, samplemeta.sample)))
# transform!(fmp_alltp, :sid_old_etoh => ByRow(id-> (ismissing(id) || !haskey(etoh_map, id)) ? missing : etoh_map[id]) => :sample_etoh)

sort!(fmp_alltp, [:subject, :timepoint])

# ## Getting values for presence of data types
#
# We want to plot set intersections for having various kinds of data.
# In some cases, additional wrangling is necessary



fmp_alltp.race = map(fmp_alltp."Merge_Dem_Child_Race") do r
    ismissing(r) && return missing
    r == "Unknown" && return missing
    r == "Decline to Answer" && return missing
    contains(r, "\n") && return "Mixed"
    r ∈ ("Mixed", "Mixed Race") && return "Mixed"
    r ∈ ("Other Asian", "Asian ") && return "Asian"
    return r
end

# unique!(fmp_alltp, [:subject, :timepoint])

fmp_alltp.has_race = .!ismissing.(fmp_alltp."race")

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

# fmp_alltp = DataFrames.transform(subj, AsTable([:timepoint, :sample]) => (s -> begin
#     has_stool = .!ismissing.(s.sample)
#     has_prevstool = fill(false, length(s[1]))
#     for i in 1:length(s[1])
#         i == 1 && continue
#         any(!ismissing, s.sample[1:i-1]) && (has_prevstool[i] = true)
#     end
#     (; has_stool, has_prevstool)
# end)=> [:has_stool, :has_prevstool])

fmp_alltp.has_everbreast = .!ismissing.(fmp_alltp."everBreastFed")
fmp_alltp.has_bfperc = .!ismissing.(fmp_alltp."breastFedPercent")


# ## Dealing with ages

fmp_alltp.ageMonths = map(eachrow(fmp_alltp)) do row
    if ismissing(row.scanAgeMonths)
        (row.assessmentAgeDays ./ 365 .* 12) .+ row.assessmentAgeMonths
    else
        (row.scanAgeDays ./ 365 .* 12) .+ row.scanAgeMonths
    end
end
count(row-> !ismissing(row."Blood Pressure::Diastolic") && !ismissing(row.ageMonths),
                eachrow(fmp_alltp))

fmp_alltp.age0to3mo   = map(a-> !ismissing(a) && a < 3,       fmp_alltp.ageMonths)
fmp_alltp.age3to6mo   = map(a-> !ismissing(a) && 3 <= a < 6,  fmp_alltp.ageMonths)
fmp_alltp.age6to12mo  = map(a-> !ismissing(a) && 6 <= a < 12, fmp_alltp.ageMonths)
fmp_alltp.age12moplus = map(a-> !ismissing(a) && 12 <= a,     fmp_alltp.ageMonths)

for n in names(fmp_alltp)
    fmp_alltp[!, n] = map(fmp_alltp[!, n]) do v
        v isa AbstractString ? replace(v, r"\n"=>"___") : v
    end
end

CSV.write("data/wrangled/timepoints.csv", fmp_alltp)
CSV.write("data/wrangled/omnisamples.csv", @rsubset(samplemeta, :Fecal_EtOH == "F"))
CSV.write("data/wrangled/etohsamples.csv", @rsubset(samplemeta, :Fecal_EtOH == "E"))
CSV.write("data/wrangled/covid.csv", fmp_covid)

##

using Airtable

base = AirBase("appSWOVVdqAi5aT5u")
tab = AirTable("Samples", base)
joinedsamples = leftjoin(select(samplemeta, [:airtable_id, :subject, :timepoint, :sample, :childAgeMonths]), 
                         select(fmp_alltp, [:subject, :timepoint, :ageMonths, :has_segmentation]), 
                         on=[:subject, :timepoint]
)

agepatches = AirRecord[]
brainpatches = AirRecord[]

for row in eachrow(joinedsamples)
    if !ismissing(row[:childAgeMonths]) && ismissing(row[:ageMonths])
        push!(agepatches, AirRecord(row[:airtable_id], tab, (; childAgeMonths=row[:ageMonths])))
    end
    if coalesce(row[:has_segmentation], false)
        push!(brainpatches, AirRecord(row[:airtable_id], tab, (; concurrent_scan=row[:has_segmentation])))
    end
end

!isempty(agepatches) && Airtable.patch!(tab, agepatches)

unique!(brainpatches)
!isempty(brainpatches) && Airtable.patch!(tab, brainpatches)