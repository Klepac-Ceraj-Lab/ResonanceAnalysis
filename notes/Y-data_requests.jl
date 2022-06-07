# Vanja + Aditi CBCL

using Resonance
omni, etoh, tps, complete_brain, metabolites, species = startup()
genes = Resonance.read_gfs_arrow()

aditi_sids = readlines("data/aditi_samples.txt")

CSV.write(joinpath(ENV["SCRATCH_SPACE"], "aditi_taxa.csv"), species[:, aditi_sids])
CSV.write(joinpath(ENV["SCRATCH_SPACE"], "aditi_genefamilies.csv"), genes[:, aditi_sids])

 

## Korean group, blood pressure

using Resonance

fecalsamples = CSV.read("data/wrangled/samples.csv", DataFrame)
@rsubset! fecalsamples :Fecal_EtOH == "F"

timepoints = CSV.read("data/wrangled/timepoints.csv", DataFrame)

function mapbp(label)
    bp = Dict(
        "E2a" => "Systolic",
        "E2b" => "Diastolic",
        "E2c" => "Heart rate"
    )
    for pair in pairs(bp)
        bp = replace(label, pair)
    end

    # bp = replace(bp)
end


@assert size(subset(timepoints, "ChildAnthropometryOverTwo::E2a1"=> ByRow(!ismissing), "ChildAnthropometryUnderTwo::E2a1" => ByRow(!ismissing)), 1) == 0

timepoints.Systolic_1 = map(eachrow(timepoints)) do row
    coalesce(row."ChildAnthropometryOverTwo::E2a1", row."ChildAnthropometryUnderTwo::E2a1")
end
timepoints.Systolic_2 = map(eachrow(timepoints)) do row
    coalesce(row."ChildAnthropometryOverTwo::E2a2", row."ChildAnthropometryUnderTwo::E2a2")
end
timepoints.Systolic_3 = map(eachrow(timepoints)) do row
    coalesce(row."ChildAnthropometryOverTwo::E2a3", row."ChildAnthropometryUnderTwo::E2a3")
end
timepoints.Diastolic_1 = map(eachrow(timepoints)) do row
    coalesce(row."ChildAnthropometryOverTwo::E2b1", row."ChildAnthropometryUnderTwo::E2b1")
end
timepoints.Diastolic_2 = map(eachrow(timepoints)) do row
    coalesce(row."ChildAnthropometryOverTwo::E2b2", row."ChildAnthropometryUnderTwo::E2b2")
end
timepoints.Diastolic_3 = map(eachrow(timepoints)) do row
    coalesce(row."ChildAnthropometryOverTwo::E2b3", row."ChildAnthropometryUnderTwo::E2b3")
end
timepoints.HeartRate_1 = map(eachrow(timepoints)) do row
    coalesce(row."ChildAnthropometryOverTwo::E2c1", row."ChildAnthropometryUnderTwo::E2c1")
end
timepoints.HeartRate_2 = map(eachrow(timepoints)) do row
    coalesce(row."ChildAnthropometryOverTwo::E2c2", row."ChildAnthropometryUnderTwo::E2c2")
end
timepoints.HeartRate_3 = map(eachrow(timepoints)) do row
    coalesce(row."ChildAnthropometryOverTwo::E2c3", row."ChildAnthropometryUnderTwo::E2c3")
end

df = select(timepoints, ["subject", "timepoint",  "childGender", "ageMonths", "Systolic_1", "Systolic_2", "Systolic_3", "Diastolic_1", "Diastolic_2", "Diastolic_3", "HeartRate_1", "HeartRate_2", "HeartRate_3"])
df = leftjoin(df, select(fecalsamples, 
                    ["subject", "timepoint", "sample", "sid_old", "Mgx_batch"]
              ),
              on=["subject", "timepoint"])

count(groupby(df, :subject)) do grp
    !any(row-> (!ismissing(row.sample) && !ismissing(row.ageMonths) && (6 < row.ageMonths < 18)), eachrow(grp)) && return false
    any(row-> (!ismissing(row.Systolic_1) && !ismissing(row.ageMonths) && (36 < row.ageMonths < 84)), eachrow(grp))
end

count(groupby(df, :subject)) do grp
    !any(row-> (!ismissing(row.sample) && !ismissing(row.ageMonths) && (6 < row.ageMonths < 18)), eachrow(grp)) && return false
    any(row-> !ismissing(row.Systolic_1), eachrow(grp))
end

count(groupby(df, :subject)) do grp
    !any(row-> !ismissing(row.sample), eachrow(grp)) && return false
    any(row-> !ismissing(row.Systolic_1), eachrow(grp))
end

count(groupby(df, :subject)) do grp
    !any(row-> (!ismissing(row.sample) && !ismissing(row.ageMonths) && row.ageMonths < 24), eachrow(grp)) && return false
    any(row-> !ismissing(row.Systolic_1), eachrow(grp))
end

CSV.write("/home/kevin/Desktop/RESONANCE_bp.csv", sort(df, ["subject", "timepoint"]))

## Kirsty + 3rd quarter

using Resonance
omni, etoh, tps, complete_brain, metabolites, species = startup()

mo3 = @subset tps begin
    .!ismissing.(:ageMonths)
    2 .< :ageMonths .< 5
end

@rsubset! omni :sample in samplenames(species)
unique!(omni, [:subject, :timepoint])
mo3 = leftjoin(mo3, select(omni, [:subject, :timepoint, :sample]), on=[:subject, :timepoint])
@rsubset! mo3 !ismissing(:sample)
mo3spec = species[:, in.(samplenames(species), Ref(mo3.sample))]

set!(mo3spec, mo3)

#-

ginisimpson!(mo3spec)
prev = prevalence(mo3spec) |> vec
count(prev .> 0.56)
#- 

using CairoMakie
using AlgebraOfGraphics

#-

fig = Figure(; resolution = (1100, 850))
ax1 = Axis(fig[1,1]; title="Species prevalence", xlabel="prevalence", ylabel="number of taxa")
prevhist = hist!(ax1, prev)

ax2 = Axis(fig[2,1]; title="Age and diversity", xlabel="Age (months)", ylabel="Gini-Simpson Diversity")
scatter!(ax2, get(mo3spec, :ageMonths), get(mo3spec, :ginisimpson))

topprev = featurenames(mo3spec)[prev .> 0.55]


lcol = GridLayout()
fig.layout[1:2, 2] = lcol

df = vcat(collect(
            DataFrame("abundance" => vec(abundances(mo3spec[Regex(topprev[i]), :])),
                      "ageMonths" => get(mo3spec, :ageMonths),
                      "species"   => fill(topprev[i], nsamples(mo3spec))
                      ) for i in eachindex(topprev)
            )...
)


pl = data(df) * mapping(:ageMonths=> "Age (months)", :abundance=> "Relative abundance"; layout=:species)
draw!(lcol, pl; facet = (; linkyaxes = :none))

supertitle = Label(fig[0, 2], "Most prevalent taxa", textsize = 20)
# Legend(lcol[2,1], [sp1, sp2, sp3], topprev; tellheight=true, tellwidth=false, )

colsize!(fig.layout, 2, Relative(2/3))

fig

#- 

mo3pcoa = pcoa(mo3spec)

#- 

fig = Figure()

ax1 = Axis(fig[1,1], xlabel="MDS1 ($(round(varexplained(mo3pcoa)[1] * 100, digits=2))%)",
                     ylabel="MDS2 ($(round(varexplained(mo3pcoa)[2] * 100, digits=2))%)")

sc = scatter!(ax1, Resonance.loadings(mo3pcoa)[:,1], Resonance.loadings(mo3pcoa)[:,2],
        color=get(mo3spec, :ginisimpson)
)


ax2 = Axis(fig[1,2], xlabel="Age (months)",
ylabel="MDS1 ($(round(varexplained(mo3pcoa)[1] * 100, digits=2))%)")

sc2 = scatter!(ax2, get(mo3spec, :ageMonths), Resonance.loadings(mo3pcoa)[:,1],
    color=get(mo3spec, :ginisimpson)
)

Colorbar(fig[2,1:2], sc; label="Gini-Simpson diversity", vertical = false, flipaxis=false)

fig

#-


## Genetics info

using Resonance
using Statistics

omni, etoh, tps, complete_brain, metabolites, species = startup()

echoidmap = CSV.read("data/echoids.csv", DataFrame)
snps = CSV.read("data/85_snps_20220401.csv", DataFrame)
rename!(snps, :Column1=> :ECHOProtocolDChild)

snps = leftjoin(snps, select(echoidmap, ["subject", "ECHOProtocolDChild"]), on="ECHOProtocolDChild"; matchmissing=:notequal)
subset!(snps, :subject => ByRow(!ismissing))
snps = leftjoin(snps, combine(groupby(select(tps, [:subject, :childGender]), :subject),
                        :childGender=> (g-> coalesce(g...) => :childGender)
                        );
                on=:subject)



gentps = subset(tps, :subject=> ByRow(s-> !ismissing(s) && s ∈ snps.subject); skipmissing=true)

scores = combine(groupby(gentps, :subject), :cogScore => mean ∘ skipmissing => :mean_cog, :childGender=> first=> :childGender)

snps = select(leftjoin(snps, scores, on=:subject), Cols("ECHOProtocolDChild", "mean_cog", "childGender", :))
subset!(snps, :mean_cog => ByRow(!ismissing))

CSV.write("data/snps_cog.csv", snps)
