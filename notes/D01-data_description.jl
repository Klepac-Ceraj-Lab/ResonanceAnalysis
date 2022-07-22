# #Machine Learning Models

# ## Loading packages

using Resonance
using CSV
using DataFrames
using MLJ

# ## Parsing tidy CSV data

rawDf = CSV.read(
    joinpath("data", "tidy_timepoints_with_brain.csv"),
    DataFrame;
    delim = ',',
    types = [Int64, Int64, String, String],
    header = ["subject", "timepoint", "variable", "value"],
    skipto=2
)

# ## Removing duplicates

# ### Removing actual duplicates
let
    n_unique_rows = size(unique(rawDf[:, 1:3]),1)

    if (n_unique_rows != size(rawDf, 1))

        n_nonunique = size(rawDf, 1) - n_unique_rows

        @warn "Long Dataframe contains non-unique rows! **Removing duplicated data**"
        @warn "After removal, $(n_unique_rows) will remain. $(n_nonunique) rows were rmoved from the original $(size(rawDf, 1))"

    end

end

rawDf = unique!(rawDf)

# ### Removing biological replicates with same metadata

# replicate_lines = [141436, 141445]
# delete!(rawDf, replicate_lines)

# ## Unstacking dataFrame

retainfirst = true

if retainfirst

    cgl = combine(groupby(rawDf,[:subject,:timepoint, :variable]),:value=>first)
    unstackedDf = unstack(cgl, :variable, :value_first)

else
    
    unstackedDf = unstack(rawDf, :variable, :value; allowduplicates=true)

end

colnames = names(unstackedDf)

# ## Description of the Dataset

# ### Group 00: Large-scale summaries

# #### All samples and subjects
g00_df = copy(unstackedDf)
println("Beginning data description:\n")
println("$(size(g00_df, 1)) unique samples from $(length(unique(g00_df.subject))) unique subjects.\n")

# #### Child samples and subjects
g00C_df = g00_df[.!(ismissing.(g00_df.ageMonths)), :]
println("Children:\n")
println("$(size(g00C_df, 1)) unique samples from $(length(unique(g00C_df.subject))) subjects.")
println("of which")
# #### Children samples with no information available (Empty)
g00CE_df = copy(g00C_df)
g00CE_df = g00CE_df[ismissing.(g00CE_df.Absiella_dolichum), :]
g00CE_df = g00CE_df[ismissing.(g00CE_df.cogScore), :]
g00CE_df = g00CE_df[ismissing.(g00CE_df.Left_Thalamus), :]
println("$(size(g00CE_df, 1)) unique samples from $(length(unique(g00CE_df.subject))) subjects have no data collected.\n")

# #### Mother samples and subjects
g00M_df = g00_df[ismissing.(g00_df.ageMonths), :]
println("Mothers:\n")
println("$(size(g00M_df, 1)) unique samples from $(length(unique(g00M_df.subject))) subjects.")
println("of which")
# #### Mother samples with no information available (Empty)
g00ME_df = copy(g00M_df)
g00ME_df = g00ME_df[ismissing.(g00ME_df.Absiella_dolichum), :]
g00ME_df = g00ME_df[ismissing.(g00ME_df.cogScore), :]
g00ME_df = g00ME_df[ismissing.(g00ME_df.Left_Thalamus), :]
println("$(size(g00ME_df, 1)) unique samples from $(length(unique(g00ME_df.subject))) subjects have no data collected.\n")

n_subjects_motheronly = sum( .!( unique(g00M_df.subject) .∈ Ref(unique(g00C_df.subject)) ) )
n_subjects_childonly = sum( .!( unique(g00C_df.subject) .∈ Ref(unique(g00M_df.subject)) ) )
n_subjects_motherchild = sum( unique(g00M_df.subject) .∈ Ref(unique(g00C_df.subject)) )

println("Of all $(length(unique(g00_df.subject))) unique subjects:
    $(n_subjects_motheronly) have only Mother samples;
    $(n_subjects_childonly) have only Child samples;
    $(n_subjects_motherchild) have both mother and child samples.\n")


g00C_Scb_df = g00C_df[.!(ismissing.(g00C_df.Absiella_dolichum)),:]
g00C_sCb_df = g00C_df[.!(ismissing.(g00C_df.cogScore)),:]
g00C_scB_df = g00C_df[.!(ismissing.(g00C_df.Left_Thalamus)),:]
  
g00C_SCb_df = g00C_Scb_df[.!(ismissing.(g00C_Scb_df.cogScore)),:]
g00C_sCB_df = g00C_sCb_df[.!(ismissing.(g00C_sCb_df.Left_Thalamus)),:]
g00C_ScB_df = g00C_scB_df[.!(ismissing.(g00C_scB_df.Absiella_dolichum)),:]
    
println("Breakdown of children samples:
       $(size(g00C_Scb_df, 1)) child samples have at least stool;
       $(size(g00C_sCb_df, 1)) child samples have at least cogScore;
       $(size(g00C_scB_df, 1)) child samples have at least brain;
       $(size(g00C_SCb_df, 1)) child samples have at least stool AND cogScore;
       $(size(g00C_sCB_df, 1)) child samples have at least cogScore AND brain;
       $(size(g00C_ScB_df, 1)) child samples have at least stool AND brain.\n")
    
# ### Details of the Children samples

# #### Group 01: Child samples for which we have stool sample but **not** brain or cogscore

g01_df = copy(unstackedDf)
g01_df = g01_df[.!(ismissing.(g01_df.ageMonths)), :]
g01_df = g01_df[.!(ismissing.(g01_df.Absiella_dolichum)), :]
g01_df = g01_df[ismissing.(g01_df.cogScore), :]
g01_df = g01_df[ismissing.(g01_df.Left_Thalamus), :]

println("Group 01: Child samples for which we have stool sample but **not** brain or cogscore
    $(size(g01_df, 1)) unique samples from $(length(unique(g01_df.subject))) unique subjects.\n")

# #### Group 02: Child samples for which we have cogScores but **not** brain or stool

g02_df = copy(unstackedDf)
g02_df = g02_df[.!(ismissing.(g02_df.ageMonths)), :]
g02_df = g02_df[ismissing.(g02_df.Absiella_dolichum), :]
g02_df = g02_df[.!(ismissing.(g02_df.cogScore)), :]
g02_df = g02_df[ismissing.(g02_df.Left_Thalamus), :]

println("Group 02: Child samples for which we have cogScores but **not** brain or stool
    $(size(g02_df, 1)) unique samples from $(length(unique(g02_df.subject))) unique subjects.\n")

# #### Group 03: Child samples for which we have brain scans but **not** stool or cogscore

g03_df = copy(unstackedDf)
g03_df = g03_df[.!(ismissing.(g03_df.ageMonths)), :]
g03_df = g03_df[ismissing.(g03_df.Absiella_dolichum), :]
g03_df = g03_df[ismissing.(g03_df.cogScore), :]
g03_df = g03_df[.!(ismissing.(g03_df.Left_Thalamus)), :]

println("Group 03: Child samples for which we have brain scans but **not** stool or cogscore
    $(size(g03_df, 1)) unique samples from $(length(unique(g03_df.subject))) unique subjects.\n")

# #### Group 04: Child samples for which we have concurrent stool + cogScore but **not** brain

g04_df = copy(unstackedDf)
g04_df = g04_df[.!(ismissing.(g04_df.ageMonths)), :]
g04_df = g04_df[.!(ismissing.(g04_df.Absiella_dolichum)), :]
g04_df = g04_df[.!(ismissing.(g04_df.cogScore)), :]
g04_df = g04_df[ismissing.(g04_df.Left_Thalamus), :]

println("Group 04: Child samples for which we have concurrent stool + cogScore but **not** brain
    $(size(g04_df, 1)) unique samples from $(length(unique(g04_df.subject))) unique subjects.\n")

# #### Group 05: Child samples for which we have concurrent stool + brain but **not** cogscore

g05_df = copy(unstackedDf)
g05_df = g05_df[.!(ismissing.(g05_df.ageMonths)), :]
g05_df = g05_df[.!(ismissing.(g05_df.Absiella_dolichum)), :]
g05_df = g05_df[ismissing.(g05_df.cogScore), :]
g05_df = g05_df[.!(ismissing.(g05_df.Left_Thalamus)), :]

println("Group 05: Child samples for which we have concurrent stool + brain but **not** cogscore
    $(size(g05_df, 1)) unique samples from $(length(unique(g05_df.subject))) unique subjects.\n")

# #### Group 06: Child samples for which we have concurrent cogScore + brain but **not** stool

g06_df = copy(unstackedDf)
g06_df = g06_df[.!(ismissing.(g06_df.ageMonths)), :]
g06_df = g06_df[ismissing.(g06_df.Absiella_dolichum), :]
g06_df = g06_df[.!(ismissing.(g06_df.cogScore)), :]
g06_df = g06_df[.!(ismissing.(g06_df.Left_Thalamus)), :]

println("Group 06: Child samples for which we have concurrent cogScore + brain but **not** stool
    $(size(g06_df, 1)) unique samples from $(length(unique(g06_df.subject))) unique subjects.\n")

# #### Group 7: Child samples for which we have concrrent stool + cogScore + brain

g07_df = copy(unstackedDf)
g07_df = g07_df[.!(ismissing.(g07_df.ageMonths)), :]
g07_df = g07_df[.!(ismissing.(g07_df.Absiella_dolichum)), :]
g07_df = g07_df[.!(ismissing.(g07_df.cogScore)), :]
g07_df = g07_df[.!(ismissing.(g07_df.Left_Thalamus)), :]

println("Group 07: Child samples for which we have concurrent stool + cogScore + brain
    $(size(g07_df, 1)) unique samples from $(length(unique(g07_df.subject))) unique subjects.\n")

# #### Group 08: Child samples for which we have **stool and any future cogScore** (for prediction)

g08_df = build_future_df(g00C_df, :cogScore)

println("Group 08: Child samples for which we have **stool and any future cogScore**
    $(size(g08_df, 1)) datapoints from
    $(size(unique(g08_df[:, [:subject, :timepoint]]), 1)) unique samples from
    $(size(unique(g08_df[:, :subject]), 1)) unique subjects.\n")

# #### Group 09: Child samples with **ageMonths < 6** for which we have **stool and any future cogScore** (for prediction)

g09_df = g08_df[parse.(Float64, g08_df.ageMonths) .<= 6.0, :]

println("Group 09: Child samples with **ageMonths < 6** for which we have **stool and any future cogScore**
    $(size(g09_df, 1)) datapoints from
    $(size(unique(g09_df[:, [:subject, :timepoint]]), 1)) unique samples from
    $(size(unique(g09_df[:, :subject]), 1)) unique subjects.\n")

# #### Group 09: Child samples with **ageMonths < 6** for which we have **stool and any future cogScore** (for prediction)

g09_df = g08_df[parse.(Float64, g08_df.ageMonths) .<= 6.0, :]

println("Group 09: Child samples with **ageMonths < 6** for which we have **stool and any future cogScore**
    $(size(g09_df, 1)) datapoints from
    $(size(unique(g09_df[:, [:subject, :timepoint]]), 1)) unique samples from
    $(size(unique(g09_df[:, :subject]), 1)) unique subjects.\n")

# #### Group 10: Child samples with **ageMonths < 6** for which we have **stool and at least one future cogScore between 1-2 years later** (for prediction)

g10_df = g09_df[(g09_df.ageMonthsDelta .<= 24.0) .& (g09_df.ageMonthsDelta .>= 12.0), :]

println("Group 10: Child samples with **ageMonths < 6** for which we have **stool and at least one future cogScore between 1-2 years later**
    $(size(g10_df, 1)) datapoints from
    $(size(unique(g10_df[:, [:subject, :timepoint]]), 1)) unique samples from
    $(size(unique(g10_df[:, :subject]), 1)) unique subjects.\n")

# #### Group 11: Child samples with **ageMonths < 9** for which we have **stool and any future cogScore** (for prediction)

g11_df = g08_df[parse.(Float64, g08_df.ageMonths) .<= 9.0, :]

println("Group 11: Child samples with **ageMonths < 9** for which we have **stool and any future cogScore**
    $(size(g11_df, 1)) datapoints from
    $(size(unique(g11_df[:, [:subject, :timepoint]]), 1)) unique samples from
    $(size(unique(g11_df[:, :subject]), 1)) unique subjects.\n")

# #### Group 12: Child samples with **ageMonths < 9** for which we have **stool and at least one future cogScore between 1-2 years later** (for prediction)

g12_df = g11_df[(g11_df.ageMonthsDelta .<= 24.0) .& (g11_df.ageMonthsDelta .>= 12.0), :]

println("Group 12: Child samples with **ageMonths < 9** for which we have **stool and at least one future cogScore between 1-2 years later**
    $(size(g12_df, 1)) datapoints from
    $(size(unique(g12_df[:, [:subject, :timepoint]]), 1)) unique samples from
    $(size(unique(g12_df[:, :subject]), 1)) unique subjects.\n")


# #### Group 13: Child samples with **ageMonths < 12** for which we have **stool and any future cogScore** (for prediction)

g13_df = g08_df[parse.(Float64, g08_df.ageMonths) .<= 12.0, :]

println("Group 13: Child samples with **ageMonths < 12** for which we have **stool and any future cogScore**
    $(size(g13_df, 1)) datapoints from
    $(size(unique(g13_df[:, [:subject, :timepoint]]), 1)) unique samples from
    $(size(unique(g13_df[:, :subject]), 1)) unique subjects.\n")

# #### Group 14: Child samples with **ageMonths < 12** for which we have **stool and at least one future cogScore between 1-2 years later** (for prediction)

g14_df = g13_df[(g13_df.ageMonthsDelta .<= 24.0) .& (g13_df.ageMonthsDelta .>= 12.0), :]

println("Group 14: Child samples with **ageMonths < 12** for which we have **stool and at least one future cogScore between 1-2 years later**
    $(size(g14_df, 1)) datapoints from
    $(size(unique(g14_df[:, [:subject, :timepoint]]), 1)) unique samples from
    $(size(unique(g14_df[:, :subject]), 1)) unique subjects.\n")

# #### Group 14: Child samples with **ageMonths < 12** for which we have **stool and at least one future cogScore between 1-2 years later** (for prediction)

g15_df = build_future_df(g00C_df, :Left_Thalamus)

println("Group 15: Child samples for which we have **stool and any future brain scan**
    $(size(g15_df, 1)) datapoints from
    $(size(unique(g15_df[:, [:subject, :timepoint]]), 1)) unique samples from
    $(size(unique(g15_df[:, :subject]), 1)) unique subjects.\n")