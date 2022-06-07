module Resonance

export airtable_metadata,
       brain_ingest,
       findprevstool,
       count_set,
       upset_dots!,
       load_metabolites,
       pull_row,
       loadings,
       varexplained,
       mdsaxis,
       commonname,
       countmap,
       codebreastfeeding!,
       stp_overlap,
       startup

using Reexport
using Airtable
using CairoMakie
using AlgebraOfGraphics
using MultivariateStats
using Dictionaries
using CodecZlib
using Statistics
using Arrow
using Tables
using SparseArrays
using CodecZlib
using FilePaths
using ProgressLogging
using CategoricalArrays

@reexport using BiobakeryUtils
@reexport using DataFrames
@reexport using CSV
@reexport using XLSX
@reexport using DataFramesMeta

include("airtable.jl")
include("wrangle.jl")
include("metabolites.jl")
include("genefamilies.jl")
include("plotting.jl")
include("brain.jl")
include("startup.jl")
include("kneaddata.jl")

end
