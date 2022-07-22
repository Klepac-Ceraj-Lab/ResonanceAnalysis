# #Machine Learning Models

# ## Loading packages

using Resonance
using CSV
using DataFrames
using MLJ
using StableRNGs
ml_rng = StableRNG(102528)

# ## Parsing tidy CSV data

rawDf = CSV.read(
    joinpath("data", "tidy_timepoints_with_brain.csv"),
    DataFrame;
    delim = ',',
    types = [Int64, Int64, String, String],
    header = ["subject", "timepoint", "variable", "value"],
    skipto = 2
)

# ### Fixing some metadata names (see source for details)

#fix_metadata_colnames!(rawDf)

# ## Removing duplicates

# ### Removing actual duplicated lines

check_longdata_metaduplicates!(rawDf; remove_duplicates=true)

# ### Removing technical/biological replicates

retainfirst = true # if true, retain the first technical/biological replicate; if false, retain the last technical/biological replicate

if retainfirst

    unstackedDf = @chain rawDf begin
        groupby( [:subject, :timepoint, :variable] )
        combine(:value=>first)
        unstack(:variable, :value_first)
    end

else 
    unstackedDf = unstack(rawDf, :variable, :value; allowduplicates=true)
end

## Prediction workflows

### Problem 1: predict current cogScore from current taxonomic profile

#### Building predictionDf

select_columns = [ "cogScore", "ageMonths", microbiome_predictors... ]

predictionDf = @chain unstackedDf begin
    select(select_columns)
    subset(:ageMonths => x -> .!(ismissing.(x)))
    subset(:cogScore => x -> .!(ismissing.(x)))
    subset(microbiome_predictors[1] => x -> .!(ismissing.(x)))
    select(Not(:ageMonths))
end

mapcols!(col -> parse.(Float64, col), predictionDf)

#### Checking the schema and MLJ scitypes

schema(predictionDf)

#### Separating the predictors and target variable, vewing matching machine models

y, X = unpack(predictionDf, ==(:cogScore); rng=ml_rng);

for m in models(matching(X,y)) println(m) end

#### Building the MLJ machine for target problem.

RandomForestRegressor = @load RandomForestRegressor pkg=DecisionTree
rf_model = RandomForestRegressor()
rf_machine = machine(rf_model, X, y)

#### Splitting training data between train and test

train, test = partition(eachindex(y), 0.7, shuffle=true, rng=ml_rng)

#### Training the Machine on `train` rows of `X`

fit!(rf_machine, rows=train)

#### Evaluating the Machine on `test` rows of `X`

y_hat = predict(rf_machine, rows=test)

abs_errors = MLJ.MLJBase.l1(y_hat, y[test]) # Absolute ErrorS
mean_abs_error = mean(abs_errors) # Mean Absolute Error
sq_errors = MLJ.MLJBase.l2(y_hat, y[test]) # Square ErrorS
mean_sq_error = mean(sq_errors) # Mean Square Error

# It didn't perform very well. Let's try to tune some hyperparameters...
#
# Hyperparameters of DecisionTree.RandomForestRegressor:
# ```
# max_depth = -1, 
# min_samples_leaf = 1, 
# min_samples_split = 2, 
# min_purity_increase = 0.0, 
# n_subfeatures = -1, 
# n_trees = 10, 
# sampling_fraction = 0.7
# ````

rf_model = RandomForestRegressor()

range_max_depth = range(rf_model, :min_samples_leaf; lower=1, upper=5)

tuned_rf_model = TunedModel(
    model = rf_model,
    range = range_max_depth,
    resampling = CV(nfolds=3),
    tuning = Grid(resolution=10),
    measure = mae
)

rf_machine = machine(tuned_rf_model, X, y)

#### Training the Machine on `train` rows of `X`

fit!(rf_machine, rows=train)