
function load_knead(dir=ENV["ANALYSIS_FILES"])
    run(Cmd([
        "kneaddata_read_count_table",
        "--input", joinpath(dir, "kneaddata/"),
        "--output", joinpath(dir, "kneaddata", "kneaddata_read_counts.tsv")
    ]))

    df = CSV.read("/grace/sequencing/processed/mgx/kneaddata/kneaddata_read_counts.tsv", DataFrame; missingstring="NA")
    for col in names(df, r"Homo_sapiens")
        df[!, col] = map(p-> coalesce(p...), zip(df[!, col], df[!, replace(col, "Homo_sapiens"=>"hg37dec_v0.1")]))
    end
    df."decontaminated Homo_sapiens pair1" = map(x-> ismissing(x) ? missing :  x isa Real ? x : parse(Float64, x), df."decontaminated Homo_sapiens pair1")
    df."decontaminated Homo_sapiens orphan2" = map(x-> ismissing(x) ? missing :  x isa Real ? x : parse(Float64, x), df."decontaminated Homo_sapiens orphan2")

    df = select(df, Not(r"hg37"))
    CSV.write("data/read_counts.csv", df)
    return df
end
