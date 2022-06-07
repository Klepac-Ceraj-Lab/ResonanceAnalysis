const _keep_fields = [
    :sample,
    :subject,
    :timepoint,
    :ECHOTPCoded,
    :Mother_Child,
    :MaternaID,
    :Fecal_EtOH,
    :CollectionRep,
    :sid_old,
    :CovidCollectionNumber,
    :childAgeMonths
]

"""
    airtable_metadata(key=ENV["AIRTABLE_KEY"])
Get fecal sample metadata table from airtable.
The API `key` comes from https://airtable.com/account.
This is unlikely to work if you're not in the VKC lab,
but published sample metadata is available from OSF.io
using `datadep"sample metadata"`.
"""
function airtable_metadata()
    base = AirBase("appSWOVVdqAi5aT5u")

    samples      = Airtable.query(AirTable("Samples", base))
    mgxbatches   = Airtable.query(AirTable("MGX Batches", base))
    metabbatches = Airtable.query(AirTable("Metabolomics Batches", base))

    df = DataFrame()
    
    @progress "Adding samples" for samp in samples
        # skip anything not in the "ECHO" project
        "rec8zdkF1pchFNiJC" ∈ get(samp, :Project, []) || continue

        mgxids = string.(get(samp, Symbol("MGX Batches"), []))
        mgxbatch = findall(batch-> Airtable.id(batch) ∈ mgxids, mgxbatches)
        mgxbatch = isempty(mgxbatch) ? missing : mgxbatches[first(mgxbatch)][:Name]

        metabids = string.(get(samp, Symbol("Metabolomics Batches"), []))
        metabbatch = findall(batch-> Airtable.id(batch) ∈ metabids, metabbatches)
        metabbatch = isempty(metabbatch) ? missing : metabbatches[first(metabbatch)][:Name]

        rec =  (; :airtable_id        => Airtable.id(samp),
                  (kf=> get(samp, kf, missing) for kf in _keep_fields)...,
                  :Mgx_batch          => mgxbatch, 
                  :Metabolomics_batch => metabbatch, 
                  )
        push!(df, rec, cols=:union)
    end

    subset!(df, AsTable(["subject", "timepoint"])=> ByRow(r-> !any(ismissing, r)))

    transform!(df, "subject"   => ByRow(s-> parse(Int, s)) => "subject",
                   "timepoint" => ByRow(tp-> parse(Int, tp)) => "timepoint")
    return select(df, Cols(:sample, :subject, :timepoint, :))
end