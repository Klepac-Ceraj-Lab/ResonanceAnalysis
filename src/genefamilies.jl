function write_gfs_arrow(dir=ENV["ANALYSIS_FILES"])
    allfiles = String[]

    @info "getting files"
    for (root, dirs, files) in walkdir(dir)
        filter!(f-> contains(f, "genefamilies") && !contains(f, "relab") && !contains(f, "zymo") && !contains(f, r"^FE\d+"), files)
        append!(allfiles, joinpath.(Ref(root), files))
    end

    unique!(allfiles) do f
        first(split(basename(f), '_')) |> String
    end

    samples = Set(Iterators.map(f-> first(split(basename(f), '_')) |> String, allfiles))
    sdic = Dictionary((s for s in samples), 1:length(samples))
    features = Set(String[])
    
    scratch = get(ENV, "SCRATCH_SPACE", "./scratch")
    isdir(scratch) || mkpath(scratch)

    tmp = tempname(scratch)

    @info "writing temporary arrow file"
    open(tmp, "w") do io
        tbls = Tables.partitioner(allfiles) do f
            samplename = replace(splitext(basename(f))[1], r"_S\d+_genefamilies" => "")

            sdf = CSV.read(f, DataFrame; header=["feature", "value"], skipto=2)
            subset!(sdf, "feature"=> ByRow(f-> !contains(f, '|'))) # skip stratified features
            sdf.sample .= samplename
            sdf.sidx .= sdic[samplename]

            union!(features, sdf.feature)

            sdf
        end

        Arrow.write(io, tbls)
    end

    @info "Making feature dictionary"
    fdic = Dictionary((f for f in features), 1:length(features))

    @info "building new table"
    df = DataFrame(Arrow.Table(tmp))
    df.fidx = [fdic[f] for f in df.feature]

    @info "Writing table"
    Arrow.write(joinpath(scratch, "genefamilies.arrow"), df)
    open(joinpath(scratch, "features.txt"), "w") do io
        for f in keys(fdic)
            println(io, f)
        end
    end
    open(joinpath(scratch, "samples.txt"), "w") do io
        for s in keys(sdic)
            println(io, s)
        end
    end
    return nothing
end

function read_gfs_arrow()
    scratch = get(ENV, "SCRATCH_SPACE", "./scratch")
    @info "reading table"
    tbl = Arrow.Table(joinpath(scratch, "genefamilies.arrow"))
    @info "building sparse mat"
    mat = sparse(tbl.fidx, tbl.sidx, tbl.value)
    @info "getting features"
    fs = [GeneFunction(line) for line in eachline(joinpath(scratch, "features.txt"))]
    @info "getting samples"
    ss = [MicrobiomeSample(line) for line in eachline(joinpath(scratch, "samples.txt"))]
    return CommunityProfile(mat, fs, ss)
end


function get_neuroactive_kos(neuroactivepath="data/gbm.txt")
    neuroactive = Dictionary{String, Vector{String}}()
    desc = ""
    for line in eachline(neuroactivepath)
       line = split(line, r"[\t,]")
       if startswith(line[1], "MGB")
           (mgb, desc) = line
           desc = rstrip(replace(desc, r"\bI+\b.*$"=>""))
           desc = replace(desc, r" \([\w\s]+\)$"=>"")
           desc = replace(desc, r"^.+ \(([\w\-]+)\) (.+)$"=>s"\1 \2")
           @info "getting unirefs for $desc"
           !in(desc, keys(neuroactive)) && insert!(neuroactive, desc, String[])
       else
           filter!(l-> occursin(r"^K\d+$", l), line)
           append!(neuroactive[desc], String.(line))
       end
   end
   return neuroactive
end

function getneuroactive(features; neuroactivepath="data/gbm.txt", map_ko_uniref_path="data/map_ko_uniref90.txt.gz")
    neuroactivekos = get_neuroactive_kos(neuroactivepath)

    kos2uniref = Dictionary{String, Vector{String}}()
    for line in eachline(GzipDecompressorStream(open(map_ko_uniref_path)))
        line = split(line, '\t')
        insert!(kos2uniref, line[1], map(x-> String(match(r"UniRef90_(\w+)", x).captures[1]), line[2:end]))
    end

    neuroactive_index = Dictionary{String, Vector{Int}}()
    for na in keys(neuroactivekos)
        searchfor = Iterators.flatten([kos2uniref[ko] for ko in neuroactivekos[na] if ko in keys(kos2uniref)]) |> Set
        pos = findall(u-> u in searchfor, features)
        insert!(neuroactive_index, na, pos)
    end
    for k in keys(neuroactive_index)
        unique!(neuroactive_index[k])
    end
    return neuroactive_index
end

fsea(cors, pos) = (cors, pos, MannWhitneyUTest(cors[pos], cors[Not(pos)]))

function fsea(cors, allfeatures::AbstractVector, searchset::Set)
    pos = findall(x-> x in searchset, allfeatures)

    return fsea(cors, pos)
end

function fsea(occ::AbstractMatrix, metadatum::AbstractVector, pos::AbstractVector{<:Int})
    cors = let notmissing = map(!ismissing, metadatum)
        occ = occ[:, notmissing]
        metadatum = metadatum[notmissing]
        cor(metadatum, occ, dims=2)'
    end

    return fsea(cors, pos)
end