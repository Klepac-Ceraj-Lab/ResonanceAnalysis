using Resonance

profiles = String[]
for (root, dirs, files) in walkdir("/grace/sequencing/processed/mgx")
    contains(root, "metaphlan") || continue
    contains(root, "links") && continue
    filter!(f-> contains(f, r"^FG.+profile"), files)
    append!(profiles, joinpath.(Ref(root), files))
end

species = metaphlan_profiles(profiles, :species)
                                        
sns = map(samplenames(species)) do s
    replace(s, r"_S\d+_profile"=>"")
end

idx = unique([findfirst(==(x), sns) for x in sns])

species = CommunityProfile(abundances(species)[:, idx], features(species), map(samplenames(species)[idx]) do s
    MicrobiomeSample(replace(s, r"_S\d+_profile"=>""))
end)

CSV.write("data/species.csv", species)