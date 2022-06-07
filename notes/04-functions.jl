using Resonance

# - Takes ~ 6min for 1500 samples:
# - Writes to ENV["SCRATCH_SPACE"], which is "./scratch/" unless otherwise set
Resonance.write_gfs_arrow()
