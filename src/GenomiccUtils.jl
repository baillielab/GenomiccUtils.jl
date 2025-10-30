module GenomiccUtils

using CSV
using DataFrames
using DelimitedFiles
using Base.Threads

include("read_write.jl")
include("ancestry.jl")

export read_bim, read_fam, read_map, read_ped
export write_map, write_ped
export estimate_ancestry

end
