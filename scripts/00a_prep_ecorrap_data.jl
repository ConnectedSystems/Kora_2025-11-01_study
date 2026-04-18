using CSV
using DataFrames
using Statistics
using CairoMakie

OUTPUT_DIR = joinpath(@__DIR__, "..", "data")

df = CSV.read(
    joinpath(OUTPUT_DIR, "ecorrap_benthic", "EcoRRAP_Benthic_Data_Moore_TorresStrait_2026-01-21.csv"),
    DataFrame
)

mapping_df = CSV.read(
    joinpath(OUTPUT_DIR, "ecorrap_benthic", "EcoRRAP_Labelset_mapping_Maren_Toor_2026-01-20.csv"),
    DataFrame;
    missingstring="NA"
)

mapped_taxa = unique(mapping_df[:, [:DESCRIPTION, :MAPPING]])
mapped_taxa = mapped_taxa[.!ismissing.(mapped_taxa.MAPPING), :]
sub_groups_of_interest = String.(mapped_taxa.DESCRIPTION)

target_site = "ONMO"
considered_taxas = unique(mapped_taxa.MAPPING)

mean_cover = DataFrame(Dict(x => Float64[] for x in considered_taxas))
insertcols!(mean_cover, 1, "year" => Int64[])
tmp = Dict(String(x) => [] for x in considered_taxas)
tmp["transect"] = []
tmp["year"] = []

for year in unique(df.year)
    df_year = df[(df.year .== year) .& (df.site .== target_site), :]
    gdf = groupby(df_year, [:transect])

    for (t_id, t) in enumerate(gdf)
        push!(tmp["year"], year)
        push!(tmp["transect"], getproperty(keys(gdf)[t_id], :transect))
        for taxa in unique(mapped_taxa.MAPPING)
            matching_taxa = mapped_taxa[mapped_taxa.MAPPING .== taxa, :DESCRIPTION]
            matching_taxa = [replace(m, "." => " ") for m in matching_taxa]

            # @info taxa matching_taxa
            push!(tmp[taxa], mean(sum(eachcol(t[:, matching_taxa]))))
        end
    end
    tmp2 = DataFrame(tmp)
    @info tmp2
    if isempty(tmp2[tmp2.year .== year, :])
        @info "No records found for $(year)"
        continue
    end

    push!(
        mean_cover,
        [
            Int64(year),
            map(mean, eachcol(tmp2[tmp2.year .== year, Not(["transect", "year"])]))...
        ]
    )
end

# Reorder dataframe to be in expected group order
group_order = [
    "Tabular Acropora",
    "Corymbose Acropora",
    "branching non-Acropora",
    "Small massive",
    "Large massive"
]
mean_cover = mean_cover[:, ["year", group_order...]]

CSV.write(joinpath(OUTPUT_DIR, "ecorrap_benthic", "moore_estimate.csv"), mean_cover)

## Repeat to extract data for Masig Reef!
target_site = "TSMA"

mean_cover = DataFrame(Dict(x => Float64[] for x in considered_taxas))
insertcols!(mean_cover, 1, "year" => Int64[])
tmp = Dict(String(x) => [] for x in considered_taxas)
tmp["transect"] = []
tmp["year"] = []

for year in unique(df.year)
    df_year = df[(df.year .== year) .& (df.site .== target_site), :]
    gdf = groupby(df_year, [:transect])

    for (t_id, t) in enumerate(gdf)
        push!(tmp["year"], year)
        push!(tmp["transect"], getproperty(keys(gdf)[t_id], :transect))
        for taxa in unique(mapped_taxa.MAPPING)
            matching_taxa = mapped_taxa[mapped_taxa.MAPPING .== taxa, :DESCRIPTION]
            matching_taxa = [replace(m, "." => " ") for m in matching_taxa]

            # Determine mean for taxa at this transect
            push!(tmp[taxa], mean(sum(eachcol(t[:, matching_taxa]))))
        end
    end

    tmp2 = DataFrame(tmp)
    if isempty(tmp2[tmp2.year .== year, :])
        @info "No records found for $(year)"
        continue
    end

    # Determine mean for taxa across this site
    push!(
        mean_cover,
        [
            Int64(year),
            map(mean, eachcol(tmp2[tmp2.year .== year, Not(["transect", "year"])]))...
        ]
    )
end

# Reorder so that groups are in consistent order
mean_cover = mean_cover[:, ["year", group_order...]]

CSV.write(joinpath(OUTPUT_DIR, "ecorrap_benthic", "masig_estimate.csv"), mean_cover)

benthic_cover = DataFrame(;
    report_year=mean_cover.year,
    mean=map(sum, eachrow(mean_cover[:, 2:end])),
    lower=map(sum, eachrow(mean_cover[:, 2:end])),
    upper=map(sum, eachrow(mean_cover[:, 2:end]))
)

CSV.write(joinpath(OUTPUT_DIR, "Masig_Reef_EcoRRAP_estimate.csv"), benthic_cover)
