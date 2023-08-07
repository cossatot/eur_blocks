using Revise

using Oiler

using CSV
using JSON
using DataFrames:eachrow
using DataFrames, DataFramesMeta
using Setfield

using PyPlot

# options
fault_weight = 2.
save_results = true


# load data
eur_block_file = "../block_data/eur_blocks.geojson"
eur_fault_file = "../block_data/eur_faults.geojson"


ana_block_file = "../../anatolia/block_data/anatolia_blocks.geojson"
ana_fault_file = "../../anatolia/block_data/anatolia_faults.geojson"
weiss_vel_field_file = "../../anatolia/geod_data/weiss_et_al_2020_vels_down_100.geojson"

daug_vels_file = "../strain_data/daugostino_vels.geojson"
gsrm_vels_file = "../../c_asia_blocks/gnss_data/gsrm_c_asia_vels.geojson"

#boundary_file = "../block_data/cea_gnss_block_domain.geojson"
#boundary_file = "../block_data/cea_hazard_boundary.geojson"
boundary_file = "../block_data/emme23_bounds.geojson"


@info "joining blocks"
eur_block_df = Oiler.IO.gis_vec_file_to_df(eur_block_file)
ana_block_df = Oiler.IO.gis_vec_file_to_df(ana_block_file)
block_df = vcat(eur_block_df, 
                ana_block_df,
                cols=:union)

#@info "culling blocks"
#println("n blocks before ", size(block_df, 1))
#bound_df = Oiler.IO.gis_vec_file_to_df(boundary_file)
#block_df = Oiler.IO.get_blocks_in_bounds!(block_df, bound_df)
#println("n blocks after ", size(block_df, 1))

@info "doing faults"
fault_df, faults, fault_vels = Oiler.IO.process_faults_from_gis_files(
                                                        ana_fault_file,
                                                        eur_fault_file,
                                                        block_df=block_df,
                                                        subset_in_bounds=true,
                                                        check_blocks=true,
                                                        )
fault_df[:,:fid] = string.(fault_df[:,:fid])
println("n faults: ", length(faults))
println("n fault vels: ", length(fault_vels))

@info "doing non-fault block boundaries"
@time non_fault_bounds = Oiler.IO.get_non_fault_block_bounds(block_df, faults)
bound_vels = vcat(map(b->Oiler.Boundaries.boundary_to_vels(b, ee=1.0, en=1.0), 
                      non_fault_bounds)...)
println("n non-fault-bound vels: ", length(bound_vels))

@info "doing GNSS"
gsrm_vel_df = Oiler.IO.gis_vec_file_to_df(gsrm_vels_file)
daug_vel_df = Oiler.IO.gis_vec_file_to_df(daug_vels_file)

@time gsrm_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gsrm_vel_df, block_df;
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:station,
    fix="1111"
)

@time daug_vels = Oiler.IO.make_vels_from_gnss_and_blocks(daug_vel_df, block_df;
    name=:site, fix="1111",
   )



@info "doing weiss vels"
weiss_vel_field_df = Oiler.IO.gis_vec_file_to_df(weiss_vel_field_file)
weiss_vel_field_df[!,"station"] = map(x->join(["weiss_", x]), 
                                      string.(weiss_vel_field_df[!,:fid]))

weiss_vel_field_vels = Oiler.IO.make_vels_from_gnss_and_blocks(
    weiss_vel_field_df, block_df;
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:station,
    fix="1112"
)
    
gnss_vels = vcat(daug_vels,
                 gsrm_vels,
                 weiss_vel_field_vels,
                 )

println("n gnss vels: ", length(gnss_vels))
#println("n vel field vels: ", length(vel_field_vels))


@info "doing geol slip rates"
#chn_slip_rate_df = Oiler.IO.gis_vec_file_to_df(chn_slip_rate_file)
#geol_slip_rate_df = vcat(cea_slip_rate_df, 
#                         chn_slip_rate_df, 
                         #nea_slip_rate_df
#                         )

#geol_slip_rate_df, geol_slip_rate_vels = Oiler.IO.make_geol_slip_rate_vels!(
#                                                geol_slip_rate_df,
#                                                fault_df)

#println("n geol slip rates: ", length(geol_slip_rate_vels))

# tris

vels = vcat(fault_vels,
            bound_vels,
            gnss_vels, 
            #geol_slip_rate_vels, 
            )

vel_groups = Oiler.group_vels_by_fix_mov(vels);


# solve
results = Oiler.solve_block_invs_from_vel_groups(vel_groups; faults=faults,
                                                 tris=[],
                                                sparse_lhs=true,
                                                weighted=true,
                                                elastic_floor=1e-2,
                                                regularize_tris=true,
                                                tri_priors=false,
                                                tri_distance_weight=20.,
                                                predict_vels=true,
                                                check_closures=true,
                                                pred_se=false,
                                                constraint_method="kkt_sym",
                                                factorization="lu")

#Oiler.ResultsAnalysis.compare_data_results(results=results,
#                                           vel_groups=vel_groups,
#                                           #geol_slip_rate_df=geol_slip_rate_df,
#                                           #geol_slip_rate_vels=geol_slip_rate_vels,
#                                           fault_df=fault_df,
#                                           )

#block_bound_df = Oiler.IO.gis_vec_file_to_df("../block_data/c_asia_block_bounds.geojson")

function make_bound_fault(row)
    trace = Oiler.IO.get_coords_from_geom(row[:geometry])

    Oiler.Fault(trace=trace, dip_dir=row[:dip_dir], dip=89., hw=row[:hw],
                fw=row[:fw], fid=row[:fid])
end

#@info "getting block rates"
#bound_faults = []
#for i in 1:size(block_bound_df, 1)
#    push!(bound_faults, make_bound_fault(block_bound_df[i,:]))
#end
#
#block_bound_rates = Oiler.Utils.get_fault_slip_rates_from_poles(bound_faults,
#                                                                results["poles"],
#                                                                use_path=true)
if save_results == true
    Oiler.IO.write_fault_results_to_gj(results, 
        "../results/eur_fault_results.geojson",
        name="Europe and Mediterranean fault results",
        calc_rake=true,
        calc_slip_rate=true)

    Oiler.IO.write_gnss_vel_results_to_csv(results, vel_groups;
        name="../results/eur_gnss_results.csv")

end

map_fig = Oiler.Plots.plot_results_map(results, vel_groups, faults)

#slip_rate_fig = Oiler.Plots.plot_slip_rate_fig(geol_slip_rate_df,
#    geol_slip_rate_vels, fault_df, results)

Oiler.WebViewer.write_web_viewer(results=results, block_df=block_df,
                                 ref_pole="1111", directory="../web_viewer")

show()


