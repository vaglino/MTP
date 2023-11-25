using Revise
using Plots: display
using Images, ImageView, ImageIO
using ImageSegmentation
using Statistics, Random, Plots
using PyCall
using TiffImages

src = "C:/Users/stravaglino3/Documents/MTP_jl-main/src/"
includet(src*"segmentation.jl")
includet(src*"background.jl")
includet(src*"image_processing.jl")
includet(src*"QC.jl")
includet(src*"analysis_functions.jl")
includet(src*"scripting_functions.jl")

using .IlluminationCorrection, .GaussianBackgroundSubtraction#, .RollingBallSubtraction

# SET PARAMETERS
px_size = 0.1083333 #um
px_area = px_size^2 #um²

# DATA DIRECTORY (set up data directory as follows)
# data_dir
#       |-> Condition 1
#          |-> image_1.tif
#          ...
#          ...
#          |-> image_N.tif
#       ...
#       ...
#       |-> Condition N

res_dir = "Y:/zhu-lab/Stefano/MTP/20231013_MTP_mNIP vs IgM and IgD after sort_4.7 pN_locker/analyzed_cells/"

#=_______________________________________________________________________________________________________________________
# IgM on monomeric NIP gradient (20231013)
_______________________________________________________________________________________________________________________=#
data_dir = "Y:/zhu-lab/Stefano/MTP/20231013_MTP_mNIP vs IgM and IgD after sort_4.7 pN_locker/IgM/"

# check conditions for each subdirectory
conds = get_conditions(data_dir)
print(keys(conds))

# # Find MTP background illumination by averaging the BSA surfaces
illumination_MTP = average_illumination(conds["IgM_mNIP_3"][2];channel=1, directory=conds["IgM_mNIP_3"][1])
imshow(illumination_MTP)
# Find RICM background illumination by averaging the BSA surfaces
illumination_RICM = average_illumination(conds["IgM_mNIP_3"][2];channel=2, directory=conds["IgM_mNIP_3"][1])
imshow(illumination_RICM)

results = process_data(data_dir; illumination = illumination_MTP)
QC_results, ΣFIs_all, MFIs_all, SAs_all = perform_QC(results)

plotly()
# plot_results(ΣFIs_all, MFIs_all, SAs_all)
plot_results(QC_results)

#-----------------------------------------------------------------------
## SAVE ANALYSIS
using JLD2
save_conditions_to_csv(QC_results, res_dir*"20231013_IgM_mNIP")
save_object(res_dir*"20231013_IgM_mNIP_cells.jld2", results)
save_object(res_dir*"illumination_MTP_IgM.jld2", illumination_MTP)
save_object(res_dir*"illumination_RICM_IgM.jld2", illumination_RICM)
# illumination_MTP = load_object(res_dir*"illumination_MTP.jld2")
# illumination_RICM = load_object(res_dir*"illumination_RICM.jld2")
#-----------------------------------------------------------------------



#=_______________________________________________________________________________________________________________________
# IgD on monomeric NIP gradient (20231014)
_______________________________________________________________________________________________________________________=#

data_dir = "Y:/zhu-lab/Stefano/MTP/20231013_MTP_mNIP vs IgM and IgD after sort_4.7 pN_locker/IgD/"
# check conditions for each subdirectory
conds = get_conditions(data_dir)
print(keys(conds))

# # Find MTP background illumination by averaging the BSA surfaces
illumination_MTP = average_illumination(conds["IgD_mNIP_18"][2];channel=1, directory=conds["IgD_mNIP_18"][1])
imshow(illumination_MTP)
# Find RICM background illumination by averaging the BSA surfaces
illumination_RICM = average_illumination(conds["IgD_mNIP_18"][2];channel=2, directory=conds["IgD_mNIP_18"][1])
imshow(illumination_RICM)

gr()
results = process_data(data_dir; illumination = illumination_MTP)
QC_results, ΣFIs_all, MFIs_all, SAs_all = perform_QC(results)

plotly()
# plot_results(ΣFIs_all, MFIs_all, SAs_all)
plot_results(QC_results)

##-----------------------------------------------------------------------
## SAVE ANALYSIS
using JLD2
save_conditions_to_csv(QC_results, res_dir*"20231013_IgD_mNIP")
save_object(res_dir*"20231013_IgD_mNIP_cells.jld2", results)
save_object(res_dir*"illumination_MTP_IgD.jld2", illumination_MTP)
save_object(res_dir*"illumination_RICM_IgD.jld2", illumination_RICM)
# illumination_MTP = load_object(res_dir*"illumination_MTP.jld2")
# illumination_RICM = load_object(res_dir*"illumination_RICM.jld2")
##-----------------------------------------------------------------------



#=_______________________________________________________________________________________________________________________
# IgG on monomeric NIP gradient (20231023)
_______________________________________________________________________________________________________________________=#
data_dir = "Y:/zhu-lab/Stefano/MTP/20231013_MTP_mNIP vs IgM and IgD after sort_4.7 pN_locker/IgG/"

# check conditions for each subdirectory
conds = get_conditions(data_dir)
print(keys(conds))

# # Find MTP background illumination by averaging the BSA surfaces
illumination_MTP = average_illumination(conds["IgG_mNIP_3"][2];channel=1, directory=conds["IgG_mNIP_3"][1])
imshow(illumination_MTP)
# Find RICM background illumination by averaging the BSA surfaces
illumination_RICM = average_illumination(conds["IgG_mNIP_3"][2];channel=2, directory=conds["IgG_mNIP_3"][1])
imshow(illumination_RICM)

results = process_data(data_dir; illumination = illumination_MTP)
QC_results, ΣFIs_all, MFIs_all, SAs_all = perform_QC(results)

plotly()
# plot_results(ΣFIs_all, MFIs_all, SAs_all)
plot_results(QC_results)

#-----------------------------------------------------------------------
## SAVE ANALYSIS
using JLD2
save_conditions_to_csv(QC_results, res_dir*"20231023_IgG_mNIP")
save_object(res_dir*"20231023_IgG_mNIP_cells.jld2", results)
save_object(res_dir*"illumination_MTP_IgG.jld2", illumination_MTP)
save_object(res_dir*"illumination_RICM_IgG.jld2", illumination_RICM)
# illumination_MTP = load_object(res_dir*"illumination_MTP.jld2")
# illumination_RICM = load_object(res_dir*"illumination_RICM.jld2")
#-----------------------------------------------------------------------



