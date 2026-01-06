include("src/utils_V1_V2.jl")
include("src/utils_V3.jl")
include("src/visualize.jl")
include("src/positions_calculate.jl")
include("src/matrix_calculations.jl")

using Plots
using LinearAlgebra
using CairoMakie

#generating positions
test_pos = calculate_positions_hex_PBC(8, 10, 1.0)

result_matrix = calculate_matrix_OBC_2D(8, 10, 1.0, 1.0, 1.0, test_pos)

f = visualize_phonon_mode_from_positions(result_matrix, test_pos, 1)
# save("plots/phonon_mode_1_hex_big.png", f)  # Optional: Save the figure

ev = Plots.scatter(eigvals(result_matrix), framestyle = :box, legend=:none)
# save("plots/phonon_mode_1_hex_big_eigvals.png", ev) 
