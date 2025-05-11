################################################################################
# 2D system - size: N_x*N_y, OBC, +[δR_x, δR_y]
################################################################################

using Plots
using LinearAlgebra

function sum_of_distances_diag_elem(M::Float64, Ω::Float64, C_3::Float64, i_index::Int, j_index::Int, N_x::Int, N_y::Int, rand_matrix::Array{Float64})
    result = 0
    for x_pos in 1:N_x
        for y_pos in 1:N_y
            if (i_index != x_pos && j_index != y_pos)
                result += (sqrt(((i_index-1)*a + rand_matrix[i_index, j_index, 1] - (x_pos-1)*a - rand_matrix[x_pos, y_pos, 1])^2 + 
                               ((j_index-1)*a + rand_matrix[i_index, j_index, 2] - (y_pos-1)*a - rand_matrix[x_pos, y_pos, 2])^2))^-5
            end
        end
    end
    return result
end

function matrix_diag_elem(M::Float64, Ω::Float64, C_3::Float64, i_index::Int, j_index::Int, N_x::Int, N_y::Int, rand_matrix::Array{Float64})
    return M*Ω^2 + 12*C_3*sum_of_distances_diag_elem(M, Ω, C_3, i_index, j_index, N_x, N_y, rand_matrix)
end

function calculate_matrix2d(N_x::Int, N_y::Int, M::Float64, Ω::Float64, C_3::Float64, rand_rng::Float64)
	rand_matrix = 2 * rand_rng * rand(Float16, N_x, N_y, 2) .- rand_rng
    result_matrix_diag = zeros(Float64, N_x*N_y, N_x*N_y)
    diag_ind = 1
    for i in 1:N_x
        for j in 1:N_y
            result_matrix_diag[diag_ind, diag_ind] = matrix_diag_elem(M, Ω, C_3, i, j, N_x, N_y, rand_matrix)
            diag_ind += 1
        end
    end
    result_matrix_upper = zeros(N_x * N_y, N_x * N_y)
    for row_index in 1:N_x * N_y
        i, j = divrem(row_index - 1, N_x)
        i += 1
        j += 1
        for col_index in (row_index + 1):(N_x * N_y)
            k, l = divrem(col_index - 1, N_x)
            k += 1
            l += 1
            dist = sqrt(
                ((i-1)*a + rand_matrix[j, i, 1] - (k-1)*a - rand_matrix[l, k, 1])^2 # Δx^2
                +
                ((j-1)*a + rand_matrix[j, i, 2] - (l-1)*a - rand_matrix[l, k, 2])^2 # Δy^2
			)
            result_matrix_upper[row_index, col_index] = dist^(-5)
        end
    end
    return result_matrix_diag .+ result_matrix_upper .+ transpose(result_matrix_upper)
end
################################################################################

N_x = 5
N_y = 7
a = 1.1
M = 1.0
Ω = 1.0
C_3 = 1.5
rand_rng = 0.3

result_matrix = calculate_matrix2d(N_x, N_y, M, Ω, C_3, rand_rng)

scatter(eigvals(result_matrix), title="Eigenvalues: 2D, size: $N_x x $N_y, δR=$rand_rng, OBC", framestyle = :box)

################################################################################
#	N_x, N_y	- size of the system
#	a 			- lattice constant
#	M			- mass of each atom
#	Ω, C_3		- constants
#	rand_rng	- determines the range (δR_i ∈ [-rand_rng,rand_rng]) from which
#				we add random vector [δR_x, δR_y] to the position of each atom
################################################################################
