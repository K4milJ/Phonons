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
	rand_matrix = 2 * rand_rng * rand(Float16, N_x, N_y, 2) .- rand_rng  #matrix of random displacements
    result_matrix_diag = zeros(Float64, N_x*N_y, N_x*N_y)			#preparing matrix of zeros for diag elems
    diag_ind = 1
    for i in 1:N_x
        for j in 1:N_y
            #println("diagind=", diag_ind, "x=", i, "y=", j)
            result_matrix_diag[diag_ind,diag_ind] = matrix_diag_elem(M, Ω, C_3, i, j, N_x, N_y, rand_matrix) #calculating diagonal element
            diag_ind += 1
        end
    end
    result_matrix_upper = zeros(N_x * N_y, N_x * N_y)			#calculating upper half of the matrix
    # Iterate through all pairs of indices
    for row_index in 1:N_x * N_y
        # Convert row index to (i, j) in the original matrix
        i, j = divrem(row_index - 1, N_x)
        i += 1
        j += 1
        for col_index in (row_index + 1):(N_x * N_y)
            # Convert column index to (k, l) in the original matrix
            k, l = divrem(col_index - 1, N_y)
            k += 1
            l += 1
            dist = sqrt(
                ((i-1)*a + rand_matrix[i, j, 1] - (k-1)*a - rand_matrix[k, l, 1])^2
                +
                ((j-1)*a + rand_matrix[i, j, 2] - (l-1)*a - rand_matrix[k, l, 2])^2
			)
            println("[i,j] = [$i,$j]\t [k,l] = [$k, $l]\t [col,row] = [$col_index,$row_index]\tdist = $dist")
            result_matrix_upper[col_index, row_index] = dist#^-5
        end
    end
	printstyled("Half of matrix\n", color=:red, bold=true)
	display(result_matrix_upper)
	printstyled("Diagonal elements\n", color=:red, bold=true)
	display(result_matrix_diag)
    return result_matrix_diag .+ result_matrix_upper .+ transpose(result_matrix_upper)
end
################################################################################

N_x = 2
N_y = 3
a = 1.1
M = 1.0
Ω = 1.0
C_3 = 1.5
rand_rng = 0.0

result_matrix = calculate_matrix2d(N_x, N_y, M, Ω, C_3, rand_rng)

scatter(eigvals(result_matrix), title="Eigenvalues: 2D, size: $N_x x $N_y, δR=$rand_rng, OBC", framestyle = :box)
#eigvals(result_matrix)
################################################################################
#	N_x, N_y	- size of the system
#	a 			- lattice constant
#	M			- mass of each atom
#	Ω, C_3		- constants
#	rand_rng	- determines the range (δR_i ∈ [-rand_rng,rand_rng]) from which
#				we add random vector [δR_x, δR_y] to the position of each atom
################################################################################