"""
Function calculating final matrix (OBC 2D) for V1 and V2 (no anisotropy).
"""
function calculate_V1_V2_matrix_part(N_x::Int, N_y::Int, M::Float64, Ω::Float64, C_3::Float64, positions::Array{Float64})
	#creating matrix of 2x2 matrices of 0.0
	result_matrix_diag = [zeros(2, 2) for _ in 1:N_x*N_y, _ in 1:N_x*N_y]
	#calculating diagonal elements (2x2 matrices)
	diag_ind = 1
	for i in 1:N_x
		for j in 1:N_y
			result_matrix_diag[diag_ind, diag_ind] = matrix_diag_elem(M, Ω, C_3, i, j, N_x, N_y, positions)
			diag_ind += 1
		end
	end

	#making a big matrix from matrix of matrices
	result_matrix_diag = reduce(vcat, map(row -> reduce(hcat, row), eachrow(result_matrix_diag)))

	#creating matrix of 2x2 matrices of 0.0
	result_matrix_upper = [zeros(2, 2) for _ in 1:N_x*N_y, _ in 1:N_x*N_y]

	#filling the upper triangular part of a matrix
	for row_index in 1:N_x * N_y
		i, j = divrem(row_index - 1, N_x)
		i += 1
		j += 1
		for col_index in (row_index + 1):(N_x * N_y)
			k, l = divrem(col_index - 1, N_x)
			k += 1
			l += 1
			M_xx = -3*C_3*sum_of_distances(k, l, N_x, N_y, positions, "xx")
			M_xy = -3*C_3*sum_of_distances(k, l, N_x, N_y, positions, "xy")
			M_yy = -3*C_3*sum_of_distances(k, l, N_x, N_y, positions, "yy")
			result_matrix_upper[row_index, col_index] =  [M_xx M_xy; M_xy M_yy]
		end
	end

	#making a big matrix from matrix of matrices
	result_matrix_upper = reduce(vcat, map(row -> reduce(hcat, row), eachrow(result_matrix_upper)))

	#return element-wise sum of matrices (upper part is symmetric to lower)
	return (result_matrix_diag .+ result_matrix_upper .+ transpose(result_matrix_upper))
end

"""
Function calculating final matrix (OBC 2D) for V3 - anisotropic term.
"""
function calculate_V3_matrix_part(N_x::Int, N_y::Int, M::Float64, Ω::Float64, C_3::Float64, m_vec::Array{Float64}, positions::Array{Float64})
	#creating matrix of 2x2 matrices of 0.0
	result_matrix_upper = [zeros(2, 2) for _ in 1:N_x*N_y, _ in 1:N_x*N_y]

	#filling the upper triangular part of a matrix
	for row_index in 1:N_x * N_y
		i, j = divrem(row_index - 1, N_x)
		i += 1
		j += 1
		for col_index in (row_index + 1):(N_x * N_y)
			k, l = divrem(col_index - 1, N_x)
			k += 1
			l += 1
			M_xx = (3/2)*C_3*sum_of_distances_V3(k, l, N_x, N_y, m_vect, positions, "xx")
			M_xy = (3/2)*C_3*sum_of_distances_V3(k, l, N_x, N_y, m_vect, positions, "xy")
			M_yy = (3/2)*C_3*sum_of_distances_V3(k, l, N_x, N_y, m_vect, positions, "yy")
			result_matrix_upper[row_index, col_index] =  [M_xx M_xy; M_xy M_yy]
		end
	end

	#making a big matrix from matrix of matrices
	result_matrix_upper = reduce(vcat, map(row -> reduce(hcat, row), eachrow(result_matrix_upper)))

	#return element-wise sum of matrices (upper part is symmetric to lower)
	return (result_matrix_upper .+ transpose(result_matrix_upper))
end

"""
Function calculating final matrix (OBC).

Arguments
- N_x, N_y: lattice grid dimensions.
- M, Ω, C_3: physical parameters used in diagonal/off-diagonal blocks.
- m_vec - [m_x, m_y] - unit vector 
- positions: N_y × N_x × 2 array of Cartesian coordinates for every lattice site.
"""
function calculate_matrix_OBC_2D(N_x::Int, N_y::Int, M::Float64, Ω::Float64, C_3::Float64, positions::Array{Float64})
	#check if m is normalized
	if !(length(v) == 2 && isapprox(norm(v), 1.0; atol=atol))
		error("Vector m is nor normalized!");
	end

	matrix_V1_V2 = calculate_V1_V2_matrix_part(N_x, N_y, M, Ω, C_3, positions)
	matrix_V3 = calculate_V3_matrix_part(N_x, N_y, M, Ω, C_3, m_vec, positions)

	return (matrix_V1_V2 + matrix_V3)
end