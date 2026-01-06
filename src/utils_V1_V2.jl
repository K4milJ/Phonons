"""
	sum_of_distances_diag_elem(i_index, j_index, N_x, N_y, positions, direction)

Arguments
- i_index, j_index: integer indices selecting the atom whose diagonal block is computed.
- N_x, N_y: lattice dimensions.
- positions: N_y × N_x × 2 array with precomputed Cartesian positions for each site.
- direction: component selector: "xx", "xy"/"yx", or "yy".
"""
function sum_of_distances_diag_elem(i_index::Int, j_index::Int, N_x::Int, N_y::Int, positions::Array{Float64}, direction::String)
	result = 0
	for x_pos in 1:N_y
		for y_pos in 1:N_x
			if (j_index != x_pos && i_index != y_pos)
				dx = positions[j_index, i_index, 1] - positions[x_pos, y_pos, 1] #i_index - x-comp of chosen atom for which we calculate int with neighbors
				dy = positions[j_index, i_index, 2] - positions[x_pos, y_pos, 2] #j_index - y-comp of chosen atom for which we calculate int with neighbors
				distance = (sqrt(dx^2 + dy^2))
				if direction == "xx"
					result += -6 * dx^2 * distance^(-5) + distance^(-3)
				elseif direction == "xy" || direction == "yx"
					result += -6 * abs(dx * dy) * distance^(-5) + distance^(-3)
				elseif direction == "yy"
					result += -6 * dy^2 * distance^(-5) + distance^(-3)
				end
			end
		end
	end

	return result
end


"""
	matrix_diag_elem(M, Ω, C_3, i_index, j_index, N_x, N_y, positions)

Arguments
- M, Ω, C_3: physical constants (mass, frequency parameter, coupling constant).
- i_index, j_index: grid indices of the atom.
- N_x, N_y: grid dimensions.
- positions: N_y × N_x × 2 array of Cartesian coordinates.
"""
function matrix_diag_elem(M::Float64, Ω::Float64, C_3::Float64, i_index::Int, j_index::Int, N_x::Int, N_y::Int, positions::Array{Float64})
	M_xx = M*Ω^2 - 3*C_3*sum_of_distances_diag_elem(i_index, j_index, N_x, N_y, positions, "xx")
	M_xy = M*Ω^2 - 3*C_3*sum_of_distances_diag_elem(i_index, j_index, N_x, N_y, positions, "xy")
	M_yy = M*Ω^2 - 3*C_3*sum_of_distances_diag_elem(i_index, j_index, N_x, N_y, positions, "yy")
	return [M_xx M_xy; M_xy M_yy] 
end
#tutaj dodać symetryczne elementy z (d^2/dx^2)(V_3)

"""
	sum_of_distances(i_index, j_index, N_x, N_y, positions, direction)

Arguments
- i_index, j_index: indices of the site being considered.
- N_x, N_y: lattice dimensions.
- positions: array of Cartesian coordinates (N_y × N_x × 2).
- direction: "xx", "xy"/"yx", or "yy" specifying component type.
"""
function sum_of_distances(i_index::Int, j_index::Int, N_x::Int, N_y::Int, positions::Array{Float64}, direction::String)
	result = 0
	for x_pos in 1:N_y
		for y_pos in 1:N_x
			if (i_index != x_pos && j_index != y_pos)
				dx = positions[i_index, j_index, 1] - positions[x_pos, y_pos, 1]
				dy = positions[i_index, j_index, 2] - positions[x_pos, y_pos, 2]
				distance = (sqrt(dx^2 + dy^2))
				if direction == "xx"
					result += 3 * dx^2 * distance^(-5) + distance^(-3)
				elseif direction == "xy" || direction == "yx"
					result += 3 * abs(dx * dy) * distance^(-5) + distance^(-3)
				elseif direction == "yy"
					result += 3 * dy^2 * distance^(-5) + distance^(-3)
				end
			end
		end
	end
	return result
end

