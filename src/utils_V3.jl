# function mxx(m::Array{Float64}, )

"""
	sum_of_distances(i_index, j_index, N_x, N_y, positions, direction)

Arguments
- i_index, j_index: indices of the site being considered.
- N_x, N_y: lattice dimensions.
- positions: array of Cartesian coordinates (N_y × N_x × 2).
- direction: "xx", "xy"/"yx", or "yy" specifying component type.
"""
function sum_of_distances_V3(i_index::Int, j_index::Int, N_x::Int, N_y::Int, m::Array{Float64}, positions::Array{Float64}, direction::String)
	result = 0
	for x_pos in 1:N_y
		for y_pos in 1:N_x
			if (i_index != x_pos && j_index != y_pos)
				dx = positions[i_index, j_index, 1] - positions[x_pos, y_pos, 1]
				dy = positions[i_index, j_index, 2] - positions[x_pos, y_pos, 2]
				dist = (sqrt(dx^2 + dy^2))
				m_x = m[1]
				m_y = m[2]
				if direction == "xx"
					result += dist^(-7) * (
							-20 * ((m_x*dx + m_y*dy)/dist) * m_x*dx
							-5 * ((m_x*dx + m_y*dy)/dist)^2
						)
						+ dist^(-9) * (
							35 * ((m_x*dx + m_y*dy)/dist)^2 * dx^2
						)
						+ dist^(-5) * (
							2 * m_x^2
						)
				elseif direction == "xy" || direction == "yx"
					result += dist^(-7) * (
							-10 * ((m_x*dx + m_y*dy)/dist) * (m_x*dx + m_y*dy)
						)
						+ dist^(-9) * (
							35 * ((m_x*dx + m_y*dy)/dist)^2 * dx * dy
						)
						+ dist^(-5) * m_x * m_y
				elseif direction == "yy"
					result += dist^(-7) * (
							-20 * ((m_x*dx + m_y*dy)/dist) * m_y*dy
							-5 * ((m_x*dx + m_y*dy)/dist)^2
						)
						+ dist^(-9) * (
							35 * ((m_x*dx + m_y*dy)/dist)^2 * dy^2
						)
						+ dist^(-5) * (
							2 * m_y^2
						)
				end
			end
		end
	end
	return result
end

