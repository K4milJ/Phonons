"""
	calculate_positions_hex_PBC(N_x::Int, N_y::Int, a::Float64)

Arguments
- N_x, N_y: integer grid dimensions (number of sites along x and y used to arrange points).
- a: lattice constant; vertical spacing uses 0.5 * a * √3 for hex packing.
"""
function calculate_positions_hex_PBC(N_x::Int, N_y::Int, a::Float64)
	positions = zeros(N_y,N_x,2)
	# Y - component:
	for i in 1:N_x, j in 1:N_y
		positions[j, i, 2] = (j-1) * 0.5 * a * √3
	end
	# X - component:
	for j in 1:N_y
		if (mod(j,2) == 0)	#even rows
			#positions[j, :, 1] = [0.5; accumulate(+, repeat([1, 2], ceil(Int, N_x/2))[1:N_x-1]) .+ 0.5]
			positions[j, :, 1] = [0; accumulate(+, repeat([2, 1], ceil(Int, N_x/2))[1:N_x-1]) .+ 0]
		else				#odd rows
			#positions[j, :, 1] = [0; accumulate(+, repeat([2, 1], ceil(Int, N_x/2))[1:N_x-1]) .+ 0]
			positions[j, :, 1] = [0.5; accumulate(+, repeat([1, 2], ceil(Int, N_x/2))[1:N_x-1]) .+ 0.5]
		end
	end
return positions
end