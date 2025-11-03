using Plots
using LinearAlgebra
using CairoMakie

#calculating the distance in 1D lattice OBC
function distance_1D(x1::Int, x2::Int, N::Int, a::Float64)
	L = N * a
	x1 *= a
	x2 *= a
	return abs(x2 - x1)
end

#calculating the distance in 2D lattice OBC
function distance_2D(atom1::Tuple{Int, Int}, 
	atom2::Tuple{Int, Int}, 
	a::Float64, N_x::Int, N_y::Int) 
	x1, y1 = atom1 .* a
	x2, y2 = atom2 .* a
	return sqrt((x2 - x1)^2 + (y2 - y1)^2)
end

function sum_of_distances_diag_elem(M::Float64, Ω::Float64, C_3::Float64, i_index::Int, j_index::Int, N_x::Int, N_y::Int, positions::Array{Float64}, direction::String)
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

function matrix_diag_elem(M::Float64, Ω::Float64, C_3::Float64, i_index::Int, j_index::Int, N_x::Int, N_y::Int, positions::Array{Float64})
    M_xx = M*Ω^2 - 3*C_3*sum_of_distances_diag_elem(M, Ω, C_3, i_index, j_index, N_x, N_y, positions, "xx")
	M_xy = M*Ω^2 - 3*C_3*sum_of_distances_diag_elem(M, Ω, C_3, i_index, j_index, N_x, N_y, positions, "xy")
	M_yy = M*Ω^2 - 3*C_3*sum_of_distances_diag_elem(M, Ω, C_3, i_index, j_index, N_x, N_y, positions, "yy")
	return [M_xx M_xy; M_xy M_yy]
end

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

function calculate_matrix2d_hex(N_x::Int, N_y::Int, M::Float64, Ω::Float64, C_3::Float64, positions::Array{Float64})
	#rand_matrix = 2 * rand_rng * rand(Float16, N_x, N_y, 2) .- rand_rng
    result_matrix_diag = [zeros(2, 2) for _ in 1:N_x*N_y, _ in 1:N_x*N_y] #zeros(Float64, N_x*N_y, N_x*N_y)
	#creating matrix of 2x2 matrices of 0.0
    diag_ind = 1
    for i in 1:N_x
        for j in 1:N_y
            result_matrix_diag[diag_ind, diag_ind] = matrix_diag_elem(M, Ω, C_3, i, j, N_x, N_y, positions)
            diag_ind += 1
        end
    end
	result_matrix_diag = reduce(vcat, map(row -> reduce(hcat, row), eachrow(result_matrix_diag)))
	#making a big matrix from matrix of matrices
	#return result_matrix_diag #DEL

    result_matrix_upper = [zeros(2, 2) for _ in 1:N_x*N_y, _ in 1:N_x*N_y] #zeros(Float64, N_x*N_y, N_x*N_y)
    for row_index in 1:N_x * N_y
        i, j = divrem(row_index - 1, N_x)
        i += 1
        j += 1
        for col_index in (row_index + 1):(N_x * N_y)
            k, l = divrem(col_index - 1, N_x)
            k += 1
            l += 1
			dx = positions[i, j, 1] - positions[k, l, 1]
			dy = positions[i, j, 2] - positions[k, l, 2]
            dist = sqrt(dx^2 + dy^2)
            #result_matrix_upper[row_index, col_index] = dist^(-5)
			M_xx = -3*C_3*sum_of_distances(k, l, N_x, N_y, positions, "xx")
			M_xy = -3*C_3*sum_of_distances(k, l, N_x, N_y, positions, "xy")
			M_yy = -3*C_3*sum_of_distances(k, l, N_x, N_y, positions, "yy")
			# M_xx = -3 * (dx^2 * dist^(-5) + dist^(-3))
			# M_xy = -3 * (dx*dy * dist^(-5) + dist^(-3))
			# M_yy = -3 * (dy^2 * dist^(-5) + dist^(-3))
			result_matrix_upper[row_index, col_index] =  [M_xx M_xy; M_xy M_yy]
			#println("[col, row] = [$col_index, $row_index]\t[i,j] = [$i,$j]\t[k,l]=[$k,$l]")
        end
    end
	result_matrix_upper = reduce(vcat, map(row -> reduce(hcat, row), eachrow(result_matrix_upper)))
	#display(result_matrix_upper)
    return result_matrix_diag .+ result_matrix_upper .+ transpose(result_matrix_upper)
end

#########################################
function calculate_positions_hex_PBC(N_x::Int, N_y::Int, a::Float64)
	## temp ##
	#if (mod(N_x,2) != 0 || mod(N_y,2) != 0)
	#	println("ERROR - currently this shape of system is not available");
	#else
		#println("OK");
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
	#end
end

#=
function visualize_phonon_mode_hex(result_matrix::Matrix{Float64}, N_x::Int, N_y::Int, a::Float64, mode_index::Int)
	N = N_x*N_y

	# Generate lattice positions (2D square)
	positions = [Point2f(i*a, j*a) for j in 1:N_x, i in 1:N_y]
	positions = vec(positions)  # flatten to 1D array

	eig = eigen(result_matrix)
	eigvec = real.(eig.vectors[:, mode_index])
	eigvec
	# Extract displacements (dx, dy) per atom
	displacements = [Vec2f(eigvec[2i-1], eigvec[2i]) for i in 1:N]

	# Normalize displacements for visualization
	max_norm = maximum(norm.(displacements))
	disp_scaled = [d ./ max_norm * 0.3f0 for d in displacements]

	# plot
	f = Figure();
	ax = Axis(f[1, 1]; aspect=DataAspect(), xlabel="x", ylabel="y");
	# Plot atoms
	CairoMakie.scatter!(ax, positions; markersize=10, color=:black);
	# Plot displacement arrows
	arrows!(ax, positions, disp_scaled; arrowsize=10, linewidth=2, color=:blue);
	#title!(ax, "Phonon Mode $mode_index on Square Lattice");
	return f
end
=#

function visualize_phonon_mode_from_positions(result_matrix::Matrix{Float64}, positions::Array{Float64, 3}, mode_index::Int)
    N_x, N_y, _ = size(positions)
    N = N_x * N_y

    # Flatten positions into a 1D array of Point2f
    flat_positions = [Point2f(positions[i, j, 1], positions[i, j, 2]) for j in 1:N_y, i in 1:N_x]
    flat_positions = vec(flat_positions)  # i runs fastest

    # Eigen-decomposition
    eig = eigen(result_matrix)
    eigvec = real.(eig.vectors[:, mode_index])

    # Extract (dx, dy) for each atom
    displacements = [Vec2f(eigvec[2i - 1], eigvec[2i]) for i in 1:N]

    # Normalize displacements for visualization
    max_norm = maximum(norm.(displacements))
    disp_scaled = [d / max_norm * 0.3f0 for d in displacements] #normalized for better visualization
    # disp_scaled = [d for d in displacements] #no normalization

    # Plot
    f = Figure()
    ax = Axis(f[1, 1]; aspect=DataAspect(), xlabel="x", ylabel="y")
    CairoMakie.scatter!(ax, flat_positions; markersize=10, color=:black)
    arrows!(ax, flat_positions, disp_scaled; arrowsize=10, linewidth=2, color=:blue)
    return f
end

################################################################################
###  TESTS  ###
test_pos = calculate_positions_hex_PBC(8, 10, 1.0)

result_matrix = calculate_matrix2d_hex(8, 10, 1.0, 1.0, 1.0, test_pos)

f = visualize_phonon_mode_from_positions(result_matrix, test_pos, 1)
# save("plots/phonon_mode_1_hex_big.png", f)  # Optional: Save the figure

ev = Plots.scatter(eigvals(result_matrix), framestyle = :box, legend=:none)
# save("plots/phonon_mode_1_hex_big_eigvals.png", ev) 
################################################################################

# Here is rozgrzebany code for visualizing lattice - it doesn't work yet

#Let's plot our lattice to check if everything is fine
x_comp = vec(test_pos[:,:,1])
y_comp = vec(test_pos[:,:,2])
lattice_plt = Plots.scatter(x_comp, y_comp, title="", 
	aspect_ratio=1, framestyle=:box, legend=:none, color=:blue, yflip=true)

function lattice_display!(positions, plt, N_x::Int, N_y::Int, a::Float64) #to be finished
	for i in 1:N_x, j in 1:N_y-1
		println("[x,y]=[$i,$j]");
		Plots.plot!([positions[j, i, 1], positions[j+1, i, 1]], [positions[j, i, 2], positions[j+1, i, 2]]; color=:gray, lw=1, label="")
	end
	#=
	for i in 1:N_x-1
		if (mod(i,2) == 0)	#even rows
			plot!([positions[:, i, 1], positions[:, i+1, 1]], [positions[:, i, 2], positions[:, i+1, 2]]; color=:gray, lw=1, label="")
		else	#odd rows
			plot!([positions[:, i, 1], positions[:, i+1, 1]], [positions[:, i, 2], positions[:, i+1, 2]]; color=:gray, lw=1, label="")
		end
	end
	=#
	return plt
end

lattice_plt = lattice_display!(test_pos, lattice_plt, 4, 6, 1.0)
display(lattice_plt)
#lattice display, not finished
################################################################################
#now we can implement our structure from square lattice, only distances will be
#calculated from our new matrix of positions
