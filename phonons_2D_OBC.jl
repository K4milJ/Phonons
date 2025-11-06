"""
	distance_1D(x1::Int, x2::Int, N::Int, a::Float64)

Compute the scalar distance between two sites on a 1D open-boundary lattice.

Arguments
- x1, x2: integer site indices (1-based) along the 1D chain.
- N: total number of sites in the chain (keeps API symmetry; not otherwise required).
- a: lattice constant (spacing between neighbouring sites).

Returns
- Absolute distance (Float64) between positions x1 and x2 (units of the lattice spacing).

Notes
- Assumes open boundary conditions and direct scaling of indices by `a`.
"""
function distance_1D(x1::Int, x2::Int, N::Int, a::Float64)
	L = N * a
	x1 *= a
	x2 *= a
	return abs(x2 - x1)
end

"""
	distance_2D(atom1::Tuple{Int, Int}, atom2::Tuple{Int, Int}, a::Float64)

Compute the Euclidean distance between two lattice sites in a 2D rectangular lattice
with open boundary conditions.

Arguments
- atom1, atom2: tuples (i, j) of integer site indices (1-based) giving x and y grid positions.
- a: lattice constant (scales the grid indices to physical coordinates).

Returns
- Euclidean distance (Float64) between the two site positions.
"""
function distance_2D(atom1::Tuple{Int, Int}, 
	atom2::Tuple{Int, Int}, 
	a::Float64) 
	x1, y1 = atom1 .* a
	x2, y2 = atom2 .* a
	return sqrt((x2 - x1)^2 + (y2 - y1)^2)
end

"""
	sum_of_distances_diag_elem(i_index, j_index, N_x, N_y, a, rand_matrix, direction)

Compute the scalar part of a diagonal 2×2 dynamical matrix block due to interactions with all
other atoms. This function accumulates angular-dependent terms coming from the model's
distance-dependent potential (contains distance^-5 and distance^-3 contributions).

Arguments
- i_index, j_index: integer indices (grid coordinates) of the atom for which the diagonal element
  is being computed.
- N_x, N_y, a: lattice sizes and lattice constant.
- rand_matrix: an N_x × N_y × 2 array of small random offsets added to the ideal lattice positions.
- direction: string selector "xx", "xy"/"yx", or "yy" choosing which component to accumulate.

Returns
- Scalar Float64 value to be used when forming the corresponding diagonal 2×2 block.
"""
function sum_of_distances_diag_elem(i_index::Int, j_index::Int, N_x::Int, N_y::Int, a::Float64, rand_matrix::Array{Float64}, direction::String)
	result = 0
    for x_pos in 1:N_x
        for y_pos in 1:N_y
            if (i_index != x_pos && j_index != y_pos)
				dx = ((i_index-1)*a + rand_matrix[i_index, j_index, 1] - (x_pos-1)*a - rand_matrix[x_pos, y_pos, 1])
				dy = ((j_index-1)*a + rand_matrix[i_index, j_index, 2] - (y_pos-1)*a - rand_matrix[x_pos, y_pos, 2])
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
	matrix_diag_elem(M, Ω, C_3, i_index, j_index, N_x, N_y, a, rand_matrix)

Build the 2×2 diagonal block of the full dynamical/stiffness matrix for the atom at
grid coordinates (i_index, j_index).

Arguments
- M, Ω, C_3: physical constants used to form the diagonal contribution.
- i_index, j_index: integer lattice coordinates of the atom.
- N_x, N_y, a: lattice size and lattice constant.
- rand_matrix: random offsets for each atom (used to compute actual distances).

Returns
- 2×2 Float64 matrix corresponding to the diagonal block for the given atom.
"""
function matrix_diag_elem(M::Float64, Ω::Float64, C_3::Float64, i_index::Int, j_index::Int, N_x::Int, N_y::Int, a::Float64, rand_matrix::Array{Float64})
	M_xx = M*Ω^2 - 3*C_3*sum_of_distances_diag_elem(i_index, j_index, N_x, N_y, a, rand_matrix, "xx")
	M_xy = M*Ω^2 - 3*C_3*sum_of_distances_diag_elem(i_index, j_index, N_x, N_y, a, rand_matrix, "xy")
	M_yy = M*Ω^2 - 3*C_3*sum_of_distances_diag_elem(i_index, j_index, N_x, N_y, a, rand_matrix, "yy")
	return [M_xx M_xy; M_xy M_yy]
end

"""
	sum_of_distances(i_index, j_index, N_x, N_y, a, rand_matrix, direction)

Compute the off-diagonal (interaction) scalar contributions used when forming the
off-diagonal 2×2 blocks of the global dynamical matrix.

Arguments
- i_index, j_index: indices of the atom whose interaction with others is being computed.
- N_x, N_y, a: lattice size and lattice constant.
- rand_matrix: position offsets for each lattice site.
- direction: one of "xx", "xy"/"yx", or "yy" indicating which component to accumulate.

Returns
- Scalar Float64 containing the sum of angular/distance-dependent interaction terms
  for the requested direction.
"""
function sum_of_distances(i_index::Int, j_index::Int, N_x::Int, N_y::Int, a::Float64, rand_matrix::Array{Float64}, direction::String)
	result = 0
	for x_pos in 1:N_x
		for y_pos in 1:N_y
			if (i_index != x_pos && j_index != y_pos)
				dx = ((i_index-1)*a + rand_matrix[j_index, i_index, 1] - (x_pos-1)*a - rand_matrix[x_pos, y_pos, 1])
				dy = ((j_index-1)*a + rand_matrix[j_index, i_index, 2] - (y_pos-1)*a - rand_matrix[x_pos, y_pos, 2])
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

"""
	calculate_matrix2d(N_x, N_y, M, Ω, C_3, a, rand_rng)

Assemble the full (2N × 2N) dynamical/stiffness matrix for a square-lattice system with
open boundary conditions.

Arguments
- N_x, N_y: integer lattice grid dimensions (number of sites along x and y).
- M, Ω, C_3: physical parameters used to form diagonal and off-diagonal blocks.
- a: lattice constant.
- rand_rng: scalar controlling the amplitude of random displacements; 0 means perfect lattice.

Returns
- A (2*N_x*N_y) × (2*N_x*N_y) Float64 matrix representing the assembled dynamical matrix.

Notes
- The implementation builds 2×2 blocks for each pair of atoms and then constructs the
  large dense matrix by hcat/vcat reduction of the block structure.
"""
function calculate_matrix2d(N_x::Int, N_y::Int, M::Float64, Ω::Float64, C_3::Float64, a::Float64, rand_rng::Float64)
	rand_matrix = 2 * rand_rng * rand(Float16, N_x, N_y, 2) .- rand_rng
	result_matrix_diag = [zeros(2, 2) for _ in 1:N_x*N_y, _ in 1:N_x*N_y]
	#creating matrix of 2x2 matrices of 0.0
	diag_ind = 1
	for i in 1:N_x
		for j in 1:N_y
			result_matrix_diag[diag_ind, diag_ind] = matrix_diag_elem(M, Ω, C_3, i, j, N_x, N_y, a, rand_matrix)
			diag_ind += 1
		end
	end
	result_matrix_diag = reduce(vcat, map(row -> reduce(hcat, row), eachrow(result_matrix_diag)))
	#making a big matrix from matrix of matrices

	result_matrix_upper = [zeros(2, 2) for _ in 1:N_x*N_y, _ in 1:N_x*N_y]
	for row_index in 1:N_x * N_y
		i, j = divrem(row_index - 1, N_x)
		i += 1
		j += 1
		for col_index in (row_index + 1):(N_x * N_y)
			k, l = divrem(col_index - 1, N_x)
			k += 1
			l += 1
			dx = (i-1)*a + rand_matrix[j, i, 1] - (k-1)*a - rand_matrix[l, k, 1]
			dy = (j-1)*a + rand_matrix[j, i, 2] - (l-1)*a - rand_matrix[l, k, 2]
			dist = sqrt(dx^2 + dy^2)
			M_xx = -3*C_3*sum_of_distances(k, l, N_x, N_y, a, rand_matrix, "xx")
			M_xy = -3*C_3*sum_of_distances(k, l, N_x, N_y, a, rand_matrix, "xy")
			M_yy = -3*C_3*sum_of_distances(k, l, N_x, N_y, a, rand_matrix, "yy")
			result_matrix_upper[row_index, col_index] =  [M_xx M_xy; M_xy M_yy]
		end
	end
	result_matrix_upper = reduce(vcat, map(row -> reduce(hcat, row), eachrow(result_matrix_upper)))
	return result_matrix_diag .+ result_matrix_upper .+ transpose(result_matrix_upper)

end
################################################################################
#=
N_x = 4
N_y = 4
a = 1.1
M = 1.0
Ω = 1.0
C_3 = 1.5
rand_rng = 0.0
=#
################################################################
"""
	visualize_phonon_mode(result_matrix::Matrix{Float64}, N_x::Int, N_y::Int, a::Float64, mode_index::Int)

Create a visual representation of a single phonon eigenmode on a square lattice using CairoMakie.

Arguments
- result_matrix: assembled dynamical/stiffness matrix (2N × 2N) for the lattice.
- N_x, N_y: lattice dimensions used to generate the background atom positions.
- a: lattice constant used to place atoms in Cartesian coordinates.
- mode_index: index of the eigenmode (column in the eigenvector matrix) to visualize.

Returns
- A CairoMakie Figure containing scatter points for atom positions and arrows for displacement vectors.

Notes
- This normalizes mode amplitudes for display and assumes eigenvectors are ordered consistently with
  the 2-component-per-site layout (dx, dy) interleaved.
"""
function visualize_phonon_mode(result_matrix::Matrix{Float64}, N_x::Int, N_y::Int, a::Float64, mode_index::Int)
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

################################################################################
#	N_x, N_y	- size of the system
#	a 			- lattice constant
#	M			- mass of each atom
#	Ω, C_3		- constants
#	rand_rng	- determines the range (δR_i ∈ [-rand_rng,rand_rng]) from which
#				we add random vector [δR_x, δR_y] to the position of each atom
################################################################################

#result_matrix = calculate_matrix2d(N_x, N_y, M, Ω, C_3, a, rand_rng)
result_matrix = calculate_matrix2d(4, 4, 1.0, 1.0, 1.5, 1.1, 0.0)
eigvals(result_matrix)
ev = Plots.scatter(eigvals(result_matrix), framestyle = :box, legend=:none)
# save("plots/phonon_mode_1_square_OBC_4x4_eigvals.png", ev) 
#evec = eigvecs(result_matrix)

#visualize_phonon_mode(result_matrix::Matrix{Float64}, N_x::Int, N_y::Int, a::Float32, mode_index::Int)
f = visualize_phonon_mode(result_matrix, 4, 4, 1.1, 1)
save("plots/phonon_mode_1_square_OBC_4x4.png", f)  # Optional: Save the figure