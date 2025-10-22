#########################################################
#														#
#		##    ## ##     ## ####### ##    ##				#
#		##   ##  ###   ### ##      ##   ## 				#
#		##  ##   #### #### ##      ##  ##  				#
#		#####    ## ### ## ##      #####   				#
#		##  ##   ##     ## ##      ##  ##  				#
#		##   ##  ##     ## ##      ##   ## 				#
#		##    ## ##     ## ##      ##    ##   path		#
#														#
############### D(q) along the KMΓK path ################

#calculating the shortest distance in 1D periodic chain
function shortest_distance_1D(x1::Int, x2::Int, N::Int, a_x::Float64)
	L = N * a_x
	x1 *= a_x
	x2 *= a_x
	dist = mod(x2 - x1 + L / 2, L) - L / 2
	return abs(dist)
end

#calculating the shortest distance in 2D periodic square lattice
function shortest_distance_2D(atom1::Tuple{Int, Int}, 
	atom2::Tuple{Int, Int}, 
	a::Float64, N_x::Int, N_y::Int) 
	x1, y1 = atom1 .* a
	x2, y2 = atom2 .* a
	Lx, Ly = N_x*a, N_y*a 
	dist_x = mod(x2 - x1 + Lx / 2, Lx) - Lx / 2
	dist_y = mod(y2 - y1 + Ly / 2, Ly) - Ly / 2
	distance = sqrt(dist_x^2 + dist_y^2)
	return distance
end

#calculating the shortest distance in 1D periodic chain - signed
function shortest_distance_1D_sign(x1::Int, x2::Int, N::Int, a::Float64)
	L = N * a 
	x1 *= a
	x2 *= a
	dist = mod(x2 - x1 + L / 2, L) - L / 2
	return dist
end

function interaction_term(x, y, a, N_x, N_y, direction)
	r = shortest_distance_2D((1,1), (x,y), a, N_x, N_y)
	dx = shortest_distance_1D(1, x, N_x, a)
	dy = shortest_distance_1D(1, y, N_y, a)
	if direction == "xx"
		return -3 * r^(-5) * dx^2 + r^(-3)
	elseif direction == "yy"
		return -3 * r^(-5) * dy^2 + r^(-3)
	elseif direction == "xy" || direction == "yx"
		return -3 * r^(-5) * dx * dy + r^(-3)
	end
end

function sum_over_entire_lattice(N_x::Int, N_y::Int, a::Float64, direction::String)
	result = 0.0
	sum_lim_x = Int(ceil(N_x/2) - 1)
	sum_lim_y = Int(ceil(N_y/2) - 1)

	for x in 2:N_x, y in 2:N_y
		result += interaction_term(x, y, a, N_x, N_y, direction)
	end
	return result
end

function sum_over_entire_lattice_Rn(N_x::Int, N_y::Int, a::Float64, q_x::Float64, q_y::Float64, direction::String)
	result = 0.0 + 0im
	sum_lim_x = Int(ceil(N_x/2) - 1)	#upper sum limit
	sum_lim_y = Int(ceil(N_y/2) - 1)	#lower sum limit
	for x in 2:N_x, y in 2:N_y
		dx = shortest_distance_1D_sign(1, x, N_x - 1, a)
		dy = shortest_distance_1D_sign(1, y, N_y - 1, a)
		phase = exp(1im * (q_x * dx + q_y * dy))
		result += sum_over_entire_lattice(N_x, N_y, a, direction) * phase
	end
	return result
end

function generate_k_path(a::Real, nseg::Integer)
    # Define points
    K = [π / a, 0.0]
    M = [π / a, π / a]
    Γ = [0.0, 0.0]

    # Helper for linear interpolation between two points
    function segment(p1, p2, n)
        return [p1 .* (1 - t) + p2 .* t for t in range(0, 1, length=n+1)]
    end

    # Segments
    path_KM = segment(K, M, nseg)
    path_MG = segment(M, Γ, nseg)
    path_GK = segment(Γ, K, nseg)

    # Concatenate (avoid duplicate points at junctions)
    path = vcat(path_KM[1:end-1], path_MG[1:end-1], path_GK)

    # Convert to 2xN matrix
    return hcat(path...)
end

### matrix elements ###

function D_xx(q_x, q_y)
	return		(1/M) * (M*Ω^2
				- 6 * C_3 * sum_over_entire_lattice(N_x, N_y, a, "xx")
				+ 3 * C_3 * sum_over_entire_lattice_Rn(N_x, N_y, a, q_x, q_y, "xx"))
end

function D_yy(q_x, q_y)
	return		(1/M) * (M*Ω^2
				- 6 * C_3 * sum_over_entire_lattice(N_x, N_y, a, "yy")
				+ 3 * C_3 * sum_over_entire_lattice_Rn(N_x, N_y, a, q_x, q_y, "yy"))
end

function D_xy(q_x, q_y)
	return		(1/M) * (- 6 * C_3 * sum_over_entire_lattice(N_x, N_y, a, "xy")
				+ 3 * C_3 * sum_over_entire_lattice_Rn(N_x, N_y, a, q_x, q_y, "xy"))
end

function D_eigvals(q_x, q_y) 
	return eigvals([D_xx(q_x, q_y) D_xy(q_x, q_y); D_xy(q_x, q_y) D_yy(q_x, q_y)])
end

##########

M = 1
Ω = 2*π
N_x = 23
N_y = 32
a = 1.0
C_3 = 1.2

##########

using LinearAlgebra
using Plots

#####

nseg = 50
path = generate_k_path(a, nseg)
#checking path
scatter(path[1, :], path[2, :], aspect_ratio=:equal)

#checking results, Re, Im and Abs
eigvals_list = D_eigvals.(path[1, :], path[2, :])
imag_parts = imag.(eigvals_list)
max_imag = maximum(imag_parts)
min_imag = minimum(imag_parts)
eigvals_real = [real.(subvec) for subvec in eigvals_list]
eigvals_imag = [imag.(subvec) for subvec in eigvals_list]
eigvals_abs = [abs.(subvec) for subvec in eigvals_list]


bands = hcat(eigvals_real...)
bandsi = hcat(eigvals_imag...)
bandabs = hcat(eigvals_abs...)
# Compute distances along path
dq = [0.0; cumsum(sqrt.(sum(diff(path, dims=2).^2, dims=1))[:])]


xticks_pos = [dq[1], dq[nseg+1], dq[2nseg+1], dq[end]]
xticks_labels = ["K", "M", "Γ", "K"]

plot(dq, bands', xlabel="k-path", ylabel="Eigenvalues Re()", label="", lw=2,
     xticks=(xticks_pos, xticks_labels), title="Band structure along high-symmetry path", color=:blue)

plot(dq, bandsi', xlabel="k-path", ylabel="Eigenvalues Im()", label="", lw=2,
     xticks=(xticks_pos, xticks_labels), title="Band structure along high-symmetry path", color=:orange)

plot(dq, bandabs', xlabel="k-path", ylabel="Eigenvalues abs()", label="", lw=2,
     xticks=(xticks_pos, xticks_labels), title="Band structure along high-symmetry path", color=:green)

#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
# Periodic Boundary Conditions (PBC)

using Plots
using LinearAlgebra
using CairoMakie

function shortest_distance_1D_PBC(x1::Int, x2::Int, N::Int, a_x::Float64)
	L = N * a_x
	x1 *= a_x
	x2 *= a_x
	dist = mod(x2 - x1 + L / 2, L) - L / 2
	return abs(dist)
end

#calculating the shortest distance in 2D periodic square lattice
function shortest_distance_2D_PBC(atom1::Tuple{Int, Int}, 
	atom2::Tuple{Int, Int}, 
	a::Float64, N_x::Int, N_y::Int) 
	x1, y1 = atom1 .* a
	x2, y2 = atom2 .* a
	Lx, Ly = N_x*a, N_y*a 
	dist_x = mod(x2 - x1 + Lx / 2, Lx) - Lx / 2
	dist_y = mod(y2 - y1 + Ly / 2, Ly) - Ly / 2
	distance = sqrt(dist_x^2 + dist_y^2)
	return distance
end

function sum_of_distances_diag_elem(M::Float64, Ω::Float64, C_3::Float64, i_index::Int, j_index::Int, N_x::Int, N_y::Int, a::Float64, direction::String)
    result = 0
    for x_pos in 1:N_x
        for y_pos in 1:N_y
            if (i_index != x_pos && j_index != y_pos)
				dx = shortest_distance_1D_PBC(i_index, x_pos, N_x, a)
				#((i_index-1)*a - (x_pos-1)*a)
				dy = shortest_distance_1D_PBC(j_index, y_pos, N_y, a)
				#((j_index-1)*a - (y_pos-1)*a)
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

function matrix_diag_elem(M::Float64, Ω::Float64, C_3::Float64, i_index::Int, j_index::Int, N_x::Int, N_y::Int, a::Float64)
    M_xx = M*Ω^2 - 3*C_3*sum_of_distances_diag_elem(M, Ω, C_3, i_index, j_index, N_x, N_y, a, "xx")
	M_xy = M*Ω^2 - 3*C_3*sum_of_distances_diag_elem(M, Ω, C_3, i_index, j_index, N_x, N_y, a, "xy")
	M_yy = M*Ω^2 - 3*C_3*sum_of_distances_diag_elem(M, Ω, C_3, i_index, j_index, N_x, N_y, a, "yy")
	return [M_xx M_xy; M_xy M_yy]
end

function sum_of_distances(M::Float64, Ω::Float64, C_3::Float64, i_index::Int, j_index::Int, N_x::Int, N_y::Int, a::Float64, direction::String)
    result = 0
    for x_pos in 1:N_x
        for y_pos in 1:N_y
            if (i_index != x_pos && j_index != y_pos)
				dx = shortest_distance_1D_PBC(i_index, x_pos, N_x, a)
				#((i_index-1)*a - (x_pos-1)*a)
				dy = shortest_distance_1D_PBC(j_index, y_pos, N_y, a)
				#((j_index-1)*a - (y_pos-1)*a)
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

function calculate_matrix2d_PBC(N_x::Int, N_y::Int, M::Float64, Ω::Float64, C_3::Float64, a::Float64)
    result_matrix_diag = [zeros(2, 2) for _ in 1:N_x*N_y, _ in 1:N_x*N_y] #zeros(Float64, N_x*N_y, N_x*N_y)
	#creating matrix of 2x2 matrices of 0.0
    diag_ind = 1
    for i in 1:N_x
        for j in 1:N_y
            result_matrix_diag[diag_ind, diag_ind] = matrix_diag_elem(M, Ω, C_3, i, j, N_x, N_y, a)
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
			dx = (i-1)*a - (k-1)*a
			dy = (j-1)*a - (l-1)*a
            dist = sqrt(dx^2 + dy^2)
            #result_matrix_upper[row_index, col_index] = dist^(-5)
			M_xx = -3*C_3*sum_of_distances(M, Ω, C_3, k, l, N_x, N_y, a, "xx")
			M_xy = -3*C_3*sum_of_distances(M, Ω, C_3, k, l, N_x, N_y, a, "xy")
			M_yy = -3*C_3*sum_of_distances(M, Ω, C_3, k, l, N_x, N_y, a, "yy")
			#M_xx = -3 * (dx^2 * dist^(-5) + dist^(-3))
			#M_xy = -3 * (dx*dy * dist^(-5) + dist^(-3))
			#M_yy = -3 * (dy^2 * dist^(-5) + dist^(-3))
			result_matrix_upper[row_index, col_index] =  [M_xx M_xy; M_xy M_yy]
			#println("[col, row] = [$col_index, $row_index]\t[i,j] = [$i,$j]\t[k,l]=[$k,$l]")
        end
    end
	result_matrix_upper = reduce(vcat, map(row -> reduce(hcat, row), eachrow(result_matrix_upper)))
	#display(result_matrix_upper)
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

##########################################################
#result_matrix = calculate_matrix2d(N_x, N_y, M, Ω, C_3, a)
result_matrix_PBC = calculate_matrix2d_PBC(4, 4, 1.0, 1.0, 1.5, 1.1)
eigvals(result_matrix_PBC)
ev = Plots.scatter(eigvals(result_matrix_PBC), legend=:none, framestyle = :box)#, title="Eigenvalues: 2D, size: $N_x x $N_y, δR=$rand_rng, OBC", framestyle = :box)
save("plots/phonon_mode_1_square_PBC_4x4_eigvals.png", ev)
#evec = eigvecs(result_matrix)

#visualize_phonon_mode(result_matrix::Matrix{Float64}, N_x::Int, N_y::Int, a::Float32, mode_index::Int)
f = visualize_phonon_mode(result_matrix_PBC, 10, 10, 1.1, 1)
save("plots/phonon_mode_1_square_PBC_10x10.png", f)  # Optional: Save the figure


#################
eig = eigen(result_matrix_PBC)
	eigvec = real.(eig.vectors[:, 1])
####################### TESTS ############################
#=
#shortest_distance_1D_sign(7, 1, 7, 1.1)
qx = path[1, 12]
qy = path[2, 12]
#path
D = [D_xx(qx, qy) D_xy(qx, qy); D_xy(qx, qy) D_yy(qx, qy)]

## REST OF THE CODE ##############################################

# for D_xx
function sum_over_entire_lattice_xx(N_x::Int, N_y::Int, a::Float64)
	result = 0
	sum_lim_x = Int(ceil(N_x/2) - 1)	#górny limit sumy
	sum_lim_y = Int(ceil(N_y/2) - 1)	#górny limit sumy
	for x in 2:sum_lim_x, y in 2:sum_lim_y
		result += -3 * shortest_distance_2D((1,1), (x,y), a, N_x, N_y)^(-5) * shortest_distance_1D(1, x, N_x, a)^2 + shortest_distance_2D((1,1), (x,y), a, N_x, N_y)^(-3)
	end
	return result
end

# for D_yy
function sum_over_entire_lattice_yy(N_x::Int, N_y::Int, a::Float64)
	result = 0
	sum_lim_x = Int(ceil(N_x/2) - 1)	#górny limit sumy
	sum_lim_y = Int(ceil(N_y/2) - 1)	#górny limit sumy
	for x in 2:sum_lim_x, y in 2:sum_lim_y
		result += -3 * shortest_distance_2D((1,1), (x,y), a, N_x, N_y)^(-5) * shortest_distance_1D(1, y, N_y, a)^2 + shortest_distance_2D((1,1), (x,y), a, N_x, N_y)^(-3)
	end
	return result
end

function sum_over_entire_lattice_xy(N_x::Int, N_y::Int, a::Float64)
	result = 0
	sum_lim_x = Int(ceil(N_x/2) - 1)	#górny limit sumy
	sum_lim_y = Int(ceil(N_y/2) - 1)	#górny limit sumy
	for x in 2:sum_lim_x, y in 2:sum_lim_y
		result += -3 * shortest_distance_2D((1,1), (x,y), a, N_x, N_y)^(-5) * shortest_distance_1D(1, y, N_y, a) * shortest_distance_1D(1, x, N_x, a) + shortest_distance_2D((1,1), (x,y), a, N_x, N_y)^(-3)
	end
	return result
end

function sum_over_entire_lattice_Rn_xy(N_x::Int, N_y::Int, a::Float64, q_x::Float64, q_y::Float64)
	#q_x, q_y = q #unpack coords
	result = 0
	sum_lim_x = Int(ceil(N_x/2) - 1)	#górny limit sumy
	sum_lim_y = Int(ceil(N_y/2) - 1)	#górny limit sumy
	for x in 2:sum_lim_x, y in 2:sum_lim_y
		result += sum_over_entire_lattice_xy(N_x, N_y, a) * exp(1im * (q_x * shortest_distance_1D_sign(1, x, N_x, a) + q_y * shortest_distance_1D_sign(1, y, N_y, a)))
					#assuming function 'shortest_distance_1D_sign()' is correct
	end
	return result
end

function sum_over_entire_lattice_Rn_xx(N_x::Int, N_y::Int, a::Float64, q_x::Float64, q_y::Float64)
	#q_x, q_y = q #unpack coords
	sum_lim_x = Int(ceil(N_x/2) - 1)	#górny limit sumy
	sum_lim_y = Int(ceil(N_y/2) - 1)	#górny limit sumy
	result = 0
	for x in 2:sum_lim_x, y in 2:sum_lim_y
		result += sum_over_entire_lattice_xx(N_x, N_y, a) * exp(1im * (q_x * shortest_distance_1D_sign(1, x, N_x, a) + q_y * shortest_distance_1D_sign(1, y, N_y, a)))
					#assuming function 'shortest_distance_1D_sign()' is correct
	end
	return result
end

function sum_over_entire_lattice_Rn_yy(N_x::Int, N_y::Int, a::Float64, q_x::Float64, q_y::Float64)
	#q_x, q_y = q #unpack coords
	sum_lim_x = Int(ceil(N_x/2) - 1)	#górny limit sumy
	sum_lim_y = Int(ceil(N_y/2) - 1)	#górny limit sumy
	result = 0
	for x in 2:sum_lim_x, y in 2:sum_lim_y
		result += sum_over_entire_lattice_xx(N_x, N_y, a) * exp(1im * (q_x * shortest_distance_1D_sign(1, x, N_x, a) + q_y * shortest_distance_1D_sign(1, y, N_y, a)))
					#assuming function 'shortest_distance_1D_sign()' is correct
	end
	return result
end


###################3
function shortest_distance_vec(atom1::Tuple{Int, Int}, 
	atom2::Tuple{Int, Int}, 
	a::Float64, N_x::Int, N_y::Int) 
	# Unpack coordinates
	x1, y1 = atom1 .* a
	x2, y2 = atom2 .* a

	Lx, Ly = N_x*a, N_y*a 

	# Calculate the minimum distance considering periodic boundary conditions
	dist_x = mod(x2 - x1 + Lx / 2, Lx) - Lx / 2
	dist_y = mod(y2 - y1 + Ly / 2, Ly) - Ly / 2
	#println("dx = ", dist_x)
	#println("dy = ", dist_y)
	# Calculate the Euclidean distance
	#distance = sqrt(dist_x^2 + dist_y^2)
	return [dist_x, dist_y]
end

for i in 1:N_x
	for j in 1:N_y
		#println("diagind=", diag_ind, "x=", i, "y=", j)
		result_matrix_diag[diag_ind,diag_ind] = matrix_diag_elem(a, M, Ω, C_3, i, j, N_x, N_y, rand_matrix)
		diag_ind += 1
	end
end

function sum_of_interactions_with_lattice_PBC(a, M, Ω, C_3, current_i, current_j, N_x, N_y) #ta wewnętrzna suma
	result = 0
	for x_pos in 1:N_x
		for y_pos in 1:N_y
			if (current_i != x_pos && current_j != y_pos)
				atom1 = (x_pos, y_pos)
				atom2 = (current_i, current_j)
				result += shortest_distance(atom1, atom2, a, N_x, N_y)^-5
			end
		end
	end
	return result
end

function sum_over_all_atoms_PBC(a, N_x, N_y, M, Ω, C_3) #ta zewnętrzna suma
	result = 0
	for i in 1:N_x
		for j in 1:N_y
			R = shortest_distance_vec([i, j], [current_i, current_j], a, N_x, N_y)
			result .+= M*Ω^2 + 12*C_3 * sum_of_interactions_with_lattice_PBC(a, M, Ω, C_3, current_i, current_j, N_x, N_y) * exp(1im.*g.*R)
		end
	end
	return result
end

D_xx = 1/M * (sum_over_all_atoms_PBC(a, N_x, N_y, M, Ω, C_3)) #ta zewnętrzna suma
D_yy = 1/M * (sum_over_all_atoms_PBC(a, N_x, N_y, M, Ω, C_3))

D_xy = 1/M * ()

function sum_of_interactions_with_lattice_PBC_nondiag(a, M, Ω, C_3, current_i, current_j, N_x, N_y) #ta wewnętrzna suma
	result = 0
	for x_pos in 1:N_x
		for y_pos in 1:N_y
			if (current_i != x_pos && current_j != y_pos)
				atom1 = (x_pos, y_pos)
				atom2 = (current_i, current_j)
				result += shortest_distance(atom1, atom2, a, N_x, N_y)^-5
			end
		end
	end
	return result
end

=#