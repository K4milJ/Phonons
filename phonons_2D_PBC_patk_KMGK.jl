################################################################################
# 2D system - size: N_x*N_y, PBC
################################################################################

using Plots
using LinearAlgebra

function shortest_distance(atom1::Tuple{Int, Int}, 
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
	distance = sqrt(dist_x^2 + dist_y^2)
	return distance
end
#function calculates the shortest distance between two atoms in 2D lattice with PBC

function sum_of_distances_diag_elem(a, M, Ω, C_3, i_index, j_index, N_x, N_y, rand_matrix)
    result = 0
    for x_pos in 1:N_x
        for y_pos in 1:N_y
            if (i_index != x_pos && j_index != y_pos)
				atom1 = (x_pos, y_pos)
				atom2 = (i_index, j_index)
                result += shortest_distance(atom1, atom2, a, N_x, N_y)^-5
							#(sqrt(((i_index-1)*a + rand_matrix[j_index, i_index, 1] - (x_pos-1)*a - rand_matrix[y_pos, x_pos, 1])^2 + 
                            #   ((j_index-1)*a + rand_matrix[j_index, i_index, 2] - (y_pos-1)*a - rand_matrix[y_pos, x_pos, 2])^2))^-5
            end
        end
    end
    return result
end
#function calculates the sum of |distance|^-5 according to the equation

function matrix_diag_elem(a, M, Ω, C_3, i_index, j_index, N_x, N_y, rand_matrix)
    return M*Ω^2 + 12*C_3*sum_of_distances_diag_elem(a, M, Ω, C_3, i_index, j_index, N_x, N_y, rand_matrix)
end
#function calculates the diagonal element

function calculate_matrix2d(a, N_x, N_y, M, Ω, C_3, rand_rng)
	rand_matrix = 2 * rand_rng * rand(Float16, N_y, N_x, 2) .- rand_rng
    result_matrix_diag = zeros(Float64, N_x*N_y, N_x*N_y)
    diag_ind = 1
    for i in 1:N_x
        for j in 1:N_y
            #println("diagind=", diag_ind, "x=", i, "y=", j)
            result_matrix_diag[diag_ind,diag_ind] = matrix_diag_elem(a, M, Ω, C_3, i, j, N_x, N_y, rand_matrix)
            diag_ind += 1
        end
    end
    result_matrix_upper = zeros(N_x * N_y, N_x * N_y)
    # Iterate through all pairs of indices
    for row_index in 1:N_x * N_y
        # Convert row index to (i, j) in the original matrix
        i, j = divrem(row_index - 1, N_y)
        i += 1
        j += 1
        for col_index in row_index+1:N_x * N_y
            # Convert column index to (k, l) in the original matrix
            k, l = divrem(col_index - 1, N_y)
            k += 1
            l += 1
			atom1 = (i, j)
			atom2 = (k, l)
            dist = shortest_distance(atom1, atom2, a, N_x, N_y)
			#sqrt(
            #    ((i-1)*a + rand_matrix[i, j, 1] - (k-1)*a - rand_matrix[k, l, 1])^2
            #    +
            #    ((j-1)*a + rand_matrix[i, j, 2] - (l-1)*a - rand_matrix[k, l, 2])^2
            #)
            result_matrix_upper[row_index, col_index] = dist^-5
        end
    end
    return result_matrix_diag .+ result_matrix_upper .+ transpose(result_matrix_upper)
end

################################################################################

N_x = 8
N_y = 10
a = 1.0
M = 1
Ω = 1
C_3 = 1
#rand_rng = 0.0

result_matrix = calculate_matrix2d(a, N_x, N_y, M, Ω, C_3, 0.0)

scatter(eigvals(result_matrix), title="Eigenvalues: 2D, size: $N_x x $N_y, PBC", framestyle = :box)

################################################################################
#	N_x, N_y	- size of the system
#	a 			- lattice constant
#	M			- mass of each atom
#	Ω, C_3		- constants
#	rand_rng	- determines the range (δR_i ∈ [-rand_rng,rand_rng]) from which
#				we add random vector [δR_x, δR_y] to the position of each atom
################################################################################


################## Dynamical matrix - D(q) #####################################

#calculating the shortest distance in 1D periodic chain
function shortest_distance_1D(x1::Int, x2::Int, N::Int, a::Float64)
	L = N * a 
	x1 *= a
	x2 *= a
	dist = mod(x2 - x1 + L / 2, L) - L / 2
	return abs(dist)
end

#calculating the shortest distance in 2D periodic square lattice
function shortest_distance_2D(atom1::Tuple{Int, Int}, 
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
	distance = sqrt(dist_x^2 + dist_y^2)
	return distance
end

#shortest_distance_1D(1, 5, 5, 1.2)

#This calculates our
function D_xx_sum_elem(ξ, q_x, a_x, x_1, x_2)
	result = 2 * cos(ξ * q_x * a_x)
	#deriv = (-15 * C_3 * ξ) / (shortest_distance_1D(x1, x2, N_x, a_x)^(7/2)) + (6 * C_3) / (shortest_distance_1D(x1, x2, N_x, a_x)^(5/2))
	deriv = (-9 * C_3) / (shortest_distance_1D(x_1, x_2, N_x, a_x)^(5/2))
	return result
end

function main_sum_xx(q_x::Float64, N_x::Int, a_x::Float64)		#, atom1::Tuple{Int, Int}, atom2::Tuple{Int, Int})
	# N_x - system size in x-direction
	upper_sum_limit = ceil(N_x/2) - 1	#górny limit sumy
	result = 0
	for i in 1:upper_sum_limit
		ξ = i * a_x
		# x_1 = 1
		# x_2 = 1 + i
		result += D_xx_sum_elem(ξ, q_x, a_x, 1, Int(i+1))
	end
	return result
end

function diag_sum_elem_xx(a, N_x, N_y)
	sum = 0
	x_1 = 1	#i chose to be 1, because of PBC
	x_2 = 1
	for i in 2:N_x		#sum over all neighbors
		for j in 2:N_y
			sum += 15 * C_3 * shortest_distance_1D(i, x_1, N_x, a) / shortest_distance_2D((i, j), (x_1, x_2), a, N_x, N_y)^(7/2)
			sum += -6 * C_3 / shortest_distance_2D((i, j), (x_1, x_2), a, N_x, N_y)^(5/2)
		end
	end
	return sum
end

#########################
function D_yy_sum_elem(ξ, q_y, a, x_1, x_2)
	result = 2 * cos(ξ * q_y * a)
	#deriv = (-15 * C_3 * ξ) / (shortest_distance_1D(x1, x2, N_x, a_x)^(7/2)) + (6 * C_3) / (shortest_distance_1D(x1, x2, N_x, a_x)^(5/2))
	deriv = (-9 * C_3) / (shortest_distance_1D(x_1, x_2, N_y, a)^(5/2))
	return result
end

function main_sum_yy(q_y::Float64, N_y::Int, a::Float64)		#, atom1::Tuple{Int, Int}, atom2::Tuple{Int, Int})
	# N_x - system size in x-direction
	upper_sum_limit = ceil(N_y/2) - 1	#górny limit sumy
	result = 0
	for i in 1:upper_sum_limit
		ξ = i * a
		# x_1 = 1
		# x_2 = 1 + i
		result += D_xx_sum_elem(ξ, q_y, a, 1, Int(i + 1))
	end
	return result
end

function diag_sum_elem_yy(a, N_x, N_y)
	sum = 0
	x_1 = 1	#i chose to be 1, because of PBC
	x_2 = 1
	for i in 2:N_x		#sum over all neighbors
		for j in 2:N_y
			sum += 15 * C_3 * shortest_distance_1D(i, x_1, N_y, a) / shortest_distance_2D((i, j), (x_1, x_2), a, N_x, N_y)^(7/2)
			sum += -6 * C_3 / shortest_distance_2D((i, j), (x_1, x_2), a, N_x, N_y)^(5/2)
		end
	end
	return sum
end

q_x = 0:0.01:1.2
q_y = 1.0

function D_xx(q_x)
	D_xx_kk = M*Ω^2 + diag_sum_elem_xx(a, N_x, N_y)
	return (main_sum_xx(q_x, N_x, a) + D_xx_kk) / M
end

function D_yy(q_y)
	D_yy_kk = M*Ω^2 + diag_sum_elem_yy(a, N_x, N_y)
	return (main_sum_yy(q_y, N_y, a) + D_yy_kk) / M
end

#using Plots
#D_yy(q_y)
#plot(D_xx.(q_x))

### non-diagonal elements

function shortest_distance_1D_sign(x1::Int, x2::Int, N::Int, a::Float64)
	L = N * a 
	x1 *= a
	x2 *= a
	dist = mod(x2 - x1 + L / 2, L) - L / 2
	return dist
end


function sum_non_diag_xy(a, N_x, N_y, q_x, q_y)
	sum = 0
	x_1 = 1	#i chose to be 1, because of PBC
	x_2 = 1
	upper_sum_limit_x = ceil(N_x/2) - 1	#górny limit sumy
	upper_sum_limit_y = ceil(N_y/2) - 1	#górny limit sumy
	for i in 2:upper_sum_limit_x		#sum over all neighbors
		for j in 2:upper_sum_limit_y
			# sum += 15 * C_3 * shortest_distance_1D(i, x_1, N_y, a) / shortest_distance_2D((i, j), (x_1, x_2), a, N_x, N_y)^(7/2)
			# sum += -6 * C_3 / shortest_distance_2D((i, j), (x_1, x_2), a, N_x, N_y)^(5/2)
			ξ_x = shortest_distance_1D_sign(Int(i), x_1, N_y, a)
			ξ_y = shortest_distance_1D_sign(Int(j), x_1, N_y, a)
			sum += (shortest_distance_1D(Int(j), x_1, N_y, a) / shortest_distance_2D((Int(i), Int(j)), (x_1, x_2), a, N_x, N_y)^(7/2)) * 2 * cos(q_x*ξ_x + q_y * ξ_y)
		end
	end
	return sum
end

function D_xy(q_x, q_y)
	return (-15*C_3/M) * sum_non_diag_xy(a, N_x, N_y, q_x, q_y)
end

#

function sum_non_diag_yx(a, N_x, N_y, q_x, q_y)
	sum = 0
	x_1 = 1	#i chose to be 1, because of PBC
	x_2 = 1
	upper_sum_limit_x = ceil(N_x/2) - 1	#górny limit sumy
	upper_sum_limit_y = ceil(N_y/2) - 1	#górny limit sumy
	for i in 2:upper_sum_limit_x		#sum over all neighbors
		for j in 2:upper_sum_limit_y
			# sum += 15 * C_3 * shortest_distance_1D(i, x_1, N_y, a) / shortest_distance_2D((i, j), (x_1, x_2), a, N_x, N_y)^(7/2)
			# sum += -6 * C_3 / shortest_distance_2D((i, j), (x_1, x_2), a, N_x, N_y)^(5/2)
			ξ_x = shortest_distance_1D_sign(Int(i), x_1, N_y, a)
			ξ_y = shortest_distance_1D_sign(Int(j), x_1, N_y, a)
			sum += (shortest_distance_1D(Int(i), x_1, N_y, a) / shortest_distance_2D((Int(i), Int(j)), (x_1, x_2), a, N_x, N_y)^(7/2)) * 2 * cos(q_x*ξ_x + q_y * ξ_y)
		end
	end
	return sum
end

function D_yx(q_x, q_y)
	return (-15*C_3/M) * sum_non_diag_yx(a, N_x, N_y, q_x, q_y)
end


#D_xy(1.0, 1.2)
#D_yx(1.0, 1.2)

#####################################################
using LinearAlgebra

function generate_k_path(a::Real, nseg::Integer)
    # Define points
    K = [4π / (3a), 0.0]
    M = [π / a, π / (√3 * a)]
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

nseg = 1000
path = generate_k_path(a, nseg)

function matrix_D(q_x, q_y)
	mx = [D_xx(q_x) D_xy(q_x, q_y); D_xy(q_y, q_x) D_yy(q_y)]
	return eigvals(mx)
end

#plot(matrix_D.(path[1,:], path[2,:]))

eigvals_list = matrix_D.(path[1, :], path[2, :])
bands = hcat(eigvals_list...)  # Now each row is one band, columns are along k-path

# Compute distances along path
dq = [0.0; cumsum(sqrt.(sum(diff(path, dims=2).^2, dims=1))[:])]

# Plot each band
using Plots
plot(dq, bands', xlabel="k-path", ylabel="Eigenvalues", label="", lw=2)

#

xticks_pos = [dq[1], dq[nseg+1], dq[2nseg+1], dq[end]]
xticks_labels = ["K", "M", "Γ", "K"]

plot(dq, bands', xlabel="k-path", ylabel="Eigenvalues", label="", lw=2,
     xticks=(xticks_pos, xticks_labels), title="Band structure along high-symmetry path")

##################### END OF INCORRECT CODE #############
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

function shortest_distance_1D_sign(x1::Int, x2::Int, N::Int, a::Float64)
	L = N * a 
	x1 *= a
	x2 *= a
	dist = mod(x2 - x1 + L / 2, L) - L / 2
	return dist
end

function interaction_term(x, y, a, N_x, N_y, direction)
	r = shortest_distance_2D((1,1), (x,y), a, N_x, N_y)
	dx = shortest_distance_1D(1, x, N_x, a)		#unSIGNED!
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
	sum_lim_x = Int(ceil(N_x/2) - 1)	#górny limit sumy
	sum_lim_y = Int(ceil(N_y/2) - 1)	#górny limit sumy

	for x in 2:N_x, y in 2:N_y
		result += interaction_term(x, y, a, N_x, N_y, direction)
	end
	return result
end

function sum_over_entire_lattice_Rn(N_x::Int, N_y::Int, a::Float64, q_x::Float64, q_y::Float64, direction::String)
	#q_x, q_y = q #unpack coords
	result = 0.0 + 0im
	sum_lim_x = Int(ceil(N_x/2) - 1)	#górny limit sumy
	sum_lim_y = Int(ceil(N_y/2) - 1)	#górny limit sumy
	#println("sum_over_entire_lattice_Rn - ", direction)
	for x in 2:N_x, y in 2:N_y
		dx = shortest_distance_1D_sign(1, x, N_x - 1, a)
		dy = shortest_distance_1D_sign(1, y, N_y - 1, a)
		phase = exp(1im * (q_x * dx + q_y * dy))
		result += sum_over_entire_lattice(N_x, N_y, a, direction) * phase
		#term = interaction_term(x, y, a, N_x, N_y, direction)
		#result += term * phase
	end
	return result
end

##########

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
scatter(path[1, :], path[2, :], aspect_ratio=:equal)

eigvals_list = D_eigvals.(path[1, :], path[2, :])
imag_parts = imag.(eigvals_list)
max_imag = maximum(imag_parts)
min_imag = minimum(imag_parts)
eigvals_real = [real.(subvec) for subvec in eigvals_list]
eigvals_imag = [imag.(subvec) for subvec in eigvals_list]
eigvals_abs = [abs.(subvec) for subvec in eigvals_list]


bands = hcat(eigvals_real...)  # Now each row is one band, columns are along k-path
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


####################### TESTY ############################
#shortest_distance_1D_sign(7, 1, 7, 1.1)
qx = path[1, 12]
qy = path[2, 12]
#path
D = [D_xx(qx, qy) D_xy(qx, qy); D_xy(qx, qy) D_yy(qx, qy)]

D - D'
D - transpose(D)
## RESZTKI ##############################################

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