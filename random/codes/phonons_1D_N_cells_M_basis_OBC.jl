################################################################################
# 1D system - N cells with M atoms each, OBC, +δR
################################################################################

using Plots
using LinearAlgebra

function calculate_positions_of_atoms(N_cells::Int, N_basis::Int, interatomic_dist::Vector{Float64})
    position_vec = zeros(Float64, N_basis*N_cells)  #vector of zeros of size N
    for i in 2:N_basis*N_cells
        position_vec[i] = position_vec[i-1] + interatomic_dist[mod1(i-1, N_basis)]
    end
    return position_vec
end

function calculate_random_displacement(N_cells::Int, N_basis::Int, rand_rng::Float64)
    return 2 * rand_rng * rand(Float16, N_basis*N_cells) .- rand_rng 
end

function sum_of_distances_diag_elem(N_cells::Int, N_basis::Int, i_index::Int, rand_rng::Float64, rand_vec::Vector{Float64})
    result = 0
    position_vec = calculate_positions_of_atoms(N_cells, N_basis, interatomic_dist)
    #rand_vec = calculate_random_displacement(N_cells, N_basis, rand_rng)
    for i in 1:N_basis*N_cells
        if (i != i_index) #1/0=Inf
            result += abs(position_vec[i_index] + rand_vec[i_index] - position_vec[i] - rand_vec[i])^-5
        end
    end
    return result
end

function diagonal_element(N_cells::Int, N_basis::Int, M::Float64, Ω::Float64, C_3::Float64, i_index::Int, rand_rng::Float64, rand_vec::Vector{Float64})
    return M*Ω^2 + 12*C_3*sum_of_distances_diag_elem(N_cells, N_basis, i_index, rand_rng, rand_vec)    
end

function calculating_matrix(N_cells::Int, N_basis::Int, M::Float64, Ω::Float64, C_3::Float64, rand_rng::Float64, interatomic_dist::Vector{Float64})
    position_vec = calculate_positions_of_atoms(N_cells, N_basis, interatomic_dist)
    rand_vec = calculate_random_displacement(N_cells, N_basis, rand_rng)
    result_matrix_diagonal = zeros(Float64, N_basis*N_cells, N_basis*N_cells) #initiating matrix NxN
    for i in 1:N_basis*N_cells
        result_matrix_diagonal[i,i] = diagonal_element(N_cells, N_basis, M, Ω, C_3, i, rand_rng, rand_vec)
    end 
    result_matrix_upper_part = zeros(Float64, N_basis*N_cells, N_basis*N_cells) #initializing another matrix
    for i in 1:N_basis*N_cells
        for j in i+1:N_basis*N_cells
            result_matrix_upper_part[i,j] = -12*C_3*abs(position_vec[i] + rand_vec[i] - position_vec[j] - rand_vec[j])^-5
        end
    end
    return result_matrix_upper_part .+ transpose(result_matrix_upper_part) .+ result_matrix_diagonal
end

################################################################################

N_cells = 5
N_basis = 2
interatomic_dist = [1.1, 1.6]#, 1.4]#, 1.0]#, 2.0, 3.0, 2.1, 1.9, 1.7, 2.5]
M = 1.0
Ω = 1.0
C_3 = 1.0
rand_rng = 0.2

result_matrix = calculating_matrix(N_cells, N_basis, M, Ω, C_3, rand_rng, interatomic_dist)

scatter(eigvals(result_matrix), title="Eigenvalues: $N_cells cells, $N_basis-ion basis, δR=$rand_rng", framestyle = :box)

eigvecs(result_matrix)

##########

# Example: 1D chain of N atoms
N = 10
positions = collect(1:N)  # atoms at positions 1, 2, ..., N

# Simulate a dummy dynamical matrix (symmetric, real)
D = randn(N, N)
D = 0.5 * (D + D')  # Make it Hermitian for real eigenvalues

# Compute eigenvalues and eigenvectors
evals = eigvals(D)
evecs = eigvecs(D)

# Pick one mode to visualize (e.g. mode 2)
mode_index = 2
mode_vec = evecs[:, mode_index]

# Normalize for display
disp_magnitudes = real.(mode_vec) ./ maximum(abs.(real.(mode_vec))) * 0.4

# Plot atoms and arrows
scatter(positions, zeros(N); label="Atoms", legend=false, markersize=6, xlabel="Atom index", ylabel="", aspect_ratio=1)
quiver!(positions, zeros(N), quiver=(disp_magnitudes, zeros(N)), arrow=true, linecolor=:red, lw=2)

# Optional: add title
title!("Phonon mode $mode_index")



#######################################
## 2D case - test - random positions ##
#######################################
using LinearAlgebra
using Plots

# Number of atoms
N = 6

# Atom positions in 2D (could be a lattice)
positions = [rand(2) .* 10 for _ in 1:N]  # random positions in 2D
X = [p[1] for p in positions]
Y = [p[2] for p in positions]

# Dynamical matrix: 2N x 2N (real symmetric)
D = randn(2N, 2N)
D = 0.5 * (D + D')  # make symmetric

# Compute eigen decomposition
eig = eigen(D)
evecs = eig.vectors

# Select mode to plot
mode_index = 4
mode_vec = real.(evecs[:, mode_index])

# Split eigenvector into (ux, uy) for each atom
displacements = [mode_vec[2i-1:2i] for i in 1:N]

# Normalize displacements for better visuals
max_disp = maximum(norm.(displacements))
scaled_disp = [d ./ max_disp * 0.6 for d in displacements]

# Create dx, dy arrays
DX = [d[1] for d in scaled_disp]
DY = [d[2] for d in scaled_disp]

# Plot atoms and displacement arrows
scatter(X, Y; label="Atoms", markersize=6, legend=false, aspect_ratio=1, xlabel="x", ylabel="y")
quiver!(X, Y, quiver=(DX, DY), arrow=true, lw=2, color=:red)

title!("Phonon mode $mode_index (2D)")





#######################################
## 2D case                           ##
#######################################
using LinearAlgebra
using CairoMakie

# ========== Lattice Parameters ==========
L = 10                    # lattice size (LxL)
a = 1.0                  # lattice spacing
N = L^2                  # total atoms

# Generate lattice positions (2D square)
positions = [Point2f(i*a, j*a) for j in 0:L-1, i in 0:L-1]
positions = vec(positions)  # flatten to 1D array

# ========== Dynamical Matrix ==========
D = randn(2N, 2N)
D = 0.5 * (D + D')  # Symmetrize

eig = eigen(D)
mode_index = 3
eigvec = real.(eig.vectors[:, mode_index])

# Extract displacements (dx, dy) per atom
displacements = [Vec2f(eigvec[2i-1], eigvec[2i]) for i in 1:N]

# Normalize displacements for visualization
max_norm = maximum(norm.(displacements))
disp_scaled = [d ./ max_norm * 0.3f0 for d in displacements]

# ========== Plot ==========
f = Figure()
ax = Axis(f[1, 1]; aspect=DataAspect(), xlabel="x", ylabel="y")

# Plot atoms
CairoMakie.scatter!(ax, positions; markersize=10, color=:black);

# Plot displacement arrows
arrows!(ax, positions, disp_scaled; arrowsize=10, linewidth=2, color=:blue);

CairoMakie.title!(ax, "Phonon Mode $mode_index on Square Lattice");

f


#########################################################################################################
#	variable			Type				Description													#
#-------------------------------------------------------------------------------------------------------#
#	N_cells				Int					- number of unit cells										#
#	N_basis				Int					- number of atoms in each cells								#
#	interatomic_dist	Vector{Float64}		- vector of lattice constants between atoms,				#
#											Ai - ith atom,												#
#											A1-A2 - lattice constant between atom A1 and A2,			#
#											interatomic_dist = [A1-A2, A2-A3, A3-A4, ..., An-A1]		#
#	M					Float64				- mass of atoms (for now all of them have the same mass)	#
#	Ω, C_3				Float64				- constants													#
#	rand_rng			Float64				- determines the range [-rand_rng,rand_rng] from which		#
#											we add random value to the position of each atom			#
#########################################################################################################
