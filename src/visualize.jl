"""
	visualize_phonon_mode_from_positions(result_matrix, positions, mode_index)

Arguments
- result_matrix: assembled dynamical/stiffness matrix (2N × 2N).
- positions: N_x × N_y × 2 array of Cartesian coordinates used to place atoms.
- mode_index: index of the eigenmode to visualize.
"""
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
# function visualize_phonon_mode_from_positions(result_matrix::Matrix{Float64}, positions::Array{Float64, 3}, mode_index::Int)
#	 N_x, N_y, _ = size(positions)
#	 N = N_x * N_y

#	 # Flatten positions into a 1D array of Point2f
#	 flat_positions = [Point2f(positions[i, j, 1], positions[i, j, 2]) for j in 1:N_y, i in 1:N_x]
#	 flat_positions = vec(flat_positions)  # i runs fastest

#	 # Eigen-decomposition
#	 eig = eigen(result_matrix)
#	 eigvals = real.(eig.values)  # Eigenvalues
#	 eigvec = real.(eig.vectors[:, mode_index])

#	 # Extract (dx, dy) for each atom
#	 displacements = [Vec2f(eigvec[2i - 1], eigvec[2i]) for i in 1:N]

#	 # Normalize displacements for visualization
#	 max_norm = maximum(norm.(displacements))
#	 disp_scaled = [d / max_norm * 0.3f0 for d in displacements] # normalized for better visualization

#	 # Plot
#	 f = Figure()
	
#	 # Left plot: Density of States (DOS)
#	 ax_dos = Axis(f[1, 1]; xlabel="DOS", ylabel="Density", title="Density of States", aspect=DataAspect())
	
#	 # Compute DOS by creating a histogram of eigenvalues
#	 hist = histogram(eigvals, bins=30, density=true, color=:blue, alpha=0.6)
	
#	 # Right plot: Phonon mode visualization
#	 ax_mode = Axis(f[1, 2]; xlabel="x", ylabel="y", title="Phonon Mode", aspect=DataAspect())
#	 CairoMakie.scatter!(ax_mode, flat_positions; markersize=10, color=:black)
#	 arrows!(ax_mode, flat_positions, disp_scaled; arrowsize=10, linewidth=2, color=:blue)
	
#	 # Adjust layout: Make the figure wide enough to show both plots side by side
#	 f.layout = (1, 2)
	
#	 return f
# end
# function visualize_phonon_mode_from_positions(
#		 result_matrix::Matrix{Float64},
#		 positions::Array{Float64,3},
#		 mode_index::Int
#	 )

#	 N_x, N_y, _ = size(positions)
#	 N = N_x * N_y

#	 # Flatten positions into Point2f
#	 flat_positions = [Point2f(positions[i, j, 1], positions[i, j, 2])
#					   for j in 1:N_y, i in 1:N_x]
#	 flat_positions = vec(flat_positions)

#	 # Eigen-decomposition
#	 eig = eigen(result_matrix)
#	 eigvals = real.(eig.values)
#	 eigvec  = real.(eig.vectors[:, mode_index])

#	 # Extract displacements
#	 displacements = [Vec2f(eigvec[2i - 1], eigvec[2i]) for i in 1:N]

#	 # Normalize displacements
#	 max_norm = maximum(norm.(displacements))
#	 disp_scaled = [d / max_norm * 0.3f0 for d in displacements]

#	 # -----------------------------
#	 # Create a 1×2 layout Figure
#	 # -----------------------------
#	 f = Figure(size = (950, 450))

#	 # ===============================
#	 # Left: DOS histogram
#	 # ===============================
#	 ax_dos = Axis(f[1, 1], xlabel="Density of States", ylabel="ω")

#	 hist!(
#		 ax_dos,
#		 eigvals;
#		 bins=30,
#		 normalization=:pdf,
#		 color=RGBAf0(0.2, 0.4, 0.8, 0.6)  # <- includes alpha!
#	 )

#	 # Flip DOS horizontally so it sits "on the left"
#	 ax_dos.xreversed[] = true

#	 # ===============================
#	 # Right: phonon mode visualization
#	 # ===============================
#	 ax_mode = Axis(f[1, 2], xlabel="x", ylabel="y", aspect=DataAspect())

#	 scatter!(ax_mode, flat_positions; markersize=10, color=:black)
#	 arrows!(ax_mode, flat_positions, disp_scaled;
#			 arrowsize=10, linewidth=2, color=:blue)

#	 return f
# end