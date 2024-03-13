# # # import numpy as np
# # # import matplotlib.pyplot as plt

# # # # Example 2D matrix (replace this with your own data)
# # # matrix = ['AAAATAAA', 'AAAATAAA', 'AAAATAAA', 'AAAAAAAA', 'AAAAGCAA', 'AAAATAAA', 'AAAATCAA', 'AAAATAAA', 'AAAATAAA', 'AAAATCAA'] 


# # # # Plotting the matrix as an image
# # plt.imshow(matrix, cmap='gray', interpolation='nearest')
# # plt.title('Consensus Image')
# # plt.colorbar(label='Intensity')
# # plt.show()
# import numpy as np
# import matplotlib.pyplot as plt

# def plot_motif_logo(motif_matrix, figsize=(10, 4), title="Motif Logo"):
#     """
#     Plot a motif logo from a motif matrix.
    
#     Parameters:
#         motif_matrix (list of lists): The motif matrix.
#         figsize (tuple, optional): The size of the figure (width, height) in inches. Default is (10, 4).
#         title (str, optional): The title of the plot. Default is "Motif Logo".
#     """
#     # Convert motif matrix to numpy array
#     motif_matrix = np.array(motif_matrix)
    
#     # Calculate the information content for each position
#     information_content = np.sum(-motif_matrix * np.log2(motif_matrix + 1e-10), axis=1)
    
#     # Plot the motif logo
#     fig, ax = plt.subplots(figsize=figsize)
#     ax.set_title(title)
#     ax.set_ylabel("Information Content")
#     ax.set_xlabel("Position")
    
#     for i in range(len(motif_matrix)):
#         for j in range(4):
#             if motif_matrix[i, j] > 0:
#                 ax.text(i+0.5, information_content[i], 
#                         ['A', 'C', 'G', 'T'][j], 
#                         color='blue', fontsize=12, ha='center', va='bottom')
    
#     ax.set_xticks(np.arange(len(motif_matrix)) + 0.5)
#     ax.set_xticklabels(np.arange(1, len(motif_matrix)+1))
#     ax.set_xlim(0, len(motif_matrix))
#     ax.set_ylim(0, max(information_content) + 1)
#     ax.tick_params(axis='x', length=0)
#     ax.spines['top'].set_visible(False)
#     ax.spines['right'].set_visible(False)
#     ax.spines['bottom'].set_visible(False)
#     ax.spines['left'].set_visible(False)
#     ax.yaxis.set_ticks_position('left')
#     ax.xaxis.set_ticks_position('bottom')
#     ax.yaxis.set_tick_params(width=0)
#     plt.show()

# # # Example motif matrix (replace this with your own data)
# # motif_matrix = [
# #     [0.1, 0.2, 0.3, 0.4],
# #     [0.4, 0.3, 0.2, 0.1],
# #     [0.25, 0.25, 0.25, 0.25],
# #     [0.35, 0.15, 0.25, 0.25]
# # ]

# # # Plot the motif logo
# # plot_motif_logo(motif_matrix, title="Motif Logo")

# def build_motif_matrix(matrix):
#     """
#     Build a motif matrix from a list of sequences.
    
#     Parameters:
#         matrix (list of str): The list of sequences.
    
#     Returns:
#         motif_matrix (list of lists): The motif matrix.
#     """
#     motif_matrix = []
#     for i in range(len(matrix[0])):
#         counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
#         for seq in matrix:
#             counts[seq[i]] += 1
#         motif_matrix.append([counts['A']/len(matrix),
#                              counts['C']/len(matrix),
#                              counts['G']/len(matrix),
#                              counts['T']/len(matrix)])
#     return motif_matrix

# # Input matrix
# matrix = ['AAAATAAA', 'AAAATAAA', 'AAAATAAA', 'AAAAAAAA', 'AAAAGCAA',
#           'AAAATAAA', 'AAAATCAA', 'AAAATAAA', 'AAAATAAA', 'AAAATCAA']

# # Build motif matrix
# motif_matrix = build_motif_matrix(matrix)
# print(motif_matrix)

# # Plot motif logo
# plot_motif_logo(motif_matrix, title="Motif Logo")
import numpy as np
import logomaker

# Input motif matrix
# matrix = ['AAAATAAA', 'AAAATAAA', 'AAAATAAA', 'AAAAAAAA', 'AAAAGCAA',
#           'AAAATAAA', 'AAAATCAA', 'AAAATAAA', 'AAAATAAA', 'AAAATCAA']
matrix = ["CGGGCCCC","TGGGCCGC","CGGGCCGC","CCGGCCGC","CGGGCCGC","CGGGCCCC",
             "CGGGCCGC","CGGGCCGC","TGGGCCCC","CGGGCCCC"]

# Create a counts matrix
counts_matrix = logomaker.alignment_to_matrix(matrix)

# Convert counts to probabilities
print(counts_matrix)
# prob_matrix = counts_matrix / counts_matrix.sum(axis=1)
row_sums = counts_matrix.sum(axis=0)
prob_matrix = counts_matrix / row_sums

# Create a logo object
logo = logomaker.Logo(counts_matrix, color_scheme='classic')

# Customize the appearance (optional)
logo.ax.set_title('Motif Logo')
logo.style_spines(visible=False)  # Hide spines
logo.style_xticks(fmt='%d', anchor = 0,spacing=50)  # Adjust x-axis ticks
logo.ax.set_ylabel("Probability")

# Show the logo
logo.draw()
print(prob_matrix)

