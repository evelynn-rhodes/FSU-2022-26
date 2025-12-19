import argparse
import numpy as np
from phylip import Phylip
from tree import Tree, TreeNode

class AncestralReconstructionML:
    def __init__(self, phylip_file, matrix_file):
        self.phylip_file = phylip_file
        self.matrix = self.load_matrix(matrix_file)
        self.phylip_parser = Phylip(phylip_file)
        self.tree = Tree()
        self.amino_acids = list(self.matrix.keys())
        sequences = self.phylip_parser.read_phylip()
        self.stationary_freq = self.calculate_stationary_frequencies(sequences)

    def load_matrix(self, file_path):
        """Load substitution matrix from file."""
        matrix = {}
        with open(file_path, 'r') as file:
            lines = file.readlines()
            headers = lines[0].split()
            for line in lines[1:]:
                row = line.split()
                amino_acid = row[0]
                scores = list(map(float, row[1:]))
                matrix[amino_acid] = dict(zip(headers, scores))
        return matrix

    def calculate_stationary_frequencies(self, sequences):
        """Calculate amino acid frequencies from input sequences."""
        counts = {aa: 0 for aa in self.amino_acids}
        total = 0
        for seq in sequences.values():
            for aa in seq:
                if aa in counts:
                    counts[aa] += 1
                    total += 1
        return {aa: count / total for aa, count in counts.items()}

    def compute_distance(self, seq1, seq2):
        """Compute distance using substitution matrix."""
        if seq1 is None or seq2 is None:
            return float('inf')
        total_score = 0
        for a, b in zip(seq1, seq2):
            if a in self.matrix and b in self.matrix[a]:
                total_score += self.matrix[a][b]
            else:
                total_score += -8  # Penalty for undefined pairs or gaps
        return abs(total_score)

    def construct_tree_stepwise(self, sequences):
        """Construct a tree step-by-step using a distance matrix."""
        taxa = list(sequences.keys())
        clusters = {taxon: TreeNode(name=taxon, sequence=sequences[taxon]) for taxon in taxa}
        dist_matrix = self.compute_pairwise_distances(sequences)

        while len(clusters) > 1:
            # Find the closest pair of clusters
            closest_pair = None
            min_dist = float('inf')
            for taxon1 in clusters:
                for taxon2 in clusters:
                    if taxon1 != taxon2 and dist_matrix[taxon1][taxon2] < min_dist:
                        min_dist = dist_matrix[taxon1][taxon2]
                        closest_pair = (taxon1, taxon2)

            if closest_pair is None:
                raise ValueError("No valid clusters to merge. Check the distance matrix for inconsistencies.")

            taxon1, taxon2 = closest_pair

            # Merge the closest clusters
            new_name = f"{taxon1}_{taxon2}"
            new_node = TreeNode(name=new_name)
            new_node.add_child(clusters[taxon1])
            new_node.add_child(clusters[taxon2])

            # Set branch lengths
            new_node.children[0].branch_length = min_dist / 2
            new_node.children[1].branch_length = min_dist / 2

            # Update the distance matrix
            dist_matrix[new_name] = {}
            for other_taxon in clusters:
                if other_taxon not in (taxon1, taxon2):
                    dist1 = dist_matrix[taxon1][other_taxon]
                    dist2 = dist_matrix[taxon2][other_taxon]
                    dist_matrix[new_name][other_taxon] = (dist1 + dist2) / 2
                    dist_matrix[other_taxon][new_name] = (dist1 + dist2) / 2

            # Remove old taxa from the matrix and cluster
            del dist_matrix[taxon1]
            del dist_matrix[taxon2]
            for taxon in dist_matrix:
                dist_matrix[taxon].pop(taxon1, None)
                dist_matrix[taxon].pop(taxon2, None)

            del clusters[taxon1]
            del clusters[taxon2]
            clusters[new_name] = new_node

        self.tree.root = next(iter(clusters.values()))

    def compute_pairwise_distances(self, sequences):
        """Compute pairwise distances between sequences using substitution matrix."""
        taxa = list(sequences.keys())
        dist_matrix = {taxon: {} for taxon in taxa}
        for i, taxon1 in enumerate(taxa):
            for j, taxon2 in enumerate(taxa):
                if i < j:
                    seq1 = sequences[taxon1]
                    seq2 = sequences[taxon2]
                    distance = self.compute_distance(seq1, seq2)
                    dist_matrix[taxon1][taxon2] = distance
                    dist_matrix[taxon2][taxon1] = distance
        return dist_matrix

    def reconstruct(self, num_trees=5):
        """Construct multiple trees for robustness testing."""
        sequences = self.phylip_parser.read_phylip()
        self.construct_tree_stepwise(sequences)
        likelihood = self._prune(self.tree.root)
        total_likelihood = self._compute_total_likelihood(likelihood)
        return self.tree.to_newick(), total_likelihood

    def _prune(self, node):
        """Recursive likelihood calculation with scaling."""
        if node.is_leaf():
            likelihoods = np.zeros((len(node.sequence), len(self.amino_acids)))
            for i, observed in enumerate(node.sequence):
                for j, aa in enumerate(self.amino_acids):
                    likelihoods[i][j] = 1.0 if observed == aa else 0.0
            return likelihoods

        child_likelihoods = [self._prune(child) for child in node.children]
        likelihoods = np.ones((len(child_likelihoods[0]), len(self.amino_acids)))

        for site in range(len(likelihoods)):
            for i, aa in enumerate(self.amino_acids):
                product = 1.0
                for child, child_like in zip(node.children, child_likelihoods):
                    branch_length = child.branch_length
                    transition_probs = np.exp(
                        np.array([self.matrix[aa][child_aa] for child_aa in self.amino_acids]) * branch_length
                    )
                    product *= np.dot(child_like[site], transition_probs)
                likelihoods[site][i] = product

            max_likelihood = np.max(likelihoods[site])
            likelihoods[site] /= max_likelihood

        node.sequence = ''.join(
            self.amino_acids[np.argmax(likelihoods[site])] for site in range(len(likelihoods))
        )
        return likelihoods


    def _compute_total_likelihood(self, likelihoods):
        """Compute the total likelihood at the root."""
        total_likelihood = 0.0
        epsilon = 1e-10  # Small value to prevent log(0)
        for i, aa in enumerate(self.amino_acids):
             total_likelihood += self.stationary_freq[aa] * max(likelihoods[:, i].prod(), epsilon)
        return total_likelihood


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ancestral sequence reconstruction using ML and substitution matrix.")
    parser.add_argument("phylip_file", type=str, help="Path to the PHYLIP file containing sequences.")
    parser.add_argument("matrix_file", type=str, help="Path to the substitution matrix file (e.g., PAM250).")

    args = parser.parse_args()

    reconstructor = AncestralReconstructionML(args.phylip_file, args.matrix_file)
    newick, likelihood = reconstructor.reconstruct()
    print("\nAncestral Protein Sequences:")
    def traverse(node):
        if node.sequence:
            print(f"{node.name}: {node.sequence}")
        for child in node.children:
            traverse(child)
    traverse(reconstructor.tree.root)
    print("\nTree in Newick Format:")
    print(newick)
    print("\nTotal Log Likelihood:")
    print(np.log(likelihood))