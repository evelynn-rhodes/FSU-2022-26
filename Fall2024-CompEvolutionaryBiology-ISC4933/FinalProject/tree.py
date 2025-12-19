import numpy as np
import sys

class TreeNode:
    def __init__(self, name: str = '', parent=None, branch_length: float = 0.0, sequence=None):
        self.name = name
        self.parent = parent
        self.children = []
        self.branch_length = branch_length
        self.sequence = sequence


    def add_child(self, child_node):
        """
        Add a child node to this node.

        Parameters:
        -----------
        child_node : TreeNode
            The child node to add.
        """
        child_node.parent = self
        self.children.append(child_node)

    def is_tip(self):
        return self.is_leaf()
    
    def is_leaf(self):
        """
        Check if this node is a leaf (i.e., has no children).

        Returns:
        --------
        bool
            True if the node is a leaf, False otherwise.
        """
        return len(self.children) == 0

    def to_newick(self):
        """
        Convert the subtree rooted at this node to a Newick string.

        Returns:
        --------
        str
            The Newick string representing the subtree.
        """
        if self.is_leaf():
            return f"{self.name}:{self.branch_length}"
        else:
            children_str = ",".join(child.to_newick() for child in self.children)
            return f"({children_str}):{self.branch_length}"


class Tree:
    """
    A class to represent a phylogenetic tree.

    Attributes:
    -----------
    root : TreeNode
        The root of the tree.
    """

    def __init__(self, root_name: str = None):
        """
        Initialize a Tree instance.

        Parameters:
        -----------
        root_name : str, optional
            The name of the root node (default is None).
        """
        if root_name:
            self.root = TreeNode(name=root_name)
        else:
            self.root = None

    def from_newick(self, newick: str):
        """
        Parse a Newick string and build the tree.

        Parameters:
        -----------
        newick : str
            The Newick string representing the tree.
        """
        stack = []
        current_node = TreeNode(name='Root')
        i = 0

        while i < len(newick):
            char = newick[i]
            
            if char == '(':
                stack.append(current_node)
                current_node = TreeNode(name='')
            elif char == ',':
                parent = stack[-1]
                parent.add_child(current_node)
                current_node = TreeNode(name='')
            elif char == ')':
                parent = stack.pop()
                parent.add_child(current_node)
                current_node = parent
            elif char == ':':
                j = i + 1
                while j < len(newick) and (newick[j].isdigit() or newick[j] == '.'):
                    j += 1
                branch_length = float(newick[i+1:j])
                current_node.branch_length = branch_length
                i = j - 1
            elif char == ';':
                break
            else:
                j = i
                while j < len(newick) and newick[j] not in '(),:;':
                    j += 1
                current_node.name = newick[i:j]
                i = j - 1

            i += 1

        self.root = current_node

    def to_newick(self):
        """
        Convert the tree to a Newick string.

        Returns:
        --------
        str
            The Newick string representing the tree.
        """
        if self.root:
            return self.root.to_newick() + ';'
        else:
            return ''

    def read_from_file(self, file_path: str):
        """
        Read a Newick string from a file and build the tree.

        Parameters:
        -----------
        file_path : str
            The path to the file containing the Newick string.
        """
        with open(file_path, 'r') as file:
            newick = file.read().strip()
            self.from_newick(newick)

    def write_to_file(self, file_path: str):
        """
        Write the tree to a file in Newick format.

        Parameters:
        -----------
        file_path : str
            The path to the file where the Newick string will be written.
        """
        with open(file_path, 'w') as file:
            newick = self.to_newick()
            file.write(newick)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python tree.py <input_tree_file>")
        sys.exit(1)

    input_file = sys.argv[1]

    # Read the tree from the input file
    my_tree = Tree()
    my_tree.read_from_file(input_file)
    
    # Output the Newick string to standard output
    print(my_tree.to_newick())
