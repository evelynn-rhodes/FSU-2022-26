class Phylip:
    def __init__(self, file_path, sequence_type="amino_acid"):
        """
        Initialize the PHYLIP parser.

        Parameters:
        -----------
        file_path : str
            Path to the PHYLIP file.
        sequence_type : str
            Type of sequences: "amino_acid" or "nucleotide".
        """
        self.file_path = file_path
        self.sequence_type = sequence_type
        self.valid_amino_acids = set("ARNDCQEGHILKMFPSTWYV")
        self.valid_nucleotides = set("ACGT")

    def read_phylip(self):
        """
        Reads a PHYLIP format file and returns a dictionary where the keys are taxon names
        and the values are their sequences.

        Returns:
            dict: A dictionary with taxon names as keys and sequences as values.
        """
        sequences = {}
        with open(self.file_path, 'r') as file:
            # Read the first line to determine number of sequences and sequence length
            first_line = file.readline().strip()
            num_sequences, sequence_length = map(int, first_line.split())
            
            for line in file:
                # Split line into taxon name and sequence based on whitespace
                parts = line.strip().split(maxsplit=1)
                if len(parts) != 2:
                    raise ValueError(f"Invalid PHYLIP format in line: {line}")
                
                taxon, sequence = parts[0], parts[1].replace(":", "")
                sequence = sequence.upper()
                
                # Validate sequence
                valid_chars = (
                    self.valid_amino_acids if self.sequence_type == "amino_acid"
                    else self.valid_nucleotides
                )
                if not set(sequence).issubset(valid_chars):
                    raise ValueError(f"Invalid characters in sequence for {taxon}: {sequence}")
                
                if taxon in sequences:
                    sequences[taxon] += sequence
                else:
                    sequences[taxon] = sequence
            
            # Validate sequence counts and lengths
            if len(sequences) != num_sequences:
                raise ValueError(f"Expected {num_sequences} sequences, found {len(sequences)}.")
            for taxon, seq in sequences.items():
                if len(seq) != sequence_length:
                    raise ValueError(f"Sequence length mismatch for {taxon}: {len(seq)} != {sequence_length}.")
        
        return sequences