module CodeGenerator

# Include my files -
include("Include.jl")

# export methods -
export read_sequence_file
export build_translation_reactions_for_protein_sequence
export build_transcription_reactions_for_gene_sequence
export build_stoichiometric_matrix_from_vff_file

end # module
