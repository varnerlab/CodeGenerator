module CodeGenerator

# Include my files -
include("Include.jl")

# export methods -
export read_sequence_file
export build_translation_reactions_for_protein_sequence
export build_transcription_reactions_for_gene_sequence
export build_stoichiometric_matrix_from_vff_file
export build_reaction_bounds_matrix_from_vff_file
export write_debug_report
export generate_net_reaction_string

export parse_vff_reaction_string


end # module
