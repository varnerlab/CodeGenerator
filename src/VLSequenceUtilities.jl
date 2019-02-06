function read_sequence_file(path_to_seq_file::String)

    # TODO: Is the path a legit path?
    # TODO: handle the error -

    # Load the file, read the sequence into a string and return
    file = open(path_to_seq_file)
    seq = read(file, String)
    close(file)

    # return the sequence string -
    return seq
end

function build_translation_reactions_for_protein_sequence(protein_name::String,protein_seq::String)

  # Load the AA symbol map -
  path_to_package = dirname(pathof(CodeGenerator))
  path_to_mapping_file = joinpath(path_to_package,"config/AAMap.csv")

  map_array = readdlm(path_to_mapping_file,','); #metabolite 1, one letter 2
  protein_aa_dictionary = Dict();

  # Create a mapping dictionary -
  symbol_metabolite_map = Dict();
  for map_index in collect(1:20)

    one_letter_aa_symbol = map_array[map_index,2];
    metabolite_symbol = map_array[map_index,1];
    symbol_metabolite_map[one_letter_aa_symbol] = metabolite_symbol;
    protein_aa_dictionary[metabolite_symbol] = 0.0;
  end


  # Parse the protein seq -
  number_aa_residues = length(protein_seq);
  local_counter = 0;
  for aa_index in collect(1:number_aa_residues)

    # What AA do we have?
    aa_value = string(protein_seq[aa_index]);
    if (aa_value != "\n" && aa_value != " ")

      key = symbol_metabolite_map[aa_value];

      # Update the dictionary -
      quantity_aa = protein_aa_dictionary[key];
      protein_aa_dictionary[key] = quantity_aa + 1;
      local_counter+=1;
    end
  end

  # Ok, we have the protein sequence , build the reaction string buffer -
  buffer="";
  buffer*="translation_initiation_$(protein_name),mRNA_$(protein_name)+RIBOSOME,RIBOSOME_START_$(protein_name),0,inf;\n"
  buffer*="translation_$(protein_name),RIBOSOME_START_$(protein_name)+$(2*local_counter)*M_gtp_c+$(2*local_counter)*M_h2o_c";
  for aa_index in collect(1:20)

    # Get charged tRNA -
    metabolite_symbol = map_array[aa_index,1];

    # number of this AA -
    value = protein_aa_dictionary[metabolite_symbol];

    # Add charged tRNA species to buffer -
    buffer*="+$(value)*$(metabolite_symbol)_tRNA";
  end

  # products -
  buffer*=",RIBOSOME+mRNA_$(protein_name)+PROTEIN_$(protein_name)+$(2*local_counter)*M_gdp_c+$(2*local_counter)*M_pi_c+$(local_counter)*tRNA,0,inf;\n"

  # Write the reactions for charing the tRNA -
  for aa_index in collect(1:20)

    # Get charged tRNA -
    metabolite_symbol = map_array[aa_index,1];

    # number of this AA -
    value = protein_aa_dictionary[metabolite_symbol];

    # Add charged tRNA species to buffer -
    buffer*="tRNA_charging_$(metabolite_symbol)_$(protein_name),$(value)*$(metabolite_symbol)+$(value)*M_atp_c+$(value)*tRNA+$(value)*M_h2o_c,";
    buffer*="$(value)*$(metabolite_symbol)_tRNA+$(value)*M_amp_c+$(value)*M_ppi_c,0,inf;\n";
  end

  # outfile = open("./translation_$(protein_name).txt", "w")
  # write(outfile,buffer);
  # close(outfile);

  return buffer;
end
