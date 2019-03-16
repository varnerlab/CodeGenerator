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

function build_transcription_reactions_for_gene_sequence(gene_name::String, gene_seq::String)

  # function variables -
  buffer= "";
  total_ntp = 0;

  # generate the sequence dictionary -
  nucleotide_dictionary = parse_gene_seq(gene_seq);

  # write the RNAP binding step -
  buffer*="transcriptional_initiation_$(gene_name),GENE_$(gene_name)+RNAP,OPEN_GENE_$(gene_name),0,inf;\n";
  buffer*="transcription_$(gene_name),OPEN_GENE_$(gene_name)";


  # go through by dictionary, and get the base count -
  for (key,value) in nucleotide_dictionary

    if (key == "a")

      # write the M_atp_c line -
      buffer*="+$(value)*M_atp_c"

      # How many a's do we have?
      total_ntp += value;

    elseif (key == "t")

      # write the M_utp_c line -
      buffer*="+$(value)*M_utp_c"

      # How many u's do we have?
      total_ntp += value;

    elseif (key == "g")

      # write the M_gtp_c line -
      buffer*="+$(value)*M_gtp_c"

      # How many g's do we have?
      total_ntp += value;

    else

      # write the M_gtp_c line -
      buffer*="+$(value)*M_ctp_c"

      # How many c's do we have?
      total_ntp += value;
    end

  end

  # mRNA+GENE+RNAP+1320*M_pi_c,0,inf;
  buffer*="+$(total_ntp)*M_h2o_c,mRNA_$(gene_name)+GENE_$(gene_name)+RNAP+$(total_ntp)*M_ppi_c,0,inf;\n"

  # mRNA_decay degradation reaction -
  # mRNA_decay,[],mRNA,144*M_cmp_c+151*M_gmp_c+189*M_ump_c+176*M_amp_c,0,inf;
  buffer*="mRNA_degradation_$(gene_name),mRNA_$(gene_name),"
  local_buffer = "";
  for (key,value) in nucleotide_dictionary

    if (key == "a")

      # write the M_atp_c line -
      local_buffer*="+$(value)*M_amp_c"

    elseif (key == "t")

      # write the M_utp_c line -
      local_buffer*="+$(value)*M_ump_c"

    elseif (key == "g")

      # write the M_gtp_c line -
      local_buffer*="+$(value)*M_gmp_c"

    else

      # write the M_gtp_c line -
      local_buffer*="+$(value)*M_cmp_c"

    end
  end

  buffer*="$(local_buffer[2:end]),0,inf;\n"

  # dump -
  # outfile = open("./transcription_$(gene_name).txt", "w")
  # write(outfile,buffer);
  # close(outfile);

  # return the buffer -
  return buffer;
end

function parse_gene_seq(gene_seq::AbstractString)

  # We will return a dictionary w/nucleotide keys and the number of nucleotides as values -
  nucleotide_dictionary = Dict();
  nucleotide_dictionary["a"] = 0;
  nucleotide_dictionary["t"] = 0;
  nucleotide_dictionary["g"] = 0;
  nucleotide_dictionary["c"] = 0;

  # What is the length of the gene sequence -
  number_of_nucleotides = length(gene_seq);

  for nucleotide_index = collect(1:number_of_nucleotides)

    # get the test nucleotide -
    test_nucleotide = gene_seq[nucleotide_index];
    if (test_nucleotide == 'a')
      value = nucleotide_dictionary["a"];
      nucleotide_dictionary["a"] = value + 1;
    end

    if (test_nucleotide == 't')
      value = nucleotide_dictionary["t"];
      nucleotide_dictionary["t"] = value + 1;
    end

    if (test_nucleotide == 'g')
      value = nucleotide_dictionary["g"];
      nucleotide_dictionary["g"] = value + 1;
    end

    if (test_nucleotide == 'c')
      value = nucleotide_dictionary["c"];
      nucleotide_dictionary["c"] = value + 1;
    end
  end

  # return -
  return nucleotide_dictionary;
end

function build_translation_reactions_for_protein_sequence(protein_name::String,protein_seq::String)

  # Load the AA symbol map -
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
