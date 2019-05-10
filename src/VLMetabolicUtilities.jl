function contains(string,token)
    return occursin(token,string)
end

function build_flux_bounds_array(reaction_array::Array{VLReaction,1})

    # dimensions -
    number_of_reactions = length(reaction_array)

    # initialize -
    default_bounds_array = zeros(number_of_reactions,2)

    # build the flux bounds array -
    for reaction_index = 1:number_of_reactions

        # get bounds strings -
        lower_bound_string = reaction_array[reaction_index].reverse
        upper_bound_string = reaction_array[reaction_index].forward

        # build -
        default_bounds_array[reaction_index,1] = parse(Float64,lower_bound_string)
        default_bounds_array[reaction_index,2] = parse(Float64,upper_bound_string)
    end

    # return -
    return default_bounds_array
end

# Function to build the stoichiometric_matrix from reaction list -
function build_stoichiometric_matrix(species_symbol_array::Array{VLSpeciesSymbol,1},reaction_array::Array{VLReaction,1})

  # Method variables -
  number_of_species = length(species_symbol_array)
  number_of_reactions = length(reaction_array);
  stoichiometric_matrix = zeros(number_of_species,number_of_reactions);

  function _parse_reaction_phrase(lexeme,reaction_phrase)

    # Split around + -
    coefficient = 0.0;
    fragment_array = split(reaction_phrase,"+")
    for fragment in fragment_array

      if (contains(fragment,"*"))
          local_fragment_array = split(fragment,"*");
          test_lexeme = local_fragment_array[2];

          if (lexeme == test_lexeme)
            coefficient = parse(Float64,local_fragment_array[1]);
            break
          end

      else

        # Build -
        test_lexeme = fragment;
        if (lexeme == test_lexeme)
          coefficient = 1.0;
          break
        end
      end
    end

    return coefficient;
  end

  function _find_stoichiometric_coefficient(species_model::VLSpeciesSymbol,reaction::VLReaction)

    # Method variables -
    stoichiometric_coefficient = 0.0

    # Check the left and right phrase -
    stoichiometric_coefficient += -1.0*(_parse_reaction_phrase(species_model.lexeme,reaction.left_phrase))
    stoichiometric_coefficient += _parse_reaction_phrase(species_model.lexeme,reaction.right_phrase)
    return stoichiometric_coefficient;
  end


  # setup counters -
  for (row_index,species_symbol) in enumerate(species_symbol_array)
    for (col_index,reaction) in enumerate(reaction_array)

      # Is this species involved in this reaction?
      stoichiometric_matrix[row_index,col_index] = _find_stoichiometric_coefficient(species_symbol,reaction);

    end
  end

  # return -
  return stoichiometric_matrix
end


# Function to build symbol list -
function build_symbol_array(reaction_array::Array{VLReaction,1})

  # Method variables -
  species_symbol_array::Array{VLSpeciesSymbol,1} = []

  # Helper function to parse the reaction phrases, split out the symbols
  function _parse_phrase(reaction_phrase::String)

    # Method variables -
    local_species_array::Array{VLSpeciesSymbol} = []

    # Split around + -
    fragment_array = split(reaction_phrase,"+")
    for fragment in fragment_array

      if (contains(fragment,"*"))

          local_fragment_array = split(fragment,"*");
          species_symbol = VLSpeciesSymbol();
          species_symbol.lexeme = local_fragment_array[2];

      else

        # Build -
        species_symbol = VLSpeciesSymbol()
        species_symbol.lexeme = fragment
      end

      # grab -
      push!(local_species_array,species_symbol)
    end

    # return -
    return local_species_array
  end

  function _isequal(species_model_1::VLSpeciesSymbol,species_model_2::VLSpeciesSymbol)
    if (species_model_1.lexeme == species_model_2.lexeme)
      return true
    end
    return false
  end

  function _add_symbol!(species_symbol_array,species_symbol)

    contains_species_already = false
    for cached_species_model in species_symbol_array

      if (_isequal(cached_species_model,species_symbol))
        contains_species_already = true
        break
      end
    end

    if (contains_species_already == false)
      push!(species_symbol_array,species_symbol)
    end
  end

  # iterate through and get the symbols -
  for reaction in reaction_array

    tmp_species_array_left = _parse_phrase(reaction.left_phrase)
    tmp_species_array_right = _parse_phrase(reaction.right_phrase)
    append!(tmp_species_array_left,tmp_species_array_right)

    for species_model in tmp_species_array_left
      _add_symbol!(species_symbol_array,species_model)
    end
  end

  # ok, so the species symbol array is *not* sorted)
  # let's sort the species -
  tmp_array = String[]
  for species_symbol::VLSpeciesSymbol in species_symbol_array
      push!(tmp_array,species_symbol.lexeme)
  end

  # generate permutation array -
  idxa_sorted = sortperm(tmp_array)
  sorted_symbol_array = species_symbol_array[idxa_sorted]

  # ok, if a species contains [], then put it at the end -
  partitioned_symbol_array = VLSpeciesSymbol[]
  external_symbol_array = VLSpeciesSymbol[]
  for species_symbol::VLSpeciesSymbol in sorted_symbol_array

      if (occursin(r"\[*\]",species_symbol.lexeme) == false)
          push!(partitioned_symbol_array, species_symbol)
      else
          push!(external_symbol_array,species_symbol)
      end
  end

  # add tmp's back to the bottom -
  for species_symbol in external_symbol_array
      push!(partitioned_symbol_array, species_symbol)
  end

  # return -
  return partitioned_symbol_array
end


# function to parse a specific vff reaction -
function parse_vff_reaction_string(reaction_string::String, reversible::Bool; delimiter::String=",")::VLReaction

    # initialize -
    reaction_wrapper = VLReaction()

    # spilt into fragments -
    fragment_array = split(reaction_string,delimiter)

    # Grab the fields -
    local_name = fragment_array[1]
    left_phrase = fragment_array[2]
    right_phrase = fragment_array[3]

    # add stuff to the reaction warpper -
    reaction_wrapper.reaction = reaction_string
    reaction_wrapper.reaction_name = local_name
    reaction_wrapper.left_phrase = left_phrase
    reaction_wrapper.right_phrase = right_phrase

    # handle reversible -
    if reversible == true
        reaction_wrapper.reverse = "-inf"
        reaction_wrapper.forward = "inf"
    else
        reaction_wrapper.reverse = "0.0"
        reaction_wrapper.forward = "inf"
    end

    # return -
    return reaction_wrapper
end

# Function to load VFF reactions from file -
function parse_vff_reaction_file(path_to_reaction_file::String)

  # Method variables -
  reaction_array::Array{VLReaction,1} = []

  try

    # Open the file, read it line by line, build reaction wrapper along the way -
    file = open(path_to_reaction_file)
    for reaction_line in eachline(file)

        # check do we have a comment?
        if (occursin("//",reaction_line) == false && isempty(reaction_line) == false)

            # let the user know ...
            println("Processing: $(reaction_line)")

            # spilt into fragments -
            fragment_array = split(reaction_line,",")

            # Grab the fields -
            local_name = fragment_array[1]
            left_phrase = fragment_array[3]
            right_phrase = fragment_array[4]
            reverse_flag = fragment_array[5]
            forward_flag = fragment_array[6]

            # Make a new reaction type, store in array -
            reaction_wrapper = VLReaction()
            reaction_wrapper.reaction = reaction_line
            reaction_wrapper.reaction_name = local_name
            reaction_wrapper.left_phrase = left_phrase
            reaction_wrapper.right_phrase = right_phrase
            reaction_wrapper.reverse = reverse_flag
            reaction_wrapper.forward = forward_flag
            push!(reaction_array,reaction_wrapper)
        end
    end

    # Make sure to close the file!
    close(file)

  catch error
    showerror(stdout, error, backtrace()); println()
  end

  # return my reaction array to caller -
  return reaction_array
end

function generate_net_reaction_string(uptake_array::Array{Float64,1},epsilon::Float64,data_dictionary::Dict{AbstractString,Any})

  # get list of metabolite symbols -
  list_of_metabolite_symbols = data_dictionary["list_of_metabolite_symbols"]

  # check for smalls -
  idx_small = find(abs(uptake_array).<epsilon)
  uptake_array[idx_small] = 0.0

  # which elememts are positive (products)?
  idx_product_array = find(uptake_array.>0)

  # which elements are negative (reactants?)
  idx_reactant_array = find(uptake_array.<0)

  # build the string ...
  net_reaction_buffer = ""
  for idx_reactant in idx_reactant_array

    metabolite_symbol = list_of_metabolite_symbols[idx_reactant]
    st_coeff = round(abs(uptake_array[idx_reactant]),2)

    if (st_coeff != 1.0)
      net_reaction_buffer *= "$(st_coeff)*$(metabolite_symbol) + "
    else
      net_reaction_buffer *= "$(metabolite_symbol) + "
    end
  end

  # cutoff trailing * -
  net_reaction_buffer = net_reaction_buffer[1:end-3]

  # add the arrow -
  net_reaction_buffer *= " --> "

  # write the trailing stuff -
  for idx_product in idx_product_array

    metabolite_symbol = list_of_metabolite_symbols[idx_product]
    st_coeff = round(abs(uptake_array[idx_product]),2)

    if (st_coeff != 1.0)
      net_reaction_buffer *= "$(st_coeff)*$(metabolite_symbol) + "
    else
      net_reaction_buffer *= "$(metabolite_symbol) + "
    end
  end

  # cutoff trailing * -
  net_reaction_buffer = net_reaction_buffer[1:end-3]

  # return -
  return net_reaction_buffer
end

function build_reaction_bounds_matrix_from_vff_file(path_to_reaction_file::String)

    # TODO - is the reaction file path legit?
    if (isempty(path_to_reaction_file) == true || length(path_to_reaction_file) == 1)
        error_message = "ERROR: Path to reaction file is empty"
        throw(error(error_message))
    end

    # get reaction array -
    reaction_array = parse_vff_reaction_file(path_to_reaction_file)

    # build the flux bounds -
    return build_flux_bounds_array(reaction_array)
end

function build_stoichiometric_matrix_from_vff_file(path_to_reaction_file::String)

    # TODO - is the reaction file path legit?
    if (isempty(path_to_reaction_file) == true || length(path_to_reaction_file) == 1)
        error_message = "ERROR: Path to reaction file is empty"
        throw(error(error_message))
    end

    # get reaction array -
    reaction_array = parse_vff_reaction_file(path_to_reaction_file)

    # get the species list -
    sorted_symbol_array = build_symbol_array(reaction_array)

    # build the stoichiometric_matrix -
    return build_stoichiometric_matrix(sorted_symbol_array, reaction_array)
end
