function contains(string::String,token::String)
    return occursin(token,string)
end

function build_flux_bounds_array(reaction_array::Array{VLReaction,1})
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
            coefficient = float(local_fragment_array[1]);
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

  # return -
  return species_symbol_array
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
            reaction_wrapper.reaction_name = local_name
            reaction_wrapper.left_phrase = left_phrase
            reaction_wrapper.right_phrase = right_phrase
            reaction_wrapper.reverse = reverse_flag
            reaction_wrapper.forward = chop(chomp(forward_flag))
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

function build_stoichiometric_matrix_from_vff_file(path_to_reaction_file::String)

    # TODO - is the reaction file path legit?
    if (isempty(path_to_reaction_file) == true || length(path_to_reaction_file) == 1)
        error_message = "ERROR: Path to reaction file is empty"
        throw(error(error_message))
    end

    # get reaction array -
    reaction_array = parse_vff_reaction_file(path_to_reaction_file)

    # get the species list -
    symbol_array = build_symbol_array(reaction_array)

    # build the stoichiometric_matrix -
    return build_stoichiometric_matrix(symbol_array, reaction_array)
end
