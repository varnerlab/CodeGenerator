mutable struct VLReaction

    reaction::String
    reaction_name::String
    left_phrase::String
    right_phrase::String
    reverse::String
    forward::String

    function VLReaction()
		this = new()
	end
end

mutable struct VLSpeciesSymbol

  lexeme::String

  function VLSpeciesSymbol()
		this = new()
	end
end
