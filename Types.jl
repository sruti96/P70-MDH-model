mutable struct TranslationParameters

    vmax::Float64
    tau_factor::Float64
    KL::Float64

    function TranslationParameters()
		this = new()
	end
end
