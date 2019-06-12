module SobolSearches

export SobolSearch, setup_sobol_search, sobol_search_next!, sobol_search!, sobol_search

using ArgCheck: @argcheck
using DocStringExtensions: SIGNATURES
using Parameters: @unpack
import Sobol

struct SobolSearch{TV,TF,TS,TM}
    lower_bounds::TV
    upper_bounds::TV
    minimand::TF
    sobol_seq::TS
    minima::TM
end

function _transform_position(lower_bounds, upper_bounds, s)
    map((l, u, α) -> l*(1-α) + u*α, lower_bounds, upper_bounds, s)
end

function _evaluated_position(minimand, lower_bounds, upper_bounds, s)
    position = _transform_position(lower_bounds, upper_bounds, s)
    (value = minimand(position), position = position)
end

function setup_sobol_search(minimand, lower_bounds, upper_bounds, keep;
                            skip = 0)
    @argcheck length(lower_bounds) == length(upper_bounds)
    @argcheck all(lower_bounds .< upper_bounds)
    @argcheck keep ≥ 1
    sobol_seq = Sobol.SobolSeq(length(lower_bounds))
    skip > 0 && Sobol.skip(sobol_seq, skip)
    minima = [_evaluated_position(minimand, lower_bounds, upper_bounds, Sobol.next!(sobol_seq))
              for _ in 1:keep]
    sort!(minima, by = x -> x.value)
    SobolSearch(lower_bounds, upper_bounds, minimand, sobol_seq, minima)
end

function sobol_search_next!(search::SobolSearch)
    @unpack lower_bounds, upper_bounds, minimand, sobol_seq, minima = search
    candidate = _evaluated_position(minimand, lower_bounds, upper_bounds,
                                    Sobol.next!(sobol_seq))
    len = length(minima)
    if candidate.value < minima[len].value
        n = searchsortedfirst(minima, candidate, by = x -> x.value)
        for i in (len - 1):-1:n
            minima[i + 1] = minima[i]
        end
        minima[n] = candidate
        true
    else
        false
    end
end

function sobol_search!(search::SobolSearch, steps)
    for _ in 1:steps
        sobol_search_next!(search)
    end
    search
end

"""
$(SIGNATURES)

Find the `keep` lowers minima of `miniman` along a Sobol sequence of length `steps` in the
hypercube between `lower_bounds` and `upper_bounds`.

Return a [`SobolSearch`](@ref) object.

```julia
```
"""
function sobol_search(minimand, lower_bounds, upper_bounds, keep, steps; skip = 0)
    search = setup_sobol_search(minimand, lower_bounds, upper_bounds, keep;
                                skip = skip)
    sobol_search!(search, steps)
end

Base.IndexStyle(::Type{<:SobolSearch}) = Base.IndexStyle(search.minima)

Base.size(search::SobolSearch) = size(search.minima)

Base.getindex(search::SobolSearch, I) = getindex(search.minima, I)

end # module
