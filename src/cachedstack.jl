"""
    CachedStack(stack; loaded=Dict{Symbol,Any}()) 

A stack that wraps another Stack and stores any layers loaded 
from it in a `Dict`. This prevents loading them multiple times
from disk.

It also means the are locked from the garbage collector if they are 
still in the Dict, which can cause problems with memory use.  Use 
`clear!(stack)` to return to the unloaded state, freeing memory.
"""
mutable struct CachedStack{S} <: AbstractGeoStack{S}
    stack::S
    cache::Dict{Symbol,Any}
end
CachedStack(stack; loaded=Dict{Symbol,Any}()) = 
    CachedStack(stack, loaded)

stack(s::CachedStack) = s.stack
cache(s::CachedStack) = s.cache

# GeoData methods

DD.dims(s::CachedStack) = DD.dims(stack(s))
DD.refdims(s::CachedStack) = DD.refdims(stack(s))
DD.metadata(s::CachedStack) = DD.metadata(stack(s))
DD.layerdims(s::CachedStack) = DD.layerdims(stack(s))
DD.layermetadata(s::CachedStack) = DD.layermetadata(stack(s))
GD.layermissingval(s::CachedStack) = GD.layermissingval(stack(s))

DD.layers(s::CachedStack) = map(k -> s[k], keys(s))
function DD.rebuild(s::CachedStack, args...; kw...)
    inner = rebuild(stack(s), args...; kw...)
    CachedStack(inner, cache(s))
end

# Base methods

Base.keys(s::CachedStack) = keys(stack(s))
Base.values(s::CachedStack) = values(stack(s))
Base.names(s::CachedStack) = names(stack(s))
Base.length(s::CachedStack) = length(stack(s))

function Base.getindex(cstack::CachedStack, key::Symbol)
    c = cache(cstack)
    if haskey(c, key)
        return c[key]
    else
        s = read(stack(cstack)[key])
        c[key] = s
        return s
    end
end

# Clear

clear!(stack::AbstractGeoStack) = clear(data(stack))
clear!(stack) = nothing
function clear!(stack::CachedStack) 
    for k in keys(cache(stack))
        delete!(cache(cache), k)
    end
end
