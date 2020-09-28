# Interpolation methods for cyclic min/max values like temperature

"""
    Cosine()

Similar to a `Linear` interpolator, but using `sin` to calculate
the weights, as if it is cyclic between two values.
"""
struct Cosine <: Interpolations.Degree{1} end

Interpolations.positions(::Cosine, ax::AbstractUnitRange{<:Integer}, x) =
    Interpolations.positions(Linear(), ax, x)

Interpolations.value_weights(::Cosine, δx) = begin
    # calculate the fraction from a sin
    cosfrac = -cos(π * δx)
    # normalise between 0 and 1
    normed = (cosfrac + 1) / 2
    (1-normed, normed)
end
Interpolations.gradient_weights(::Cosine, δx) = (-oneunit(δx), oneunit(δx))
Interpolations.hessian_weights(::Cosine, δx) = (zero(δx), zero(δx))
Interpolations.padded_axis(ax::AbstractUnitRange, ::BSpline{Cosine}) = ax

"""
    HypTan{S}()

Similar to a `Linear` interpolator, but using hyperbolic tangent to 
calculate the weights. It has a multiplier `X` that controls how long the
curve is near each index. This is automatically scaled so the max and min 
values are always 0 and 1 for all `X`. At the limits, HypTan is essentially 
`Linear` and `NoInterp`. 

`HypTan{0}` will error, instead uses `Linear` - which is identical.
"""
struct HypTan{X,S} <: Interpolations.Degree{1} end
HypTan() = HypTan{2}()
HypTan{X}() where X = HypTan{X,1 / tanh(π * X / 2)}()
HypTan{0}() = Linear()

Interpolations.positions(::HypTan, ax::AbstractUnitRange{<:Integer}, x) =
    Interpolations.positions(Linear(), ax, x)

Interpolations.value_weights(ht::HypTan{X,S}, δx) where {X,S} = begin
    # calculate the fraction from a sin
    hypfrac = tanh(X * π * (δx - 0.5)) * S
    # normalise between 0 and 1
    normed = (hypfrac + 1) / 2
    (1-normed, normed)
end
Interpolations.gradient_weights(::HypTan, δx) = (-oneunit(δx), oneunit(δx))
Interpolations.hessian_weights(::HypTan, δx) = (zero(δx), zero(δx))
Interpolations.padded_axis(ax::AbstractUnitRange, ::BSpline{HypTan}) = ax
