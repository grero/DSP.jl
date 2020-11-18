# Filter response functions

#
# Frequency response of a digital filter
#

using ..DSP: xcorr

"""
    H, w = freqresp(filter)

Frequency response `H` of a `filter` at (normalized) frequencies `w` in
radians/sample for a digital filter or radians/second for an analog filter
chosen as a reasonable default.
"""
function freqresp(filter::FilterCoefficients)
    w = _freqrange(filter)
    return (freqresp(filter, w), w)
end

"""
    freqresp(filter, w)

Frequency response of `filter` at (normalized) frequency or frequencies `w` in
radians/sample for a digital filter or radians/second for an analog filter.
"""
freqresp(filter::FilterCoefficients, w::AbstractVector) = [freqresp(filter, i) for i = w]

"""
    freqresp(filter, hz, fs)

Frequency response of a `filter` at frequency or frequencies `hz` with sampling rate `fs`
for a digital filter or frequencies `hz/fs` for an analog filter.
"""
freqresp(filter::FilterCoefficients, hz::Union{Number, AbstractVector}, fs::Number) =
    freqresp(filter,  hz * ((2 * pi) / fs))


freqresp(filter::FilterCoefficients{:z}, w::Number) = _freq(filter, exp(im * w))

freqresp(filter::FilterCoefficients{:s}, w::Number) = _freq(filter, im * w)

function _freq(filter::FilterCoefficients, x::Number)
    filter = convert(PolynomialRatio, filter)
    return filter.b(x) ./ filter.a(x)
end

_freq(filter::ZeroPoleGain, x::Number) =
    filter.k * prod([x - z for z in filter.z]) / prod([x - p for p in filter.p])

function _freq(filter::Biquad, x::Number)
    x2 = x*x
    return (filter.b0*x2 + filter.b1*x + filter.b2) / (x2 + filter.a1*x  + filter.a2)
end

_freq(filter::SecondOrderSections, x::Number) =
    filter.g * prod([_freq(b, x) for b in filter.biquads])


"""
    phi, w = phaseresp(filter)

Phase response `phi` of a `filter` at (normalized) frequencies `w` in
radians/sample for a digital filter or radians/second for an analog filter
chosen as a reasonable default.
"""
function phaseresp(filter::FilterCoefficients)
    w = _freqrange(filter)
    return (phaseresp(filter, w), w)
end

"""
    phaseresp(filter, w)

Phase response of a `filter` at (normalized) frequency or frequencies `w` in
radians/sample for a digital filter or radians/second for an analog filter.
"""
function phaseresp(filter::FilterCoefficients, w)
    h = freqresp(filter, w)
    unwrap(angle.(h); dims=ndims(h))
end

"""
    tau, w = grpdelay(filter)

Group delay `tau` of a `filter` at (normalized) frequencies `w` in
radians/sample for a digital filter or radians/second for an analog filter
chosen as a reasonable default.
"""
function grpdelay(filter::FilterCoefficients)
    w = _freqrange(filter)
    return (grpdelay(filter, w), w)
end

"""
    grpdelay(fliter, w)

Group delay of a digital 'filter' at normalized frequency
or frequencies 'w' in radians/sample.
"""
function grpdelay(filter::FilterCoefficients{:z}, w)
    filter = convert(PolynomialRatio, filter)
    b, a = coefb(filter), coefa(filter)

    # Linear Phase FIR
    if (length(a) == 1) & (_is_sym(b) | _is_anti_sym(b))
        return fill((length(b)-1)/2, length(w))
    end

    c = xcorr(b, a; padmode = :none)
    cr = range(0, stop=length(c)-1) .* c
    ejw = exp.(-im .* w)
    num = Polynomial(cr).(ejw)
    den = Polynomial(c).(ejw)
    return real.(num ./ den) .- (length(a) - 1)
end

function grpdelay(filter::FilterCoefficients{:s}, w)
    filter = convert(PolynomialRatio, filter)
    b, a = filter.b, filter.a
    bd = derivative(b)
    ad = derivative(a)
    s = im .* w
    return real.((bd*a - ad*b).(s) ./ (a * b).(s))
end

"""
    impresp(filter, n=100)

Impulse response of a digital `filter` with `n` points.
"""
function impresp(filter::FilterCoefficients{:z}, n=100)
  i = [1; zeros(n-1)]
  filt(filter, i)
end

"""
    stepresp(filter, n=100)

Step response of a digital `filter` with `n` points.
"""
function stepresp(filter::FilterCoefficients{:z}, n=100)
  i = ones(n)
  filt(filter, i)
end



#
# Helper functions
#

function _is_sym(x::AbstractArray)
    n = length(x) ÷ 2
    return all(x[1+i] == x[end-i] for i in 0:n-1)
end

function _is_anti_sym(x::AbstractArray)
    n = length(x) ÷ 2
    return all(x[1+i] == -x[end-i] for i in 0:n)
end

_freqrange(::FilterCoefficients{:z}) = range(0, stop=π, length=257)
function _freqrange(filter::FilterCoefficients{:s})
    filter = convert(ZeroPoleGain, filter)
    w_interesting = sort!(Float64.(abs.([filter.p; filter.z])))
    include_zero = !isempty(w_interesting) && iszero(w_interesting[1])
    w_interesting = collect(Iterators.dropwhile(iszero, w_interesting))
    if isempty(w_interesting) # no non-zero poles or zeros
        if !include_zero || !isfinite(1/filter.k)
            return [0; 10 .^ (0:6)] # fallback
        end
        # include the point where |H|=1 (if any) and go further by factor 10
        return range(0.0, stop=10 * Float64(max(filter.k, 1/filter.k)), length=200)
    end
    # normal case: go from smalles to largest pole/zero, extended by factor 10
    w_min, w_max = w_interesting[[1,end]]
    w = 10 .^ range(log10(w_min)-1, stop=log10(w_max)+1, length=200)
    return include_zero ? [0.0; w] : w
end
