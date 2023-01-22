# window functions

"""
    raised_cosine(N; symmetric=false)
Raised Cosine Window
Hann window corresponds to Raised Cosine Window with α = 1/2, β = 1/4
    1/2 + 1/2 \cos(2*\pi*n/N)
"""
function raised_cosine(N; symmetric=false)
    if symmetric
        # symmetric
        w = zeros(Float64, N)
        for k in 1:N
            n = (k-1) - N/2 + 1/2
            w[k] = 1/2 + 1/2*cos(2π*n/N)
        end
    else
        # cyclic
        w = zeros(Float64, N)
        for k in 1:N
            n = (k-1) - div(N,2)
            w[k] = 1/2 + 1/2*cos(2π*n/N)
        end
    end
    w
end

root_raised_cosine(N) = sqrt.(raised_cosine(N))