

# binary set partitioning
function bin_set_partition(symbol_size::Integer)
    bitlen = ndigits(symbol_size-1, base=2) # max bit length
    set = Vector{Set{Int}}(undef, bitlen) # set vector
    for i in 1:bitlen
        inds = Set{Int}()
        base = 2^(bitlen-i)
        for x in 0:symbol_size-1
            b = (x&base) >>> (bitlen-i)
            b > 0 && push!(inds, x)
        end
        set[i] = inds
    end
    set
end # function
