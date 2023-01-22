
"""
AWGN通信路(実数)
"""
function awgn!(input::AbstractArray{T}, N0) where T<:Real
    input .+= sqrt(N0/2)*randn(size(input))
end

"""
AWGN通信路(複素数)
"""
function awgn!(input::AbstractArray{T}, N0) where T<:Complex
    input .+= sqrt(N0)*randn(ComplexF64, size(input))
end

awgn(input, N0) = awgn!(copy(input), N0)


"""
マルチパスフェージング通信路
"""
function multipath_fading!(output, input, chresp; fdTs=0.0, sps=1, istail=false)
    Nrx, Ntx = size(output,2), size(input,2)
    # 畳み込み
    for (m,n) in Iterators.product(1:Nrx, 1:Ntx)
        if ndims(chresp) <= 3
            @views transversal_filter!(output[:,m], input[:,n], chresp[:,m,n]) # 時不変応答
        else
            @views transversal_filter!(output[:,m], input[:,n], chresp[:,:,m,n]) # 時変応答
        end
    end
    output
end

"""
    multipath_fading(input, pdp::AbstractVector, ndim, N0; fdTs, sps, istail)

マルチパスフェージング通信路

    Arguments: 
        input  : transmit signal
        pdp    : power delay profile
        ndim   : Tx / Rx antenna size
        N0     : noise density
        fdTs   : Normalized Doppler Frequency
        sps    : Sample / Symbol
        istail : tail biting
    
    Returns:
        output : 
        chresp :
"""
function multipath_fading(input, pdp::AbstractVector, ndim, N0=0.0; fdTs=0.0, sps=1, istail=false)
    Ntx, Nrx = ndim
    N, L = size(input,1), length(pdp)
    M = N + Int(istail)*(L - 1)
    output = zeros(eltype(input), M, Nrx)

    # 通信路応答生成
    chresp = gen_fading(parent(pdp), M, Nrx, Ntx; fdTs=fdTs, sps=sps, istail=istail, by=Jakes())

    # チャネル応答取得
    output = multipath_fading!(output, input, chresp; fdTs=fdTs, sps=sps, istail=istail)
    N0 > 0.0 && awgn!(output, N0)
    output, chresp
end


function multipath(input, pdp, ndim; ch="fading", fdTs=0.0, sps=1, istail=false)
    Ntx, Nrx = ndim
    N, L = size(input,1), length(pdp)
    M = N + Int(istail)*(L - 1)
    output = zeros(eltype(input), M, Nrx) # 出力配列

    if ch=="static"
        path_coefs = repeat(sqrt.(pdp), Nrx, Ntx) .* exp.(im*(2*pi*rand(L, Nrx, Ntx)))
        static_channel!(output, input, path_coefs)
    elseif ch=="fading"

    end
    output, path_coefs
end

function static_channel!(output, input, path_coefs)
    # 畳み込み
    Nrx, Ntx = size(output,2), size(input,2)
    for (m,n) in Iterators.product(1:Nrx, 1:Ntx)
        @views transversal_filter!(output[:,m], input[:,n], path_coefs[:,m,n], conj_flag=false) # 時不変応答
    end
end

"""
Rayleigh Fading
"""
function rayfading(pdp, nrx=1, ntx=1)
    L = length(pdp)
    h = pdp .* randn(ComplexF64, L, nrx, ntx)
    h
end
