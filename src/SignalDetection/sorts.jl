using LinearAlgebra: norm

struct SIR end
struct SNR end
struct SINR end
struct Norm end
struct PostSNR end
struct PostSINR end

# columnsort:
function columnsort!(ord, H, criterion=:SIR; N0=nothing, rev=true)
    Nrx, Ntx = size(H)
    evals = Array{Float64}(undef, Ntx) # 評価値配列
    if criterion==:SIR # SIR基準
        hnorm = zeros(Ntx); sumnorm = 0.0
        for j = axes(H,2)
            for i = axes(H,1)
                hnorm[j] += abs2(H[i,j])
            end
            hnorm[j] = sqrt(hnorm[j])
            sumnorm += hnorm[j]
        end
        for i in axes(H,2)
            evals[i] = hnorm[i] / (sumnorm - hnorm[i])
        end

    elseif criterion==:SNR && !isnothing(W) # SNR基準
        for j in 1:Ntx
            w = conj(view(W,j,:)) # ウェイトベクトル
            eqnvar = real(w' * (nvar .* w)) # 二次形式
            evals[j] = 1/eqnvar
        end

    elseif criterion==:SINR && !isnothing(W) # SINR基準
        for j in 1:Ntx
            w = transpose(view(W,j,:)) # ウェイトベクトル
            mu = real(w * view(H,:,j)) # 利得
            evals[j] = 1/(1 - mu)
        end
    else
        error("{:SIR, :SNR, :SINR} のいずれかを選択してください.")
    end
    sortperm!(ord, evals, rev=rev) # 降順
end

function columnsort(H, criterion; N0=nothing, rev=true)
    ord = Array{Int}(undef, size(H,2))
    columnsort!(ord, H, criterion; N0=N0, rev=rev)
    return ord
end

# 大きい順に並べる
function columnsort!(::SIR, order::AbstractVector, A::AbstractMatrix; rev=true)
    squared_norms = Vector{Float64}(undef, size(A,2))
    SIR = Vector{Float64}(undef, size(A,2))
    for a in eachcol(A)
        squared_norms[i] = dot(a,a) |> real # 二乗ノルム
    end
    sum_norms = sum(squared_norms)
    for i in eachindex(evals)
        SIR[i] /= (sum_norms - squared_norms[i])
    end
    sort!(order, evals, rev=rev)
end


function columnsort!(::PostSNR, order::AbstractVector, W_zf::AbstractMatrix, N0; rev=true)
    Pn = [norm(w)^2 for w in eachcol(W)] .* N0
    sort!(order, Pn, rev=rev)
    order
end

function columnsort!(::PostSINR, order::AbstractVector, W_mmse::AbstractMatrix, H::AbstractMatrix; rev=true)
    Pin = Vector{Float64}(undef, size(H,2))
    for i in axes(H,2)
        μi = dot(W_mmse[:,i], H[:,i]) # 等価利得
        Pin[i] = μi - μi^2 # 等価雑音分散
    end
    sort!(order, Pin, rev=rev)
    order
end

function columnsort!(::Norm, order, H; rev=true)
    norms = [norm(h) for h in eachcol(H)]
    sort!(order, norms, rev=rev)
    order
end
