
# パンクチャの実行
function punct(input, P)
    period = length(P) # パンクチャの周期
    len = Int(ceil(length(input) * sum(P) / period)) # 推定出力長
    output = eltype(input)[]; sizehint!(output, len) # 出力配列

    for n in eachindex(input)
        if P[ind] > 0
            push!(output, input[n])
        end
        ind = (ind+1)%period
    end
    output
end

# デパンクチャ
function depunct(out::AbstractMatrix, inp::AbstractVector, P) where T
    p = length(P);
    len = 1; n = 1; l = 1;
    while l <= length(inp)
        if P[((n-1)%p)+1]
            out[n] = inp[l]
            l += 1
        end
        n += 1
    end
    out
end

# パックチャ行列取得
function get_punc_mtx(inout::Tuple{<:Integer, <:Integer}, rate::Rational{<:Integer})
    if inout[2]==2
    elseif inout[2]==3
        if rate == 1//2
            P = Bool.([1 1;
                       1 0;
                       0 1])
        else
            @error "現在サポートされていません"
        end
    end
    P
end
