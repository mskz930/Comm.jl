
# ２乗距離
function _squared_error(z, R, x)
    e = z - R * x |> r -> dot(r,r)
    return e
end

# 
function _back_tracking(x_tmp, cand_lists, refs, layer)
    n = length(cand_lists)

    # 候補点を変更する
    while layer > 0
        idx = n - layer + 1
        cand_list = cand_lists[idx]
        if isempty(cand_list)
            layer -= 1
        else
            x_tmp[idx] = refs[popfirst!(cand_list)]
            layer += 1
            break
        end
    end
    return layer
end

# 境界条件を計算する
function _calc_bounds(z, R, x_tmp, srad, layer)
    n = length(x_tmp)
    α = β = 0.0
    i = n - layer + 1 # current index

    # calc alpha
    α_tmp = 0.0
    for k in n:-1:i+1
        temp = z[k]
        for j in i:n
            temp -= R[k,j] * x_tmp[j]
        end
        α_tmp += abs2(temp)
    end
    α = sqrt(srad - α_tmp)

    # calc beta
    β = z[i]
    for j in i+1:n
        β -= R[i,j] * x_tmp[j]
    end
    
    sign_r = sign(R[i,i])
    upper_bound  =  (α + sign_r*β) / abs(R[i,i])
    lower_bound  = (-α + sign_r*β) / abs(R[i,i])

    lower_bound, upper_bound
end