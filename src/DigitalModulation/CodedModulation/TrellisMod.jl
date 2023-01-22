# Trellis.jl
# トレリス符号化変調
module TCM

using Commun.DigiMod

struct TrellisMod

end

struct GenPoly

end

function get_genpoly()

end

struct APPDemod()
end

# softmax関数
function softmax!(x::AbstractVector{T}) where T <: Real
    c = max(x)
    exp.(x .- c) / sum(exp.(x .- c))
end

# BCJR復号
function _bcjr_calc(pos, pri, alpha, beta, gamma, operator)
    lp_vec = zeros();
    # forward calc
    for k in 1:trellis_size-1
        alpha_tmp = 0.0
        for j in 1:n_states
            for (i,l) in fwd_paths[j] # i 状態遷移における出力
                log_sum = alpha[i,k] + gamma[l,k] + pri[i,l]
                alpha_tmp = operator(alpha_tmp, log_sum)
            end
        end
    end
    # backward calc
    for k in trellis_size:-1:2
        for i in 1:n_states
            for (j,l) in bkwd_paths[j] # 状態iへのパス
                log_sum = beta[j,k] + gamma[l,k] + pri[i,l]
                beta_tmp = operator(beta_temp, log_sum)
            end
            beta[i,k-1] = beta_tmp;
        end
        # A Posteriori Probability
        if target == "input" # 入力シンボルに対する事後確率計算
            for j in 1:n_inputs
                for
                    lp_temp = operator(lp_temp, alpha+gamma+beta)
                end
                lp_vec[]
            end
            # scaling
            pos[:,k] .= softmax(lp_vec);
        else # 出力シンボルに対する事後確率計算

        end
    end
end





end # module
