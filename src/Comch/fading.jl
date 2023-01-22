# 定数
const PLANE_WAVES = 34 # 素波数
const PHASE_POINTS = let N = PLANE_WAVES
    div(div(N,2)-1,2) # points(位相点の数)
end
const jakes_params = let N=PLANE_WAVES, N0=PHASE_POINTS
    (
        cosθ = [cos(2π*i/N) for i in 1:N0],
        cosϕ = [cos(π*i/(N0+1)) for i in 1:N0],
        sinϕ = [sin(π*i/(N0+1)) for i in 1:N0],
    )
end

mutable struct Jakes end

"""
jakesモデルによるフェージング系列の生成
"""
function gen_fading!(output, P0, fdTs=0.0, sps=1, st::Symbol=:rand, by::Jakes=Jakes())
    N    = PLANE_WAVES
    N0   = PHASE_POINTS
    cosθ = jakes_params[:cosθ]    # cosθ
    cosϕ = jakes_params[:cosϕ]    # cosϕ
    sinϕ = jakes_params[:sinϕ]    # sinϕ
    df = fdTs/sps                 # １サンプルあたりの周波数偏移
    wd = 2π*df                    # ドップラー各周波数

    # 系列の生成
    ncoef = size(output, 1)
    t0 = st==:rand ? 1000*rand(ncoef) : zeros(ncoef) # 初期時刻
    t = cos_wn_t = 0.0
    val = complex(0.0)
    amp = @. sqrt(P0)/sqrt(2*N0+1)
    @inbounds for n in axes(output,2)
        for i in axes(output,1)
            t = t0[i] + (n-1)*wd
            val = sqrt(2.0)*cos(t) + 0.0im
            for n in 1:N0
                cos_wn_t = cos(cosθ[n]*t)
                val += 2.0*cos_wn_t*(cosϕ[n] + im*sinϕ[n])
            end
            output[i,n] = amp[i] * val
        end
    end
    output
end

"""
フェージング系列の生成
    gen_fading(pdp, len, Nr, Nt)
"""
function gen_fading(pdp, len, Nr, Nt; fdTs=0.001, sps=1, istail=true, by=Jakes())
    if fdTs == 0.0
        chresp = sqrt.(pdp) .* randn(ComplexF64, length(pdp), Nr, Nt)
    else
        # 出力配列確保
        inds = PDP.nzinds(pdp)
        chresp = zeros(ComplexF64, length(pdp), len, Nr, Nt)

        # フェージング生成
        for (m,n) in Iterators.product(1:Nt, 1:Nr)
            @views gen_fading!(chresp[inds,:,n,m], pdp[inds], fdTs, sps, :rand)
        end
    end
    chresp
end







#==================== Jakes型の試験 =====================#
# Jakes型
#=
struct Jakes
    N::Int64                # 素波数
    N0::Int64               # 位相点の数
    fdTs::Float64           # 正規化ドップラー周波数
    sps::Int64              # 単位時間あたりのサンプル数
    df::Float64             # 1サンプルあたりの周波数偏移
    wd::Float64             # 1サンプルあたりの各周波数偏移
end
=#

# Jakesモデルによるフェージング乱数系列の生成
function (jakes::Jakes)(fading_out, P0, st=:rand)
    N, N0 = jakes.N, jakes.N0
    cosθ = [cos(2π*i/N) for i in 1:N0]
    cosϕ = [cos(π*i/(N0+1)) for i in 1:N0]
    sinϕ = [sin(π*i/(N0+1)) for i in 1:N0]
    t = cos_wn_t = 0.0
    val = complex(0.0)
    amp = sqrt(P0)/sqrt(2*N0+1)        # 絶対振幅の平均
    t0 = ifelse(st==:rand, 1000*rand(), 0.0) # 初期時刻

    for m in eachindex(fading_out)
        t = t0 + (m-1)*jakes.wd
        cos_wn_t = cos(t)
        val = sqrt(2.0)*cos_wn_t + 0.0im
        for n in 1:N0
            cos_wn_t = cos(cosθ[n]*t)
            val += 2.0*cos_wn_t*(cosϕ[n] + im*sinϕ[n])
        end
        fading_out[m] = amp * val
    end
end
