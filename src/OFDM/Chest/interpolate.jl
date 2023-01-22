# interpolate.jl

# method type
struct Linear end
struct Square end
struct Spline{N} end
struct Flat end

# オブジェクトの生成
function interpolator(name)
    if name == :linear
        Linear()
    elseif name == :square
        Square()
    elseif name == :flat
        Flat()
    else
    end
end


""" interp!(::Linear, y::AbstractArray, xq::Union{AbstractVector,AbstractRange}, exterp=true)

線形内挿補完

    arguments:

    returns:
"""
function interp!(y::AbstractArray, xq::Union{AbstractRange,AbstractVector}, ::Linear; exterp=true)
    @assert length(xq) > 1
    if typeof(xq) <:AbstractVector && !issorted(xq)
        sort!(xq)
    end
    @assert 1<= xq[1] <= length(y)
    @assert length(y) == size(y,1)
    length(xq)==1 && return

    for i = 1:length(xq)-1
        x1, x2 = xq[i], xq[i+1]
        y1, y2 = y[x1], y[x2]
        dy = y2 - y1
        dx = x2 - x1
        a = dy/dx     # 傾き
        b = y1 - a*x1 # 切片

        for x = x1+1:x2-1
            y[x] = a * x + b
        end
        if exterp && i == 1
            for x = 1:x1-1
                y[x] = a * x + b
            end
        elseif exterp && i == length(xq)-1
            for x = x2+1:length(y)
                y[x] = a * x + b
            end
        end
    end
    return y
end

interp(y, xq, l::Linear; exterp=true) = interp!(l, copy(y), xq, exterp=exterp)

# 二次補間
# xqはクエリ点の位置
function interp!(y::AbstractVector{T}, xq, ::Square; exterp=true) where T
    @assert length(xq) >= 3
    if typeof(xq) <:AbstractVector && !issorted(xq)
        sort!(xq)
    end

    for i = 1:length(xq)-2
        x1, x2, x3 = xq[i], xq[i+1], xq[i+2]
        y1, y2, y3 = y[x1], y[x2], y[x3]
        den1 = (x1-x2)*(x1-x3)
        den2 = (x2-x1)*(x2-x3)
        den3 = (x3-x1)*(x3-x2)
        a = 1/den1*y1 + 1/den2*y2 + 1/den3*y3
        b = (-(x2+x3)/den1)*y1 + (-(x1+x3)/den2)*y2 + (-(x1+x2)/den3)*y3
        c = (x2*x3)/den1*y1 + (x1*x3)/den2*y2 + (x1*x2)/den3*y3

        # 2次補間
        for x = x1+1:x2-1
            y[x] = a*x^2 + b*x + c
        end
        if i == 1
            for x = 1:x1
                @show x
                y[x] = a*x^2 + b*x + c
            end
        elseif i == length(xq)-2
            for x = x2+1:x3
                @show x
                y[x] = a*x^2 + b*x + c
            end
        end
    end
    return y
end
interp(y, xq, s::Square; exterp=true) = interp!(s, copy(y), xq, exterp=exterp)

# 3次スプライン補間
function interp!(y, xq, ::Spline{3}; exterp=true)
    length(xq) < 3 && error("補間には少なくとも3点必要です．")
    s = zeros(T, 3)
    A = zeros(T, 4)
    B = zeros(T, 4)

    for i in 1:length(xq)-2
        x_tmp = @views xq[i:i+2]
        y_tmp = @views y[x_tmp]
        h_tmp =
        for i in axes(A,1)
            for j in 1:i
                A[i,j] = h[i]
            end
        end
    end
end

# 端の値をコピーする外挿
function exterp!(y, idxs, ::Flat)
    l, r = idxs[1], idxs[end]
    for k = 1:l
        y[k] = y[l]
    end
    for k = r+1:length(y)
        y[k] = y[r]
    end
    y
end


# パイロット補間(時間方向)
function pinterp!(pilot::Pilot, CFRs, ::Time; method, exterp=method)
    size(CFRs,2)==1 && return CFRs
    order = pilot.order
    frame_len = size(CFRs,2)
    interpolator = select_interp_method(method)
    time_inds = get_time_idxs(pilot, size(CFRs,2))

    # time domain interpolation
    step = get_pilot_period(pilot)
    for q in axes(CFRs,3)
        for p in axes(CFRs,4)
            for id in 1:n_group
                st = @views t0 + findfirst(x->x==id, ord[p,:])
                knots = st : step : frame_len
                for i in pilot.gidxs[id]
                    @views interp!(CFRs[i,:,p,q], knots, interpolator, exterp=true) # 外挿ON
                end
            end
        end
    end
    CFRs
end

# 周波数方向
function pinterp!(pilot::Pilot, CFRs::AbstractArray{T,3}, ::Freq; method, exterp=method) where T
    fidxs = pilot.idxs
    for q in axes(CFRs,2)
        for p in axes(CFRs,3)
            flag = interp == exterp
            y = @views CFRs[:,p,q]
            @views interp!(y, fidxs, interpolator(method), exterp=flag)
            if !flag # 別の種法で外挿する
                @views exterp!(y, fidxs, interpolator(exterp))
            end
        end
    end
    CFRs
end
function pinterp!(pilot::Pilot, CFRs::AbstractArray{T,4}, ::Freq; method, exterp=method) where T
    table = make_time_table(pilot, size(CFRs,2))
    for k in axes(CFRs,2)
        for q in axes(CFRs,3)
            for p in axes(CFRs,4)
                gid = table[q,k] # group id
                fidxs = pilot.gidxs[gid]
                @views interp!(CFRs[:,p,q], fidxs, interpolator(method))
            end
        end
    end
end


# 3次スプライン補完
# S(x) = \sum_{j=0}^{N-1} Sj(x) = a_j(x-x_j)^3 + b_j(x-x_j)^2 + c_j(x-x_j) + d_j
function spline(xq, yq, x; exterp=true)
    a, b, c, d = _solve_equation(xq, yq) # 区分多項式の係数ベクトル

    y = Array{eltype(yq)}(undef, length(x))
    N = length(yq)-1
    i = 1; j = 1
    while i < N+1 && j <= length(x)
        if x[j] >= xq[i+1] && i < N # 区間を更新する条件
            i += 1; continue
        end
        if !exterp && (x[j] < xq[i] || x[j] > xq[i+1]) # 外挿しない場合
            j += 1; continue
        end
        y[j] = a[i]*(x[j]-xq[i])^3 + b[i]*(x[j]-xq[i])^2 + c[i]*(x[j]-xq[i]) + d[i]

        j += 1
    end
    x, y
end

# 連立式から各クエリ点の2次微分値を求め、係数a,b,c,dを求める
function _solve_equation(xq, yq)
    N = length(xq)-1
    A, b, h = _make_equation(xq,yq)
    u = A \ b # 各点の2階微分値
    pushfirst!(u, zero(eltype(yq))); push!(u, zero(eltype(yq)))
    a = [(u[i+1]-u[i])/(6*(xq[i+1]-xq[i])) for i in 1:N]
    b = u ./ 2
    c = [(yq[i+1]-yq[i])/(xq[i+1]-xq[i]) - 1/6 * (xq[i+1]-xq[i])*(u[i+1]+2u[i])  for i in 1:N]
    d = yq
    a, b, c, d
end

# 連立式を立てる
function _make_equation(xq, yq)
    N = length(xq)-1
    h = zeros(N)
    b = zeros(N-1)
    h = [xq[i+1]-xq[i] for i in 1:N]
    b = [6*((yq[i+1]-yq[i])/h[i] - (yq[i]-yq[i-1])/h[i-1])  for i in 2:N]
    A = zeros(eltype(xq), N-1, N-1)
    for i in axes(A,1)
        A[i,i] = 2*(h[i]+h[i+1]) # 対角成分
        if i > 1
            A[i,i-1] = h[i]
        end
        if i < size(A,1)
            A[i,i+1] = h[i+1]
        end
    end
    A, b, h
end


# インターフェース
function interpolate!(pilot::Pilot{Comb}, CFR; method=:linear)
    interpolate!(CFR, pilot, Time(), method)
    interpolate!(CFR, pilot, Freq(), method)
    CFR
end
function interpolate!(pilot::Pilot{Block}, CFR; method=:linear)
    interpolate!(CFR, pilot, Freq(), method)
    interpolate!(CFR, pilot, Time(), method)
    CFR
end
function interpolate!(pilot::Pilot{LTE}, CFR; method=:linear)
    interpolate!(CFR, pilot, Time(), method)
    interpolate!(CFR, pilot, Freq(), method)
    CFR
end
