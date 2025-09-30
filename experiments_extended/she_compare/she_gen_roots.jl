using HomotopyContinuation
using MAT

function gray_code(N::Int)
    # https://cn.julialang.org/LeetCode.jl/dev/democards/problems/problems/89.gray-code/#:~:text=The%20gray%20code%20is%20a,sequence%20must%20begin%20with%200.
    powN = 1 << N
    res = Array{Int}(undef, powN)
    for i in 0:powN-1
        res[i + 1] = i ⊻ (i >> 1);
    end
    res
end


function get_path_qw(d::Int, Lmax=2)    
    #generate all qw-symmetric paths
    G = gray_code(d)
    mask =  digits.(G, base=2, pad=d)
    du = [[0; (2 .* g .- 1)] for g in mask]
    uall = cumsum.(du)

    u = [u for u in uall if all(u .>= 0) && all(u .<= Lmax)]
end

struct she_roots
    u
    M
    sols
    real_sols
    candidate
    alpha
end

function she_track(u, M)
    #solve the selective harmonics elimination problem
    # u: path from get_path_qw
    # M: modulation level for first fundamental frequency

    d = size(u)[1] - 1
    @var c[1:d]

    #form polynomial system (chebyshev polynomials)
    T = c*ones(1, 2*d+1)

    T[:, 1] = ones(d, 1)

    
    for i = 3:(2*d+1)
        T[:, i] = 2 .* c .* T[:, i-1, :] - T[:, i-2]
    end

    T_active = T[:, 2:2:end]
    du = diff(u)
    
    sumopt = (du' * T_active)'

    lhs = sumopt
    rhs = zeros(d, 1)
    rhs[1] = M * pi/4


    parameter_sampler = 10 .* randn(ComplexF64, d)
    ps = parameter_sampler
    F = System(lhs .- vec(rhs))
    result = solve(F)



    sols = solutions(result)
    real_sols = real_solutions(result)

    # filter_flag =  [all(diff(r) .< 0) && all(r .> 0) for r in real_sols]

    nonneg_flag =  [ all(diff(r) .< 0) for r in real_sols]
    one_flag =  [ all(r .< 1) for r in real_sols]
    decr_flag =  [ all(r .> 0) for r in real_sols]
    filter_flag = nonneg_flag .& decr_flag .& one_flag

    candidate = real_sols[filter_flag]

    if size(candidate)[1]!=0
        cv = candidate[1]
        sv = sqrt.(1 .- cv.^2)

        alpha = atan.(sv, cv)
    else
        alpha = []
    end

    she_roots(u, M, sols, real_sols, candidate, alpha)

end

# path_max = 12
# path_max = 2
# path_max = 3
# path_max = 3
path_max = 7
paths = [0.5 .* get_path_qw(i, 2) for i in 1:path_max]

# path_curr = paths[4]
M = 0.8
# M = 1

# root = she_track(path_curr[1], M)

root_list = []
root_list_all = []
for i = 1:path_max
    for j = 1:size(paths[i],1)
        print("u: ", paths[i][j], "\n")
        root_curr = she_track(paths[i][j], M)
        push!(root_list_all, root_curr)
        if size(root_curr.alpha,1)!=0
            push!(root_list, root_curr)
            print(root_curr.alpha, '\n')
        end
        
    end    
end

file = matopen("she_gen_08_1.mat", "w")
write(file, "root_list_all", root_list_all)
close(file)