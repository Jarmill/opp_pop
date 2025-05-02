L =1;

va = kron(kron(-L:L, ones(1, 2*L+1)), ones(1, 2*L+1));
vb = kron(kron(ones(1, 2*L+1), -L:L), ones(1, 2*L+1));
vc = kron(kron(ones(1, 2*L+1), ones(1, 2*L+1)), -L:L);

V = [va; vb; vc];

P3 = [0 0 1; 1 0 0; 0 1 0];
            
L0 = V;
L1 = P3*L0;
L2 = P3*L1;

N = length(va);

G = speye(N, N);
ind_i = 1:N;
ind_j = 1:N;

for i = 1:N
    l1 = find(all(L1 == L0(:, i)));
    l2 = find(all(L2 == L0(:, i)));
    % G([i, i], [l1, l2]) = 1;
    % G([l1, l2], [i, i]) = 1;
    ind_i = [ind_i, i, i, l1, l2];
    ind_j = [ind_j, l1, l2, i, i];
end

G = sparse(ind_i, ind_j, ones(1, length(ind_i)), N, N);

sG = graph(G);

c = conncomp(sG);
figure(4)
clf
plot(sG)
