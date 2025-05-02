function orbit = get_orbits(L)
            %find the orbits associated with the cyclic symmetry

            %L: a vector of size 3 x N
            %TODO: half-wave symmetry
            %start with full-wave symmetry
            
            P3 = [0 0 1; 1 0 0; 0 1 0];
            
            

            [k, N] = size(L);
            
            P = kron(eye(k/3), P3);
            
            L0 = L;
            L1 = P*L0;
            L2 = P*L1;

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

            Go = sparse(ind_i, ind_j, ones(1, length(ind_i)), N, N);

            sG = graph(Go);
            
            orbit = conncomp(sG);

        end