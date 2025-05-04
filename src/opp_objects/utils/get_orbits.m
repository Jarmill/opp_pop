function orbits = get_orbits(L)
            %find the orbits associated with the cyclic symmetry

            %L: a vector of size 3d x N for some integer d>=1
            %
            %
            %output:
            %   orbits:             structure storing the cyclic orbits
            %       id:         the identifier/index of the orbit
            %       ind_std:    the index in L of the orbit
            %       ind_std:    the index in -L of the orbit (for HW/QW)
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

            %collect together the indices describing the orbits
            ind_std = zeros(3, N);
            ind_flip = zeros(3, N);
            ind_std(1, :) = 1:N;
            
            for i = 1:N
                %standard indices
                l1 = find(all(L1 == L0(:, i)));
                l2 = find(all(L2 == L0(:, i)));
                
                if isempty(l1)
                    l1 = NaN;
                end
                if isempty(l2)
                    l2 = NaN;
                end

                ind_std(2:3, i) = [l1; l2];
                % G([i, i], [l1, l2]) = 1;
                % G([l1, l2], [i, i]) = 1;
                % ind_i = [ind_i, i, i, l1, l2];
                % ind_j = [ind_j, l1, l2, i, i];

                %flipped indices
                %slightly inefficient to do this flipping, but whatever
                %this isn't the bottleneck
                lf0 = find(all(L0 == -L0(:, i)));
                lf1 = find(all(L1 == -L0(:, i)));
                lf2 = find(all(L2 == -L0(:, i)));
                if isempty(lf0)
                    lf0 = NaN;
                end
                if isempty(lf1)
                    lf1 = NaN;
                end
                if isempty(lf2)
                    lf2 = NaN;
                end
                ind_flip(:, i) = [lf0; lf1; lf2];
            end

            
            % Go = sparse(ind_i, ind_j, ones(1, length(ind_i)), N, N);
            % 
            % sG = graph(Go);

            orbits = struct;
            
           
            % orbits.id = conncomp(sG);
            orbits.ind_std = ind_std;
            orbits.ind_flip = ind_flip;

        end