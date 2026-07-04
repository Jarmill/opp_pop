function [n, bounds, types] = harm_screen(opts)

            %screen the harmonics constraints for nontrivial terms (given
            %the symmetry configuration)

            Sym = opts.Symmetry;
            if Sym==1 && opts.quarter_match
                Sym = 2;
            end

            % Sym_Scale = double(2^(-Sym));
            Sym_Scale = 1;
            n = [];
            bounds = [];
            types = [];

            harm = opts.harmonics;

            %process the cosine harmonics
            for i = 1:length(harm.index_cos)
                
                ncurr = harm.index_cos(i);
                bndcurr = harm.bound_cos(i, :);
                
                keep = 1;
                if Sym==1
                    %half-wave
                    keep = mod(ncurr, 2)==1;

                elseif Sym==2
                    %quarter-wave
                    keep = 0;
                end

                if keep
                    n = [n; ncurr];
                    types = [types; 0];
                    bounds = [bounds; Sym_Scale*bndcurr]; %TODO: check dimensions here
                end
            end
            
            %process the sine harmonics
            for i = 1:length(harm.index_sin)

                ncurr = harm.index_sin(i);
                bndcurr = harm.bound_sin(i, :);

                keep = 1;
                if Sym>0
                    %half-wave and quarter-wave
                    keep = mod(ncurr, 2)==1;
                end

                if keep
                    n = [n; ncurr];
                    types = [types; 1];
                    bounds = [bounds; Sym_Scale*bndcurr]; %TODO: check dimensions here
                end
            end

        end