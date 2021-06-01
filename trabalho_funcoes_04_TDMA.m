function [ T, ite, res_global] = TDMA( A, alpha, beta, D C_dot, tol, div_x, div_y, n_vol

    while (abs(res_global) > tol) % CODIGO DE RESOLUÇÂO PELO TDMA

        T_last = T;

        for k = 0:(div_y-1)% k REPRESENTA CADA LINHA, OU SEJA, POSSUI div_y linhas

            A(1 + k*div_x) = alpha(1 + k*div_x)/D(1 + k*div_x);
            C_dot(1 + k*div_x) = C(1 + k*div_x)/D(1 + k*div_x);

            for j = 2:div_x % FORWARD ELIMINATION

                A(j + k*div_x)     = alpha(j + k*div_x )/(D(j + k*div_x) - beta(j + k*div_x)*A((j-1) + k*div_x));
                C_dot(j + k*div_x) = (C(j + k*div_x) + beta(j + k*div_x)*C_dot((j-1)+ k*div_x))/(D(j+ k*div_x) - beta(j+ k*div_x)*A((j-1)+ k*div_x));

            end

            T(div_x + k*div_x) = C_dot(div_x + k*div_x);
           %T(k+1,div_x) = C_dot(div_x + k*div_x) % CASO T SEJA MATRIZ E NÃO VETOR

            for i = 1:(div_x - 1)% BACK SUBSTITUTION

                T((div_x - i)+ k*div_x) = A((div_x - i)+ k*div_x)*T((div_x + 1-i)+ k*div_x) + C_dot((div_x - i)+ k*div_x);
               %T(k+1, div_x - i ) = A((div_x - i)+ k*div_x)*T(k+1,div_x+1-i) + C_dot((div_x - i)+ k*div_x);
            end

            if ite == 1 % ESSE IF ATUALIZA OS VALOR DE C COM OS VIZINHOS DE CIMA E
                        % DE BAIXO PARA A PRIMEIRA VIZINHANÇA POR NA PRIMEIRA VEZ
                        % NÂO PRECISO ATUALIZAR A PRIMEIRA LINHA.

                for i = 1:(div_x) % AQUI ELE VARRE CADA LINHA 

                    if k < (div_x-1) % k se refere a linha  
                        if k == (div_x-2)
                            C( i + div_y + k*div_x) = a_n( i + div_y + k*div_x)*T( i + k*div_y) + a_s( i + div_y + k*div_y)*T_S_vec(i + div_y + k*div_x) + s_u( i + div_y + k*div_y);

                        else
                            C( i + div_y + k*div_x) = a_n( i + div_y + k*div_x)*T( i + k*div_y) + a_s( i + div_y + k*div_y)*T( i + div_y*2 + k*div_y) + s_u( i + div_y + k*div_y);

                        end
                    end
                end

            else
                for i = 1:(div_x) % ATUALIZA C PARA VALORES MAIORES QUE ITE > 1 

                    if k < (div_x) % ISSO AQUI É Q EU CHAMO DE CHUCHO kLKKKKKKKK
                        if k == 0
                            % AQUI TEM T_N MAS NÂO INFLUENCIA EM NADA POIS A_N = 0
                            C( i + k*div_x) = a_n( i + k*div_y)*T_N_vec( i + k*div_x) + a_s( i + k*div_y)*T( i + div_x +k*div_x) + s_u( i +  k*div_x);

                        elseif k == div_x-1

                            C( i + k*div_x) = a_n( i + k*div_x)*T( i - div_x + k*div_x) + a_s( i + k*div_x)*T_S_vec(i + k*div_x) + s_u( i + k*div_x);
                        else

                            C( i + k*div_x) = a_n( i + k*div_x)*T( i - div_x + k*div_x) + a_s( i + k*div_x)*T( i + div_x + k*div_x) + s_u( i + k*div_x);

                        end
                    end
                 end 
            end
        end

        % RESÍDUO PARA TDMA
        for i = 1:n_vol % MONTA O CADA VETOR A_N, A_S, A_W, A_E, A_P, S_P e S_ para cada elemento


            if mod(i,div_x) == 0 % MONTANDO VETORES PARA RESIDUO

                if i ~= div_x && i ~= n_vol
                res(i) = abs(a_n(i)*T( i - div_x) + a_s(i)*T(i + div_x) + a_w(i)*T(i - 1) + a_e(i)*T_E_vec(i) + s_u(i) - a_p(i)*T(i));
                end
                res(div_x) = abs(a_n(div_x)*T_N_vec(div_x) + a_s(div_x)*T(div_x + div_x) + a_w(div_x)*T(div_x - 1)+ a_e(div_x)*T_E_vec(div_x) + s_u(div_x) - a_p(div_x)*T(div_x));
                res(n_vol) = abs(a_n(n_vol)*T(n_vol - div_x) + a_s(n_vol)*T_S_vec(n_vol) + a_w(n_vol)*T(n_vol - 1)+ a_e(n_vol)*T_E_vec(n_vol) + s_u(n_vol) - a_p(n_vol)*T(n_vol));

            elseif mod(i,div_x) == 1 % MONTANDO VETORES PARA PAREDE ESQUERDA

                if i ~= 1 && i ~= (n_vol-div_x+1)
                res(i) = abs(a_n(i)*T(i - div_x) + a_s(i)*T(i + div_x) + a_w(i)*T_W_vec(i) + a_e(i)*T( i + 1) + s_u(i) - a_p(i)*T(i));
                end
                res(1) = abs(a_n(1)*T_N_vec(1) + a_s(1)*T(1 + div_x) + a_w(1)*T_W_vec(1) + a_e(1)*T(1+1)+ s_u(1) - a_p(1)*T(1));
                res(n_vol-div_x+1) = abs(a_n(n_vol-div_x+1)*T(n_vol - div_x - div_x + 1) + a_s(n_vol-div_x+1)*T_S_vec(n_vol-div_x+1) + a_w(n_vol-div_x+1)*T_W_vec(n_vol-div_x+1) + a_e(n_vol-div_x+1)*T(n_vol-div_x+1+1) + s_u(n_vol-div_x+1) - a_p(n_vol-div_x+1)*T(n_vol-div_x+1));

            elseif (i > 1) && (i < div_x)%SUPERIOR SEM AS LATERAIS

                res(i) = abs(a_n(i)*T_N_vec(i) + a_s(i)*T(i + div_x) + a_w(i)*T(i - 1) + a_e(i)*T( i + 1) + s_u(i) - a_p(i)*T(i));

            elseif (i > ( n_vol - div_x + 1 )) && (i < n_vol )%INFERIOR SEM AS LATERAIS

                res(i) = abs(a_n(i)*T(i - div_x) + a_s(i)*T_S_vec(i) + a_w(i)*T(i - 1) + a_e(i)*T( i + 1) + s_u(i) - a_p(i)*T(i));

            else % PARA OS NÓS CENTRAIS

                res(i) = abs(a_n(i)*T(i - div_x) + a_s(i)*T(i + div_x) + a_w(i)*T(i - 1) + a_e(i)*T( i + 1) + s_u(i) - a_p(i)*T(i));

            end
        end

        res_global = max(res);

        ite = ite + 1;         % ITERAÇÔES    

    end

end
