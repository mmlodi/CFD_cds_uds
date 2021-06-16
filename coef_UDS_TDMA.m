function [a_n, a_s, a_e, a_w, a_p, s_u, s_p, max_peclet_u_UDS, max_peclet_v_UDS] = coef_UDS_TDMA( dx, dy, div_x, n_vol, c_p, k_ele, rho, T_N_vec, T_S_vec, T_E_vec, T_W_vec, angle, r_mat , omega)

    s_u = zeros( 1, n_vol);
    s_p = zeros( 1, n_vol); % SEM CONSIDERAR AS CONDIÇÔES DE CONTROLE
    a_w = zeros( 1, n_vol);      %
    a_e = zeros( 1, n_vol);      %
    a_n = zeros( 1, n_vol);
    a_s = zeros( 1, n_vol); 
    a_p = zeros( 1, n_vol);

    V_r = omega.*r_mat;

    V_u = -V_r.*sind(angle);

    V_v = V_r.*cosd(angle);
    
    Dx = (k_ele)/(c_p*dx);
    Dy = (k_ele)/(c_p*dy);
    Fx = rho*V_u;
    Fy = rho*V_v;
    
    peclet_u = Fx./Dx;
    peclet_v = Fy./Dy;
    
    max_peclet_u_UDS = max(max(abs(peclet_u)));
    max_peclet_v_UDS = max(max(abs(peclet_v)));


    for i = 1:n_vol % MONTA O CADA VETOR A_N, A_S, A_W, A_E, A_P, S_P e S_U PARA CADA ELEMENTO


        if mod(i,div_x) == 0 % MONTANDO VETORES PARA PAREDE DIREITA

            a_e(i) = 0;
            a_w(i) = + Dx;
            a_n(i) = + Dy;
            a_s(i) = + Fy(i) + Dy;

            s_u(i) =  + Fx(i)*T_E_vec(i) + (2*Dx)*T_E_vec(i); % SEM GERAÇÂO DE CALOR
            s_p(i) = - ( -Fx(i) +  2*Dx);

            % DEFINIÇÔES DE CONTORNO
            
            a_n(div_x) = 0; % CORREÇÂO PARA O PRIMEIRO ELEMENTO DA LINHA
            a_s(n_vol) = 0; % CORREÇÂO PARA O ULTIMO ELEMENTO 

            s_u(div_x) =  + (2*Dy)*T_N_vec(div_x) + Fx(div_x)*T_E_vec(div_x) + (2*Dx)*T_E_vec(div_x);% CORREÇÂO PARA O ULTIMO ELEMENTO
            s_u(n_vol) =  + Fy(n_vol)*T_S_vec(n_vol) + (2*Dy)*T_S_vec(n_vol) - Fx(n_vol)*T_E_vec(n_vol) + (2*Dx)*T_E_vec(n_vol);

            s_p(div_x) = -(  - Fx(div_x) + (2*Dy) + (2*Dx));% CORREÇÂO PARA O ULTIMO ELEMENTO
            s_p(n_vol) = -( - Fx(div_x) + Fy(n_vol)  + (2*Dy) + (2*Dx));

            a_p(i) =  a_e(i) + a_w(i) + a_n(i) + a_s(i) - s_p(i) ; 


        elseif mod(i,div_x) == 1 % MONTANDO VETORES PARA PAREDE ESQUERDA
            a_w(i) = 0;
            a_e(i) = - Fx(i) + Dx;
            a_n(i) = + Dy;
            a_s(i) = + Fy(i) + Dy;

            s_u(i) =  + Fx(i)*T_W_vec(i) + (2*Dx)*T_W_vec(i); 
            s_p(i) =  -( + Fx(i) + 2*Dx); 

            % DEFINIÇÔES DE CONTORNO

            a_n(1) = 0; % CORREÇÂO PARA O PRIMEIRO ELEMENTO
            a_s(n_vol-div_x+1) = 0; % CORREÇÂO PARA O ULTIMO ELEMENTO DA ULTIMA LINHA

            s_u(1) =    - Fy(1)*T_N_vec(1) + (2*Dy)*T_N_vec(1) + Fx(1)*T_W_vec(1) + (2*Dx)*T_W_vec(1);% CORREÇÂO PARA O ULTIMO ELEMENTO
            s_u(n_vol-div_x+1) =  + Fy(n_vol-div_x+1)*T_S_vec(n_vol-div_x+1) + (2*Dy)*T_S_vec(n_vol-div_x+1) + Fx(n_vol-div_x+1)*T_W_vec(n_vol-div_x+1) + (2*Dx)*T_W_vec(n_vol-div_x+1);

            s_p(1) = -( + Fx(1) - Fy(1) + (2*Dy) + (2*Dx)) ;% CORREÇÂO PARA O PRIMEIRO ELEMENT
            s_p(n_vol - div_x+1) = -( + Fx(n_vol-div_x+1) + Fy(n_vol-div_x+1) + (2*Dx) + (2*Dy) );% CORREÇÂO PARA O PRIMEIRO ELEMENTO DA ULTIMA LINHA

            a_p(i) =  a_e(i) + a_w(i) + a_n(i) + a_s(i) - s_p(i); 


        elseif (i > 1) && (i < div_x)%SUPERIOR SEM AS LATERAIS

            a_w(i) = + Dx; % F = rho V_u ou V_v
            a_e(i) = - Fx(i) + Dx; 
            a_n(i) = 0;
            a_s(i) = + Fy(i) + Dy;

            s_u(i) =  - Fy(i)*T_N_vec(i) + (2*Dy)*T_N_vec(i); % SEM GERAÇÂO DE CALOR
            s_p(i) = -( - Fy(i) + 2*Dy);

            a_p(i) =  a_e(i) + a_w(i) + a_n(i) + a_s(i) - s_p(i); 

        elseif (i > ( n_vol - div_x + 1 )) && (i < n_vol )%INFERIOR SEM AS LATERAIS

            a_w(i) = + Dx; % F = rho V_u ou V_v
            a_e(i) = - Fx(i) + Dx; 
            a_n(i) = + Dy;
            a_s(i) = 0;

            s_u(i) = + Fy(i)*T_S_vec(i) + (2*Dy)*T_S_vec(i); % SEM GERAÇÂO DE CALOR
            s_p(i) = -( + Fy(i) + 2*Dy);

            a_p(i) =  a_e(i) + a_w(i) + a_n(i) + a_s(i) - s_p(i);

        else

            a_w(i) = + Dx; % F = rho V_u ou V_v
            a_e(i) = - Fx(i) + Dx; 
            a_n(i) = + Dy;
            a_s(i) = + Fy(i) + Dy;

            s_u(i) = 0; % SEM GERAÇÂO DE CALOR
            s_p(i) = 0;

            a_p(i) =  a_e(i) + a_w(i) + a_n(i) + a_s(i) - s_p(i); 



         end
    end


end