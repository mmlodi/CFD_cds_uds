%UFSC - MECÂNICA DOS FLUIDOS COMPUTACIONAL
%Prof. Dr. Ernane
%Autor: Mateus Mischel Lodi
%Matrícula: 15159424
%%
% - - - - - - - - - HIPÓTESES - - - - - - - - - -

% REGIME PERMANENTE
% BIDIMENSIONAL
% SEM GERAÇÃO DE CALOR
% PAREDES COM TEMPERATURA CONSTANTE 
% ESQUEMA DE INTERPOLAÇÂO EM UDS 
% ESQUEMA DE INTERPOLAÇÂO EM CDS 

clear
%clc
%close all
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - - - - - - - - - - AQUISIÇÂO DE DADOS - - - - - - - - - - - -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 prompt = {'Raio parede externa [m]: ',...
           'Raio parede interna [m]: : ',...
           'Temperatura na parede externa [°C]: ',...
           'Temperatura na parede interna [°C]: ',... % NÂO FUNCIONA MUITO BEM PRA TEMPERATURA MAIORES EM T1
           'Número de volumes em x ou y [malha quadrada]: ',...
           'Cond. Térmica, k [W/mK]: ',...
           'Calor específ., cp [[J/kgK]]: ',...
           'densidade, rho [[kg/m^3]]: ',...
           'velocidade angular 1 [rad/s]: ',...
           'velocidade angular 2 [rad/s]: ',...
           'velocidade angular 3 [rad/s]: ',...
           'velocidade angular 4 [rad/s]: ',...
           'Tolerância p/ resíduo: ',...
           'Número máximo de iterações: '};
      
dlgtitle = 'Definindo as condições do problema ';

definput = {'0.05','0.02','300','100','5','60.5','434','7854','0','0.06','0.1','0.9','1e-5','20000'};
answer = inputdlg(prompt,dlgtitle,[1 40;1 40;1 40;1 40;1 40;1 40;1 40;1 40;1 40;1 40;1 40;1 40;1 40;1 40],definput);
r_e = str2num(answer{1});
r_i = str2num(answer{2});
T_e = str2num(answer{3});
T_i = str2num(answer{4});
n_vol = str2num(answer{5});
k_ele = str2num(answer{6});
c_p = str2num(answer{7});
rho = str2num(answer{8});
w1 = str2num(answer{9});
w2 = str2num(answer{10});
w3 = str2num(answer{11});
w4 = str2num(answer{12});
tol = str2num(answer{13});
max_ite = str2num(answer{14});

n_vol = n_vol^2;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - - - - - - - - - - DADOS INICIAIS - - - - - - - - - - - - - -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% T_i = 100;                % [°C]
% T_e = 300;                % [°C]
% 
% r_i = 0.02;               % [m]
% r_e = 0.05;               % [m]

% n_vol = 25;
L = r_e - r_i;
Lx = L/(sqrt(2));
Ly = Lx;

omega = [ w1 w2 w3 w4];     % [ rad/s]

div_x = sqrt(n_vol);
div_y = sqrt(n_vol);
dx = Lx/div_x;              % [m] 
dy = Ly/div_y;              % [m]

% k_ele = 60.5;             % [W/mK]
% c_p = 434;                % [J/kgK]
% rho = 7854;               % [kg/m^3]

%ite = 1;                   % É IMPORTANTE QUE COMECE DO 1 
res = zeros( 1, n_vol);     % resíduo começando em um pra entrar no laço while
res_global = 1;
% tol = 1e-5;
% max_ite = 50;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - - - - - - - INICIALIZANDO VETORES E MATRIZES - - - - - - - - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
angle = zeros(div_y,div_x);

x = (dx/2):dx:(Lx-dx/2);
y = (dy/2):dy:(Ly-dy/2);

T_ana = zeros(sqrt(n_vol),sqrt(n_vol));

T_N = zeros(1,div_x);
T_S = T_N;
T_E = T_N;
T_W = T_N;

T_W_mat = zeros(div_x,div_x);
T_N_mat = zeros(div_x,div_x);
T_E_mat = zeros(div_x,div_x);
T_S_mat = zeros(div_x,div_x);

T_W_vec = zeros(1,n_vol);
T_N_vec = zeros(1,n_vol);
T_E_vec = zeros(1,n_vol);
T_S_vec = zeros(1,n_vol);
 
T     = zeros( 1, n_vol);
T_plot = [];

r_mat = zeros( div_x, div_y);

r_mat_N = zeros( div_x, div_y);
r_mat_E = zeros( div_x, div_y);
r_mat_W = zeros( div_x, div_y);
r_mat_S = zeros( div_x, div_y);

T_r_ana = zeros(1,n_vol);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - - - - - - - ENCONTRANDO OS RAIOS PARA CADA ETAPA - - - - - - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SOFRI, MAS TERMINEI 
% MANO 

h =   0:(dx*sqrt(2)/2):((div_x-1)*(dx*sqrt(2)/2)) ;% ((div_x-2)*(dx*sqrt(2)/2)):-(dx*sqrt(2)/2):0 ]
l =   (r_i + dx*sqrt(2)/2) : (dx*sqrt(2)/2) :(r_e - dx*sqrt(2)/2) ;

for i = 0:(div_x-1)
    for j = (i+1):((2*(length(diag(r_mat,i))-1))+1*(i+1)) % ATENÇÂO PODE CAUSAR BUG DEVIDO A EXCESSIVA TEIMOSIA MINHA     
        if mod(i,2) == 0               
            r_mat( [(i+1):(div_x+1):(n_vol - i*div_x)] ) = sqrt( h(i+1)^2 + l( (i+1):2:((2*(length(diag(r_mat,i))-1))+1*(i+1)) ).^2 ) ; %(r_i - dx*sqrt(2)/2)  + j*dx*sqrt(2)
        else
            r_mat( [(i+1):(div_x+1):(n_vol - i*div_x)] ) = sqrt( h(i+1)^2 + l( (i+1):2:((2*(length(diag(r_mat,i))-1))+1*(i+1)) ).^2 ) ; %(r_i - dx*sqrt(2)/2)  + j*dx*sqrt(2)
        end
    end
end

% DEIXANDO A MATRIZ SIMETRICA

r_mat = r_mat + tril(r_mat,-1).';

% ROTACIONANDO A MATRIZ

r_mat = rot90(r_mat);
r = zeros(1,n_vol);
% CONVERTENDO MATRIZ PARA VETOR
for i = 0:(div_y-1) % TRANSFORMA O VETOR RESPOSTA EM UMA MATRIZ
    for j = 1:div_x
         r( j + div_x*i) = r_mat(i+1,j);
    end 
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - - - - - - - - - CÁLCULO DO VALOR ANALÍTICO - - - - - - - - -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:n_vol % variação dos valores de "x"

    T_r_ana(i) = T_e + (T_i - T_e)/log(r_i/r_e)*log(r(i)/r_e);

end

for i = 0:(div_y-1) % TRANSFORMA O VETOR RESPOSTA EM UMA MATRIZ
    for j = 1:div_x
        T_r_plot(i+1,j) = T_r_ana( j + div_x*i);
    end 
end

% APENAS ORGANIZANDO A MATRIZ PRA O PLOT
T_r_plot = rot90(rot90(rot90(T_r_plot)));
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - - - - - -  LOCALIZAÇÂO DOS RAIO DAS FRONTEIRAS - - - - - - -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PARA LATERAL ESQUERDA E DE BAIXO
h =   [0:(dx*sqrt(2)/2):((div_x-1)*(dx*sqrt(2)/2))] +  (dx*sqrt(2)/4) ;% ((div_x-2)*(dx*sqrt(2)/2)):-(dx*sqrt(2)/2):0 ]
ht =  [0:-(dx*sqrt(2)/2):-((div_x-1)*(dx*sqrt(2)/2))] +  (dx*sqrt(2)/4); 
l =   [(r_i + dx*sqrt(2)/2) : (dx*sqrt(2)/2) :(r_e - dx*sqrt(2)/2)] - (dx*sqrt(2)/4);

for i = 0:(div_x-1)
    
    for j = (i+1):((2*(length(diag(r_mat,i))-1))+1*(i+1)) 
       
        
        if mod(i,2) == 0               
            r_mat_W([(i+1):(div_x+1):(n_vol - i*div_x)] ) = sqrt( h(i+1)^2 + l( (i+1):2:((2*(length(diag(r_mat,i))-1))+1*(i+1)) ).^2 ) ; %(r_i - dx*sqrt(2)/2)  + j*dx*sqrt(2)
            if i > 0
                r_mat_W([(i*div_x+1):(div_x+1):(n_vol - i )]) = sqrt( ht(i)^2 + l( (i+1):2:((2*(length(diag(r_mat,i))-1))+1*(i+1)) ).^2 ) ; %(r_i - dx*sqrt(2)/2)  + j*dx*sqrt(2)
            end
     else
            r_mat_W([(i+1):(div_x+1):(n_vol - i*div_x)] ) = sqrt( h(i+1)^2 + l( (i+1):2:((2*(length(diag(r_mat,i))-1))+1*(i+1)) ).^2 ) ; %(r_i - dx*sqrt(2)/2)  + j*dx*sqrt(2)
            if i > 0
                r_mat_W([(i*div_x+1):(div_x+1):(n_vol - i )]) = sqrt( ht(i)^2 + l( (i+1):2:((2*(length(diag(r_mat,i))-1))+1*(i+1)) ).^2 ) ; %(r_i - dx*sqrt(2)/2)  + j*dx*sqrt(2)
            end
        end
    end
end
r_mat_S = r_mat_W';

r_left =  r_mat_W(:,1)'; 
r_bottom = r_left; % SIMETRIA

% PARA PARTE DE CIMA E DIREITA
h =   [0:(dx*sqrt(2)/2):((div_x-1)*(dx*sqrt(2)/2))] +  (dx*sqrt(2)/4) ;% ((div_x-2)*(dx*sqrt(2)/2)):-(dx*sqrt(2)/2):0 ]
ht =  [0:-(dx*sqrt(2)/2):-((div_x-1)*(dx*sqrt(2)/2))] +  (dx*sqrt(2)/4);
l =   [(r_i + dx*sqrt(2)/2) : (dx*sqrt(2)/2) :(r_e - dx*sqrt(2)/2)] + (dx*sqrt(2)/4);

for i = 0:(div_x-1)
    
    for j = (i+1):((2*(length(diag(r_mat,i))-1))+1*(i+1)) 
       
        if mod(i,2) == 0               
            r_mat_E( [(i+1):(div_x+1):(n_vol - i*div_x)] ) = sqrt( h(i+1)^2 + l( (i+1):2:((2*(length(diag(r_mat,i))-1))+1*(i+1)) ).^2 ) ; %(r_i - dx*sqrt(2)/2)  + j*dx*sqrt(2)
            if i > 0
                r_mat_E([(i*div_x+1):(div_x+1):(n_vol - i )]) = sqrt( ht(i)^2 + l( (i+1):2:((2*(length(diag(r_mat,i))-1))+1*(i+1)) ).^2 ) ; %(r_i - dx*sqrt(2)/2)  + j*dx*sqrt(2)
            end 
        else
            r_mat_E( [(i+1):(div_x+1):(n_vol - i*div_x)] ) = sqrt( h(i+1)^2 + l( (i+1):2:((2*(length(diag(r_mat,i))-1))+1*(i+1)) ).^2 ) ; %(r_i - dx*sqrt(2)/2)  + j*dx*sqrt(2)
            if i > 0
                r_mat_E([(i*div_x+1):(div_x+1):(n_vol - i )]) = sqrt( ht(i)^2 + l( (i+1):2:((2*(length(diag(r_mat,i))-1))+1*(i+1)) ).^2 ) ; %(r_i - dx*sqrt(2)/2)  + j*dx*sqrt(2)
            end        
        end
    end
end

r_mat_N = r_mat_E';

r_right =  r_mat_E(div_x,:); 
r_top = r_right; % SIMETRIA
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - - - ENCONTRANDO AS TEMP. PARA COND. DE CONTORNO - - - - - - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALCULANDO A TEMPERATURA NAS LATERAIS

for i = 1:div_x % variação dos valores de "x"

    T_N(i) = T_e + (T_i - T_e)/log(r_i/r_e)*log(r_top(i)/r_e);
    T_S(i) = T_e + (T_i - T_e)/log(r_i/r_e)*log(r_bottom(i)/r_e);
    T_W(i) = T_e + (T_i - T_e)/log(r_i/r_e)*log(r_left(i)/r_e);
    T_E(i) = T_e + (T_i - T_e)/log(r_i/r_e)*log(r_right(i)/r_e);
    
end

% INVERTENDO O VETOR T_E PARA OS NUMEROS FICAREM EM ORDEM
T_E = fliplr(T_E);
T_W = fliplr(T_W);

T_W_mat(:,1) = T_W
T_E_mat(:,div_x) = T_E
T_N_mat(1,:) = T_N
T_S_mat(div_y,:) = T_S

for i = 0:(div_y-1) % TRANSFORMA A MATRIZ RESPOSTA EM UM VETOR RESPOSTA
    for j = 1:div_x
         T_W_vec( j + div_x*i) = T_W_mat(i+1,j);
         T_E_vec( j + div_x*i) = T_E_mat(i+1,j);
         T_S_vec( j + div_x*i) = T_S_mat(i+1,j);
         T_N_vec( j + div_x*i) = T_N_mat(i+1,j);
    end 
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - - ENCONTRANDO ANGULOS PARA ENCONTRAR GRADIENTE DE V em x e y 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R_u = diag(rot90(rot90(rot90(r_mat)))); % Raio na diagonal princial da matriz

L_x = R_u/sqrt(2); % Projeção da diagonal prinpal dos raios em x

for i = 1:div_y
    for j = 1:div_x
    
    angle(i,j) = acosd(L_x(j)/r_mat(i,j)); % EM "graus" BASEADO -> alpha= acosd( CA/HIP) 
    
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - - ENCONTRANDO OS COEF. AN, AS, AE, AW, SU, SP ETC - - - - -  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CALCULA OS COEFICIENTES DO CDS PARA TDMA 4 RESULTADOS;

[a_n_uds_1, a_s_uds_1, a_e_uds_1, a_w_uds_1, a_p_uds_1, s_u_uds_1, s_p_uds_1, max_peclet_u_UDS_1, max_peclet_v_UDS_1] = coef_UDS_TDMA( dx, dy, div_x, n_vol, c_p, k_ele, rho, T_N_vec, T_S_vec,T_E_vec,T_W_vec, angle, r_mat, omega(1));
[a_n_uds_2, a_s_uds_2, a_e_uds_2, a_w_uds_2, a_p_uds_2, s_u_uds_2, s_p_uds_2, max_peclet_u_UDS_2, max_peclet_v_UDS_2] = coef_UDS_TDMA( dx, dy, div_x, n_vol, c_p, k_ele, rho, T_N_vec, T_S_vec,T_E_vec,T_W_vec, angle, r_mat, omega(2));
[a_n_uds_3, a_s_uds_3, a_e_uds_3, a_w_uds_3, a_p_uds_3, s_u_uds_3, s_p_uds_3, max_peclet_u_UDS_3, max_peclet_v_UDS_3] = coef_UDS_TDMA( dx, dy, div_x, n_vol, c_p, k_ele, rho, T_N_vec, T_S_vec,T_E_vec,T_W_vec, angle, r_mat, omega(3));
[a_n_uds_4, a_s_uds_4, a_e_uds_4, a_w_uds_4, a_p_uds_4, s_u_uds_4, s_p_uds_4, max_peclet_u_UDS_4, max_peclet_v_UDS_4] = coef_UDS_TDMA( dx, dy, div_x, n_vol, c_p, k_ele, rho, T_N_vec, T_S_vec,T_E_vec,T_W_vec, angle, r_mat, omega(4));

%CALCULA OS COEFICIENTES DO CDS PARA TDMA 4 RESULTADOS;
                                                                                               
[a_n_cds_1, a_s_cds_1, a_e_cds_1, a_w_cds_1, a_p_cds_1, s_u_cds_1, s_p_cds_1, max_peclet_u_CDS_1, max_peclet_v_CDS_1] = coef_CDS_TDMA( dx, dy, div_x, n_vol, c_p, k_ele, rho, T_N_vec, T_S_vec, T_E_vec, T_W_vec, angle, r_mat, omega(1));
[a_n_cds_2, a_s_cds_2, a_e_cds_2, a_w_cds_2, a_p_cds_2, s_u_cds_2, s_p_cds_2, max_peclet_u_CDS_2, max_peclet_v_CDS_2] = coef_CDS_TDMA( dx, dy, div_x, n_vol, c_p, k_ele, rho, T_N_vec, T_S_vec, T_E_vec, T_W_vec, angle, r_mat, omega(2));
[a_n_cds_3, a_s_cds_3, a_e_cds_3, a_w_cds_3, a_p_cds_3, s_u_cds_3, s_p_cds_3, max_peclet_u_CDS_3, max_peclet_v_CDS_3] = coef_CDS_TDMA( dx, dy, div_x, n_vol, c_p, k_ele, rho, T_N_vec, T_S_vec, T_E_vec, T_W_vec, angle, r_mat, omega(3));
[a_n_cds_4, a_s_cds_4, a_e_cds_4, a_w_cds_4, a_p_cds_4, s_u_cds_4, s_p_cds_4, max_peclet_u_CDS_4, max_peclet_v_CDS_4] = coef_CDS_TDMA( dx, dy, div_x, n_vol, c_p, k_ele, rho, T_N_vec, T_S_vec, T_E_vec, T_W_vec, angle, r_mat, omega(4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - - - - - - - - - - APLICÃO DO MÉTODO TDMA - - - - - - - - - -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEFINIÇÂO DOS VALORES PARA O METODO TDMA PARA O UDS

[ T_uds_1, ite_uds_1, res_global_uds_1] = TDMA( a_n_uds_1, a_s_uds_1, a_e_uds_1, a_w_uds_1, a_p_uds_1, s_u_uds_1, tol, div_x, div_y, n_vol, T, T_N_vec, T_S_vec,T_E_vec,T_W_vec, max_ite);
[ T_uds_2, ite_uds_2, res_global_uds_2] = TDMA( a_n_uds_2, a_s_uds_2, a_e_uds_2, a_w_uds_2, a_p_uds_2, s_u_uds_2, tol, div_x, div_y, n_vol, T, T_N_vec, T_S_vec,T_E_vec,T_W_vec, max_ite);
[ T_uds_3, ite_uds_3, res_global_uds_3] = TDMA( a_n_uds_3, a_s_uds_3, a_e_uds_3, a_w_uds_3, a_p_uds_3, s_u_uds_3, tol, div_x, div_y, n_vol, T, T_N_vec, T_S_vec,T_E_vec,T_W_vec, max_ite);
[ T_uds_4, ite_uds_4, res_global_uds_4] = TDMA( a_n_uds_4, a_s_uds_4, a_e_uds_4, a_w_uds_4, a_p_uds_4, s_u_uds_4, tol, div_x, div_y, n_vol, T, T_N_vec, T_S_vec,T_E_vec,T_W_vec, max_ite);

% DEFINIÇÂO DOS VALORES PARA O METODO TDMA PARA O CDS

[ T_cds_1, ite_cds_1, res_global_cds_1] = TDMA(a_n_cds_1, a_s_cds_1, a_e_cds_1, a_w_cds_1, a_p_cds_1, s_u_cds_1, tol, div_x, div_y, n_vol, T, T_N_vec, T_S_vec,T_E_vec,T_W_vec, max_ite);
[ T_cds_2, ite_cds_2, res_global_cds_2] = TDMA(a_n_cds_2, a_s_cds_2, a_e_cds_2, a_w_cds_2, a_p_cds_2, s_u_cds_2, tol, div_x, div_y, n_vol, T, T_N_vec, T_S_vec,T_E_vec,T_W_vec, max_ite);
[ T_cds_3, ite_cds_3, res_global_cds_3] = TDMA(a_n_cds_3, a_s_cds_3, a_e_cds_3, a_w_cds_3, a_p_cds_3, s_u_cds_3, tol, div_x, div_y, n_vol, T, T_N_vec, T_S_vec,T_E_vec,T_W_vec, max_ite);
[ T_cds_4, ite_cds_4, res_global_cds_4] = TDMA(a_n_cds_4, a_s_cds_4, a_e_cds_4, a_w_cds_4, a_p_cds_4, s_u_cds_4, tol, div_x, div_y, n_vol, T, T_N_vec, T_S_vec,T_E_vec,T_W_vec, max_ite);

% APENAS ORGANIZANDO A MATRIZ PRA O PLOT
T_uds_1 = rot90(rot90(rot90(T_uds_1)));
T_uds_2 = rot90(rot90(rot90(T_uds_2)));
T_uds_3 = rot90(rot90(rot90(T_uds_3)));
T_uds_4 = rot90(rot90(rot90(T_uds_4)));

T_cds_1 = rot90(rot90(rot90(T_cds_1)));
T_cds_2 = rot90(rot90(rot90(T_cds_2)));
T_cds_3 = rot90(rot90(rot90(T_cds_3)));
T_cds_4 = rot90(rot90(rot90(T_cds_4)));

% % % - - - - - - - -  CALCULANDO O VALOR ANALÍTICO  - - - - - - - - - - - -



T_uds = zeros(div_x,div_y,length(omega));
T_uds(:,:,1) =  T_uds_1;
T_uds(:,:,2) =  T_uds_2;
T_uds(:,:,3) =  T_uds_3;
T_uds(:,:,4) =  T_uds_4;

T_cds = zeros(div_x,div_y,length(omega));
T_cds(:,:,1) =  T_cds_1;
T_cds(:,:,2) =  T_cds_2;
T_cds(:,:,3) =  T_cds_3;
T_cds(:,:,4) =  T_cds_4;

%PECLES SÃO IGUAIS TANTO PARA O UDS QUANTO CDS

max_peclet_UDS = [ max_peclet_u_UDS_1, max_peclet_v_UDS_1;
                    max_peclet_u_UDS_2, max_peclet_v_UDS_2;
                    max_peclet_u_UDS_3, max_peclet_v_UDS_3;
                    max_peclet_u_UDS_4, max_peclet_v_UDS_4];
                
max_peclet_CDS      = [ max_peclet_u_CDS_1, max_peclet_v_CDS_1;
                    max_peclet_u_CDS_2, max_peclet_v_CDS_2;
                    max_peclet_u_CDS_3, max_peclet_v_CDS_3;
                    max_peclet_u_CDS_4, max_peclet_v_CDS_4];
                
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - - - - - - - - - - PLOTS DE GRÁFICOS - - - - - - - - - - - - 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%------------------------------------ PLOT ANALÍTICA
figure('Name','Resultado analítico','NumberTitle','off');
contourf(T_r_plot,70)
legend('ANALÍTICO');
hold on
contour(T_r_plot,70)
hold off
axis equal
ylabel('Volumes em y')
xlabel('Volumes em x')
c.LineWidth = 1;
colormap(hot)
colorbar
caxis([100 300])

figure('Name','Gráficos para UDS','NumberTitle','off');

%------------------------------------ PLOT UDS
for i = 1:length(omega)

    subplot(2, 2, i);
    contourf(T_uds(:,:,i),70)
    hold on
    contour(T_uds(:,:,i),70)
    hold off
    axis equal
    title(['w = ',num2str(omega(i)),' [rad/s]']);
    xlabel('Volumes em X');
    ylabel('Volumes em Y');
    colormap(hot)
    colorbar
    caxis([100 300])
      
end

figure('Name','Gráficos para CDS','NumberTitle','off');

%------------------------------------ PLOT CDS
for i = 1:length(omega)

    subplot(2, 2, i);
    contourf(T_cds(:,:,i),70)
    hold on
    contour(T_cds(:,:,i),70)
    hold off
    axis equal
    title(['w = ',num2str(omega(i)),' [rad/s]']);
    xlabel('Volumes em X');
    ylabel('Volumes em Y');
    colormap(hot)
    colorbar
    caxis([100 300])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - - - - - - IMPRESSÃO E FORMATAÇÃO NO COMMAND WINDOW - - - - -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

fprintf('\n')
disp('-------------- RESULTADO ---------------------------------')
fprintf('                               Velocidade angular (rad/s)    \n');
fprintf('                                %.2f | %.2f | %.2f | %.2f    \n', omega);
%fprintf('\n')
disp('-----------UDS--------------------------------------------')
fprintf('Resíduo máximo  [10^(-5)]:      %.2f | %.2f | %.2f | %.2f     \n', res_global_uds_1*10^5, res_global_uds_2*10^5, res_global_uds_3*10^5,res_global_uds_4*10^5);
fprintf('Nº de Iterações :                %.f  |  %.f  |  %.f  |  %.f    \n', ite_uds_1 , ite_uds_2, ite_uds_3, ite_uds_4  );
fprintf('Máximo Peclet abs. entre u e v: %.2f | %.2f | %.2f | %.2f     \n', max(max_peclet_UDS(1,:)), max(max_peclet_UDS(2,:)), max(max_peclet_UDS(3,:), max(max_peclet_UDS(4,:) )));
fprintf('\n')
disp('-----------CDS--------------------------------------------')
fprintf('Resíduo máximo  [10^(-5)]:      %.2f | %.2f | %.2f | %.2f     \n', res_global_cds_1*10^5, res_global_cds_2*10^5, res_global_cds_3*10^5,res_global_cds_4*10^5);
fprintf('Nº de Iterações :                %.f  |  %.f  |  %.f  |  %.f    \n', ite_cds_1 , ite_cds_2, ite_cds_3, ite_cds_4  );
fprintf('Máximo Peclet abs. entre u e v: %.2f | %.2f | %.2f | %.2f     \n', max(max_peclet_CDS(1,:)), max(max_peclet_CDS(2,:)), max(max_peclet_CDS(3,:)), max(max_peclet_CDS(4,:)));
disp('----------------------------------------------------------')