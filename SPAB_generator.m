function u_spab = SPAB_generator(N, p, Te, u0, delta_u)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                          % Generator SPAB %%
                      % Last modified: 21.10.2025
                      % Author: Teme Marius
    %____________________________________________________________%


%% INPUT DATA  %%%%%%
clc;
disp('Generare SPAB')

% N = 11;
% p = 1;
% Te = 0.175;
% u0 = round((60*255)/100);
% delta_u = round((15*255)/100);
% u_spab = SPAB_generator(11, 1, 0.175, 60, 15);

u0 = round((u0*255)/100);
delta_u = round((delta_u*255)/100);

disp(' ')
disp(['N    ' num2str(N)])
disp(['p    ' num2str(p)])
disp(['Te   ' num2str(Te)])
disp(['u0   ' num2str(u0)])
disp(['du   ' num2str(delta_u)])

disp(' ')
disp(['u0 - delta_u     '  num2str(u0-delta_u)]);
disp(['u0 + delta_u     '  num2str(u0+delta_u)]); 


%% CALCUL SPAB  %%%%%%
L_SPAB = p*(2^N-1);
spab_init = idinput(L_SPAB,'prbs',[0 1/p],[u0-delta_u u0+delta_u]);

% Modificare spab_init pe front crescator
v1 = spab_init(1:N*p);
v2 = spab_init(N*p+1: end);
spab = [v2 ; v1];

spab_half = spab(1:floor(length(spab)/2));
u0_init = u0*ones(2*N*p,1);

% SPAB FINAL pt instalatie
u_spab = [u0_init ; spab ; spab_half];

t_SPAB = length(u_spab)*Te;
% total=length(u0_init)+length(spab)+length(spab_half);
total=length(u_spab);


disp(' ')
disp('Dimensiune u_spab')
disp('u0        L_SPAB         L_SPAB/2         Total'); 
disp([num2str(length(u0_init)) '        ' num2str(length(spab)) '                ' num2str(length(spab_half)) '                     '  num2str(total)]) 


disp(' ')
disp(' ')
disp(['Timp experiment     '  num2str(t_SPAB) '  sec']) 
disp(['Timp experiment     '  num2str(t_SPAB/60) '  min']) 

% plot(u_spab)
end