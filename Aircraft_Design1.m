clc
clear
close all;
tic
%% Finding A & B , ploting Regression Line
W_E_e = xlsread('AD', 1, 'D2:D27');
W_TO_to = xlsread('AD', 1, 'E2:E27');
num = xlsread('AD', 1, 'A2:A27');
show = [num, W_E_e, W_TO_to];
%%
W_E = log10(W_E_e)';
W_TO = log10(W_TO_to)';
%%
[R] = regression(W_E, W_TO);  
format long
A = 0.108637;
B = 1.050745;
fprintf("A: %d\n",A)
fprintf("B: %d\n",B)
fprintf("R: %d\n",R)
figure, p = plotregression(W_E, W_TO);
xlim([3.5 4.5]);
ylim([3.7 4.53]);
legend off
title("R^2: ",R);
xlabel("Take off Weight");
ylabel("Empty Weight");
%% Error Estimating
Error = 1;
W_TO = 0;
for W_TO_guess= 1000:0.01:60000
W_E1 = regAns(W_TO_guess);
W_E2 = berguetAns(W_TO_guess);
delta = abs(W_E1 - W_E2);
if W_E1 > W_E2
    x = delta / W_E1;
else
    x = delta / W_E2;
end
if x < 0.005
    if x < Error
        Error = x;
        W_TO = W_TO_guess;
    end
end
end
fprintf('\nError: %d\n',Error) ;
fprintf('Max Take-off Weight: %d\n',W_TO) ;
fprintf('Emplty Weight: %d\n',berguetAns(W_TO));
%% Sensitivities
% Defining Variables
syms W_TO y W_PL dW_TO dy Sen_W_PL Sen_Range Sen_Endurance Sen_Speed Sen_SFC Sen_L_Over_D
M_ff = 0.825028803329556;
M_res = 0;
M_tfo = 0.005;
W_PL = 4175;
W_crew = 525;
R = 805.54561361645;
E = 0.8;
C_p_cr = 0.5;
C_p_ltr = 0.5;
L_over_D_cr = 11; 
L_over_D_ltr = 10;
V = 207.1403;
etha_p_cr = 0.82;
etha_p_ltr = 0.72;
A = 0.108637;
B = 1.050745;
C_1 = (1 - (1 + M_res) * (1 - M_ff) - M_tfo);
C_2 = (M_ff * (1 + M_res) - M_tfo -M_res);
D = (W_PL + W_crew);
W_TO = vpasolve(A + B * log10(C_1 * W_TO - D) - log10(W_TO)==0,W_TO);
F = - B * (W_TO ^ 2) * (1 / (C_2 * W_TO * (1 - B) - D)) * (1 + M_res) * M_ff;
%%
% Calculating Sensitivities
Sen_W_PL = (B * W_TO)/((D - C_1 *(1 - B) * W_TO));
Sen_W_E = (B * W_TO) / regAns(W_TO);
Sen_Range = (F * C_p_cr) / (375 * etha_p_cr * L_over_D_cr);
Sen_Endurance = F * V * C_p_ltr * (1 / (375 * etha_p_ltr * L_over_D_ltr));
Sen_C_p_cr = F * R * (1 / (375 * etha_p_cr * L_over_D_cr));
Sen_C_p_ltr = (F * E * V) / (375 * etha_p_ltr * L_over_D_ltr);
Sen_Speed = (F * E * C_p_ltr) / (375 * etha_p_ltr * L_over_D_cr);
Sen_etha_p_cr = F * -R * C_p_cr * (1 / (375 * (etha_p_cr ^ 2) * L_over_D_cr));
Sen_etha_p_ltr = F * -E * V * C_p_ltr * (1 / (375 * (etha_p_ltr ^ 2) * L_over_D_ltr));
Sen_L_Over_D_cr = F * -R * C_p_cr * (1 / (375 * etha_p_cr * (L_over_D_cr ^ 2)));
Sen_L_Over_D_ltr = F * -E * V * C_p_ltr * (1 / (375 * etha_p_ltr * (L_over_D_ltr ^ 2)));
%% 
% Tabeling The Results
Name = ["Sen_W_PL", "Sen_W_E","Sen_Range","Sen_Endurance","Sen_C_p_cr","Sen_C_p_ltr","Sen_Speed","Sen_etha_p_cr","Sen_etha_p_ltr","Sen_L_Over_D_cr","Sen_L_Over_D_ltr"]';
Value = [Sen_W_PL, Sen_W_E, Sen_Range, Sen_Endurance, Sen_C_p_cr, Sen_C_p_ltr, Sen_Speed, Sen_etha_p_cr, Sen_etha_p_ltr, Sen_L_Over_D_cr, Sen_L_Over_D_ltr]';
Value = double(Value);
format long
Sensitivities = table(Name, Value);
disp(Sensitivities);
%% Stall Sizing
V_S_clean = 124.8979294254888; % Unit: ft/s
V_S_TOandL = 101.26859142607199; % Unit: ft/s
rho_cruise = 0.00175549; % Unit: slug/ft^3
rho_TOandL = 0.00200335; % Unit: slug/ft^3
C_L_max_clean = [1.2, 1.4, 1.6, 1.8];
C_L_max_takeoff = [1.4, 1.6, 1.8, 2.0];
C_L_max_landing = [1.6, 1.8, 2.0, 2.2, 2.5];
disp("                           ""For Cruise Phase""         ");

subplot(3,1,1), hold on
for C_L = C_L_max_clean
W_over_S = 0.5 * rho_cruise * (V_S_clean ^ 2) * C_L;
fprintf('For C_L_max_clean = %d    ---->    W_over_S = %d\n',C_L,W_over_S);
xline(W_over_S, "Color",rand(1,3),"LineWidth",3,"Label",C_L);
end
title("\fontsize{14}\fontname{Times}Stall Speed Sizing for V_s Clean (Flaps up)");
xlim([10, 40]);
xlabel("\fontsize{14}\fontname{Times}W/S [psf]");
ylabel("\fontsize{14}\fontname{Times}W/P [{shp}/{lb}]");
disp("                           ""For Take off Phase""         ");

subplot(3,1,2), hold on
for C_L = C_L_max_takeoff
W_over_S = 0.5 * rho_TOandL * (V_S_clean ^ 2) * C_L;
fprintf('For C_L_max_Takeoff = %d    ---->    W_over_S = %d\n',C_L,W_over_S);
xline(W_over_S, "Color",rand(1,3),"LineWidth",3,"Label",C_L);
end
title("\fontsize{14}\fontname{Times}Stall Speed Sizing for V_s Take off (Flaps Down)");
xlim([10, 40]);
xlabel("\fontsize{14}\fontname{Times}W/S [psf]");
ylabel("\fontsize{14}\fontname{Times}W/P [{shp}/{lb}]");
disp("                           ""For Landing Phase""         ");

subplot(3,1,3), hold on
for C_L = C_L_max_landing
W_over_S = 0.5 * rho_TOandL * (V_S_clean ^ 2) * C_L;
fprintf('For C_L_max_Landing = %d    ---->    W_over_S = %d\n',C_L,W_over_S);
xline(W_over_S, "Color",rand(1,3),"LineWidth",3,"Label",C_L);
end
title("\fontsize{14}\fontname{Times}Stall Speed Sizing for V_s Landing (Flaps Down)");
xlim([10, 40]);
xlabel("\fontsize{14}\fontname{Times}W/S [psf]");
ylabel("\fontsize{14}\fontname{Times}W/P [{shp}/{lb}]");
%% Sizing of Take-off Field length (1)
C_L_max_takeoff = [1.4, 1.6, 1.8, 2.0];
s_tofl = 1200;
sigma = 0.842744;
W_over_S = 1:0.01:500;
figure, hold on
for C_L_max = C_L_max_takeoff
    T_over_W = W_over_S./((s_tofl/37.5) .* sigma .* C_L_max);
    plot(W_over_S, T_over_W,"LineWidth",3);
    grid on
end
legend('(Cl_m_a_x)_T_O = 1.4', '(Cl_m_a_x)_T_O = 1.6', '(Cl_m_a_x)_T_O = 1.8', '(Cl_m_a_x)_T_O = 2.0');
title("\fontsize{14}\fontname{Times}Thrust Loading to Wing Loading");
xlabel("W/S");
ylabel("T/W");
hold off
%% T_P Diagram
T_ = 0:0.01:700;
P_ = T_./(20/7);
figure, plot(P_, T_, 'Color', 'red', 'LineWidth', 3);
title('\fontsize{16}\fontname{Times}T-P Diagram');
xlabel('\fontsize{14}\fontname{Times}Take-off Thrust [T]');
ylabel('\fontsize{14}\fontname{Times}Take-off Shaft Horse Power [P]');
grid on
%% Sizing of Take-off Field length (2)
C_L_max_takeoff = [1.4, 1.6, 1.8, 2.0];
s_tofl = 1200;
sigma = 0.842744;
W_TO = 14512.05;
W_over_S = 1:0.01:500;
figure, hold on
for C_L_max = C_L_max_takeoff
    T_over_W = W_over_S./((s_tofl/37.5).* sigma .* C_L_max);
    T = T_over_W .* W_TO;
    P = T./(20/7);
    W_over_P = W_TO./P;
    plot(W_over_S, W_over_P,"LineWidth",3);
    grid on
end
legend('(Cl_m_a_x)_T_O = 1.4', '(Cl_m_a_x)_T_O = 1.6', '(Cl_m_a_x)_T_O = 1.8', '(Cl_m_a_x)_T_O = 2.0');
title("\fontsize{14}\fontname{Times}Power Loading to Wing Loading");
xlim([0 30]);
ylim([0 40]);
xlabel("W/S");
ylabel("W/P");
hold off
%% Sizing of Landing Field length (1)
rho_TOandL = 0.00200335; % Unit: slug/ft^3
W_TO = 14512.05;
W_L = [0.88 0.99 1.00] * W_TO;
fprintf('W_L = %d\n', W_L);
C_L_max_landing = [1.6, 1.8, 2.0, 2.2, 2.5];
S_L = 1020;
S_FL = S_L / 0.6;
V_A = sqrt(S_FL/0.3);
V_S_L = V_A/1.3;
sva = 0:0.01:10; % Square of V approach
sfl = (6/20) * sva;   % Square of V approach to FAR25 Landing Field Length
figure, plot(sva, sfl, ':','LineWidth', 3);
title("\fontsize{16}\fontname{Times}Square of V approach to FAR25 Landing Field Length");
xlabel("\fontsize{14}\fontname{Times}(V_A)^2  [(KN^2)10^-^3]");
ylabel("\fontsize{14}\fontname{Times}FAR25 Landing Field Length - S_F_L  [ft(10^-^3)]");
grid on
hold on
v_a_p = (V_A ^ 2) * (10 ^ (-3));
sfl_p = S_FL * (10 ^ (-3));
plot(v_a_p, sfl_p,'O','MarkerSize',7, 'MarkerFaceColor', 'r');
legend('\fontsize{14}\fontname{Times}Average Line', '\fontsize{14}\fontname{Times}Our Airplane');
%% Sizing of Landing Field length (2)
disp("                           ""For Landing Wing Loading""         ");
V_S_L = V_S_L * 1.6878098571;   % Changing the unit from KN to ft/s 
figure, hold on
for C_L = C_L_max_landing
    W_over_S_Landing = 0.5 * rho_TOandL * (V_S_L ^ 2) * C_L;
    fprintf('For C_L_max_Landing = %d    ---->    W_over_S = %d\n',C_L,W_over_S_Landing);
    xline(W_over_S_Landing, "Color",rand(1,3),"LineWidth",3,"Label",C_L);
end
title("\fontsize{16}\fontname{Times}Wing Loading in Landing");
xlabel("\fontsize{16}\fontname{Times}(W/S)_L_a_n_d_i_n_g");
ylabel("\fontsize{16}\fontname{Times}(W/P)_L_a_n_d_i_n_g");
xlim([10 25]);
%% Sizing of Landing Field length (3)
rho_TOandL = 0.00200335; % Unit: slug/ft^3
C_L_max_landing = [1.6, 1.8, 2.0, 2.2, 2.5];
figure,
disp("                           ""Medium Aircraft Take-off Wing Loading""         ");
subplot(3,1,1), hold on
for C_L = C_L_max_landing
W_over_S_TakeOff = ((0.5 * rho_TOandL * (V_S_L ^ 2))/0.88) * C_L;
fprintf('For C_L_max_clean = %d    ---->    W_over_S_Take-off = %d\n',C_L,W_over_S_TakeOff);
xline(W_over_S_TakeOff, "Color",rand(1,3),"LineWidth",3,"Label",C_L);
end
title("\fontsize{16}\fontname{Times}Landing Sizing for Medium Aircraft");
% xlabel("\fontsize{14}\fontname{Times}W/S [psf]");
ylabel("\fontsize{14}\fontname{Times}W/P");
xlim([13, 30]);

disp("                           ""Average Aircraft Take-off Wing Loading""         ");
subplot(3,1,2), hold on
for C_L = C_L_max_landing
W_over_S_TakeOff = ((0.5 * rho_TOandL * (V_S_L ^ 2))/0.99) * C_L;
fprintf('For C_L_max_Takeoff = %d    ---->    W_over_S_Take-off = %d\n',C_L,W_over_S_TakeOff);
xline(W_over_S_TakeOff, "Color",rand(1,3),"LineWidth",3,"Label",C_L);
end
title("\fontsize{16}\fontname{Times}Landing Sizing for Average Aircraft");
% xlabel("\fontsize{14}\fontname{Times}W/S [psf]");
ylabel("\fontsize{14}\fontname{Times}W/P");
xlim([13, 30]);

disp("                           ""Maximum Aircraft Take-off Wing Loading""         ");
subplot(3,1,3), hold on
for C_L = C_L_max_landing
W_over_S_TakeOff = ((0.5 * rho_TOandL * (V_S_L ^ 2))/1) * C_L;
fprintf('For C_L_max_Landing = %d    ---->    W_over_S_Take-off = %d\n',C_L,W_over_S_TakeOff);
xline(W_over_S_TakeOff, "Color",rand(1,3),"LineWidth",3,"Label",C_L);
end
title("\fontsize{16}\fontname{Times}Landing Sizing for Maximum Aircraft");
xlabel("\fontsize{14}\fontname{Times}W/S [psf]");
ylabel("\fontsize{14}\fontname{Times}W/P");
xlim([13, 30]);
%% ُClimb Sizing
S_wet = 1810.185;
m = 0.0100;
x = 100:0.1:10000;
x_1 = 0;
x_2 = -6.4;
x_3 = -12.7;
x_4 = -19;
y_1 = m * x + x_1;
y_2 = m * x + x_2;
y_3 = m * x + x_3;
y_4 = m * x + x_4;
p_1 = m * S_wet + x_1;
p_2 = m * S_wet + x_2;
p_3 = m * S_wet + x_3;
p_4 = m * S_wet + x_4;
figure, hold on
plot(x, y_1, 'LineWidth',3);
plot(x, y_2, 'LineWidth',3);
plot(x, y_3, 'LineWidth',3);
plot(x, y_4, 'LineWidth',3);
ylim([1 100]);
xlim([100 10000]);
xline(S_wet, 'LineWidth',3,'Label','\fontsize{16}\fontname{Times}S_W_e_t','Color','r');
plot(S_wet, [p_1, p_2, p_3, p_4],'O','MarkerSize',9, 'MarkerFaceColor', 'r','MarkerEdgeColor','black');
text(S_wet,p_3,cellstr(append('    ',num2str(p_3))));
text(S_wet,p_2,cellstr(append('    ',num2str(p_2))));
text(S_wet,p_1,cellstr(append('    ',num2str(p_1))));
legend('\fontsize{16}\fontname{Times}C_f = 0.0020', '\fontsize{16}\fontname{Times}C_f = 0.0030','\fontsize{16}\fontname{Times}C_f = 0.0040', '\fontsize{16}\fontname{Times}C_f = 0.0050');
xlabel('\fontsize{16}\fontname{Times}Wetted Area S_W_e_t ft^2');
ylabel('\fontsize{16}\fontname{Times}(Equivalent Parasite Area ft^2)*10^2');
grid on
%% Drag Polar
syms C_L
AR = 9;
e = [0.85 0.8 0.75 0.8 0.75];
deltaC_D_0 = [0 0.020 0.075 0.025 0.025];
C_D_0 = 0.01645775;
C_D = (C_D_0 + deltaC_D_0);
CC = 1./(pi.* AR .* e);
fprintf('For Clean Phase: %d+(%d)C_L^2\n', C_D(1), CC(1));
fprintf('For Take-ff (gear up) Phase: %d+(%d)C_L^2\n', C_D(2), CC(2));
fprintf('For Landing (gear up) Phase: %d+(%d)C_L^2\n', C_D(3), CC(3));
fprintf('For Take-off (gear down) Phase: %d+(%d)C_L^2\n', C_D(4), CC(4));
fprintf('For Landing (gear down) Phase: %d+(%d)C_L^2\n', C_D(5), CC(5));
disp('');
% Rate of Climb
C_L_max_takeoff = [1.4, 1.6, 1.8, 2.0];
C_L_max_landing = [1.6, 1.8, 2.0, 2.2, 2.5];
V_S_TO = 101.268591;
rho_airport = 0.00200335;
sigma = 0.842744;
RCP = 2/45;
C_D = (C_D + (CC .* (C_L.^2)))';
figure, hold on
% RC FAR25.111 Take-off
disp('');
etha_p = 0.72;
C_L = C_L_max_takeoff';
V_2 = 1.2 * V_S_TO;
W_over_S = 0.5 * rho_airport * (V_2 ^ 2) * C_L;
temp = sqrt(W_over_S)./(19 .* ((C_L .^ (3/2))./subs(C_D(2),C_L)) .* sqrt(sigma));
W_over_P = double(((1 ./ (RCP + temp))) .* etha_p);
fprintf('V_2:%d fps\n',V_2);
disp(table(C_L, W_over_S, W_over_P, 'VariableNames',{'C_L', 'W_over_S psf', 'W_over_P lb/hp'}))
disp(' ')
W_over_S = 15:0.1:70;
temp = sqrt(W_over_S)./(19 .* ((C_L .^ (3/2))./subs(C_D(2),C_L)) .* sqrt(sigma));
W_over_P = double(((1 ./ (RCP + temp))) .* etha_p);
plot(W_over_S, W_over_P);

% RC FAR25.121 Take-off
disp('');
C_L = C_L_max_takeoff';
V_2 = (V_2 + V_S_TO) / 2;
W_over_S = 0.5 * rho_airport * (V_2 ^ 2) * C_L;
temp = sqrt(W_over_S)./(19 .* ((C_L .^ (3/2))./subs(C_D(4),C_L)) .* sqrt(sigma));
W_over_P = double(((1 ./ (RCP + temp))) .* etha_p);
fprintf('V_2:%d fps\n',V_2);
disp(table(C_L, W_over_S, W_over_P, 'VariableNames',{'C_L', 'W_over_S psf', 'W_over_P lb/hp'}))
disp(' ')
W_over_S = 15:0.1:70;
temp = sqrt(W_over_S)./(19 .* ((C_L .^ (3/2))./subs(C_D(4),C_L)) .* sqrt(sigma));
W_over_P = double(((1 ./ (RCP + temp))) .* etha_p);
plot(W_over_S, W_over_P);
% RC FAR25.119 Landing
disp('');
C_L = C_L_max_landing';
V_2 = 1.3 * V_S_TO;
W_over_S = 0.5 * rho_airport * (V_2 ^ 2) * C_L;
temp = sqrt(W_over_S)./(19 .* ((C_L .^ (3/2))./subs(C_D(3),C_L)) .* sqrt(sigma));
W_over_P = double(((1 ./ (RCP + temp))) .* etha_p);
fprintf('V_2:%d fps\n',V_2);
disp(table(C_L, W_over_S, W_over_P, 'VariableNames',{'C_L', 'W_over_S psf', 'W_over_P lb/hp'}))
disp(' ')
W_over_S = 15:0.1:70;
temp = sqrt(W_over_S)./(19 .* ((C_L .^ (3/2))./subs(C_D(3),C_L)) .* sqrt(sigma));
W_over_P = double(((1 ./ (RCP + temp))) .* etha_p);
plot(W_over_S, W_over_P);
% RC FAR25.121 Landing
disp('');
C_L = C_L_max_landing';
V_2 = 1.65 * V_S_TO;
W_over_S = 0.5 * rho_airport * (V_2 ^ 2) * C_L;
temp = sqrt(W_over_S)./(19 .* ((C_L .^ (3/2))./subs(C_D(5),C_L)) .* sqrt(sigma));
W_over_P = double(((1 ./ (RCP + temp))) .* etha_p);
fprintf('V_2:%d fps\n',V_2);
disp(table(C_L, W_over_S, W_over_P, 'VariableNames',{'C_L', 'W_over_S psf', 'W_over_P lb/hp'}))
disp(' ')
W_over_S = 15:0.1:70;
temp = sqrt(W_over_S)./(19 .* ((C_L .^ (3/2))./subs(C_D(5),C_L)) .* sqrt(sigma));
W_over_P = double(((1 ./ (RCP + temp))) .* etha_p);
plot(W_over_S, W_over_P);
xlim([20 70]);
xlabel("\fontsize{16}\fontname{Times}Wing Loading {lb}/{ft^2}")
ylabel("\fontsize{16}\fontname{Times}Power Loading {lb}/{hp}")
title("\fontsize{18}\fontname{Times}Effect of FAR25 Climb Requirements on Power Loading and Wing Loading")
hold off
legend('FAR25.111 C_L_{max} = 1.4', 'FAR25.111 C_L_{max} = 1.6', 'FAR25.111 C_L_{max} = 1.8', 'FAR25.111 C_L_{max} = 2', 'FAR25.121 C_L_{max} = 1.4', 'FAR25.121 C_L_{max} = 1.6', 'FAR25.121 C_L_{max} = 1.8', 'FAR25.121 C_L_{max} = 2','FAR25.119 C_L_{max} = 1.6', 'FAR25.119 C_L_{max} = 1.8', 'FAR25.119 C_L_{max} = 2', 'FAR25.119 C_L_{max} = 2.2', 'FAR25.119 C_L_{max} = 2.5', 'FAR25.121 C_L_{max} = 1.6', 'FAR25.121 C_L_{max} = 1.8', 'FAR25.121 C_L_{max} = 2', 'FAR25.121 C_L_{max} = 2.2', 'FAR25.121 C_L_{max} = 2.5');
% Climb Gradient
sigma = 0.842744;
etha_p = 0.72;
L_over_D = 11;
C_L_max_takeoff = [1.4, 1.6, 1.8, 2.0];
C_L_max_landing = [1.6, 1.8, 2.0, 2.2, 2.5];
V_S_TO = 101.268591;
rho_airport = 0.00200335;
figure, hold on
% CGR FAR25.111 Take-off
disp('');
C_L = C_L_max_takeoff';
V_2 = 1.2 * V_S_TO;
CGR = 0.012;
W_over_S = 0.5 .* rho_airport .* (V_2 ^ 2) .* C_L;
temp = ((18.97 .* etha_p .* sqrt(sigma)) ./ sqrt(W_over_S));
W_over_P = (sqrt(C_L) ./ (CGR + (1 ./ L_over_D))) .* temp;
fprintf('V_2:%d fps\n',V_2);
disp(table(C_L, W_over_S, W_over_P, 'VariableNames',{'C_L', 'W_over_S psf', 'W_over_P lb/hp'}))
disp(' ')
W_over_S = 15:0.1:70;
temp = (18.97 .* etha_p .* sqrt(sigma)) ./ sqrt(W_over_S);
W_over_P = (sqrt(C_L) ./ (CGR + (1 ./ L_over_D))) .* temp;
plot(W_over_S, W_over_P);
% CGR FAR25.121 (1) Take-off
disp('');
C_L = C_L_max_takeoff';
V_2 =  (V_2 * V_S_TO)/2;
CGR = 0;
W_over_S = 0.5 .* rho_airport .* (V_2 ^ 2) .* C_L;
temp = (18.97 .* etha_p .* sqrt(sigma)) ./ sqrt(W_over_S);
W_over_P = (sqrt(C_L) ./ (CGR + (1 ./ L_over_D))) .* temp;
fprintf('V_2:%d fps\n',V_2);
disp(table(C_L, W_over_S, W_over_P, 'VariableNames',{'C_L', 'W_over_S psf', 'W_over_P lb/hp'}))
disp(' ')
W_over_S = 15:0.1:70;
temp = (18.97 .* etha_p .* sqrt(sigma)) ./ sqrt(W_over_S);
W_over_P = (sqrt(C_L) ./ (CGR + (1 ./ L_over_D))) .* temp;
plot(W_over_S, W_over_P);
% CGR FAR25.121 (2) Take-off
disp('');
C_L = C_L_max_takeoff';
V_2 =  1.2 * V_S_TO;
CGR = 0.024;
W_over_S = 0.5 .* rho_airport .* (V_2 ^ 2) .* C_L;
temp = (18.97 .* etha_p .* sqrt(sigma)) ./ sqrt(W_over_S);
W_over_P = (sqrt(C_L) ./ (CGR + (1 ./ L_over_D))) .* temp;
fprintf('V_2:%d fps\n',V_2);
disp(table(C_L, W_over_S, W_over_P, 'VariableNames',{'C_L', 'W_over_S psf', 'W_over_P lb/hp'}))
disp(' ')
W_over_S = 15:0.1:70;
temp = (18.97 .* etha_p .* sqrt(sigma)) ./ sqrt(W_over_S);
W_over_P = (sqrt(C_L) ./ (CGR + (1 ./ L_over_D))) .* temp;
plot(W_over_S, W_over_P);
% CGR FAR25.119 Landing
disp('');
C_L = C_L_max_landing';
V_2 =  1.3 * V_S_TO;
CGR = 0.032;
W_over_S = 0.5 .* rho_airport .* (V_2 ^ 2) .* C_L;
temp = (18.97 .* etha_p .* sqrt(sigma)) ./ sqrt(W_over_S);
W_over_P = (sqrt(C_L) ./ (CGR + (1 ./ L_over_D))) .* temp;
fprintf('V_2:%d fps\n',V_2);
disp(table(C_L, W_over_S, W_over_P, 'VariableNames',{'C_L', 'W_over_S psf', 'W_over_P lb/hp'}))
disp(' ')
W_over_S = 15:0.1:70;
temp = (18.97 .* etha_p .* sqrt(sigma)) ./ sqrt(W_over_S);
W_over_P = (sqrt(C_L) ./ (CGR + (1 ./ L_over_D))) .* temp;
plot(W_over_S, W_over_P);
% CGR FAR25.121 Landing
disp('');
C_L = C_L_max_landing';
V_2 =  1.65 * V_S_TO;
CGR = 0.021;
W_over_S = 0.5 .* rho_airport .* (V_2 ^ 2) .* C_L;
temp = (18.97 .* etha_p .* sqrt(sigma)) ./ sqrt(W_over_S);
W_over_P = (sqrt(C_L) ./ (CGR + (1 ./ L_over_D))) .* temp;
fprintf('V_2:%d fps\n',V_2);
disp(table(C_L, W_over_S, W_over_P, 'VariableNames',{'C_L', 'W_over_S psf', 'W_over_P lb/hp'}))
disp(' ')
W_over_S = 15:0.1:70;
temp = (18.97 .* etha_p .* sqrt(sigma)) ./ sqrt(W_over_S);
W_over_P = (sqrt(C_L) ./ (CGR + (1 ./ L_over_D))) .* temp;
plot(W_over_S, W_over_P);
xlim([15 70]);
xlabel("\fontsize{16}\fontname{Times}Wing Loading {lb}/{ft^2}")
ylabel("\fontsize{16}\fontname{Times}Power Loading {lb}/{hp}")
title("\fontsize{18}\fontname{Times}Effect of FAR25 Climb Requirements on Power Loading and Wing Loading")
hold off
legend('FAR25.111 C_L_{max} = 1.4', 'FAR25.111 C_L_{max} = 1.6', 'FAR25.111 C_L_{max} = 1.8', 'FAR25.111 C_L_{max} = 2', 'FAR25.121 C_L_{max} = 1.4', 'FAR25.121 C_L_{max} = 1.6', 'FAR25.121 C_L_{max} = 1.8', 'FAR25.121 C_L_{max} = 2','FAR25.119 C_L_{max} = 1.6', 'FAR25.119 C_L_{max} = 1.8', 'FAR25.119 C_L_{max} = 2', 'FAR25.119 C_L_{max} = 2.2', 'FAR25.119 C_L_{max} = 2.5', 'FAR25.121 C_L_{max} = 1.6', 'FAR25.121 C_L_{max} = 1.8', 'FAR25.121 C_L_{max} = 2', 'FAR25.121 C_L_{max} = 2.2', 'FAR25.121 C_L_{max} = 2.5');
%% Linear Relation Between RC and h
h = 22000; 
t_cl = 15; % Unit min
h_abs = 45000; % Unit ft
RC_0 = (h_abs ./ t_cl) .* log(1 ./ (1 - (h ./ h_abs)));
h = 0:1:45000;
RC_h = RC_0 .* (1 - (h ./ h_abs));
sigma = 0.842744;
figure, plot(RC_h, h, 'LineWidth', 3, 'Color', 'red');
title('\fontsize{18}\fontname{Times}Linearized Rate-of-Climb With Altitude');
xlabel('\fontsize{16}\fontname{Times}Rate-of-Climb [fpm]');
ylabel('\fontsize{16}\fontname{Times}Altitude [ft]');
%% Time-to-Climb Sizing
C_L_max_clean = [1.2, 1.4, 1.6, 1.8];
h = 22000;
RC_0 = (h_abs ./ t_cl) .* log(1 ./ (1 - (h ./ h_abs)));
RC = RC_0 .* (1 - (h ./ h_abs));
RCP = RC / 33000;
W_over_S = 0:0.1:100;
etha_p = 0.72;
figure, hold on;
for C_L = C_L_max_clean
    C_D_clean = 1.645775e-02 + (4.160914e-02) .* C_L .^ 2;
    temp = sqrt(W_over_S)./(19 .* ((C_L .^ (3/2))./ C_D_clean) .* sqrt(sigma));
    W_over_P = double(((1 ./ (RCP + temp))) .* etha_p);
    plot(W_over_S, W_over_P, 'LineWidth', 2);
end
xlim([10 100])
legend('\fontsize{14}\fontname{Times}C_L_{max} = 1.2' , '\fontsize{14}\fontname{Times}C_L_{max} = 1.4' , '\fontsize{14}\fontname{Times}C_L_{max} = 1.6' , '\fontsize{14}\fontname{Times}C_L_{max} = 1.8')
title('\fontsize{18}\fontname{Times} Time to Climb Sizing');
xlabel('\fontsize{16}\fontname{Times}W/S [psf]');
ylabel('\fontsize{16}\fontname{Times}W/P [lb/hsp]');
hold off
%% Ceiling
C_L_max_clean = [1.2, 1.4, 1.6, 1.8];
RC = 100;
RCP = RC / 33000;
W_over_S = 0:0.1:100;
etha_p = 0.72;
figure, hold on;
for C_L = C_L_max_clean
    C_D_clean = 1.645775e-02 + (4.160914e-02) .* C_L .^ 2;
    temp = sqrt(W_over_S)./(19 .* ((C_L .^ (3/2))./ C_D_clean) .* sqrt(sigma));
    W_over_P = double(((1 ./ (RCP + temp))) .* etha_p);
    plot(W_over_S, W_over_P, 'LineWidth', 2);
end
xlim([10 100])
legend('\fontsize{14}\fontname{Times}C_L_{max} = 1.2' , '\fontsize{14}\fontname{Times}C_L_{max} = 1.4' , '\fontsize{14}\fontname{Times}C_L_{max} = 1.6' , '\fontsize{14}\fontname{Times}C_L_{max} = 1.8')
title('\fontsize{18}\fontname{Times} Ceiling Sizing');
xlabel('\fontsize{16}\fontname{Times}W/S [psf]');
ylabel('\fontsize{16}\fontname{Times}W/P [lb/hsp]');
hold off
%% Manuvering
W_over_S = 10:0.1:100;
C_L_max = 1.8;
AR = 9;
e = 0.85;
M = 0.2764;
delta = 10000;
q = 109.692;
C_D_0 = 0.01645775;
n_max = (1.482 .* (M .^ 2) .* delta .* C_L_max) ./ W_over_S;
W_over_P =(22 ./ 7) ./ ((C_D_0 .* q) ./ W_over_S) + ((W_over_S .* (n_max .^ 2)) ./ (pi .* AR .* e .* q));
figure, plot(W_over_S, W_over_P, 'Color', 'red', 'LineWidth', 3);
title('\fontsize{18}\fontname{Times}Sizing to Maneuvering');
xlabel('\fontsize{16}\fontname{Times} W/S [psf]');
ylabel('\fontsize{16}\fontname{Times} W/P [lb/hsp]');
%% Power Index
x = 0:0.01:3.5;
y = (500 ./ 2.92) .* x;
figure, hold on;
plot(x, y);
plot(1.2097, 207.14, 'Marker','O', 'MarkerFaceColor', 'red', 'MarkerSize', 7);
text(1.2097, 207.14, '   \fontsize{16}\fontname{Times}Power Index = 1.2097');
xlabel('\fontsize{16}\fontname{Times}Power Index');
ylabel('\fontsize{16}\fontname{Times}Speed [mph]');
title('\fontsize{16}\fontname{Times}Correlation of airplane speed with power index for retractable gear, cantileverd wing configurations')
%% Cruise
sigma = 0.738479;
Ip = 1.2097;
W_over_S = 0:0.01:100;
W_over_P = W_over_S ./ ((Ip .^ 3) * sigma);
figure, plot(W_over_S, W_over_P, 'LineWidth', 3, "Color", 'red');
title('\fontsize{18}\fontname{Times}Cruise Sizing');
xlabel('\fontsize{16}\fontname{Times} W/S [psf]');
ylabel('\fontsize{16}\fontname{Times} W/P [lb/hsp]');
%% Matching Diagram
% rho_airport = 0.00177847;
C_L_max_takeoff = 2;
C_L_max_landing = 2.8;
W_over_S = 0:0.01:200;
% V_S_TO = 101.268591;
figure, hold on;

% Take off
s_tofl = 1200;
sigma = 0.842744;
W_TO = 14512.05;
T_over_W = W_over_S./((s_tofl/37.5).* sigma .* C_L_max_takeoff);
T = T_over_W .* W_TO;
P = T./(20/7);
W_over_P = W_TO./P;
plot(W_over_S, W_over_P,"LineWidth",3);
name_1 = '\fontname{Times}Take-off C_L = 2'; 
% Landing
% S_L = 1050;
S_FL = 1050;
V_A = sqrt(S_FL/0.3);
V_S_L = V_A/1.3;
V_S_L = V_S_L * 1.6878098571;   % Changing the unit from KN to ft/s 
rho_TOandL = 0.00200335; % Unit: slug/ft^3
for C_L = C_L_max_landing
W_over_S_TakeOff = ((0.5 .* rho_TOandL .* (V_S_L .^ 2))/0.88) .* C_L;
xline(W_over_S_TakeOff, "LineWidth", 3);
end
name_2 = '\fontname{Times}Landing C_L = 2.5';
% Time to climb
t_cl = 15;
C_L_max_clean = 1.8;
h = 22000;
h_abs = 45000;
RC_0 = (h_abs ./ t_cl) .* log(1 ./ (1 - (h ./ h_abs)));
RC = RC_0 .* (1 - (h ./ h_abs));
RCP = RC / 33000;
etha_p = 0.72;
C_D_clean = 1.645775e-02 + (4.160914e-02) .* C_L_max_clean .^ 2;
temp = sqrt(W_over_S)./(19 .* ((C_L_max_clean .^ (3/2))./ C_D_clean) .* sqrt(sigma));
W_over_P = double(((1 ./ (RCP + temp))) .* etha_p);
plot(W_over_S, W_over_P, 'LineWidth', 2);
name_3 = '\fontname{Times}Time to Climb C_L = 1.8';
% Cruise
sigma = 0.738479;
Ip = 1.2097;
W_over_S = 0:0.01:100;
W_over_P = W_over_S ./ ((Ip .^ 3) * sigma);
plot(W_over_S, W_over_P, 'LineWidth', 3);
name_4 = '\fontname{Times}Cruise';
% Ceiling
C_L_max_clean = 1.8;
RC = 100;
RCP = RC / 33000;
etha_p = 0.72;
C_D_clean = 1.645775e-02 + (4.160914e-02) .* C_L_max_clean .^ 2;
temp = sqrt(W_over_S)./(19 .* ((C_L_max_clean .^ (3/2))./ C_D_clean) .* sqrt(sigma));
W_over_P = double(((1 ./ (RCP + temp))) .* etha_p);
plot(W_over_S, W_over_P, 'LineWidth', 2);
name_5 = '\fontname{Times}Ceiling C_L = 1.8';
% Climb
syms C_L
AR = 9;
e = [0.85 0.8 0.75 0.8 0.75];
deltaC_D_0 = [0 0.020 0.075 0.025 0.025];
C_D_0 = 0.01645775;
C_D = (C_D_0 + deltaC_D_0);
CC = 1./(pi.* AR .* e);
% V_S_TO = 101.268591;
sigma = 0.842744;
RCP = 2/45;
C_D = (C_D + (CC .* (C_L.^2)))';
% RC FAR25.111 Take-off
etha_p = 0.72;
C_L = C_L_max_takeoff';
% V_2 = 1.2 * V_S_TO;
temp = sqrt(W_over_S)./(19 .* ((C_L .^ (3/2))./subs(C_D(2),C_L)) .* sqrt(sigma));
W_over_P = double(((1 ./ (RCP + temp))) .* etha_p);
plot(W_over_S, W_over_P, 'LineWidth', 3, 'LineStyle', ':');
name_6 = '\fontname{Times}RC FAR25.111 Take-off';
% RC FAR25.121 Take-off
C_L = C_L_max_takeoff';
% V_2 = (V_2 + V_S_TO) / 2;
temp = sqrt(W_over_S)./(19 .* ((C_L .^ (3/2))./subs(C_D(4),C_L)) .* sqrt(sigma));
W_over_P = double(((1 ./ (RCP + temp))) .* etha_p);
plot(W_over_S, W_over_P, 'LineWidth', 3, 'LineStyle', ':');
name_7 = '\fontname{Times}RC FAR25.121 Take-off';
% RC FAR25.119 Landing
C_L = C_L_max_landing';
% V_2 = 1.3 * V_S_TO;
temp = sqrt(W_over_S)./(19 .* ((C_L .^ (3/2))./subs(C_D(3),C_L)) .* sqrt(sigma));
W_over_P = double(((1 ./ (RCP + temp))) .* etha_p);
plot(W_over_S, W_over_P, 'LineWidth', 3, 'LineStyle', ':');
name_8 = '\fontname{Times}RC FAR25.119 Landing';
% RC FAR25.121 Landing
C_L = C_L_max_landing';
% V_2 = 1.65 * V_S_TO;
temp = sqrt(W_over_S)./(19 .* ((C_L .^ (3/2))./subs(C_D(5),C_L)) .* sqrt(sigma));
W_over_P = double(((1 ./ (RCP + temp))) .* etha_p);
plot(W_over_S, W_over_P, 'LineWidth', 3, 'LineStyle', ':');
name_9 = '\fontname{Times}RC FAR25.121 Landing';
% Climb Gradient
sigma = 0.842744;
etha_p = 0.72;
L_over_D = 11;
V_S_TO = 101.268591;
rho_airport = 0.00200335;
% CGR FAR25.111 Take-off
C_L = C_L_max_takeoff';
% V_2 = 1.2 * V_S_TO;
CGR = 0.012;
temp = (18.97 .* etha_p .* sqrt(sigma)) ./ sqrt(W_over_S);
W_over_P = (sqrt(C_L) ./ (CGR + (1 ./ L_over_D))) .* temp;
plot(W_over_S, W_over_P, 'LineWidth', 3, 'LineStyle', '--');
name_10 = '\fontname{Times}CGR FAR25.111 Take-off';
% CGR FAR25.121 (1) Take-off
C_L = C_L_max_takeoff';
% V_2 =  (V_2 * V_S_TO)/2;
CGR = 0;
temp = (18.97 .* etha_p .* sqrt(sigma)) ./ sqrt(W_over_S);
W_over_P = (sqrt(C_L) ./ (CGR + (1 ./ L_over_D))) .* temp;
plot(W_over_S, W_over_P, 'LineWidth', 3, 'LineStyle', '--');
name_11 = '\fontname{Times}CGR FAR25.121 (1) Take-off';
% CGR FAR25.121 (2) Take-off
C_L = C_L_max_takeoff';
% V_2 =  1.2 * V_S_TO;
CGR = 0.024;
temp = (18.97 .* etha_p .* sqrt(sigma)) ./ sqrt(W_over_S);
W_over_P = (sqrt(C_L) ./ (CGR + (1 ./ L_over_D))) .* temp;
plot(W_over_S, W_over_P, 'LineWidth', 3, 'LineStyle', '--');
name_12 = '\fontname{Times}CGR FAR25.121 (2) Take-off';
% CGR FAR25.119 Landing
C_L = C_L_max_landing';
% V_2 =  1.3 * V_S_TO;
CGR = 0.032;
temp = (18.97 .* etha_p .* sqrt(sigma)) ./ sqrt(W_over_S);
W_over_P = (sqrt(C_L) ./ (CGR + (1 ./ L_over_D))) .* temp;
plot(W_over_S, W_over_P, 'LineWidth', 3, 'LineStyle', '--');
name_13 = '\fontname{Times}CGR FAR25.119 Landing';
% CGR FAR25.121 Landing
C_L = C_L_max_landing';
V_2 =  1.65 * V_S_TO;
CGR = 0.021;
temp = (18.97 .* etha_p .* sqrt(sigma)) ./ sqrt(W_over_S);
W_over_P = (sqrt(C_L) ./ (CGR + (1 ./ L_over_D))) .* temp;
plot(W_over_S, W_over_P, 'LineWidth', 3, 'LineStyle', '--');
name_14 = '\fontname{Times}CGR FAR25.121 Landing';
% Stall Take off
% W_over_S_Stall = 0.5 * rho_TOandL * (V_S_TO ^ 2) * C_L_max_takeoff;
% xline(W_over_S_Stall, "LineWidth",3);
% name_15 = '\fontname{Times}Stall Speed for Take-off C_L = 2';
% Design Point
plot(28.05, 5.49382, 'Marker','O', 'MarkerFaceColor', 'red', 'MarkerSize', 10);
text(28.05, 5.49382, '\fontname{Times}   Design Point');
fprintf('Design Point:\n');
fprintf('       Wing Loading: 28.05\n');
fprintf('       Power Loading: 5.49382\n');
xlabel('\fontsize{16}\fontname{Times}Wing Loading [psf]');
ylabel('\fontsize{16}\fontname{Times}Power Loading [lb/hsp]');
title('\fontsize{18}\fontname{Times} Matching Diagram');
legend(name_1, name_2, name_3, name_4, name_5, name_6, name_7, name_8, name_9, name_10, name_11, name_12, name_13, name_14);
ylim([0 70]);
xlim([0 100]);
grid on;

toc