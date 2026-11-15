%% =====================================================================
%  PUNTO 2 - Diseño de un controlador en espacio de estados
%  Motor DC con reductora y encoder JGB37-555B (12 V, 1:30)
%  Control de posición angular del eje de salida (rueda)
%
%  Autor : Santiago Alexander Bravo Rodriguez - 2224734
%  Curso : Sistemas de Control II - UIS
%
%  Contenido:
%   1. Parámetros del motor
%   2. Modelo en espacio de estados
%   3. Análisis en lazo abierto (polos, controlabilidad, observabilidad)
%   4. Validación del modelo contra la hoja de datos
%   5. Especificaciones y sistema aumentado con acción integral
%   6. Diseño 1: ubicación de polos (place)
%   7. Diseño 2: regulador óptimo LQR (lqr)
%   8. Observador de Luenberger con estimación del torque de carga
%   9. Simulaciones (lineal, con observador, saturación y encoder)
%  10. Tabla de resultados y figuras
%
%  Requiere: Control System Toolbox.
%  Las figuras se guardan en la carpeta IMAGENES/ y los resultados
%  numéricos en resultados_punto2.txt
%  =====================================================================
clear; close all; clc;

carpeta = fullfile(pwd,'IMAGENES');
if ~exist(carpeta,'dir'); mkdir(carpeta); end
fid = fopen('resultados_punto2.txt','w');
imprime = @(varargin) [fprintf(varargin{:}); fprintf(fid,varargin{:})];

%% 1. PARÁMETROS DEL MOTOR -------------------------------------------------
% Valores nominales del JGB37-555B 12 V 1:30 (hoja de datos del fabricante)
V_nom   = 12;       % [V]     tensión nominal
n0_out  = 333;      % [rpm]   velocidad en vacío en el eje de salida
I0      = 0.15;     % [A]     corriente en vacío
N       = 30;       % [-]     relación de reducción
PPR     = 11;       % [pulsos/rev] encoder Hall en el eje del motor

% Parámetros eléctricos y mecánicos (estimados, ver documento)
R   = 4.0;          % [Ohm]   resistencia de armadura
L   = 2.5e-3;       % [H]     inductancia de armadura
Ke  = (V_nom - R*I0)/(n0_out*N*2*pi/60);   % [V·s/rad] constante FCEM
Kt  = Ke;           % [N·m/A] constante de torque (SI: Kt = Ke)
Jm  = 5e-6;         % [kg·m^2] inercia del rotor
mw  = 0.05;         % [kg]    masa de la rueda
rw  = 0.0325;       % [m]     radio de la rueda (65 mm de diámetro)
Jw  = 0.5*mw*rw^2;  % [kg·m^2] inercia de la rueda (disco)
Jeq = Jm + Jw/N^2;  % [kg·m^2] inercia equivalente vista desde el motor
w0m = n0_out*N*2*pi/60;                    % [rad/s] vel. vacío del motor
b   = Kt*I0/w0m;    % [N·m·s/rad] fricción viscosa (equilibrio en vacío)

imprime('==== 1. PARAMETROS ====\n');
imprime('Ke = Kt = %.5f V.s/rad\n',Ke);
imprime('b       = %.4e N.m.s/rad\n',b);
imprime('Jw      = %.4e kg.m^2\n',Jw);
imprime('Jeq     = %.4e kg.m^2\n',Jeq);
imprime('tau_e = L/R = %.3f ms\n',1e3*L/R);
imprime('tau_m = Jeq*R/(Kt*Ke) = %.3f s\n\n',Jeq*R/(Kt*Ke));

%% 2. MODELO EN ESPACIO DE ESTADOS ------------------------------------------
% Estados : x1 = theta (posición eje de salida) [rad]
%           x2 = omega (velocidad eje de salida) [rad/s]
%           x3 = i     (corriente de armadura)   [A]
% Entrada : u  = voltaje de armadura [V]
% Perturb.: TL = torque de carga en el eje de salida [N·m]
% Salida  : y  = theta (medida por el encoder)
A = [0      1              0;
     0     -b/Jeq          Kt/(N*Jeq);
     0     -Ke*N/L        -R/L];
B = [0; 0; 1/L];
E = [0; -1/(N^2*Jeq); 0];
C = [1 0 0];
D = 0;

sys_ol = ss(A,B,C,D);
sys_ol.StateName  = {'theta','omega','i'};
sys_ol.InputName  = {'u'};
sys_ol.OutputName = {'theta'};

imprime('==== 2. MODELO ====\n');
imprime('A =\n'); imprime('  %12.4f %12.4f %12.4f\n',A.');
imprime('B = [%g; %g; %g]\n',B);
imprime('E = [%g; %g; %g]\n',E);
imprime('C = [1 0 0], D = 0\n\n');

G = tf(sys_ol);
[numG,denG] = tfdata(G,'v');
imprime('G(s) = theta(s)/U(s): num = %s\n',mat2str(numG,6));
imprime('                      den = %s\n\n',mat2str(denG,6));

%% 3. ANÁLISIS EN LAZO ABIERTO ----------------------------------------------
p_ol = eig(A);
Mc   = ctrb(A,B);
Mo   = obsv(A,C);
imprime('==== 3. ANALISIS LAZO ABIERTO ====\n');
imprime('Polos de lazo abierto: %s\n',mat2str(p_ol.',6));
imprime('rango(ctrb) = %d de %d  -> controlable\n',rank(Mc),size(A,1));
imprime('rango(obsv) = %d de %d  -> observable\n',rank(Mo),size(A,1));
imprime('det(ctrb) = %.4e, det(obsv) = %.4e\n\n',det(Mc),det(Mo));

%% 4. VALIDACIÓN DEL MODELO CONTRA LA HOJA DE DATOS -------------------------
% Escalón de 12 V en lazo abierto: debe llegar a ~333 rpm y ~0.15 A
sys_vel = ss(A,B,eye(3),zeros(3,1));
t_ol = linspace(0,1,2000);
[x_ol,~] = lsim(sys_vel, V_nom*ones(size(t_ol)), t_ol);
w_rpm = x_ol(:,2)*60/(2*pi);
x_ss  = -A(2:3,2:3)\B(2:3)*V_nom;
imprime('==== 4. VALIDACION (escalon 12 V) ====\n');
imprime('Velocidad final modelo = %.1f rpm (hoja de datos: %d rpm)\n',x_ss(1)*60/(2*pi),n0_out);
imprime('Corriente final modelo = %.3f A  (hoja de datos: %.2f A)\n',x_ss(2),I0);
imprime('Corriente de arranque  = %.2f A  (V/R = %.2f A)\n\n',max(x_ol(:,3)),V_nom/R);

nuevaFigura('Name','Lazo abierto','Position',[100 100 700 450]);
subplot(2,1,1); plot(t_ol,w_rpm,'LineWidth',1.5); grid on;
yline(n0_out,'--r','333 rpm (hoja de datos)');
ylabel('\omega [rpm]'); title('Respuesta en lazo abierto a un escalón de 12 V');
subplot(2,1,2); plot(t_ol,x_ol(:,3),'LineWidth',1.5); grid on;
ylabel('i [A]'); xlabel('t [s]');
exportgraphics(gcf,fullfile(carpeta,'fig1_lazo_abierto.png'),'Resolution',200);

%% 5. ESPECIFICACIONES Y SISTEMA AUMENTADO ----------------------------------
% Especificaciones (referencia escalón de pi rad = 180°):
%   ts(2%) <= 0.5 s, Mp <= 5 %, ess = 0 (también ante torque de carga),
%   |u| <= 12 V
ts_d = 0.5; Mp_d = 5; Umax = 12; r_ref = pi; TL = 0.3; t_TL = 1.5;

% Acción integral: xi' = r - theta  ->  xa = [theta omega i xi]'
Aa = [A zeros(3,1); -C 0];
Ba = [B; 0];
Ea = [E; 0];
Br = [0;0;0;1];
imprime('==== 5. SISTEMA AUMENTADO ====\n');
imprime('rango(ctrb aumentado) = %d de %d\n\n',rank(ctrb(Aa,Ba)),4);

%% 6. DISEÑO 1: UBICACIÓN DE POLOS ------------------------------------------
zeta = 0.8;                 % Mp teórico = 1.5 %
wn   = 10;                  % ts ~ 4/(zeta*wn) = 0.5 s
p_dom = roots([1 2*zeta*wn wn^2]).';
p_des = [p_dom, -40, -1600];  % polo integral y polo eléctrico rápido
Kpp = place(Aa,Ba,p_des);
imprime('==== 6. UBICACION DE POLOS ====\n');
imprime('Polos deseados: %s\n',mat2str(p_des,5));
imprime('K = [k1 k2 k3 ki] = %s\n',mat2str(Kpp,6));
imprime('Polos obtenidos: %s\n\n',mat2str(eig(Aa-Ba*Kpp).',5));

%% 7. DISEÑO 2: LQR ---------------------------------------------------------
Q  = diag([50 0.5 0 8000]);  % pesos: theta, omega, i, integral del error
Rw = 0.3;                    % peso del esfuerzo de control
[Klqr,Plqr,p_lqr] = lqr(Aa,Ba,Q,Rw);
imprime('==== 7. LQR ====\n');
imprime('Q = diag(%s), R = %g\n',mat2str(diag(Q).'),Rw);
imprime('K = [k1 k2 k3 ki] = %s\n',mat2str(Klqr,6));
imprime('Polos de lazo cerrado: %s\n\n',mat2str(p_lqr.',5));

nuevaFigura('Name','Polos','Position',[100 100 700 450]);
plot(real(p_ol),imag(p_ol),'kx','MarkerSize',10,'LineWidth',2); hold on;
pp = eig(Aa-Ba*Kpp); pl = eig(Aa-Ba*Klqr);
plot(real(pp),imag(pp),'bo','MarkerSize',8,'LineWidth',1.5);
plot(real(pl),imag(pl),'rs','MarkerSize',8,'LineWidth',1.5);
grid on; xlim([-60 5]); xlabel('Re'); ylabel('Im');
legend('Lazo abierto','Ubicación de polos','LQR','Location','best');
title('Polos dominantes (los polos eléctricos en ~-1600 quedan fuera de la vista)');
exportgraphics(gcf,fullfile(carpeta,'fig2_polos.png'),'Resolution',200);

%% 8. OBSERVADOR DE LUENBERGER CON ESTIMACIÓN DE PERTURBACIÓN ---------------
% Solo se mide theta (encoder). Se estiman omega, i y el torque de carga TL,
% modelado como constante (TL' = 0): xo = [theta omega i TL]'
Ao = [A E; zeros(1,4)];
Bo = [B; 0];
Co = [C 0];
p_obs = [-50 -60 -70 -1600];   % 5-8 veces más rápidos que los dominantes;
                               % el polo eléctrico se deja cerca del de la planta
Lobs  = place(Ao',Co',p_obs)';
imprime('==== 8. OBSERVADOR ====\n');
imprime('rango(obsv aumentado) = %d de %d\n',rank(obsv(Ao,Co)),4);
imprime('Polos del observador: %s\n',mat2str(p_obs));
imprime('L = %s\n\n',mat2str(Lobs.',6));

%% 9a. SIMULACIÓN LINEAL CON REALIMENTACIÓN DE ESTADOS COMPLETA -------------
% Entradas del lazo cerrado: [r ; TL]. Salidas: theta, u, i
t  = (0:1e-4:3)';
rr = r_ref*ones(size(t));
dd = TL*(t>=t_TL);
mk_cl = @(K) ss(Aa-Ba*K, [Br Ea], [1 0 0 0; -K; 0 0 1 0], zeros(3,2));
y_pp  = lsim(mk_cl(Kpp),  [rr dd], t);
y_lqr = lsim(mk_cl(Klqr), [rr dd], t);

nuevaFigura('Name','Comparacion','Position',[100 100 750 600]);
subplot(3,1,1);
plot(t,y_pp(:,1),'b',t,y_lqr(:,1),'r','LineWidth',1.4); hold on;
plot(t,rr,'k--'); grid on; ylabel('\theta [rad]');
legend('Ubicación de polos','LQR','Referencia','Location','southeast');
title('Escalón de \pi rad y torque de carga de 0.3 N·m en t = 1.5 s');
subplot(3,1,2);
plot(t,y_pp(:,2),'b',t,y_lqr(:,2),'r','LineWidth',1.4); grid on;
yline(Umax,'--k'); yline(-Umax,'--k'); ylabel('u [V]');
subplot(3,1,3);
plot(t,y_pp(:,3),'b',t,y_lqr(:,3),'r','LineWidth',1.4); grid on;
ylabel('i [A]'); xlabel('t [s]');
exportgraphics(gcf,fullfile(carpeta,'fig3_comparacion.png'),'Resolution',200);

% Métricas (solo el escalón, antes de la perturbación)
idx = t < t_TL;
S_pp  = stepinfo(y_pp(idx,1),  t(idx), r_ref, 'SettlingTimeThreshold',0.02);
S_lqr = stepinfo(y_lqr(idx,1), t(idx), r_ref, 'SettlingTimeThreshold',0.02);
dev_pp  = max(abs(y_pp(~idx,1)-r_ref));
dev_lqr = max(abs(y_lqr(~idx,1)-r_ref));

%% 9b. SIMULACIÓN REALISTA: OBSERVADOR + SATURACIÓN + ENCODER ---------------
% Integración de paso fijo (Euler, dt = 10 us) del motor y del observador.
% La medida del encoder se cuantiza a 2*pi/(PPR*4*N) rad (1320 cuentas/rev).
q   = 2*pi/(PPR*4*N);
dt  = 1e-5;  ts_ = (0:dt:3)';  n = numel(ts_);
x  = zeros(3,1); xh = zeros(4,1); xi = 0;
X = zeros(n,3); XH = zeros(n,4); U = zeros(n,1); TLv = zeros(n,1);
K = Klqr;
for k = 1:n
    tk = ts_(k);
    y  = q*round(x(1)/q);                       % encoder cuantizado
    u  = -K(1:3)*xh(1:3) - K(4)*xi;             % ley de control con x estimado (TL no se realimenta)
    u  = min(max(u,-Umax),Umax);                % saturación del driver
    TLk = TL*(tk>=t_TL);
    X(k,:) = x.'; XH(k,:) = xh.'; U(k) = u; TLv(k) = TLk;
    dx  = A*x  + B*u + E*TLk;                   % planta
    dxh = Ao*xh + Bo*u + Lobs*(y - Co*xh);      % observador
    x  = x  + dt*dx;
    xh = xh + dt*dxh;
    xi = xi + dt*(r_ref - y);
end
idx2 = ts_ < t_TL;
S_obs  = stepinfo(X(idx2,1), ts_(idx2), r_ref, 'SettlingTimeThreshold',0.02);
dev_obs = max(abs(X(~idx2,1)-r_ref));

nuevaFigura('Name','Observador','Position',[100 100 750 760]);
subplot(4,1,1);
plot(ts_,X(:,1),'r',ts_,XH(:,1),'b--','LineWidth',1.3); hold on;
plot(ts_,r_ref*ones(n,1),'k:'); grid on; ylabel('\theta [rad]');
legend('Real','Estimado','Referencia','Location','southeast');
title('LQR + observador, con saturación de \pm12 V y encoder de 1320 cuentas/rev');
subplot(4,1,2);
plot(ts_,X(:,2),'r',ts_,XH(:,2),'b--','LineWidth',1.3); grid on;
ylabel('\omega [rad/s]');
subplot(4,1,3);
plot(ts_,X(:,3),'r',ts_,XH(:,3),'b--','LineWidth',1.3); grid on;
ylabel('i [A]');
subplot(4,1,4);
plot(ts_,TLv,'r',ts_,XH(:,4),'b--','LineWidth',1.3); grid on;
ylabel('T_L [N·m]'); xlabel('t [s]');
exportgraphics(gcf,fullfile(carpeta,'fig4_observador.png'),'Resolution',200);

nuevaFigura('Name','Senal de control','Position',[100 100 750 350]);
plot(ts_,U,'LineWidth',1.3); grid on; hold on;
yline(Umax,'--k'); yline(-Umax,'--k');
ylabel('u [V]'); xlabel('t [s]');
title('Señal de control en la simulación con observador');
exportgraphics(gcf,fullfile(carpeta,'fig5_control_observador.png'),'Resolution',200);

% Error de estimación con condición inicial distinta (prueba del observador)
sys_err = ss(Ao-Lobs*Co, zeros(4,1), eye(4), zeros(4,1));
[e_obs,te] = initial(sys_err,[0.5; 5; 0.2; 0.1],0:1e-5:0.2);
nuevaFigura('Name','Error observador','Position',[100 100 700 350]);
plot(te*1e3,e_obs,'LineWidth',1.3); grid on;
legend('e_\theta [rad]','e_\omega [rad/s]','e_i [A]','e_{TL} [N·m]');
xlabel('t [ms]'); ylabel('error');
title('Convergencia del error de estimación, e(0) = [0.5 5 0.2 0.1]^T');
exportgraphics(gcf,fullfile(carpeta,'fig6_error_observador.png'),'Resolution',200);

%% 10. TABLA DE RESULTADOS --------------------------------------------------
imprime('==== 9. RESULTADOS (escalon de pi rad) ====\n');
imprime('%-28s %10s %10s %10s %12s\n','Metrica','Especif.','Polos','LQR','LQR+obs+sat');
imprime('%-28s %10s %10.3f %10.3f %12.3f\n','Tiempo de establec. [s]','<= 0.5',S_pp.SettlingTime,S_lqr.SettlingTime,S_obs.SettlingTime);
imprime('%-28s %10s %10.2f %10.2f %12.2f\n','Sobrepaso [%]','<= 5',S_pp.Overshoot,S_lqr.Overshoot,S_obs.Overshoot);
imprime('%-28s %10s %10.4f %10.4f %12.4f\n','Error final [rad]','0',abs(y_pp(end,1)-r_ref),abs(y_lqr(end,1)-r_ref),abs(X(end,1)-r_ref));
imprime('%-28s %10s %10.2f %10.2f %12.2f\n','|u| maximo [V]','<= 12',max(abs(y_pp(:,2))),max(abs(y_lqr(:,2))),max(abs(U)));
imprime('%-28s %10s %10.2f %10.2f %12.2f\n','|i| maxima [A]','-',max(abs(y_pp(:,3))),max(abs(y_lqr(:,3))),max(abs(X(:,3))));
imprime('%-28s %10s %10.4f %10.4f %12.4f\n','Desv. max. por TL [rad]','-',dev_pp,dev_lqr,dev_obs);
imprime('\nTL estimado al final = %.4f N.m (real = %.2f N.m)\n',XH(end,4),TL);
imprime('Error max. de estimacion de i = %.3f A\n',max(abs(X(:,3)-XH(:,3))));
fclose(fid);
disp('Listo: figuras en IMAGENES/ y resultados en resultados_punto2.txt');

%% FUNCIÓN AUXILIAR ---------------------------------------------------------
function f = nuevaFigura(varargin)
% Crea una figura con tema claro (fondo blanco) para el documento.
f = figure(varargin{:});
try
    theme(f,'light');   % disponible desde MATLAB R2025a
catch
    set(f,'Color','w');
end
end
