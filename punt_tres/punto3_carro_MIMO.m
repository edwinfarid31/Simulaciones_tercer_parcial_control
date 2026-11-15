%% =====================================================================
%  PUNTO 3 - Diseño e implementación (en simulación) de un controlador
%  MIMO en espacio de estados
%  Carro eléctrico pequeño de tracción trasera con dos motores JGB37-555B
%  independientes en el eje trasero (diferencial electrónico)
%
%  Autor : Santiago Alexander Bravo Rodriguez - 2224734
%  Curso : Sistemas de Control II - UIS
%
%  Entradas : uL, uR  (voltajes de los motores izquierdo y derecho) [V]
%  Salidas  : v       (velocidad de avance del carro) [m/s]
%             W       (velocidad de giro, yaw rate) [rad/s]
%
%  Contenido:
%   1. Parámetros (motor del punto 2 + carro)
%   2. Modelo completo (4 estados: v, W, iL, iR)
%   3. Modelo reducido de diseño (2 estados: v, W)
%   4. Análisis: polos, controlabilidad, observabilidad, matriz de transferencia
%   5. Discretización y LQR discreto MIMO con acción integral
%   6. Filtro de Kalman discreto con estimación de perturbaciones
%   7. Controlador de comparación: dos PI independientes (uno por rueda)
%   8. Respuesta lineal nominal (matriz de escalones 2x2)
%   9. Simulación no lineal: planta con corrientes, motores desiguales,
%      fricción de Coulomb, saturación, encoders cuantizados, giroscopio
%      con ruido, pendiente y freno en una rueda
%  10. Figuras y tabla de resultados
%
%  Requiere: Control System Toolbox.
%  Las figuras se guardan en IMAGENES/ (prefijo p3_) y los resultados en
%  resultados_punto3.txt
%  =====================================================================
clear; close all; clc;
rng(1);                                   % ruido reproducible

carpeta = fullfile(pwd,'IMAGENES');
if ~exist(carpeta,'dir'); mkdir(carpeta); end
fid = fopen('resultados_punto3.txt','w');
imprime = @(varargin) [fprintf(varargin{:}); fprintf(fid,varargin{:})];

%% 1. PARÁMETROS ------------------------------------------------------------
% --- Motor JGB37-555B 12 V 1:30 (mismos valores del punto 2) ---
mot.R   = 4.0;                                  % [Ohm]
mot.L   = 2.5e-3;                               % [H]
mot.N   = 30;                                   % reducción
w0m     = 333*mot.N*2*pi/60;                    % vel. vacío del motor [rad/s]
mot.Ke  = (12 - mot.R*0.15)/w0m;                % [V.s/rad]
mot.Kt  = mot.Ke;                               % [N.m/A]
mot.Jm  = 5e-6;                                 % [kg.m^2]
mot.b   = mot.Kt*0.15/w0m;                      % [N.m.s/rad] (eje motor)
mot.PPR = 11;                                   % pulsos/vuelta del motor

% --- Carro ---
car.m  = 1.5;            % [kg]    masa total
car.r  = 0.0325;         % [m]     radio de rueda (65 mm)
car.d  = 0.15;           % [m]     trocha (distancia entre ruedas traseras)
car.Lc = 0.25;           % [m]     largo del chasis
car.Jb = car.m*(car.Lc^2 + car.d^2)/12;   % [kg.m^2] inercia de giro (placa)
car.cv = 0.10;           % [N.s/m]   fricción viscosa de avance
car.cw = 0.005;          % [N.m.s/rad] fricción viscosa de giro
car.Jw = 0.5*0.05*car.r^2;                % inercia de una rueda
g = 9.81;

% Inercia de cada rueda incluyendo el rotor reflejado: J_eff = N^2 Jm + Jw
Jeff = mot.N^2*mot.Jm + car.Jw;

% Cinemática: [wL; wR] = T [v; W]
T = [1/car.r, -car.d/(2*car.r);
     1/car.r,  car.d/(2*car.r)];

% Matrices generalizadas: Mb*xdot = T'*tau - (N^2 b T'T + Cb) x
Mb = diag([car.m, car.Jb]) + Jeff*(T.'*T);
Cb = diag([car.cv, car.cw]);

imprime('==== 1. PARAMETROS ====\n');
imprime('Ke = Kt = %.5f V.s/rad,  b = %.4e N.m.s/rad\n',mot.Ke,mot.b);
imprime('Jb = %.4e kg.m^2,  Jeff (por rueda) = %.4e kg.m^2\n',car.Jb,Jeff);
imprime('Masa equivalente  M  = %.3f kg   (carro: %.2f kg)\n',Mb(1,1),car.m);
imprime('Inercia equiv.    Jz = %.4f kg.m^2 (chasis: %.4f)\n\n',Mb(2,2),car.Jb);

%% 2. MODELO COMPLETO (4 ESTADOS) -------------------------------------------
% x = [v; W; iL; iR],  u = [uL; uR]
Mi = inv(Mb);
Af = [ -Mi*(mot.N^2*mot.b*(T.'*T) + Cb),   Mi*T.'*mot.N*mot.Kt;
       -mot.Ke*mot.N/mot.L*T,              -mot.R/mot.L*eye(2) ];
Bf = [zeros(2); eye(2)/mot.L];
Cf = [eye(2) zeros(2)];
sys_full = ss(Af,Bf,Cf,0);
imprime('==== 2. MODELO COMPLETO ====\n');
imprime('Af =\n'); imprime('  %11.4f %11.4f %11.4f %11.4f\n',Af.');
imprime('Bf = diag(1/L) en las filas de corriente: 1/L = %.1f\n',1/mot.L);
imprime('Polos: %s\n\n',mat2str(eig(Af).',6));

%% 3. MODELO REDUCIDO DE DISEÑO (2 ESTADOS) ---------------------------------
% tau_e = L/R = 0.6 ms << Ts = 10 ms  ->  L*di/dt ~ 0:
%   i = (u - Ke N T x)/R
% x = [v; W], u = [uL; uR], y = [v; W]
A = -Mi*( mot.N^2*mot.Kt*mot.Ke/mot.R*(T.'*T) + mot.N^2*mot.b*(T.'*T) + Cb );
B =  Mi*T.'*mot.N*mot.Kt/mot.R;
E =  Mi;                      % entrada de perturbaciones [fuerza; torque]
C =  eye(2);
sys_red = ss(A,B,C,0);
sys_red.InputName  = {'uL','uR'};
sys_red.OutputName = {'v','W'};

imprime('==== 3. MODELO REDUCIDO ====\n');
imprime('A = [%9.4f %9.4f; %9.4f %9.4f]\n',A.');
imprime('B = [%9.4f %9.4f; %9.4f %9.4f]\n',B.');
imprime('Polos reducido: %s\n',mat2str(eig(A).',6));
xmax = -A\B*[12;12];
imprime('Velocidad maxima con 12 V en ambos motores: %.3f m/s\n\n',xmax(1));

% Comparación modelo completo vs reducido (escalón de 6 V en uR)
t_c = (0:1e-4:1.5)';
y_f = lsim(sys_full, [zeros(size(t_c)) 6*ones(size(t_c))], t_c);
y_r = lsim(sys_red,  [zeros(size(t_c)) 6*ones(size(t_c))], t_c);
nuevaFigura('Name','Completo vs reducido','Position',[100 100 750 450]);
subplot(2,1,1); plot(t_c,y_f(:,1),'b',t_c,y_r(:,1),'r--','LineWidth',1.4); grid on;
ylabel('v [m/s]'); legend('Completo (4 estados)','Reducido (2 estados)','Location','southeast');
title('Escalón de 6 V solo en el motor derecho (u_L = 0)');
subplot(2,1,2); plot(t_c,y_f(:,2),'b',t_c,y_r(:,2),'r--','LineWidth',1.4); grid on;
ylabel('\Omega [rad/s]'); xlabel('t [s]');
exportgraphics(gcf,fullfile(carpeta,'p3_fig1_modelos.png'),'Resolution',200);

%% 4. ANÁLISIS --------------------------------------------------------------
Cm = [T; 0 1];               % medidas: [wL_encoder; wR_encoder; W_giroscopio]
imprime('==== 4. ANALISIS ====\n');
imprime('rango(ctrb(A,B))  = %d de 2\n',rank(ctrb(A,B)));
imprime('rango(obsv(A,Cm)) = %d de 2\n',rank(obsv(A,Cm)));
G = tf(sys_red);
imprime('Matriz de transferencia G(s) = Y(s)/U(s):\n');
for i=1:2
    for j=1:2
        [nn,dd] = tfdata(G(i,j),'v');
        imprime('  G%d%d(s) = %s / %s\n',i,j,mat2str(nn(end),5),mat2str(dd,5));
    end
end
G0 = dcgain(sys_red);
imprime('Ganancia DC G(0) = [%8.4f %8.4f; %8.4f %8.4f]\n',G0.');
imprime('(cada voltaje afecta a v y a W: sistema acoplado 2x2)\n\n');

%% 5. DISCRETIZACIÓN Y LQR DISCRETO MIMO CON ACCIÓN INTEGRAL ----------------
Ts = 0.01;                                % periodo de muestreo [s]
sysd = c2d(ss(A,[B E],C,0),Ts,'zoh');
Ad = sysd.A; Bd = sysd.B(:,1:2); Ed = sysd.B(:,3:4);

% Integradores discretos del error: z[k+1] = z[k] + Ts (r - y)
Aa = [Ad zeros(2); -Ts*eye(2) eye(2)];
Ba = [Bd; zeros(2)];
imprime('==== 5. LQR DISCRETO MIMO ====\n');
imprime('Ts = %.3f s\n',Ts);
imprime('rango(ctrb aumentado) = %d de 4\n',rank(ctrb(Aa,Ba)));

Q  = diag([20 20 600 1000]);   % pesos: v, W, int(e_v), int(e_W)
Rl = 0.01*eye(2);              % peso de los voltajes
[K,~,pz] = dlqr(Aa,Ba,Q,Rl);
ps = log(pz)/Ts;               % polos continuos equivalentes
imprime('Q = diag(%s), R = %g*I\n',mat2str(diag(Q).'),Rl(1));
imprime('K (2x4) = [Kx | Kz]:\n');
imprime('  [%9.4f %9.4f | %10.4f %10.4f]\n',K.');
imprime('Polos z: %s\n',mat2str(pz.',5));
imprime('Polos s equivalentes: %s\n\n',mat2str(ps.',5));

%% 6. FILTRO DE KALMAN DISCRETO CON ESTIMACIÓN DE PERTURBACIONES -------------
% xe = [v; W; av; aW]: av, aW = aceleraciones debidas a perturbaciones
% (pendiente, freno), modeladas como caminata aleatoria.
Ae  = [A eye(2); zeros(2,4)];
Be  = [B; zeros(2)];
sye = c2d(ss(Ae,Be,eye(4),0),Ts,'zoh');
Ade = sye.A;  Bde = sye.B;
Ce  = [Cm zeros(3,2)];
q_enc  = 2*pi/(mot.PPR*4*mot.N)/Ts;       % resolución de velocidad del encoder
sig_g  = 0.02;                            % ruido del giroscopio [rad/s]
Rk = diag([q_enc^2/3, q_enc^2/3, sig_g^2]);
Qk = diag([1e-5 1e-5 0.05 0.05]);
[Mk,Pk] = dlqe(Ade,eye(4),Ce,Qk,Rk);      % ganancia del estimador actual
imprime('==== 6. KALMAN ====\n');
imprime('rango(obsv(Ade,Ce)) = %d de 4\n',rank(obsv(Ade,Ce)));
imprime('Resolucion de velocidad del encoder = %.4f rad/s por rueda\n',q_enc);
imprime('Polos del estimador |z|: %s\n',mat2str(abs(eig((eye(4)-Mk*Ce)*Ade)).',4));
imprime('M (4x3) =\n'); imprime('  %9.4f %9.4f %9.4f\n',Mk.');
imprime('\n');

%% 7. CONTROLADOR DE COMPARACIÓN: DOS PI INDEPENDIENTES ---------------------
% Cada rueda sigue su velocidad de referencia [wL*; wR*] = T [v*; W*]
pid.kp = 0.3;  pid.ki = 3;
imprime('==== 7. PI POR RUEDA ====\n');
imprime('Kp = %.2f V.s/rad, Ki = %.2f V/rad (iguales en ambas ruedas)\n\n',pid.kp,pid.ki);

%% 8. RESPUESTA LINEAL NOMINAL (MATRIZ DE ESCALONES 2x2) --------------------
% Lazo cerrado con realimentación de estados (modelo nominal, sin ruido)
Acl = Aa - Ba*K;
Bcl = [zeros(2); Ts*eye(2)];
syscl = ss(Acl,Bcl,[eye(2) zeros(2)],0,Ts);
syscl.InputName = {'v*','W*'}; syscl.OutputName = {'v','W'};
nuevaFigura('Name','Escalones MIMO','Position',[100 100 750 550]);
step(syscl,3); grid on;
title('Respuesta al escalón del lazo cerrado LQR (modelo nominal)');
exportgraphics(gcf,fullfile(carpeta,'p3_fig3_escalones_mimo.png'),'Resolution',200);

nuevaFigura('Name','Polos','Position',[100 100 650 450]);
plot(real(eig(A)),imag(eig(A)),'kx','MarkerSize',11,'LineWidth',2); hold on;
plot(real(ps),imag(ps),'rs','MarkerSize',9,'LineWidth',1.5);
grid on; xlabel('Re'); ylabel('Im'); xlim([-200 5]); ylim([-5 5]);
legend('Lazo abierto (reducido)','LQR (equivalente continuo)','Location','best');
title('Polos de lazo abierto y de lazo cerrado');
exportgraphics(gcf,fullfile(carpeta,'p3_fig2_polos.png'),'Resolution',200);

%% 9. SIMULACIÓN NO LINEAL --------------------------------------------------
P.mot = mot; P.car = car; P.g = g; P.T = T; P.Mi = Mi; P.Cb = Cb;
P.Ts = Ts; P.Ade = Ade; P.Bde = Bde; P.Ce = Ce; P.Mk = Mk;
P.K = K; P.pid = pid; P.sig_g = sig_g;
% Motores desiguales en la planta (tolerancia de fabricación):
P.Kt_r = 0.92; P.R_r = 1.05;     % motor derecho: 8 % menos Kt/Ke, 5 % más R

imprime('==== 8. SIMULACION NO LINEAL ====\n');
imprime('Planta: motor derecho con Kt,Ke x%.2f y R x%.2f; Coulomb; sat. +-12 V\n',P.Kt_r,P.R_r);

% --- Prueba A: escalón de velocidad (v* = 0.5 m/s, W* = 0) ---
LA_lqr = simular_carro('lqr','vstep',3,P);
LA_pid = simular_carro('pid','vstep',3,P);
% --- Prueba B: escalón de giro (W* = 1 rad/s a v* = 0.5 m/s) ---
LB_lqr = simular_carro('lqr','wstep',4,P);
LB_pid = simular_carro('pid','wstep',4,P);
% --- Prueba C: recorrido completo con perturbaciones ---
LC_lqr = simular_carro('lqr','full',11,P);
LC_pid = simular_carro('pid','full',11,P);

% Métricas
[osvL,tsvL] = metricas(LA_lqr.t,LA_lqr.v,0.5,0,0.5,0.02);
[osvP,tsvP] = metricas(LA_pid.t,LA_pid.v,0.5,0,0.5,0.02);
[oswL,tswL] = metricas(LB_lqr.t,LB_lqr.W,2,0,1,0.05);
[oswP,tswP] = metricas(LB_pid.t,LB_pid.W,2,0,1,0.05);
WpkL = max(abs(LA_lqr.W));  WpkP = max(abs(LA_pid.W));
vdL  = max(abs(LB_lqr.v(LB_lqr.t>=2)-0.5));  vdP = max(abs(LB_pid.v(LB_pid.t>=2)-0.5));
psi_ref = 2.0;   % 1 rad/s durante 2 s
iS = LC_lqr.t>=6.5 & LC_lqr.t<8.5;  iF = LC_lqr.t>=8.5;
dvsL = max(abs(LC_lqr.v(iS)-0.5));  dvsP = max(abs(LC_pid.v(iS)-0.5));
dwfL = max(abs(LC_lqr.W(iF)));      dwfP = max(abs(LC_pid.W(iF)));
umL = max(abs([LA_lqr.u(:);LB_lqr.u(:);LC_lqr.u(:)]));
umP = max(abs([LA_pid.u(:);LB_pid.u(:);LC_pid.u(:)]));

%% 10. FIGURAS --------------------------------------------------------------
% Prueba A
nuevaFigura('Name','Prueba A','Position',[100 100 750 600]);
subplot(3,1,1); plot(LA_lqr.t,LA_lqr.v,'r',LA_pid.t,LA_pid.v,'b','LineWidth',1.3); hold on;
plot(LA_lqr.t,LA_lqr.ref(:,1),'k--'); grid on; ylabel('v [m/s]');
legend('LQR MIMO','2 PI por rueda','Referencia','Location','southeast');
title('Prueba A: escalón de velocidad v* = 0.5 m/s, \Omega* = 0');
subplot(3,1,2); plot(LA_lqr.t,LA_lqr.W,'r',LA_pid.t,LA_pid.W,'b','LineWidth',1.3); grid on;
ylabel('\Omega [rad/s]');
subplot(3,1,3); plot(LA_lqr.t,LA_lqr.u,'r',LA_pid.t,LA_pid.u,'b','LineWidth',1.1); grid on;
ylabel('u_L, u_R [V]'); xlabel('t [s]');
exportgraphics(gcf,fullfile(carpeta,'p3_fig4_prueba_A.png'),'Resolution',200);

% Prueba B
nuevaFigura('Name','Prueba B','Position',[100 100 750 600]);
subplot(3,1,1); plot(LB_lqr.t,LB_lqr.W,'r',LB_pid.t,LB_pid.W,'b','LineWidth',1.3); hold on;
plot(LB_lqr.t,LB_lqr.ref(:,2),'k--'); grid on; ylabel('\Omega [rad/s]');
legend('LQR MIMO','2 PI por rueda','Referencia','Location','southeast');
title('Prueba B: escalón de giro \Omega* = 1 rad/s a v* = 0.5 m/s');
subplot(3,1,2); plot(LB_lqr.t,LB_lqr.v,'r',LB_pid.t,LB_pid.v,'b','LineWidth',1.3); grid on;
ylabel('v [m/s]');
subplot(3,1,3); plot(LB_lqr.t,LB_lqr.u,'r',LB_pid.t,LB_pid.u,'b','LineWidth',1.1); grid on;
ylabel('u_L, u_R [V]'); xlabel('t [s]');
exportgraphics(gcf,fullfile(carpeta,'p3_fig5_prueba_B.png'),'Resolution',200);

% Prueba C
nuevaFigura('Name','Prueba C','Position',[100 100 750 700]);
subplot(4,1,1); plot(LC_lqr.t,LC_lqr.v,'r',LC_pid.t,LC_pid.v,'b','LineWidth',1.2); hold on;
plot(LC_lqr.t,LC_lqr.ref(:,1),'k--'); grid on; ylabel('v [m/s]');
xline(6.5,':','Pendiente 10°'); xline(8.5,':','Freno rueda der.');
legend('LQR MIMO','2 PI por rueda','Referencia','Location','southwest');
title('Prueba C: recorrido completo con perturbaciones');
subplot(4,1,2); plot(LC_lqr.t,LC_lqr.W,'r',LC_pid.t,LC_pid.W,'b','LineWidth',1.2); hold on;
plot(LC_lqr.t,LC_lqr.ref(:,2),'k--'); grid on; ylabel('\Omega [rad/s]');
subplot(4,1,3); plot(LC_lqr.t,LC_lqr.u,'LineWidth',1.1); grid on; ylabel('u LQR [V]');
legend('u_L','u_R','Location','northwest');
subplot(4,1,4); plot(LC_pid.t,LC_pid.u,'LineWidth',1.1); grid on; ylabel('u PI [V]'); xlabel('t [s]');
exportgraphics(gcf,fullfile(carpeta,'p3_fig6_prueba_C.png'),'Resolution',200);

% Kalman
nuevaFigura('Name','Kalman','Position',[100 100 750 650]);
subplot(4,1,1); plot(LC_lqr.t,LC_lqr.v,'r',LC_lqr.t,LC_lqr.xe(:,1),'b--','LineWidth',1.1); grid on;
ylabel('v [m/s]'); legend('Real','Estimado (Kalman)','Location','southeast');
title('Filtro de Kalman en la prueba C (controlador LQR)');
subplot(4,1,2); plot(LC_lqr.t,LC_lqr.Wmeas,'Color',[0.75 0.75 0.75]); hold on;
plot(LC_lqr.t,LC_lqr.W,'r',LC_lqr.t,LC_lqr.xe(:,2),'b--','LineWidth',1.1); grid on;
ylabel('\Omega [rad/s]'); legend('Giroscopio (medida)','Real','Estimado','Location','northeast');
subplot(4,1,3); plot(LC_lqr.t,LC_lqr.wmeas(:,2),'Color',[0.75 0.75 0.75]); hold on;
plot(LC_lqr.t,LC_lqr.wR,'r','LineWidth',1.1); grid on;
ylabel('\omega_R [rad/s]'); legend('Encoder (medida)','Real','Location','southeast');
subplot(4,1,4); plot(LC_lqr.t,LC_lqr.xe(:,3),'b',LC_lqr.t,LC_lqr.xe(:,4),'m','LineWidth',1.1); grid on;
ylabel('acel. perturb.'); xlabel('t [s]');
legend('a_v estimada [m/s^2]','a_\Omega estimada [rad/s^2]','Location','southwest');
exportgraphics(gcf,fullfile(carpeta,'p3_fig7_kalman.png'),'Resolution',200);

% Trayectoria
[Xr,Yr] = trayectoria_ideal(LC_lqr.t,LC_lqr.ref,Ts);
nuevaFigura('Name','Trayectoria','Position',[100 100 650 550]);
plot(Xr,Yr,'k--','LineWidth',1.3); hold on;
plot(LC_lqr.pos(:,1),LC_lqr.pos(:,2),'r','LineWidth',1.4);
plot(LC_pid.pos(:,1),LC_pid.pos(:,2),'b','LineWidth',1.4);
axis equal; grid on; xlabel('X [m]'); ylabel('Y [m]');
legend('Ideal (referencias)','LQR MIMO','2 PI por rueda','Location','best');
title('Trayectoria del carro en la prueba C');
exportgraphics(gcf,fullfile(carpeta,'p3_fig8_trayectoria.png'),'Resolution',200);

%% TABLA DE RESULTADOS ------------------------------------------------------
imprime('\n==== 9. RESULTADOS ====\n');
imprime('%-42s %10s %10s %10s\n','Metrica','Especif.','LQR MIMO','2 PI');
imprime('%-42s %10s %10.3f %10.3f\n','A: ts de v (2%) [s]','<= 1.0',tsvL,tsvP);
imprime('%-42s %10s %10.2f %10.2f\n','A: sobrepaso de v [%]','<= 5',osvL,osvP);
imprime('%-42s %10s %10.3f %10.3f\n','A: |W| max. inducido por el escalon de v','-',WpkL,WpkP);
imprime('%-42s %10s %10.3f %10.3f\n','B: ts de W (5%) [s]','<= 1.0',tswL,tswP);
imprime('%-42s %10s %10.2f %10.2f\n','B: sobrepaso de W [%]','<= 10',oswL,oswP);
imprime('%-42s %10s %10.4f %10.4f\n','B: |v - v*| max. durante el giro [m/s]','-',vdL,vdP);
imprime('%-42s %10s %10.4f %10.4f\n','C: |v - v*| max. por la pendiente [m/s]','-',dvsL,dvsP);
imprime('%-42s %10s %10.3f %10.3f\n','C: |W| max. por el freno en una rueda','-',dwfL,dwfP);
imprime('%-42s %10s %10.3f %10.3f\n','C: rumbo final psi [rad] (ideal 2.000)','2.000',LC_lqr.pos(end,3),LC_pid.pos(end,3));
imprime('%-42s %10s %10.2f %10.2f\n','|u| maximo en todas las pruebas [V]','<= 12',umL,umP);
fclose(fid);
disp('Listo: figuras p3_*.png en IMAGENES/ y resultados en resultados_punto3.txt');

%% =========================================================================
%  FUNCIONES LOCALES
%  =========================================================================
function L = simular_carro(tipo,perfil,Tend,P)
% Simula la planta no lineal (Euler, dt = 0.1 ms) con el controlador
% discreto (Ts = 10 ms) indicado: 'lqr' o 'pid'.
mot = P.mot; car = P.car; T = P.T; Ts = P.Ts;
dt = 1e-4; nsub = round(Ts/dt); nk = round(Tend/Ts);
Kt = mot.Kt*[1 P.Kt_r]; Ke = mot.Ke*[1 P.Kt_r]; R = mot.R*[1 P.R_r];
qth = 2*pi/(mot.PPR*4*mot.N);           % resolución angular del encoder
x  = zeros(4,1);  th = zeros(2,1);  cntp = zeros(2,1);
pos = zeros(3,1); xe = zeros(4,1);  z = zeros(2,1);  integ = zeros(2,1);
L.t = zeros(nk,1); L.ref = zeros(nk,2); L.v = zeros(nk,1); L.W = zeros(nk,1);
L.u = zeros(nk,2); L.xe = zeros(nk,4); L.pos = zeros(nk,3);
L.wmeas = zeros(nk,2); L.Wmeas = zeros(nk,1); L.wR = zeros(nk,1);
for k = 1:nk
    t = (k-1)*Ts;
    switch perfil
        case 'vstep', ref = [0.5*(t>=0.5); 0];
        case 'wstep', ref = [0.5; 1.0*(t>=2)];
        case 'full',  ref = [0.5*(t>=0.5); 1.0*(t>=3 && t<5)];
    end
    % --- Sensores ---
    cnt   = floor(th/qth);                       % cuentas de los encoders
    wmeas = (cnt - cntp)*qth/Ts;  cntp = cnt;    % velocidad por diferencia
    Wg    = x(2) + P.sig_g*randn;                % giroscopio con ruido
    y     = [wmeas; Wg];
    % --- Estimador (Kalman, forma de estimador actual) ---
    xe = xe + P.Mk*(y - P.Ce*xe);
    % --- Ley de control ---
    if strcmp(tipo,'lqr')
        u = -P.K(:,1:2)*xe(1:2) - P.K(:,3:4)*z;
        z = z + Ts*(ref - xe(1:2));
    else
        e = T*ref - T*xe(1:2);                   % error de velocidad de cada rueda
        integ = integ + Ts*e;
        u = P.pid.kp*e + P.pid.ki*integ;
    end
    u = min(max(u,-12),12);                      % saturación del driver
    % registro
    L.t(k) = t; L.ref(k,:) = ref.'; L.v(k) = x(1); L.W(k) = x(2);
    L.u(k,:) = u.'; L.xe(k,:) = xe.'; L.pos(k,:) = pos.';
    L.wmeas(k,:) = wmeas.'; L.Wmeas(k) = Wg; w_ = T*x(1:2); L.wR(k) = w_(2);
    % predicción del estimador
    xe = P.Ade*xe + P.Bde*u;
    % --- Planta no lineal entre muestras ---
    for s = 1:nsub
        ts_ = t + (s-1)*dt;
        w = T*x(1:2);                            % velocidades de las ruedas
        Fd = [0;0]; tauw = [0;0];
        if strcmp(perfil,'full')
            if ts_ >= 6.5, Fd(1) = -car.m*P.g*sind(10); end   % pendiente 10°
            if ts_ >= 8.5, tauw(2) = 0.15; end               % freno rueda der.
        end
        tauc = 0.01*tanh(w/0.5);                 % fricción de Coulomb en la reductora
        Fc   = 0.02*car.m*P.g*tanh(x(1)/0.01);   % resistencia a la rodadura
        tau  = mot.N*(Kt.'.*x(3:4)) - mot.N^2*mot.b*w - tauc - tauw;
        dvw  = P.Mi*(T.'*tau - P.Cb*x(1:2) + Fd - [Fc;0]);
        di   = (u - R.'.*x(3:4) - Ke.'.*mot.N.*w)/mot.L;
        x    = x + dt*[dvw; di];
        th   = th + dt*w;
        pos  = pos + dt*[x(1)*cos(pos(3)); x(1)*sin(pos(3)); x(2)];
    end
end
end

function [os,ts] = metricas(t,y,t0,y0,yf,banda)
% Sobrepaso [%] y tiempo de establecimiento con la banda indicada.
i  = t >= t0;  tt = t(i) - t0;  yy = y(i);
os = max(0,(max(yy)-yf)/abs(yf-y0)*100);
fuera = find(abs(yy-yf) > banda*abs(yf-y0));
if isempty(fuera), ts = 0; else, ts = tt(fuera(end)); end
end

function [X,Y] = trayectoria_ideal(t,ref,Ts)
X = zeros(size(t)); Y = X; psi = 0; px = 0; py = 0;
for k = 1:numel(t)
    X(k) = px; Y(k) = py;
    px = px + Ts*ref(k,1)*cos(psi); py = py + Ts*ref(k,1)*sin(psi);
    psi = psi + Ts*ref(k,2);
end
end

function f = nuevaFigura(varargin)
f = figure(varargin{:});
try, theme(f,'light'); catch, set(f,'Color','w'); end
end
