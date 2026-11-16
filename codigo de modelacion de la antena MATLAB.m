 %% =========================================================================
%  MODELO DE ANTENA PCB 2.4 GHz (Reconstruccion visual segun foto)
%  Placa: 38 mm x 6 mm | Requiere: Antenna Toolbox (pcbStack)
%% =========================================================================
clear; clc; close all;

%% 1. PARAMETROS GENERALES Y PROPORCIONES (AJUSTAR AQUI) ------------------
P.f0   = 2.44e9;          % frecuencia central BT/Wi-Fi [Hz]
P.fmin = 2.10e9;          % barrido [Hz]
P.fmax = 2.80e9;
P.Npts = 21;              % puntos del barrido
P.Z0   = 50;              % impedancia de referencia [ohm]

% Dimensiones generales de la placa PCB [m]
P.L    = 38e-3;           % largo total (38 mm)
P.W    = 6e-3;            % ancho total (6 mm)
P.h    = 0.8e-3;          % espesor FR-4
P.epsR = 4.4;             % permitividad FR-4
P.tanD = 0.02;            % tangente de perdidas FR-4

% --- PROPORCIONES VISUALES ("A OJO" SEGUN LA FOTO) -----------------------
% Margenes al borde de la placa
P.m_bot   = 2.5e-3;       % borde inferior sin cobre (donde entra el cable)
P.m_top   = 0.5e-3;       % margen superior en la punta
P.m_lat   = 0.5e-3;       % margen lateral izquierdo y derecho

% Zona inferior de pads (donde va soldado el coaxial)
P.L_pad1  = 4.2e-3;       % largo pad inferior (malla / masa)
P.L_pad2  = 3.8e-3;       % largo pad medio (vivo / alimentacion)
P.w_pad   = 1.8e-3;       % ancho de los pads estanados derechos
P.gap_h   = 0.6e-3;       % separacion horizontal entre bloques

% Zona superior radiante (brazos y dedo central)
P.w_izq   = 2.1e-3;       % ancho brazo grueso izquierdo
P.w_cen   = 1.0e-3;       % ancho dedo vertical central
P.w_der   = 1.0e-3;       % ancho trazo vertical derecho
P.w_top   = 1.2e-3;       % grosor barra horizontal en la punta
P.gap_v   = 0.45e-3;      % ranuras verticales entre brazos
P.gap_tip = 0.65e-3;      % separacion entre punta del dedo central y barra sup.
% -------------------------------------------------------------------------

P.etaRad = 0.6;           % eficiencia de radiacion estimada

%% 2. CONSTRUCCION Y VISUALIZACION DE GEOMETRIA ---------------------------
ant = construirAntenaVisual(P);
figure('Name','Geometria PCB');
show(ant);
title('Geometria de la antena PCB (vista superior)');

%% 3. BARRIDO EN FRECUENCIA: S11, SWR, Zin y Carta de Smith ---------------
freq  = linspace(P.fmin, P.fmax, P.Npts);
s     = sparameters(ant, freq, P.Z0);
S11   = squeeze(rfparam(s,1,1));
S11dB = 20*log10(abs(S11));
SWRf  = (1+abs(S11))./(1-abs(S11));
Zf    = impedance(ant, freq);

figure('Name','S11 y SWR');
subplot(2,1,1);
plot(freq/1e9, S11dB, 'LineWidth',1.5); grid on; hold on;
yline(-9.54,'r--','SWR = 2'); xline(P.f0/1e9,'k:');
xlabel('Frecuencia [GHz]'); ylabel('S_{11} [dB]'); title('Perdida de retorno');

subplot(2,1,2);
plot(freq/1e9, SWRf, 'LineWidth',1.5); grid on; hold on;
yline(2,'r--'); xline(P.f0/1e9,'k:'); ylim([1 10]);
xlabel('Frecuencia [GHz]'); ylabel('SWR');

figure('Name','Impedancia de entrada');
plot(freq/1e9, real(Zf), 'LineWidth',1.5); hold on;
plot(freq/1e9, imag(Zf), '--', 'LineWidth',1.5); grid on;
legend('R [\Omega]','X [\Omega]'); xlabel('Frecuencia [GHz]');
title('Z_{in} = R + jX'); xline(P.f0/1e9,'k:');

figure('Name','Carta de Smith');
try
    smithplot(s);
    title('Carta de Smith (S_{11})');
catch
    polarplot(angle(S11), abs(S11), '-o', 'LineWidth', 1.5);
    rlim([0 1]);
    title('Coeficiente de Reflexion \Gamma');
end

%% 4. PARAMETROS A f0 -----------------------------------------------------
Zin  = impedance(ant, P.f0);
Gam  = (Zin - P.Z0)/(Zin + P.Z0);
SWR0 = (1+abs(Gam))/(1-abs(Gam));
RL0  = -20*log10(abs(Gam));

ok = S11dB <= -9.54;
if any(ok)
    fi = freq(ok);  BW = (max(fi)-min(fi));
    bwTxt = sprintf('%.0f MHz (%.3f - %.3f GHz)', BW/1e6, min(fi)/1e9, max(fi)/1e9);
else
    bwTxt = 'sin banda con SWR<2 en el barrido (ajustar P.gap_tip o largos)';
end

%% 5. PATRON DE RADIACION -------------------------------------------------
figure('Name','Patron 3D'); pattern(ant, P.f0);
title(sprintf('Directividad 3D a %.2f GHz [dBi]', P.f0/1e9));

figure('Name','Cortes 2D');
subplot(1,2,1); patternElevation(ant, P.f0, 0);   title('Elevacion (az = 0)');
subplot(1,2,2); patternAzimuth(ant, P.f0, 0);     title('Azimut (el = 0)');

az = 0:5:360;  el = -90:5:90;
D  = pattern(ant, P.f0, az, el);
Dmax = max(D(:));
iel0 = find(el==0);
cutAz = D(iel0, :);
iaz0  = find(az==0);  iaz180 = find(az==180);
cutEl = [D(:,iaz0); flipud(D(:,iaz180))]';
angEl = [el, fliplr(180-el)];
HPBW_az = hpbw(az, cutAz);
HPBW_el = hpbw(angEl, cutEl);

%% 6. GANANCIA, APERTURA EFECTIVA ----------------------------------------
lambda = 3e8/P.f0;
Greal  = Dmax + 10*log10(P.etaRad) + 10*log10(1-abs(Gam)^2);
G_lin  = 10^((Dmax + 10*log10(P.etaRad))/10);
Ae     = lambda^2*G_lin/(4*pi);

%% 7. RESUMEN -------------------------------------------------------------
fprintf('================ RESULTADOS a %.2f GHz ================\n', P.f0/1e9);
fprintf('Zin                     : %.1f %+.1fj ohm\n', real(Zin), imag(Zin));
fprintf('|Gamma|                 : %.3f\n', abs(Gam));
fprintf('SWR                     : %.2f   (ideal 1, aceptable <= 2)\n', SWR0);
fprintf('Return Loss             : %.1f dB\n', RL0);
fprintf('Ancho de banda SWR<2    : %s\n', bwTxt);
fprintf('Directividad max        : %.2f dBi  (dipolo lambda/2 = 2.15 dBi)\n', Dmax);
fprintf('HPBW azimut (el=0)      : %.0f deg\n', HPBW_az);
fprintf('HPBW plano xz           : %.0f deg\n', HPBW_el);
fprintf('Eficiencia rad. supuesta: %.0f %%\n', P.etaRad*100);
fprintf('Ganancia realizada      : %.2f dBi\n', Greal);
fprintf('Apertura efectiva Ae    : %.2f cm^2  (lambda^2 = %.1f cm^2)\n', Ae*1e4, lambda^2*1e4);

%% ======================= FUNCIONES LOCALES ===============================
function ant = construirAntenaVisual(P)
    rect = @(len,wid,cx,cy) antenna.Rectangle('Length',len,'Width',wid,'Center',[cx cy]);

    % Convencion de ejes:
    %   X: -L/2 (base donde entra el cable)  -->  +L/2 (punta superior)
    %   Y: -W/2 (borde derecho de la foto)   -->  +W/2 (borde izquierdo de la foto)
    y_der  = -P.W/2 + P.m_lat;
    y_izq  =  P.W/2 - P.m_lat;
    W_util =  y_izq - y_der;

    % --- 1) FILA INFERIOR DE PADS (malla coaxial / masa) ---
    x1_ini = -P.L/2 + P.m_bot;
    x1_fin =  x1_ini + P.L_pad1;
    w_pad_izq = W_util - P.w_pad - P.gap_v;

    pad1_der = rect(P.L_pad1, P.w_pad,     (x1_ini+x1_fin)/2, y_der + P.w_pad/2);
    pad1_izq = rect(P.L_pad1, w_pad_izq,   (x1_ini+x1_fin)/2, y_izq - w_pad_izq/2);

    % --- 2) FILA MEDIA DE PADS (vivo coaxial / alimentacion) ---
    x2_ini = x1_fin + P.gap_h;
    x2_fin = x2_ini + P.L_pad2;

    pad2_der = rect(P.L_pad2, P.w_pad,     (x2_ini+x2_fin)/2, y_der + P.w_pad/2);
    pad2_izq = rect(P.L_pad2, w_pad_izq,   (x2_ini+x2_fin)/2, y_izq - w_pad_izq/2);

    % --- 3) ESTRUCTURA RADIANTE SUPERIOR ---
    x_rad_ini = x2_fin + P.gap_h;          % inicio brazo grueso tras ranura
    x_tip     = P.L/2 - P.m_top;           % borde superior en la punta

    % a) Trazo vertical derecho (sale directo de pad2_der hasta la punta)
    traza_der = rect(x_tip - (x2_fin - 0.05e-3), P.w_der, ...
                     (x_tip + x2_fin - 0.05e-3)/2, y_der + P.w_der/2);

    % b) Barra horizontal superior en la punta
    barra_sup = rect(P.w_top, W_util, x_tip - P.w_top/2, (y_izq + y_der)/2);

    % c) Brazo vertical ancho izquierdo
    brazo_izq = rect(x_tip - x_rad_ini, P.w_izq, ...
                     (x_tip + x_rad_ini)/2, y_izq - P.w_izq/2);

    % d) Puente inferior y dedo vertical central abierto en el interior
    y_cen_izq = (y_izq - P.w_izq) - P.gap_v;
    y_cen_der = y_cen_izq - P.w_cen;
    x_cen_fin = (x_tip - P.w_top) - P.gap_tip;

    puente_inf = rect(1.1e-3, (y_izq - P.w_izq + 0.05e-3) - y_cen_der, ...
                      x_rad_ini + 0.55e-3, ((y_izq - P.w_izq + 0.05e-3) + y_cen_der)/2);

    dedo_cen   = rect(x_cen_fin - x_rad_ini, P.w_cen, ...
                      (x_cen_fin + x_rad_ini)/2, (y_cen_izq + y_cen_der)/2);

    % Union de toda la capa superior de cobre
    top = pad1_der + pad1_izq + pad2_der + pad2_izq + ...
          traza_der + barra_sup + brazo_izq + puente_inf + dedo_cen;

    % Plano de referencia inferior bajo la zona de los pads del cable
    Lgnd = x2_fin - (-P.L/2);
    gnd  = rect(Lgnd, P.W, -P.L/2 + Lgnd/2, 0);

    % Sustrato FR-4
    d = dielectric('FR4');
    d.EpsilonR = P.epsR;  d.LossTangent = P.tanD;  d.Thickness = P.h;

    % Ensamblaje PCB
    ant = pcbStack;
    ant.BoardShape     = antenna.Rectangle('Length',P.L,'Width',P.W);
    ant.BoardThickness = P.h;
    ant.Layers         = {top, d, gnd};

    % Puerto en pad2_der (vivo) y via a masa en pad1_der (malla)
    x_feed = (x2_ini + x2_fin)/2;
    x_gnd  = (x1_ini + x1_fin)/2;
    y_port = y_der + P.w_pad/2;

    ant.FeedLocations  = [x_feed y_port 1 3];
    ant.FeedDiameter   = 0.4e-3;
    ant.ViaLocations   = [x_gnd  y_port 1 3];
    ant.ViaDiameter    = 0.4e-3;
end

function bw = hpbw(ang, g)
    [gm, im] = max(g);  n = numel(g);
    thr = gm - 3;
    r = 0; while r < n-1 && g(mod(im+r, n)+1) >= thr, r = r+1; end
    l = 0; while l < n-1 && g(mod(im-l-2, n)+1) >= thr, l = l+1; end
    dang = mean(diff(ang));
    bw = (r + l + 1)*abs(dang);
    if bw >= 359, bw = NaN; end
end
