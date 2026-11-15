function capturar_figuras(carpeta_salida)
% CAPTURAR_FIGURAS  Guarda una captura de pantalla de cada ventana de figura
% abierta (incluida la barra de herramientas de MATLAB) como PNG.
%
%   capturar_figuras('CAPTURAS')
%
% Uso: ejecute primero el script de simulación (que deja las figuras
% abiertas) y luego esta función.

if nargin < 1, carpeta_salida = 'CAPTURAS'; end
if ~exist(carpeta_salida,'dir'), mkdir(carpeta_salida); end

figs = findall(0,'Type','figure');
[~,orden] = sort([figs.Number]);  figs = figs(orden);

pant  = get(0,'ScreenSize');                         % [1 1 ancho alto] (MATLAB)
tk    = java.awt.Toolkit.getDefaultToolkit();
jsize = tk.getScreenSize();
esc   = double(jsize.getWidth())/pant(3);            % escala MATLAB -> Java
robot = java.awt.Robot();

for k = 1:numel(figs)
    f = figs(k);
    f.Units = 'pixels';
    % Ajustar la ventana para que quepa completa en la pantalla
    op = f.OuterPosition;
    w  = min(op(3), 790);
    h  = min(op(4), pant(4) - 80);
    f.OuterPosition = [30, pant(4) - h - 40, w, h];
    figure(f); drawnow; pause(1.2);
    op = f.OuterPosition;
    % Coordenadas de pantalla para Java (origen arriba a la izquierda)
    x  = round((op(1)-1)*esc);
    y  = round((pant(4) - (op(2)-1) - op(4))*esc);
    rw = round(op(3)*esc);  rh = round(op(4)*esc);
    img = robot.createScreenCapture(java.awt.Rectangle(x, y, rw, rh));
    nombre = regexprep(lower(f.Name),'[^a-z0-9]+','_');
    if isempty(nombre), nombre = sprintf('figura_%d',f.Number); end
    archivo = fullfile(carpeta_salida, sprintf('cap%02d_%s.png',k,nombre));
    javax.imageio.ImageIO.write(img,'png',java.io.File(fullfile(pwd,archivo)));
    fprintf('Captura guardada: %s\n',archivo);
end
end
