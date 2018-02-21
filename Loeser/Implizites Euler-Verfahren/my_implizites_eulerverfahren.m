function [x y] = my_implizites_eulerverfahren(f,dfdy,a,b,y0,h)
% Beschreibung der Variablen
% Eingabe
% f: Rechte Seite der ODE, ein Function-Handle. Wichtig:
% f = @(x,y), d.h. das erste Argument ist die Koordinate.
% dfdy: Ableitung der Funktion f nach y. Wichtig:
% dfdy = @(x,y)
% a,b: Das Intervall, auf dem die ODE geloest werden soll
% y0: Startwert, y(a)
% h: Schrittweite
%
% Ausgabe:
% x,y: Vektoren, y(i) = Funktionswert der Loesung an der
% Stelle x(i)

x = a : h : b;
y = zeros(size(y0, 1), length(x));
y(:, 1) = y0;

for i = 2 : length(x)
    implK = @(K) K - f(x(i - 1) + h, y(:, i - 1) + h * K);
    dimplK = @(K) 1 - h * dfdy(x(i - 1) + h, y(:, i - 1) + h * K);
    K = my_gedaempftes_newton(implK,dimplK,f(x(i - 1), y(:, i - 1)),10^(-8),25);
    y(:, i) = y(:, i - 1) + h * K;
end

end