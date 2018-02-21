function [err,alpha] = Aufgabe_5_1_f()
% Beschreibung der Variablen
% Eingabe: nicht vorhanden
% Ausgabe:
% err: Vektor, err(i) enthaelt den Fehler zur Schrittweite h_i
% alpha: Vektor, alpha(i) enthaelt die Konvergenzordnung, die sich
% aus err(i) und err(i+1) ergibt

a = 0;
b = 3;
h = [1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256];
f = @(x,y) x * sqrt(y);
dfdy = @(x,y) x * 1 / (2 * sqrt(y));
y0 = 1;

err = zeros(1, 8);

for i = 1 : length(h)
    [x, y] = my_implizites_eulerverfahren(f,dfdy,a,b,y0,h(i));
    err(i) = abs(y(end) - 169/16);
end

alpha = diff(log(err))/log(1/2);
end