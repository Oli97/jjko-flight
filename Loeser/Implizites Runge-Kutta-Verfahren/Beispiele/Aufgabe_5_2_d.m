function [err,alpha] = Aufgabe_5_2_d()
% Beschreibung der Variablen
% Eingabe: nicht vorhanden
% Ausgabe:
% err: Vektor, err(i) enthaelt den Fehler zur Schrittweite h_i
% alpha: Vektor, alpha(i) enthaelt die Konvergenzordnung, die sich
% aus err(i) und err(i+1) ergibt

a = 0;
b = 2;
h = [1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128, 1/256];
f = @(x,y) [-21, 19, -20; 19, -21, 20;40, -40, -40] * y;
dfdy = @(x,y) [-21, 19, -20; 19, -21, 20;40, -40, -40];
y0 = [1;0;-1];

err = zeros(1, length(h));

for i = 1 : length(h)
    [x, y] = my_implizites_runge_kutta(f,dfdy,a,b,y0,h(i));
    err(i) = norm(y(:, end) - [1/2*exp(-4)+1/2*exp(-80)*(cos(80)+sin(80));1/2*exp(-4)-1/2*exp(-80)*(cos(80)+sin(80));-exp(-80)*(cos(80)-sin(80))]);
end

alpha = diff(log(err))/log(1/2);
end