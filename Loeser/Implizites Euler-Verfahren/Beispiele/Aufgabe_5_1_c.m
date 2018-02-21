function [] = Aufgabe_5_1_c()
% Beschreibung der Variablen
% Eingabe: nicht vorhanden
% Ausgabe:
% nicht Vorhanden

a = 0;
b = 10;
h = 1/40;
f = @(x,y) [10 * y(1) * (1 - y(2)); y(2) * (y(1) - 1)];
dfdy = @(x,y) [10 * (1 - y(2)), -10 * y(1); y(2), (y(1) - 1)];
y0 = [[3;1],[4;2],[3;3]];

for i = 1 : 3
    [x, y] = my_implizites_eulerverfahren(f,dfdy,a,b,y0(:, i),h);
    per = 100;
     for m = 2 : 400
        if (norm(y(:, m) - y0(:,i)) < norm(y(:, m - 1) - y0(:,i)))
            per = x(m);
            if (norm(y(:, m + 1) - y0(:,i)) > norm(y(:, m) - y0(:,i)))
                break;
            end
        end
    end
    per
end