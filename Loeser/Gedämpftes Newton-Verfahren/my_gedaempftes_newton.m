function x = my_gedaempftes_newton(f,df,x0,tol,maxiter)
% Beschreibung der Variablen
% Eingabe
%  f: Funktion, deren Nullstelle gesucht ist. Wichtig: f ist ein
%  Function-handle, welches nur von einer Variablen abhaengt!. D.h.
%  f = @(x)
%  df: Ableitung der Funktion f.  Wichtig: df ist ein Function-Handle, d.h.
%  df = @(x)
%  x0: Startwert fuer das Newton-Verfahren
%  tol: Parameter fuer die Abbruchbedingung. Das Verfahren wird
%  durchgefuehrt, bis ||f(x)||<tol gilt ...
%  maxiter: ... oder maxiter Iterationsschritte durchgefuehrt wurden.
%
% Ausgabe:
%  x: Nullstelle der Funktion f.

lambda_min = 1e-4;
x = x0;
k = 0;

while (norm(f(x)) > tol ) && (k < maxiter)
    fx = f(x);
    dfx = df(x);
    s = -dfx\fx;
    
    lambda = 4;
    xk = x + lambda*s;
    C = 1-lambda/4;
    while(norm(dfx\f(xk)) > C*norm(dfx\fx) ) && (lambda>lambda_min)
        lambda = lambda/2;
        xk = x+lambda*s;
        C = 1 - lambda/4;
    end

    if (lambda <= lambda_min)
        error('Minimales Lambda erreicht. Das gedaempfte Newton-Verfahren konvergiert nicht.')
    end
    x = x + lambda*s;

end

if (k == maxiter)
    disp('Maximale Anzahl Iterationsschritte nicht erreicht!')
    error('Keine Konvergenz mit dieser Anzahl an maximalen Iterationsschritten')
end

end 