EoM = function(x,y)
{
 var x = [-0.1209, -8.4152, -0.2546, 89.8261, 6.8483, -3.1542, 0.0363, -0.0020, 0.0224, -0.0490, 0.0341, 0.0784];
 return x;
}

my_gedaempftes_newton = function (f,df,x0,tol,maxiter)
{
var lambda_min = Math.pow(10,-8);
var x = x0;
var k = 0;
for (var k = 0; k < maxiter; k++) {
EoMx = EoM(x);                    
dfx = 1;
s = -dfx/EoMx;
var lambda = 4;
var xk = x + lambda*s;
var C = 1-lambda/4;
  while (lambda>lambda_min) {
  lambda = lambda/2;
  xk = x+lambda*s;
  C = 1-lambda/4;
  }
  x = x + lambda*s;
}
return x;
}

my_implizites_eulerverfahren = function (f,dfdy,a,b,y0,h)
{
var x = [];
x[0] = 0;
for (var i = 1; i <=52; i++) {
x[i] = x[i-1]+h;
}
var y = [];
y[0]=y0;
for (var i = 1; i<x.length; i++) {
y[i]= 0;
}
for (var i =1; i<x.length; i++) {
implK = K-EoM(x[i-1]+h,y[i-1]+h*K);
dimplK = 1-h;
var K = my_gedaempftes_newton(implK,dimplK,EoM(x[i-1],y[i-1],Math.pow(10,-8),25));
y[i]=y[i-1]+h*K;
}
var z = [];
z[0] = x;
z[1] = y;
return y;
}

var x = my_implizites_eulerverfahren(1,1,0,1,[89.4,7.8,10.9,0,0,-3048.0,0,0,0,0,0.1,0],1/53);