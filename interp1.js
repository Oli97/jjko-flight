find = function(x,y)
{
 var n = x.length;
 var ind = [];
 for (var i = 0; x[i] <= y; i++)
 {
    ind[i] = i;
 }
 return ind;
}




sort = function(x)
{
 var n = x.length;
 var k = [];
 for (var l = 0; l < n; l++)
 {
    k[l]=l;
 }
 for (var i = n-1; i>=0; i--)
 {
    for(var j = 1; j<=i; j++)
    {
      if(x[j-1]>x[j])
      {
	var temp = x[j-1];
        x[j-1] = x[j];
        x[j] = temp;

	var temp2 = j-1;
	k[j-1] = j;
	k[j] = temp2;
      }
    }
 }
 var y = [];
 y[0] = x;
 y[1] = k;
 return y;

}



interp1 = function(x,y,xi)
{
//   INTERP1Q Quick 1-D linear interpolation..
//   F=INTERP1Q(X,Y,XI) returns the value of the 1-D function Y at the
//   points XI using linear interpolation. Length(F)=length(XI).
//   The vector X to specifies the coordinates of the underlying interval.
//   
//   If Y is a matrix, then the interpolation is performed for each column
//   of Y in which case F is length(XI)-by-size(Y,2).
//
//   NaN's are returned for values of XI outside the coordinates in X.
//
//   INTERP1Q is quicker than INTERP1 on non-uniformly spaced data because
//   it does no input checking.   INTERP1(...,'*linear') is even quicker
//   for uniformly spaced data.  For INTERP1Q to work properly:
//     X must be a monotonically increasing column vector.
//     Y must be a column vector or matrix with length(x) rows.
//
//   See also INTERP1.

//   Copyright (c) 1984-98 by The MathWorks, Inc.
//   $Revision: 1.7 $  $Date: 1997/11/21 23:40:41 $

 var siz = xi.length;
 if (xi.length~=1)
 {

    alert("xi ist ein Vektor statt ein Skalar bei interp1. Bitte implementieren Sie den vektoriellen Fall!");
//     var y = sort(xi);
//     var xxi = y[0];
//     var k = y[1];
//     var y1 = sort([x;xxi]);
//     var dum = y1[0];
//     var j = y[1]; // Ab hier Copy-Paste
//     r(j)=1:j.length;
//     r=r(x.length+1:end)-(1:xxi.length);
//     r(k)=r;
//     r(xi==x(end))=length(x)-1;
//     ind=find((r>0) & (r<length(x)));
//     ind = ind(:);
//     yi=repmat(NaN,length(xxi),size(y,2));
//     rind = r(ind);
//     u = (xi(ind)-x(rind))./(x(rind+1)-x(rind));
//     yi(ind,:)=y(rind,:)+(y(rind+1,:)-y(rind,:)).*u(:,ones(1,size(y,2)));

 }
 else
 {
    // Special scalar xi case
    var r1 = find(x,xi);
    var n1 = r1.length;
    var r = r1[n1];
    var n = x.length;
    if (xi == x[n])
    {
      r = x.length-1;
    }
    /*if (isempty(r)) // Wenn r leer ist, gibt es ein Problem!
    {
      yi = NaN;  
      return
    }*/
    if ((r>0) && (r<x.length))
    {
      var u = (xi-x[r])/(x[r+1]-x[r]);
      yi=y[r,:]+(y[r+1,:]-y[r,:]).*u[:,ones[1,size(y,2)]];				// Fehlt noch!!

    }
    else
    {
      yi = NaN;
    }

 }
/* if ((min(size(yi))==1) & (prod(siz)>1))   // prod(siz) ist bei uns 1!! Sonst Fehler!!
 {
    yi = reshape(yi,siz);
 }*/
}


