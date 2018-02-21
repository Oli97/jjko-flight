numjac = function(F,t,y,Fty,thresh,vectorized,S,g,varargin){
// Initialize.
var eps = Math.pow(2.2204,-16);
var br = Math.pow(eps,0.875);
var bl = Math.pow(eps,0.75);
var bu = Math.pow(eps,0.25);
var facmin = Math.pow(eps,0.78);
var facmax = 0.1;
var ny = y.length;
var nF = Fty.length;
var fac = new Array();
if(fac.length==0){
  for(var i=0; i<ny; i++) {
    fac[i] = Math.sqrt(eps);
  }
}
// Select an increment del for a difference approximation to
// column j of dFdy.  The vector fac accounts for experience
// gained in previous calls to numjac.
var yscale = new Array();
for(var i=0;i<ny; i++) {
    yscale[i] = Math.max(Math.abs(y[i]),thresh[i]);
}
var del = new Array();

var j=[];
for(var i=0;i<ny; i++) {
  del[i] = (y[i] + fac[i] * yscale[i]) - y[i];
  if(del[i]==0){
    if (fac[i] < facmax) {
      fac[i] = Math.min(100*fac[i],facmax);
      del[i] = (y[i]+fac[i]*yscale[i])-y[i];
      if(del[i]==0){}
      else {
        del[i] = thresh[i];
      }
    }
  }
}

var s=[];
if (nF == ny) {
for(var i=0;i<ny;i++) {
  if(Math.sign(Fty)>=0){
    s[i] = 1;
  }else{
    s[i] = 0;
  }
  del[i] = ((s[i]-1+s[i])*Math.abs(del[i])
}
}

// Form a difference approximation to all columns of dFdy.
var S = new Array ();
var g = new Array ();
var ydel=[];
for(var i=0;i<ny;i++){
  for(var h=0;h<ny;h++){
    if(i==h){
      ydel[i][h] = y[i] + del[i];
    }else{
      ydel[i][h] = y[i];
    }
  }
}
var Fdel = [];
if (vectorized==1){
  Fdel = feval(F,t,ydel,varargin{:});//feval muss noch umgeschrieben werden
} //see other code
else{
  for(var i=0;i<nF;i++){
    for(var h=0;h<ny;h++){
      Fdel[i][h]=0;
    }
  }
  for (var i = 0;i<ny; i++) {
        Fdel(:,i) = feval(F,t,ydel(:,i),varargin{:});
  }
}
  var nfevals = ny;                         // stats (at least one per loop)
  var Fdiff = Fdel - Math.ones(1,ny)*Fty;
  var dFdy = Fdiff * Math.eye(1,ny)*Math.pow(del,-1) ;
  var Difmax = new Array();
  for(var i=0;i<ny; i++) {
      Difmax[i] = Math.max(Math.abs(Fdiff[i])
  }
var Rowmax = new Array();
var Rowmax = Difmax;//rowmax beinhaltet die indizes dessen eintrags der am größten ist

// If Fdel is a column vector, then index is a scalar, so indexing is okay.
var absFdelRm = new Array();
for(var i=0;i<ny; i++) {
    absFdelRm[i] = Math.abs(Fdel[i]*nF+Rowmax[i]);
}

// Adjust fac for next call to numjac.
var absFty=[];
for(var i=0;i<ny;i++) {
  absFty[i]= Math.abs(Fty[i]);
}
for(var i=0;i<ny;i++) {
if (Rowmax[i]<0) {
Rowmax[i]=-Rowmax[i];
}
}
absFtyRm = transpose(Math.absFty(Rowmax));//transpose muss noch geschrieben werden
var ydel = y;
var Fscale=[];
for(var i=0;i<ny;i++) {
  Fscale = Math.max(absFdelRm[i],absFtyRm[i]);
}
var j=[];
for(var i=0;i<ny;i++) {
  if((absFdelRm[i]!=0&&absFtyRm!=0)||Difmax==0){
    j = [j 1];
  }else{
    j = [j 0];
  }
}
  // If the difference in f values is so small that the column might be just
  // roundoff error, try a bigger increment.
var k1 = [];
for(var i=0;i<ny;i++) {
  if(Difmax <= br*Fscale){
    k1 = [k1 1];
  }else{
    k1 = [k1 0];
  }
}          // Difmax and Fscale might be zero
var k = new Array();
for(var i;i<=k1.length;i++){
  if(j[i]!=0&&k1[i]!=0){
    k = [k i]
  }
}
  for k = finde j und k1
    tmpfac = Math.min(Math.sqrt(fac[k]),facmax);
    del = (y[k] + tmpfac*yscale[k]) - y[k];
    if (tmpfac ~= fac[k] && del ~= 0) {
      if (nF == ny){
        if (Fty[k] >= 0) {                 // keep del pointing into region
          del = Math.abs(del);
          }
        else {
          del = -Math.abs(del);
      }
      }
      }
      ydel[k] = y[k] + del;
      var fdel = feval(F,t,ydel,varargin{:});
      nfevals = nfevals + 1;            // stats
      ydel[k] = y[k];
      var fdiff,temp,rowmax = [];
      var difmax=0;
      for(var i=0;i<ny; i++) {
          fdiff[i] = fdel[i] - Fty[i];
          tmp[i] = fdiff[i]*Math.pow(del,-1);
          if(difmax < Math.max(Math.abs(fdiff[i])){
            difmax = Math.max(Math.abs(fdiff[i]);
            rowmax = i;
          }
      }
      if (tmpfac * Math.norm(tmp,inf) >= Math.norm(dFdy,inf);

        // The new difference is more significant, so
        // use the column computed with this increment.
          dFdy = tmp;

        // Adjust fac for the next call to numjac.
        fscale = max(abs(fdel(rowmax)),absFty(rowmax));
        if (difmax <= bl*fscale) {

          // The difference is small, so increase the increment.
          fac[k] = Math.min(10*tmpfac, facmax);
          }
        else if (difmax > bu*fscale) {

          // The difference is large, so reduce the increment.
          fac[k] = Math.max(0.1*tmpfac, facmin);
          }
        else {
          fac[k] = tmpfac;
          }

  // If the difference is small, increase the increment.
  for(var i;i<=length(k1);i++){
    if(j[i]!=0&&k1[i]==0&&(Difmax[i] <= bl*Fscale[i])){
      k = [k i]
    }
  }
  if (k === '') {
    fac[k] = Math.min(10*fac[k], facmax);
  }

  // If the difference is large, reduce the increment.
  for(var i;i<=length(j);i++){
    if(j[i]!=0&&(Difmax[i] > bl*Fscale[i])){
      k = [k i]
    }
  }
  if (k === '') {
    fac[k] = Math.max(0.1*fac[k], facmin);
}

z=[dFdy,fac,g,nfevals];
return z;

}
