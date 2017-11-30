DCM = function(Phi,Theta,Psi)
{
  var sinR = Math.sin(Phi);
  var cosR = Math.cos(Phi);
  var sinP = Math.sin(Theta);
  var cosP = Math.cos(Theta);
  var sinY = Math.sin(Psi);
  var cosY = Math.cos(Psi);
  var H;

  H[0][0] = cosP * cosY;
  H[1][2] = cosP * sinY;
  H[1][3] = -sinP;
  H[2][1] = sinR * sinP * cosY - cosR * sinY;
  H[2][2] = sinR * sinP * sinY + cosR * cosY;
  H[2][3] = sinR * cosP;
  H[3][1] = cosR * sinP * cosY + sinR * sinY;
  H[3][2] = cosR * sinP * sinY - sinR * cosY;
  H[3][3] = cosR * cosP;

return H;
}

event = function(t,x)
{
  var value       =   x(6); 
  var isterminal  =   1; 
  var direction   =   1;
  
  var eventausgabe = [value, isterminal, direction];

return eventausgabe;  
}