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


/*	 RMQ computes the rotation matrix as a function of quaternions.
The rotation matrices transform vectors from the earth-relative frame of
reference to the body-axis frame.
Quaternion: q1, x Component of quaternion q2, y Component of quaternion
q3, z Component of quaternion q4, cos(Euler) Component of quaternion */
RMQ = function(q1,q1,q3,q4)
{
  var H;
	H[1][1] = q1^2 - q2^2 - q3^2 + q4^2;
	H[1][2] = 2*(q1*q2 + q3*q4);
	H[1][3] = 2*(q1*q3 - q2*q4);
	H[2][1] = 2*(q1*q2 - q3*q4);
	H[2][2] = -q1^2 + q2^2 - q3^2 + q4^2;
	H[2][3] = 2*(q2*q3 + q1*q4);
	H[3][1] = 2*(q1*q3 + q2*q4);
	H[3][2] = 2*(q2*q3 - q1*q4);
	H[3][3] = -q1^2 - q2^2 + q3^2 + q4^2;
  return H;
}
