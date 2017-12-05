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

/* RMQ computes the rotation matrix as a function of quaternions.
The rotation matrices transform vectors from the earth-relative frame of
reference to the body-axis frame.
Quaternion: q1, x Component of quaternion q2, y Component of quaternion
q3, z Component of quaternion q4, cos(Euler) Component of quaternion */
RMQ = function(q1,q1,q3,q4)
{
  var H;
  H[1][1] = Math.pow(q1,2) - Math.pow(q2,2) - Math.pow(q3,2) + Math.pow(q4,2);
  H[1][2] = 2*(q1*q2 + q3*q4);
  H[1][3] = 2*(q1*q3 - q2*q4);
  H[2][1] = 2*(q1*q2 - q3*q4);
  H[2][2] = -Math.pow(q1,2) + Math.pow(q2,2) - Math.pow(q3,2) + Math.pow(q4,2);
  H[2][3] = 2*(q2*q3 + q1*q4);
  H[3][1] = 2*(q1*q3 + q2*q4);
  H[3][2] = 2*(q2*q3 - q1*q4);
  H[3][3] = -Math.pow(q1,2) - Math.pow(q2,2) + Math.pow(q3,2) + Math.pow(q4,2);
  return H;
}

AeroModelMach = function(x,u,Mach,alphar,betar,V)
{
  //Typical Mass and Inertial Properties
  var m = 4536; // Mass, kg
  var Ixx = 35926.5; // Roll Moment of Inertia, kg-m^2
  var Iyy = 33940.7;// Pitch Moment of Inertia, kg-m^2
  var Izz = 67085.5;// Yaw Moment of Inertia, kg-m^2
  var Ixz = 3418.17;// Nose-high(low) Product of Inertia, kg-m^2
  //Geometric Properties
  var cBar=2.14;// Mean Aerodynamic Chord, m
  var b=10.4;// Wing Span, m
  var S=21.5;// Reference Area, m^2
  var ARw=5.02;// Wing Aspect Ratio
  var taperw=0.507;// Wing Taper Ratio
  var sweepw=13 * .01745329;// Wing 1/4-chord sweep angle, rad
  var ARh=4;// Horizontal Tail Aspect Ratio
  var sweeph=25 * .01745329;// Horiz Tail 1/4-chord sweep angle, rad
  var ARv=0.64;// Vertical Tail Aspect Ratio
  var sweepv=40 * .01745329;// Vert Tail 1/4-chord sweep angle, rad
  var lvt=4.72;// Vert Tail Length, m

  //Thrust Properties
  var StaticThrust = 26243.2;// Static Thrust @ Sea Level, N

  //Current Thrust
  var atmos = Atmos(-x[6]);
  var airDens = atmos[0];
  var airPres = atmos[1];
  var temp = atmos[2];
  var soundSpeed = atmos[3];
  var Thrust=u[4] * StaticThrust * Math.pow((airDens / 1.225),0.7)* (1 - Math.exp((-x[6] - 17000) / 2000));
  // Thrust at Altitude, N

  //Current Mach Effects, normalized to Test Condition B (Mach = 0.1734)
  var PrFac=1 / (Math.sqrt(1 - Math.pow(Mach,2) * 1.015);
  // Prandtl Factor
  var WingMach=1 / ((1 + Math.sqrt(1 + (Math.pow((ARw/(2 * Math.cos(sweepw))),2))* (1 - Math.pow(Mach,2) * Math.cos(sweepw)))) * 0.268249);
  // Modified Helmbold equation
  var HorizTailMach=1 / ((1 + Math.sqrt(1 + (Math.pow((ARh/(2 * Math.cos(sweeph))),2))* (1 - Math.pow(Mach,2) * Math.cos(sweeph)))) * 0.294539);
  // Modified Helmbold equation
  var VertTailMach=1 / ((1 + Math.sqrt(1 + (Math.pow((ARv/(2 * Math.cos(sweepv))),2))* (1 - Math.pow(Mach,2) * Math.cos(sweepv)))) * 0.480338);
  // Modified Helmbold equation

  //Current Longitudinal Characteristics
  //====================================

  //Lift Coefficient
  var CLo=0.1095;// Zero-AoA Lift Coefficient (B)
  if (GEAR >= 1)
  {
    CLo=CLo - 0.0192;// Gear-down correction
  }
  if (u[6] >= 0.65)
  {
    CLo=CLo + 0.5182;// 38 deg-flap correction
  }
  if (SPOIL >= 1)
  {
    CLo=CLo - 0.1897;// 42 deg-Symmetric Spoiler correction
  }

  var CLar=5.6575;// Lift Slope (B), per rad
  if (u[6] >= 0.65)
  {
    CLar=CLar - 0.0947;
  }

  var CLqr=4.231 * cBar / (2 * V);
  // Pitch-Rate Effect, per rad/s

  CLdSr=1.08;// Stabilator Effect, per rad
  if (u[6] >= 0.65)
  {
    CLdSr=CLdSr - 0.4802;// 38�-flap correction
  }

  var CLdEr=0.5774;// Elevator Effect, per rad
  if ( u[6] >= 0.65)
  {
    CLdEr=CLdEr - 0.2665;// 38 deg-flap correction
  }

  var CL=CLo + (CLar*alphar + CLqr*x[8] + CLdSr*u[7] + CLdEr*u[1])* WingMach;
                                    // Total Lift Coefficient, w/Mach Correction

  //Drag Coefficient
  var CDo=0.0255;// Parasite Drag Coefficient (B)
  if ( GEAR >= 1)
  {
    CDo=CDo + 0.0191;// Gear-down correction
  }
  if ( u[6] >= 0.65)
  {
    CDo=CDo + 0.0836;// 38 deg-flap correction
  }
  if ( SPOIL >= 1)
  {
    CDo=CDo + 0.0258;// 42 deg-Symmetric Spoiler correction
  }

  var epsilon=0.0718;// Induced Drag Factor
  if ( u[6] >= 0.65)
  {
    epsilon=0.079;// 38 deg-flap correction
  }

  var CD=CDo * PrFac + epsilon * Math.pow(CL,2);
                                    // Total Drag Coefficient, w/Mach Correction

  //Pitching Moment Coefficient
  var Cmo=0;// Zero-AoA Moment Coefficient (B)
  if ( GEAR >= 1)
  {
    Cmo=Cmo + 0.0255;// Gear-down correction
  }
  if ( u[6] >= 0.65)
  {
    Cmo=Cmo - 0.058;// 38 deg-flap correction
  }
  if ( SPOIL >= 1)
  {
    Cmo=Cmo - 0.0154;// 42 deg-Symmetric Spoiler correction
  }
  
  var Cmar=-1.231;// Static Stability (B), per rad
  if ( u[6] >= 0.65)
  {
    Cmar=Cmar + 0.0138;
  }

  var Cmqr    = -18.8 * cBar / (2 * V);
                                    // Pitch-Rate + Alpha-Rate Effect, per rad/s

  var CmdSr=-2.291;// Stabilator Effect, per rad
  if ( u[6] >= 0.65)
  {
    CmdSr=CmdSr + 0.121;// 38 deg-flap correction
  }

  var CmdEr=-1.398;// Elevator Effect, per rad
  if ( u[6] >= 0.65)
  {
    CmdEr=CmdEr + 0.149;// 38 deg-flap correction
  }
  
  var Cm=Cmo + (Cmar*alphar + Cmqr*x[8] + CmdSr*u[7] + CmdEr*u[1])* HorizTailMach;
                                    // Total Pitching Moment Coefficient, w/Mach Correction

  //Current Lateral-Directional Characteristics
  //===========================================

  //Side-Force Coefficient
  var CYBr=-0.7162;// Side-Force Slope (B), per rad
  if ( u[6] >= 0.65)
  {
    CYBr=CYBr + 0.0826;
  }

  var CYdAr=-0.00699;// Aileron Effect, per rad

  var CYdRr=0.1574;// Rudder Effect, per rad
  if ( u[6] >= 0.65)
  {
    CYdRr=CYdRr - 0.0093;// 38 deg-flap correction
  }

  var CYdASr=0.0264;// Asymmetric Spoiler Effect, per rad
  if ( u[6] >= 0.65)
  {
    CYdASr=CYdASr + 0.0766;
  // 38 deg-flap correction
  }

  var CY=(CYBr*betar + CYdRr*u[3]) * VertTailMach+ (CYdAr*u[2] + CYdASr*u[5]) * WingMach;
                                    // Total Side-Force Coefficient, w/Mach Correction

  //Yawing Moment Coefficient
  var CnBr=0.1194;// Directional Stability (B), per rad
  if ( u[6] >= 0.65)
  {
    CnBr=CnBr - 0.0092;
  }

  var Cnpr=CL * (1 + 3 * taperw)/(12 * (1 + taperw)) * (b / (2 * V));
  // Roll-Rate Effect, per rad/s

  var Cnrr=(-2 * (lvt / b) * CnBr * VertTailMach - 0.1 * Math.pow(CL,2))* (b / (2 * V));
  // Yaw-Rate Effect, per rad/s

  var CndAr=0;                      // Aileron Effect, per rad
  if ( u[6] >= 0.65)
  {
    CndAr=CndAr + 0.0028;
  }

  var CndRr=-0.0713;                // Rudder Effect, per rad
  if ( u[6] >= 0.65)
  {
    CndRr=CndRr - 0.0185;     // 38 deg-flap correction
  }

  var CndASr=-0.0088;                // Asymmetric Spoiler Effect, per rad
  if ( u[6] >= 0.65)
  {
    CndASr=CndASr - 0.0106;
                                        // 38 deg-flap correction
  }

  var Cn=(CnBr*betar + CndRr*u[3]) * VertTailMach+ Cnrr * x[9] + Cnpr * x[7]+ (CndAr*u[2] + CndASr*u[5]) * WingMach;
                                        // Total Yawing-Moment Coefficient, w/Mach Correction

  //Rolling Moment Coefficient
  var ClBr=-0.0918;                // Dihedral Effect (B), per rad
  if ( u[6] >= 0.65)
  {
    ClBr=ClBr - 0.0092;
  }

  var Clpr=-CLar * (1 + 3 * taperw)/(12 * (1 + taperw))* (b / (2 * V));
                                        // Roll-Rate Effect, per rad/s

  var Clrr=(CL * (1 + 3 * taperw)/(12 * (1 + taperw))* Math.pow((Mach * Math.cos(sweepw)),2) - 2) / Math.pow((Mach * Math.cos(sweepw)),2) - 1))* (b / (2 * V));
                                        // Yaw-Rate Effect, per rad/s

  var CldAr=0.1537;                 // Aileron Effect, per rad
  if ( u[6] >= 0.65)
  {
    CldAr=CldAr + 0.01178;
  }

  var CldRr=0.01208;// Rudder Effect, per rad
  if ( u[6] >= 0.65)
  {
    CldRr=CldRr + 0.01115;// 38 deg-flap correction
  }

  var CldASr=-0.01496;// Asymmetric Spoiler Effect, per rad
  if ( u[6] >= 0.65)
  {
    CldASr=CldASr - 0.02376;
  // 38 deg-flap correction
  }

  var Cl=(ClBr*betar + CldRr*u[3]) * VertTailMach+ Clrr * x[9] + Clpr * x[7]+ (CldAr*u[2] + CldASr*u[5]) * WingMach;
                                        // Total Rolling-Moment Coefficient, w/Mach Correction
}