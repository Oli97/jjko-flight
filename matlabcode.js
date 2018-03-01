/*globale Variablen*/ // Was soll das??
var Trimhist, x, u, V;
/*// Supporting Calculations for Geometric, Inertial, and Aerodynamic
// Properties of BizJet B
// June 12, 2015
// Copyright 2006-2015 by ROBERT F. STENGEL. All rights reserved.

 clear
 disp('==============================================')
 disp('Geometric and Inertial Properties for BizJet B')
 disp('==============================================')
 Date = date
 */
// Comparable Bizjet Weight Distribution, lb, based on empty weight less engine weight.
// (from Stanford AA241 Notes, Ilan Kroo)
var WingSys = 1020;
var TailSys = 288;
var BodySys = 930;
var GearSys = 425;
var NacelleSys = 241;
var PropSys = 340;
var ControlSys = 196;
var InstrSys = 76;
var HydrPneuSys = 94;
var ElecSys = 361;
var AvionSys = 321;
var FurnEquipSys = 794;
var ACSys = 188;
var AntiIceSys = 101;
var LoadHandSys = 2;
var EmptyStruc = (WingSys+TailSys+BodySys+GearSys+PropSys+ControlSys+InstrSys+HydrPneuSys+ElecSys+AvionSys+FurnEquipSys+ACSys+AntiIceSys+LoadHandSys); // Less nacelles and engines
var EngWgt = 1002;
var EmptyWgt = (EmptyStruc+EngWgt); // Less nacelles
var EmptyStrucWgt = EmptyStruc/EmptyWgt;
var WingRatio = (WingSys+GearSys+HydrPneuSys+AntiIceSys)/EmptyStruc;
var HTRatio = 0.75*(TailSys)/EmptyStruc;
var VTRatio = 0.25*(TailSys)/EmptyStruc;
var FusRatio = (BodySys+PropSys+ControlSys+InstrSys+ElecSys+AvionSys+FurnEquipSys+ACSys+LoadHandSys)/EmptyStruc;
var Total = (WingRatio+HTRatio+VTRatio+FusRatio);
var NacelleRatio = (NacelleSys)/EngWgt; // Related to engine weight rather than empty structure

// BizJet B Geometric Properties
// x measurements from nose along centerline, negative aft
// y & z measurements from centerline, positive right and down

var S = 19.51 // Reference Area, m^2
var taperw = 0.5333 // Wing Taper Ratio
var cBar = 1.56 // Mean Aerodynamic Chord, m
var sweep = 11*0.01745329 // Wing L.E. sweep angle, rad
var xcp = -5.7473 // Wing center of pressure, m
var GamWing = 3*0.01745329 // Dihedral angle of the wing, rad

// BizJet B Mass and Inertial Properties
var m = 3000 // Total mass for simulation (USER-specified), kg
var mEmpty = 2522 // Gross empty mass, kg
var mEng = (1+NacelleRatio)*240 // Mass of engines + nacelles, kg
var mStruc = mEmpty - mEng // Empty structural mass (less engines + nacelles), kg
var mWing = WingRatio*mStruc // Wing mass, kg
var mHT = HTRatio*mStruc // Horizontal tail mass, kg
var mVT = VTRatio*mStruc // Vertical tail mass, kg
var mFus = FusRatio*mStruc // Empty fuselage mass, kg
var mPay = 0.5*(m - mEmpty) // Payload mass, kg
var mFuel = 0.5*(m - mEmpty) // Fuel mass, kg
var xcm = xcp - 0.45*cBar// Center of mass from nose (USER-specified), m
var lWing = xcm - xcp // Horizontal distance between c.m and wing c.p., m
var zWing = -0.557 // Vertical distance between c.m and wing c.p., m
var b = 13.16 // Wing Span, m
var lenFus = 10.72 // Fuselage length, m
var xcpFus = -0.25*lenFus // Linear-regime fuselage center of pressure, m
var xcpFusN = -0.5*lenFus // Newtonian-regime fuselage center of pressure, m
var lFus = xcm - xcpFus; // Linear fuselage lift cp offset,m
var lFusN = xcm - xcpFusN; // Newtonian fuselage lift cp offset, m
var dFus = 1.555 // Fuselage diameter, m
var Sfus = (Math.PI/4)*lenFus*dFus // Plan or side area of fuselage, m
var Sbase = (Math.PI/4)*Math.pow(dFus,2) // Fuselage cross-sectional area, m^2
var bHT = 5.3; // Horizontal tail span, m
var cHT = 1.1; // Mean horizontal tail chord, m
var swpHT = 38*0.0174533; // Horizontal tail sweep, rad
var SHT = bHT*cHT // Horizontal tail area, m^2
var xHT = -11.3426; // Linear xcp of horizontal tail
var lHT = xcm - xHT; // Horizontal tail length, m
var zHT = 1.5; // zcp of horizontal tail, m
var xVT = -10.044; // Linear xcp of vertical tail, m
var lVT = xcm - xVT; // Vertical tail length, m
var bVT = 2.409; // Vertical tail span, m
var cVT = 1.88; // Mean vertical tail chord, m
var swpVT = 50*0.0174533; // Vertical tail sweep, rad
var SVT = bVT*cVT; // Vertical tail area, m^2
var zVT = 1.5; // zcp of vertical tail, m
var xEng = -7.735; // xcm of engine, m
var lEng = xcm - xEng; // Engine length, m
var yEng = 1.1325; // ycm of engine, m
var zEng = 0.4038; // zcm of engine, m
var xNac = -7.7252; // xcp of engine, m
var lNac = xcm - xNac; // Nacelle length, m
var bNac = 2.5; // Nacelle span, m
var cNac = 1.82; // Nacellechord, m
var dNac = 0.73; // Nacelle diameter, m
var SbaseNac = 0.25*Math.PI*Math.pow(dNac,2); // Nacelle base area, m^2
var Snac = bNac*cNac; // Nacelle plan area, m^2
var xVent = -9.94; // xcp of ventral fin, m
var lVent = xcm - xVent; // Ventral fin length, m
var zVent = 0; // zcp of ventral fin, m
var bVent = 1; // Ventral fin span, m
var cVent = 0.85; // Ventral fin chord, m
var Svent = bVent*cVent; // Ventral fin area, m^2
var swpVent = 60*0.0174533; // Ventral sweep angle, rad
var Splan = S + Sfus + Snac + SHT + Svent // Plan area of airplane
var Swet = 2*(Splan + Sfus + SVT) // Wetted area of airplane
var ARwing = (b*b) / S // Wing aspect ratio
var ARHT = (bHT*bHT) / SHT // Horizontal tail aspect ratio
var ARnac = (bNac*bNac) / Snac // Engine nacelle aspect ratio
var ARvent = (bVent*bVent)/ Svent // Ventral fin aspect ratio
var ARVT = (bVT*bVT) / SVT // Vertical tail aspect ratio

// Moments and Product of Inertia
var Ixx = (1/12)*((mWing+mFuel)*b*b + mHT*Math.pow(bHT,2) + mVT*Math.pow(bVT,2)) + (0.25*(mFus+mPay)*Math.pow(dFus,2) + mEng*Math.pow(yEng,2) + mVT*Math.pow(zVT,2));
var Iyy = (1/12)*((mFus+mPay)*Math.pow(lenFus,2) + (mWing+mFuel)*Math.pow(cBar,2) + mVT*Math.pow(cHT,2)) + (mEng*Math.pow(lEng,2) + mHT*Math.pow(lHT,2) + mVT*Math.pow(lHT,2));
var Izz = (1/12)*((mFus+mPay)*Math.pow(lenFus,2) + (mWing+mFuel)*Math.pow(b,2) + mHT*Math.pow(bHT,2)) + (mEng*Math.pow(lEng,2) + mHT*Math.pow(lHT,2) + mVT*Math.pow(lVT,2));
var Ixz = mHT*lHT*zHT + mVT*lVT*zVT + mEng*lEng*zEng
var dEmax = 20 * 0.01745329 // Maximum Elevator Deflection is ±20 deg
var dAmax = 35 * 0.01745329 // Maximum Aileron Deflection is ±35 deg
var dRmax = 35 * 0.01745329 // Maximum Rudder Deflection is ±35 deg

// BizJet B Aero Properties
var AlphaTable = [-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90];
var Points = AlphaTable.length;
var AlphaRad = 0.0174533*AlphaTable;
var SinAlpha = Math.sin(AlphaRad);
var CosAlpha = Math.cos(AlphaRad);

// Newtonian Coefficients
var CN = 2*(Splan/S)*elemmult(SinAlpha,SinAlpha);
CN = elemmult(CN,sign(AlphaRad));
var CDNewt = 2*(Splan/S)*Math.abs(elemmult(elemmult(SinAlpha,SinAlpha),SinAlpha));
var CLNewt = elemmult(CN,CosAlpha);
var cpNewt = (S*lWing + Sfus*lFusN + SHT*lHT + Svent*lVent + Snac*lNac)/ (S + Sfus + SHT + Svent + Snac);

// Longitudinal Aerodynamics
// =========================
// Lift
var CLaWing = Math.PI*ARwing / (1 + Math.sqrt(1 + Math.pow((0.5*ARwing/Math.cos(sweep)),2)));
var deda = 0; // T-tail
var CLaHT = (1 - deda)*(Math.PI*ARHT / (1 + Math.sqrt(1 + Math.pow((0.5*ARHT/Math.cos(swpHT)),2))))*SHT / S;
var CLaFus = 2*Sbase / S;
var CLaNac = 2*SbaseNac / S;
var CLaVent = (Math.PI*ARvent / (1 + Math.sqrt(1 + Math.pow((0.5*ARvent/Math.cos(swpVent)),2))))*Svent / S;
var CLaTot = CLaWing + CLaHT + CLaFus + CLaNac + CLaVent;
// Assume CL is linear to Alpha = 10 deg = 0.1745 rad
var CL10 = CLaTot*0.174533;
// Assume CL is symmetrically quadratic about CLmax
var CLmax = 1.35;
var delAlph = 4*0.0174533; // Stall occurs at 14 deg
var CLstatic = zeros(0,39);

for (var i = 1; i <= 11; i++)
{
  CLstatic[i] = CLaTot*AlphaRad[i];
}

var kStall = (CLmax - CL10)/Math.pow(delAlph,2);

for (var i = 12; i <= 21; i++)
{
  CLstatic[i] = CLmax - kStall*elemmult((AlphaRad[15] - AlphaRad[i]),(AlphaRad[15] - AlphaRad[i]));
}

CLstatic[22] = 0.73;
CLstatic[23] = 0.73;
CLstatic[24] = 0.74;
CLstatic[25] = 0.76;
CLstatic[26] = 0.78;

for (var i = 27; i <= 39; i++)
{
  CLstatic[i] = CN[i]*CosAlpha[i];
}

var CLTable = CLstatic;

// Drag
var Re = 1.225*100*lenFus / 1.725e-5; // Reynolds number at 100 m/s
var Cf = 0.46*Math.pow((Math.log10(Re)),-2.58); // Flat plate friction coefficient
var CDf = Cf*Swet / S;
var CDbase = 0.12*Sbase / S; // Base pressure drag
var CDo = CDf + CDbase;
var OEF = 1.78*(1 - 0.045*Math.pow(ARwing,0.68)) - 0.64; // Oswald Efficiency factor, Raymer (straight wing)
var CDstatic = zeros(0,39);

for (var i = 1; i <= 10; i++)
{
  CDstatic[i] = CDo + elemmult(CLstatic[i],CLstatic[i]) / (OEF*Math.PI*ARwing);
}

CDstatic[11] = CDstatic[10]*1.15;
CDstatic[12] = CDstatic[11]*1.15;
CDstatic[13] = CDstatic[12]*1.15;
CDstatic[14] = CDstatic[13]*1.15;
CDstatic[15] = CDstatic[14]*1.1;
CDstatic[16] = CDstatic[15]*1.1;
CDstatic[17] = CDstatic[16]*1.1;
CDstatic[18] = CDstatic[17]*1.1;
CDstatic[19] = CDstatic[18]*1.1;
CDstatic[20] = CDstatic[19]*1.1;
CDstatic[21] = CDstatic[20]*1.1;
CDstatic[22] = CDstatic[21]*1.1;
CDstatic[23] = CDstatic[22]*1.1;
CDstatic[24] = CDstatic[23]*1.1;
CDstatic[25] = CDstatic[24]*1.1;
CDstatic[26] = CDstatic[25]*1.1;

for (var i = 26; i <= 39; i++)
{
  CDstatic[i] = 2*(Splan/S)*Math.abs(SinAlpha[i]*SinAlpha[i]*SinAlpha[i]);
}

var CDNewt = 2*(Splan/S)*Math.abs(elemmult(elemmult(SinAlpha,SinAlpha),SinAlpha));
var CDTable = CDstatic;

// Pitching Moment (c.m. @ wing c.p.)
var CmStatic = zeros(0,39);
var SM = (CLaWing*lWing + CLaFus*lFus + CLaHT*lHT + CLaNac*lNac + CLaVent*lVent) /(cBar*(CLaWing + CLaFus + CLaHT + CLaNac + CLaVent)); // Static Margin at Alpha = 0 deg

for (var i = 1; i <= 21; i++)
{
  CmStatic[i] = -(elemmult(CLstatic[i],Math.cos(AlphaRad[i])) + elemmult(CDstatic[i],Math.sin(AlphaRad[i])))*SM;
}

CmStatic[22] = CmStatic[21]*0.9;
CmStatic[23] = CmStatic[22]*0.9;
CmStatic[24] = CmStatic[23]*0.9;
CmStatic[25] = CmStatic[24]*0.9;
CmStatic[26] = 0.9*CmStatic[25] - 0.1*CN[26]*Math.sign(AlphaRad[26])*cpNewt/cBar;
CmStatic[27] = 0.4*CmStatic[26] - 0.6*CN[27]*Math.sign(AlphaRad[27])*cpNewt/cBar;
CmStatic[28] = 0.1*CmStatic[27] - 0.9*CN[28]*Math.sign(AlphaRad[28])*cpNewt/cBar;

var CN = 2*(Splan/S)*elemmult(SinAlpha,SinAlpha);

for (var i = 29; i <= 39; i++)
{
  CmStatic[i] = -CN[i]*cpNewt/cBar;
}

var CmTable = CmStatic;

// CmdETable & CLdEo, Elevator Effect
var tauDE = 0.32; // Geometric elevator chord/horizontal tail chord
var tauCO = 0.68; // Elevator Carryover effect
var Sigmoid = [];

for (var i = 0; i <= 4; i++)
{
  Sigmoid[i] = tauCO;
}

for (var i = 5; i <= 38; i++)
{
  Sigmoid[i] = tauDE + elemdiv((tauCO - tauDE), (1 + Math.exp(-15.*(AlphaRad[26] - AlphaRad[i]))));
}

var CmdETable = zeros(0,38);
var CLdETable = zeros(0,38);
var CLdEo = tauCO*CLaHT; // CLdE at Alpha = 0
var CmdEo = -tauCO*(lHT/cBar)*CLaHT; // CmdE at Alpha = 0

for (var i = 0; i <= 38; i++)
{
  CmdETable[i] = elemmult((CmdEo*Sigmoid),Math.cos(AlphaRad[i])); // Elevator effect on moment, per rad
}

for (var i = 0; i <= 38; i++)
{
  CLdETable[i] = elemmult((CLdEo*Sigmoid),Math.cos(AlphaRad[i])); // Elevator effect on lift, per rad
}

// Longitudinal Rotary & Unsteady Derivatives
var CLqHat = 2*CLaHT*lHT/cBar;
var CmqHat = -CLqHat*lHT/cBar;

// Lateral-Directional Aerodynamics
// ================================
// CYBetaTable, Side Force Sensitivity to Sideslip Angle
var EndPlate = 1.1; // End-plate effect of T-tail
var CYBetaVT = -EndPlate*(Math.PI*ARVT / (1 + Math.sqrt(1 + Math.pow((0.5*ARVT/Math.cos(swpVT)),2))))*SVT / S;
var CYBetaFus = -2*Sbase / S;
var CDoWing = 0.005;
var CYBetaWing = -CDoWing - (Math.pow(GamWing,2))*Math.PI*ARwing / (1 + Math.sqrt(1 + Math.pow(ARwing,2)));
var CYBetaVent = -0.4*CLaVent;
var CYBetao = CYBetaVT + CYBetaFus + CYBetaWing + CYBetaVent;
var CYBetaTable = CYBetao*Math.cos(AlphaRad);

// ClBetaTable, Roll Moment Sensitivity to Sideslip Angle
var ClBetaWing = -((1 + 2*taperw)/(6*(1 + taperw)))*(GamWing*CLaWing + (elemmult(CLTable,Math.tan(sweep))));
var ClBetaWF = 1.2*Math.sqrt(ARwing)*(2*zWing*dFus/Math.pow(b,2));
var ClBetaVT = -zVT*CYBetaVT/b;
var ClBetao = (ClBetaWing + ClBetaWF + ClBetaVT);
var ClBetaTable = elemmult(ClBetao,Math.cos(AlphaRad));

// CnBetaTable, Yaw Moment Sensitivity to Sideslip Angle
var CnBetaWing = 0.075*CLTable*GamWing;
var CnBetaFus = CLaFus*lFusN / b;
var CnBetaVT = -CYBetaVT*lVT / b;
var CnBetaVent = -CYBetaVent*lVent / b;
var CnBetaTable = elemmult((CnBetaWing + CnBetaFus + CnBetaVT),Math.cos(AlphaRad));

// CldATable, Roll Moment Sensitivity to Aileron Deflection
var tauDA = 0.25;
var kDA = 0.38;
var CldAo = tauDA*(CLaWing/(1 + taperw))*((1 - Math.pow(kDA,2))/3 - (1 - Math.pow(kDA,3))*(1 - taperw)/3);
var CldATable = CldAo*Math.cos(AlphaRad);
var CYdAo = 0; // Side force due to aileron, rad

// CndATable, Yaw Moment Sensitivity to Aileron Deflection; Cessna 510 has an Aileron-Rudder Interconnect; assume CndA = 0
var CndATable = zeros(0,39);

// CldRTable, Roll Moment Sensitivity to Rudder Deflection
var tauDR = 0.5; // Geometric Rudder chord/horizontal tail chord
var tauCOR = 0.8; // Rudder Carryover effect
var CldRo = tauCOR*zVT*CYBetaVT / b;
var CldRTable = CldRo*Math.cos(AlphaRad);

// CndRTable, Yaw Moment Sensitivity to Rudder Deflection
var CndRo = -tauCOR*CnBetaVT;
var CndRTable = CndRo*Math.cos(AlphaRad);

// Lateral-Directional Rotary & Unsteady Derivatives
var CYrHat = -2*CYBetaVT*lVT/b;
var ClpHato = -(CLaWing + CLaHT*(SHT/S) - CYBetaVT*(SVT/S))*((1 + 3*taperw)/(1 + taperw))/12;
var ClpHatTable = ClpHato*Math.cos(AlphaRad);
var ClrHato = -(CLaWing + CLaHT*(SHT/S) - CYBetaVT*(SVT/S))*(1 + 3 * taperw)/(12 * (1 + taperw));
var ClrHatTable = ClrHato*Math.cos(AlphaRad);
var CnpHatTable = elemmult((- CLTable*((1 + 3*taperw)/(1 + taperw))/12), Math.cos(AlphaRad));
var CnrHatVT = -2*CnBetaVT*(lVT/b);
var CnrHatWing = -0.103*elemmult(CLTable,CLTable) - 0.4*CDoWing; // (from Seckel)
var CnrHato = CnrHatVT + CnrHatWing;
var CnrHatTable = elemmult(CnrHato, Math.cos(AlphaRad));

//läuft
sign = function(x)
{
if (x<0) { return -1}
else if (x==0) { return 0}
else { return 1}

}



//läuft
min = function(a,b)
{
  if (a<=b){
    return a;
  }
  else {
    return b;
  }
}





//läuft
max = function(a,b)
{
  if (a>=b){
    return a;
  }
  else {
    return b;
  }
}






zeros = function(i,j)
{
 var zeros=[];
 for(i;i<=j;i++)
 {
    zeros = [zeros,0]
 }

 return zeros;
}





//läuft
elemmult = function(array1, array2)
{
 var l = array1.length;
 var res = [];
 for(var i=0;i<=l-1;i++)
 {
    res[i] = array1[i]*array2[i];
 }

 return res;
}





//läuft
elemdiv = function(array1, array2)
{
 var l = array1.length;
 var res = [];
 for(var i=0;i<=l-1;i++)
 {
    res[i] = array1[i]/array2[i];
 }

 return res;
}





//läuft
elemqu = function(array1)
{
 var l = array1.length;
 var res = [];
 for(var i=0;i<=l-1;i++)
 {
    res[i] = array1[i]*array1[i];
 }

 return res;
}





//läuft
find = function(x,y) //nur fuer sortierte Vektoren
{
 var n = x.length;
 var ind = [];
 for (var i = 0; x[i] <= y; i++)
 {
    ind[i] = i;
 }
 return ind;
}



//läuft
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

 return x;

}





//läuft
interp1 = function(x,y,xi) //xi muss ein Skalar sein!!
{
 x = sort(x); // sortiere x
 var r1 = find(x,xi); // Vektor der die Indizes speichert, die kleinergleich xi sind
 var n1 = r1.length; // Laenge des Indizevektors
 var r = r1[n1-1]; // maximaler Eintrag des Indizevektors
 var n = x.length; //
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
    var u = (xi-x[r-1])/(x[r]-x[r-1]);
    yi=y[r-1]+(y[r]-y[r-1])*u;
 }
 else
 {
    yi = NaN;
 }

 /* if ((min(size(yi))==1) & (prod(siz)>1)) // prod(siz) ist bei uns 1!! Sonst Fehler!!
 {
 yi = reshape(yi,siz);
 }*/
 return yi;
}





//läuft durch aber es kommt immer der Nullvektor raus!!
Windfield=function(height,phir,thetar,psir)
{
 /*// FLIGHT Wind Field Interpolation for 3-D wind as a Function of Altitude
 // June 12, 2015
 // ===============================================================
 // Copyright 2006-2015 by ROBERT F. STENGEL. All rights reserved.
 */
 var windh = [-10, 0, 100, 200, 500, 1000, 2000, 4000, 8000, 16000]; // Height, m
 var windx = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; // Northerly wind, m/s
 var windy = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; // Easterly wind, m/s
 var windz = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]; // Vertical wind. m/s
 var winde = [interp1(windh,windx,height),interp1(windh,windy,height),interp1(windh,windz,height)]; // Earth-relative frame
 var HEB = DCM(phir,thetar,psir);
 //var windb = HEB * winde; // Body-axis frame
var windb=[];
windb[0]=HEB[0][0]*winde[0]+HEB[0][1]*winde[1]+HEB[0][2]*winde[2];
windb[1]=HEB[1][0]*winde[0]+HEB[1][1]*winde[1]+HEB[1][2]*winde[2];
windb[2]=HEB[2][0]*winde[0]+HEB[2][1]*winde[1]+HEB[2][2]*winde[2];
return windb;
}





//läuft
LinModel=function(tj,xj)
{
 /*// FLIGHT Equations of Motion for Linear Model (Jacobian) Evaluation,
 // with dummy state elements added for controls

 // June 12, 2015
 // ===============================================================
 // Copyright 2006 by ROBERT F. STENGEL. All rights reserved.
 */
 var x = [];
 var u = [];
   for (var i = 0; i <= 11; i++)
   {
      x[i] = xj[i];
   }

   for (var i = 12; i <= 18; i++)
   {
      u[i-12] = xj[i];
   }


 var xdot = EoM(tj,x);
 var xdotj = [xdot,0,0,0,0,0,0,0];

 return xdotj;
}





//läuft
Atmos = function(geomAlt)
{
 /* 1976 U.S. Standard Atmosphere Interpolation for FLIGHT

 // June 12, 2015
 // ===============================================================
 // Copyright 2006-15 by ROBERT F. STENGEL. All rights reserved.

 // Note: Function does not extrapolate outside altitude range
 // Input: Geometric Altitude, m (positive up)
 // Output: Air Density, kg/m^3
 //         Air Pressure, N/m^2
 //         Air Temperature, K
 //         Speed of Sound, m/s

 //Values Tabulated by Geometric Altitude*/
 var Z = [-1000,0,2500,5000,10000,11100,15000,20000,47400,51000];
 var H = [-1000,0,2499,4996,9984,11081,14965,19937,47049,50594];
 var ppo = [1,1,0.737,0.533,0.262,0.221,0.12,0.055,0.0011,0.0007];
 var rro = [1,1,0.781,0.601,0.338,0.293,0.159,0.073,0.0011,0.0007];
 var T = [288.15,288.15,271.906,255.676,223.252,216.65,216.65,216.65,270.65,270.65];
 var a = [340.294,340.294,330.563,320.545,299.532,295.069,295.069,295.069,329.799,329.799];
 var R = 6367435; // Mean radius of the earth, m
 var Dens = 1.225; // Air density at sea level, Kg/m^3
 var Pres = 101300; // Air pressure at sea level, N/m^2

 // Geopotential Altitude, m
 var geopAlt = R * geomAlt / (R + geomAlt);

 // Linear Interpolation in Geopotential Altitude for Temperature and Speed of Sound
 var temp  = interp1(Z,T,geopAlt);
 var soundSpeed = interp1(Z,a,geopAlt);

 // Exponential Interpolation in Geometric Altitude for Air Density and Pressure
 var betap, betar, airDens, airPres;
 for(k=1;k<=9;k++)
 {
    if (geomAlt <= Z[k])
    {
      betap = Math.log(ppo[k] / ppo[k-1]) / (Z[k] - Z[k-1]);
      betar = Math.log(rro[k] / rro[k-1]) / (Z[k] - Z[k-1]);
      airPres = Pres * ppo[k-1] * Math.exp(betap * (geomAlt - Z[k-1]));
      airDens = Dens * rro[k-1] * Math.exp(betar * (geomAlt - Z[k-1]));
      break;
    }
 }

 var atmos = [airDens,airPres,temp,soundSpeed];

 return atmos;
}





//läuft
TrimCost=function(OptParam)
{
 /* FLIGHT Cost Function for Longitudinal Trim in Steady Level Flight
 // June 12, 2015
 // ===============================================================
 // Copyright 2006-2015 by ROBERT F. STENGEL. All rights reserved.
 */
 var R = new Array(3);

 for (var i=0;i<R.length;i++){
   R[i] = new Array (3);
 }

 R[0][0] = 1;
 R[0][1] = 0;
 R[0][2] = 0;
 R[1][0] = 0;
 R[1][1] = 1;
 R[1][2] = 0;
 R[2][0] = 0;
 R[2][1] = 0;
 R[2][2] = 1;

 /*
 // Optimization Vector:
 // 1 = Stabilator, rad
 // 2 = Throttle, //
 // 3 = Pitch Angle, rad
 */

 u = [u[0], u[1], u[2], OptParam[1], u[4], u[5], OptParam[0]]; // Wird das u aus einer anderen Funktion geholt? --> this.u

 x = [V * Math.cos(OptParam[2]),  x[1],  V * Math.sin(OptParam[2]),  x[3],  x[4], x[5], x[6], x[7], x[8], x[9], OptParam[2], x[11]];  // Wird das x aus einer anderen Funktion geholt? --> this.x

 var xdot = EoM(1,x);
 var xCost = [xdot[1], xdot[3], xdot[8]];
 var J = R[0][0]*xCost[0]*xCost[0]+R[0][1]*xCost[0]*xCost[1]+R[0][2]*xCost[0]*xCost[2]+R[1][0]*xCost[0]*xCost[1]+R[1][1]*xCost[1]*xCost[1]+R[1][2]*xCost[1]*xCost[2]+R[2][0]*xCost[0]*xCost[2]+R[2][1]*xCost[1]*xCost[2]+R[2][2]*xCost[2]*xCost[2];
 var ParamCost = [];
 for (var k=0;k<OptParam.length;k++){
   ParamCost[k]=OptParam[k];
 }
 ParamCost[OptParam.length-1]=J;


 TrimHist = [TrimHist, ParamCost];

 return J;
}





//läuft
DCM = function(Phi,Theta,Psi)
{
  var sinR = Math.sin(Phi);
var cosR = Math.cos(Phi);
var sinP = Math.sin(Theta);
var cosP = Math.cos(Theta);
var sinY = Math.sin(Psi);
var cosY = Math.cos(Psi);
var H = new Array(3);

for (var i=0;i<H.length;i++){
  H[i] = new Array (3);
}

H[0][0] = cosP * cosY;
H[0][1] = cosP * sinY;
H[0][2] = -sinP;
H[1][0] = sinR * sinP * cosY - cosR * sinY;
H[1][1] = sinR * sinP * sinY + cosR * cosY;
H[1][2] = sinR * cosP;
H[2][0] = cosR * sinP * cosY + sinR * sinY;
H[2][1] = cosR * sinP * sinY - sinR * cosY;
H[2][2] = cosR * cosP;

 return H;
}





//läuft
event = function(t,x)
{
 var value = x[5];
 var isterminal = 1;
 var direction = 1;

 var eventausgabe = [value, isterminal, direction];

 return eventausgabe;
}






/* RMQ computes the rotation matrix as a function of quaternions.
The rotation matrices transform vectors from the earth-relative frame of
reference to the body-axis frame.
Quaternion: q1, x Component of quaternion q2, y Component of quaternion
q3, z Component of quaternion q4, cos(Euler) Component of quaternion */

//läuft
RMQ = function(q1,q1,q3,q4)
{
  var H = new Array(3);

  for (var i=0;i<H.length;i++){
    H[i] = new Array (3);
  }
  H[0][0] = Math.pow(q1,2) - Math.pow(q2,2) - Math.pow(q3,2) + Math.pow(q4,2);
  H[0][1] = 2*(q1*q2 + q3*q4);
  H[0][2] = 2*(q1*q3 - q2*q4);
  H[1][0] = 2*(q1*q2 - q3*q4);
  H[1][1] = -Math.pow(q1,2) + Math.pow(q2,2) - Math.pow(q3,2) + Math.pow(q4,2);
  H[1][2] = 2*(q2*q3 + q1*q4);
  H[2][0] = 2*(q1*q3 + q2*q4);
  H[2][1] = 2*(q2*q3 - q1*q4);
  H[2][2] = -Math.pow(q1,2) - Math.pow(q2,2) + Math.pow(q3,2) + Math.pow(q4,2);

 return H;
}





//läuft durch aber bis jetzt nur mit NaN
AeroModelAlpha = function(x,u,Mach,alphar,betar,V)
{

 // FLIGHT Aerodynamic Coefficients of the Aircraft, Thrust Model,
 // and Geometric and Inertial Properties

 // BizJet B
 // ========
 // June 12, 2015
 // ===============================================================
 // Copyright 2006-2015 by ROBERT F. STENGEL. All rights reserved.

 // High-Angle-of-Attack, Mach-Independent Model
 // Landing Gear: Up (GEAR = 0)
 // Flap Setting, u[6]: 0 deg
 // Symmetric Spoiler: Closed (SPOIL = 0)

 // global b cBar SMI

 // Inertial, Geometric, and Aerodynamic Properties
 // *** GeoMassAero.m must first be run to save InerGeo.mat, DataTable.mat,
 // and RotCont.mat ***

 var alphadeg = 57.2957795 * alphar;

 // Thrust Properties
 var  StaticThrust = 2*6.49*Math.pow(10,3); // Static Thrust @ Sea Level, N
 // Current Thrust
 var atmos = Atmos(-x[5]);
 var airDens = atmos[0];
 var airPres = atmos[1];
 var temp = atmos[2];
 var soundSpeed = atmos[3];
 var Thrust = u[3] * StaticThrust * Math.pow((airDens / 1.225),0.7) * (1 - Math.exp((-x[5] - 17000) / 2000)); // Thrust at Altitude, N

 // Current Longitudinal Characteristics
 // ====================================

 // Lift Coefficient
 var  CLStatic = interp1(AlphaTable,CLTable,alphadeg); // Static Lift Coefficient
 var  CLqr = CLqHat * cBar/(2*V); // Pitch-Rate Effect, per rad/s
 var  CLdEr  = interp1(AlphaTable,CLdETable,alphadeg); // Elevator Effect, per rad
 var  CLdSr = CLdEr; // Stabilator Effect, per rad
 var  CL = CLStatic + CLqr*x[7] + CLdSr*u[6] + CLdEr*u[0]; // Total Lift Coefficient

 // Drag Coefficient
 var CDStatic = interp1(AlphaTable,CDTable,alphadeg); // Static Drag Coefficient
 var CD = CDStatic; // Total Drag Coefficient
 // Pitching Moment Coefficient
 var  CmStatic = interp1(AlphaTable,CmTable,alphadeg); // Static Pitching Moment Coefficient
 var  CmdEr  = interp1(AlphaTable,CmdETable,alphadeg); // Elevator Effect, per rad
 var Cmqr = -CLqHat*(lHT/cBar) * cBar/(2*V); // Pitch-Rate + Alpha-Rate Effect, per rad/s
 var  CmdSr = CmdEr; // Stabilator Effect, per rad
 var  Cm = CmStatic - CL*SMI + Cmqr*x[7] + CmdSr*u[6] + CmdEr*u[0]; // Total Pitching Moment Coefficient

 // Current Lateral-Directional Characteristics
 // ===========================================

 // Rolling Moment Coefficient
 var  ClBr = interp1(AlphaTable,ClBetaTable,alphadeg); // Dihedral Effect, per rad
 var  ClpHat = interp1(AlphaTable,ClpHatTable,alphadeg);
 var  Clpr = ClpHat * (b / (2 * V));
 var  ClrHat = interp1(AlphaTable,ClrHatTable,alphadeg); // Roll-Rate Effect, per rad/s
 var  Clrr = ClrHat * (b / (2 * V)); // Yaw-Rate Effect, per rad/s
 var  CldAr = interp1(AlphaTable,CldATable,alphadeg); // Aileron Effect, per rad
 var  CldRr = interp1(AlphaTable,CldRTable,alphadeg); // Rudder Effect, per rad
 var  CldASr = 0; // Asymmetric Spoiler Effect, per rad
 var  Cl = (ClBr*betar + CldRr*u[2]) + Clrr * x[8] + Clpr * x[6]+ (CldAr*u[1] + CldASr*u[4]); // Total Rolling-Moment Coefficient

 // Side-Force Coefficient
 var  CYBr = interp1(AlphaTable,CYBetaTable,alphadeg); // Side-Force Slope, per rad
 var  CYdAr = CYdAo; // Aileron Effect, per rad
 var  CYdRr = 0.1574; // Rudder Effect, per rad
 var  CYdASr = 0; // Asymmetric Spoiler Effect, per rad
 var  CY = (CYBr*betar + CYdRr*u[2]) + (CYdAr*u[1] + CYdASr*u[4]); // Total Side-Force Coefficient

 // Yawing Moment Coefficient
 var  CnBr = interp1(AlphaTable,CnBetaTable,alphadeg); // Directional Stability, per rad
 var  Cnpr = CL * (1 + 3 * taperw)/(12 * (1 + taperw)) * (b / (2 * V)); // Roll-Rate Effect, per rad/s
 var  Cnrr = (-2 * (lVT / b) * CnBr - 0.1 * Math.pow(CL,2)) * (b / (2 * V));
 var  CnpHat = interp1(AlphaTable,CnpHatTable,alphadeg);
 var  Cnpr = CnpHat * (b / (2 * V));
 var  CnrHat = interp1(AlphaTable,CnrHatTable,alphadeg); // Roll-Rate Effect, per rad/s
 var  Cnrr = CnrHat * (b / (2 * V)); // Yaw-Rate Effect, per rad/s
 var  CndAr = interp1(AlphaTable,CndATable,alphadeg); // Aileron Effect, per rad
 var  CndRr = interp1(AlphaTable,CndRTable,alphadeg); // Rudder Effect, per rad
 var  CndASr = 0; // Asymmetric Spoiler Effect, per rad
 var  Cn = (CnBr*betar + CndRr*u[2]) + Cnrr * x[8] + Cnpr * x[6]+ (CndAr*u[1] + CndASr*u[4]); // Total Yawing-Moment Coefficient
 var ret = [CD,CL,CY,Cl,Cm,Cn,Thrust];

 return ret;
}





//läuft aber Werte noch überprüfen
AeroModelMach = function(x,u,Mach,alphar,betar,V)
{

 //Typical Mass and Inertial Properties
 var m = 4536; // Mass, kg
 var Ixx = 35926.5; // Roll Moment of Inertia, kg-m^2
 var Iyy = 33940.7; // Pitch Moment of Inertia, kg-m^2
 var Izz = 67085.5; // Yaw Moment of Inertia, kg-m^2
 var Ixz = 3418.17; // Nose-high(low) Product of Inertia, kg-m^2

 //Geometric Properties
 var cBar=2.14; // Mean Aerodynamic Chord, m
 var b=10.4; // Wing Span, m
 var S=21.5; // Reference Area, m^2
 var ARw=5.02; // Wing Aspect Ratio
 var taperw=0.507; // Wing Taper Ratio
 var sweepw=13 * .01745329; // Wing 1/4-chord sweep angle, rad
 var ARh=4; // Horizontal Tail Aspect Ratio
 var sweeph=25 * .01745329; // Horiz Tail 1/4-chord sweep angle, rad
 var ARv=0.64; // Vertical Tail Aspect Ratio
 var sweepv=40 * .01745329; // Vert Tail 1/4-chord sweep angle, rad
 var lvt=4.72; // Vert Tail Length, m

 //Thrust Properties
 var StaticThrust = 26243.2; // Static Thrust @ Sea Level, N

 //Current Thrust
 var atmos = [];
 atmos=Atmos(-x[5]);
 var airDens = atmos[0];
 var airPres = atmos[1];
 var temp = atmos[2];
 var soundSpeed = atmos[3];
 var Thrust=u[3] * StaticThrust * Math.pow((airDens / 1.225),0.7)* (1 - Math.exp((-x[5] - 17000) / 2000)); // Thrust at Altitude, N

 //Current Mach Effects, normalized to Test Condition B (Mach = 0.1734)
 var PrFac=1 / (Math.sqrt(1 - Math.pow(Mach,2) * 1.015)); // Prandtl Factor
 var WingMach=1 / ((1 + Math.sqrt(1 + (Math.pow((ARw/(2 * Math.cos(sweepw))),2))* (1 - Math.pow(Mach,2) * Math.cos(sweepw)))) * 0.268249); // Modified Helmbold equation
 var HorizTailMach=1 / ((1 + Math.sqrt(1 + (Math.pow((ARh/(2 * Math.cos(sweeph))),2))* (1 - Math.pow(Mach,2) * Math.cos(sweeph)))) * 0.294539); // Modified Helmbold equation
 var VertTailMach=1 / ((1 + Math.sqrt(1 + (Math.pow((ARv/(2 * Math.cos(sweepv))),2))* (1 - Math.pow(Mach,2) * Math.cos(sweepv)))) * 0.480338); // Modified Helmbold equation

 //Current Longitudinal Characteristics
 //====================================

 //Lift Coefficient
 var CLo=0.1095; // Zero-AoA Lift Coefficient (B)

 if (GEAR >= 1)
 {
    CLo=CLo - 0.0192; // Gear-down correction
 }

 if (u[5] >= 0.65)
 {
    CLo=CLo + 0.5182; // 38 deg-flap correction
 }

 if (SPOIL >= 1)
 {
    CLo=CLo - 0.1897; // 42 deg-Symmetric Spoiler correction
 }

 var CLar=5.6575; // Lift Slope (B), per rad

 if (u[5] >= 0.65)
 {
    CLar=CLar - 0.0947;
 }

 var CLqr=4.231 * cBar / (2 * V); // Pitch-Rate Effect, per rad/s
 var CLdSr=1.08; // Stabilator Effect, per rad

 if (u[5] >= 0.65)
 {
    CLdSr=CLdSr - 0.4802; // 38 deg-flap correction
 }

 var CLdEr=0.5774; // Elevator Effect, per rad

 if ( u[5] >= 0.65)
 {
    CLdEr=CLdEr - 0.2665; // 38 deg-flap correction
 }

 var CL=CLo + (CLar*alphar + CLqr*x[8] + CLdSr*u[7] + CLdEr*u[1])* WingMach; // Total Lift Coefficient, w/Mach Correction

 //Drag Coefficient
 var CDo=0.0255; // Parasite Drag Coefficient (B)

 if ( GEAR >= 1)
 {
    CDo=CDo + 0.0191; // Gear-down correction
 }

 if ( u[5] >= 0.65)
 {
    CDo=CDo + 0.0836; // 38 deg-flap correction
 }

 if ( SPOIL >= 1)
 {
    CDo=CDo + 0.0258; // 42 deg-Symmetric Spoiler correction
 }

 var epsilon=0.0718; // Induced Drag Factor

 if ( u[5] >= 0.65)
 {
    epsilon=0.079; // 38 deg-flap correction
 }

 var CD=CDo * PrFac + epsilon * Math.pow(CL,2); // Total Drag Coefficient, w/Mach Correction

 //Pitching Moment Coefficient
 var Cmo=0; // Zero-AoA Moment Coefficient (B)

 if ( GEAR >= 1)
 {
    Cmo=Cmo + 0.0255; // Gear-down correction
 }

 if ( u[5] >= 0.65)
 {
    Cmo=Cmo - 0.058; // 38 deg-flap correction
 }

 if ( SPOIL >= 1)
 {
    Cmo=Cmo - 0.0154; // 42 deg-Symmetric Spoiler correction
 }

 var Cmar=-1.231; // Static Stability (B), per rad

 if ( u[5] >= 0.65)
 {
    Cmar=Cmar + 0.0138;
 }

 var Cmqr = -18.8 * cBar / (2 * V); // Pitch-Rate + Alpha-Rate Effect, per rad/s
 var CmdSr=-2.291; // Stabilator Effect, per rad

 if ( u[5] >= 0.65)
 {
    CmdSr=CmdSr + 0.121; // 38 deg-flap correction
 }

 var CmdEr=-1.398; // Elevator Effect, per rad

 if ( u[5] >= 0.65)
 {
    CmdEr=CmdEr + 0.149; // 38 deg-flap correction
 }

 var Cm=Cmo + (Cmar*alphar + Cmqr*x[7] + CmdSr*u[6] + CmdEr*u[0])* HorizTailMach; // Total Pitching Moment Coefficient, w/Mach Correction

 //Current Lateral-Directional Characteristics
 //===========================================

 //Side-Force Coefficient
 var CYBr=-0.7162; // Side-Force Slope (B), per rad

 if ( u[5] >= 0.65)
 {
    CYBr=CYBr + 0.0826;
 }

 var CYdAr=-0.00699; // Aileron Effect, per rad
 var CYdRr=0.1574; // Rudder Effect, per rad

 if ( u[5] >= 0.65)
 {
    CYdRr=CYdRr - 0.0093; // 38 deg-flap correction
 }

 var CYdASr=0.0264; // Asymmetric Spoiler Effect, per rad

 if ( u[5] >= 0.65)
 {
    CYdASr=CYdASr + 0.0766; // 38 deg-flap correction
 }

 var CY=(CYBr*betar + CYdRr*u[2]) * VertTailMach+ (CYdAr*u[1] + CYdASr*u[4]) * WingMach; // Total Side-Force Coefficient, w/Mach Correction

 //Yawing Moment Coefficient
 var CnBr=0.1194; // Directional Stability (B), per rad

 if ( u[5] >= 0.65)
 {
    CnBr=CnBr - 0.0092;
 }

 var Cnpr=CL * (1 + 3 * taperw)/(12 * (1 + taperw)) * (b / (2 * V)); // Roll-Rate Effect, per rad/s
 var Cnrr=(-2 * (lvt / b) * CnBr * VertTailMach - 0.1 * Math.pow(CL,2))* (b / (2 * V)); // Yaw-Rate Effect, per rad/s
 var CndAr=0; // Aileron Effect, per rad

 if ( u[5] >= 0.65)
 {
    CndAr=CndAr + 0.0028;
 }

 var CndRr=-0.0713; // Rudder Effect, per rad

 if ( u[5] >= 0.65)
 {
    CndRr=CndRr - 0.0185; // 38 deg-flap correction
 }

 var CndASr=-0.0088; // Asymmetric Spoiler Effect, per rad

 if ( u[5] >= 0.65)
 {
    CndASr=CndASr - 0.0106; // 38 deg-flap correction
 }

 var Cn=(CnBr*betar + CndRr*u[3]) * VertTailMach+ Cnrr * x[9] + Cnpr * x[7]+ (CndAr*u[2] + CndASr*u[5]) * WingMach; // Total Yawing-Moment Coefficient, w/Mach Correction

 //Rolling Moment Coefficient
 var ClBr=-0.0918; // Dihedral Effect (B), per rad

 if ( u[5] >= 0.65)
 {
    ClBr=ClBr - 0.0092;
 }

 var Clpr=-CLar * (1 + 3 * taperw)/(12 * (1 + taperw))* (b / (2 * V)); // Roll-Rate Effect, per rad/s
 var Clrr=(CL * (1 + 3 * taperw)/(12 * (1 + taperw))* (Math.pow(Mach * Math.cos(sweepw),2) - 2) / (Math.pow(Mach * Math.cos(sweepw),2) - 1))* (b / (2 * V)); // Yaw-Rate Effect, per rad/s
 var CldAr=0.1537; // Aileron Effect, per rad

 if ( u[5] >= 0.65)
 {
    CldAr=CldAr + 0.01178;
 }

 var CldRr=0.01208; // Rudder Effect, per rad

 if ( u[5] >= 0.65)
 {
    CldRr=CldRr + 0.01115; // 38 deg-flap correction
 }

 var CldASr=-0.01496; // Asymmetric Spoiler Effect, per rad

 if ( u[5] >= 0.65)
 {
    CldASr=CldASr - 0.02376; // 38 deg-flap correction
 }

 var Cl=(ClBr*betar + CldRr*u[2]) * VertTailMach+ Clrr * x[8] + Clpr * x[6]+ (CldAr*u[1] + CldASr*u[4]) * WingMach; // Total Rolling-Moment Coefficient, w/Mach Correction
 var ret = [];
 ret = [CD,CL,CY,Cl,Cm,Cn,Thrust];

 return ret;
}





EoM = function(t,x)
{
 /*
 // FLIGHT Equations of Motion

 // June 12, 2015
 // ===============================================================
 // Copyright 2006-2015 by ROBERT F. STENGEL. All rights reserved.

 global m Ixx Iyy Izz Ixz S b cBar CONHIS u tuHis deluHis uInc MODEL RUNNING

 // Select Aerodynamic Model

 if MODEL == 0
 AeroModel = @AeroModelAlpha;
 else
 AeroModel = @AeroModelMach;
 end
 */
 var ev = event(t,x);
 var value = ev[0];
 var isterminal = ev[1];
 var direction = ev[2];

 // Earth-to-Body-Axis Transformation Matrix
 var HEB  = DCM(x[9],x[10],x[11]);

 // Atmospheric State
 x[5] = min(x[5],0); // Limit x[6] to <= 0 m
 var atmos = Atmos(-x[5]);
 var airDens = atmos[0];
 var airPres = atmos[1];
 var temp = atmos[2];
 var soundSpeed = atmos[3];

 // Body-Axis Wind Field
 var windb = WindField(-x[5], x[9], x[10], x[11]);

 // Body-Axis Gravity Components
 var gb = [];

  gb[0]=HEB[0][2]*9.80665;
  gb[1]=HEB[1][2]*9.80665;
  gb[2]=HEB[2][2]*9.80665;

 // Air-Relative Velocity Vector
 x[0] = max(x[0],0); // Limit axial velocity to >= 0 m/s
 var Va  = [];
 for (var i=0;i<3;i++){
   Va[i]=x[i]+windb[i];
 }
 var V  = Math.sqrt(Va[0]*Va[0]+Va[1]*Va[1]+Va[2]*Va[2]);
 var alphar = Math.atan(Va[2] / Math.abs(Va[0])); // alphar = min(alphar, (pi/2 - 1e-6)); // Limit angle of attack to <= 90 deg
 var alpha  = 57.2957795 * alphar;
 var betar =  Math.asin(Va[1] / V);
 var beta =  57.2957795 * betar;
 var Mach =  V / soundSpeed;
 var qbar = 0.5 * airDens * Math.pow(V,2);

 // Incremental Flight Control Effects
 var uTotal;

 if(CONHIS >=1 && RUNNING == 1)
 {
    [uInc] = interp1(tuHis,deluHis,t);


for(var i=0;i<lenght.uInc;i++){
    uTotal[i] = u[i] + uInc[i];
}

 }
 else
 {
   for(var i=0;i<lenght.uInc;i++){
       uTotal[i] = u[i];
   }

 }

 // Force and Moment Coefficients; Thrust
 var mod;

 if(MODEL==0)
 {
    mod = AeroModelAlpha(x,uTotal,Mach,alphar,betar,V);
 }
 else
 {
    mod = AeroModelMach(x,uTotal,Mach,alphar,betar,V);
 }

 var CD = mod[0];
 var CL = mod[1];
 var CY = mod[2];
 var Cl = mod[3];
 var Cm = mod[4];
 var Cn = mod[5];
 var Thrust = mod[6];
 var qbarS = qbar * S;
 var CX = -CD * Math.cos(alphar) + CL * Math.sin(alphar); // Body-axis X coefficient
 var CZ = -CD * Math.sin(alphar) - CL * Math.cos(alphar); // Body-axis Z coefficient

 // State Accelerations
 var Xb = (CX * qbarS + Thrust) / m;
 var Yb = CY * qbarS / m;
 var Zb = CZ * qbarS / m;
 var Lb = Cl * qbarS * b;
 var Mb = Cm * qbarS * cBar;
 var Nb = Cn * qbarS * b;
 var nz = -Zb / 9.80665; // Normal load factor

 // Dynamic Equations
 var xd1 = Xb + gb[0] + x[8] * x[1] - x[7] * x[2];
 var xd2 = Yb + gb[1] - x[8] * x[0] + x[6] * x[2];
 var xd3 = Zb + gb[2] + x[7] * x[0] - x[6] * x[1];
 var y = [];
 for(var i=0;i<3;i++){
   y[i]=HEB[i][0]*x[0]+HEB[i][1]*x[1]+HEB[i][2]*x[2];
 }
 var xd4 = y[0];
 var xd5 = y[1];
 var xd6 = y[2];
 var xd7 = (Izz * Lb + Ixz * Nb - (Ixz * (Iyy - Ixx - Izz) * x[6] +(Math.pow(Ixz,2) + Izz * (Izz - Iyy)) * x[8]) * x[7]) / (Ixx * Izz - Math.pow(Ixz,2));
 var xd8 = (Mb - (Ixx - Izz) * x[6] * x[8] - Ixz * (Math.pow(x[6],2) - Math.pow(x[8],2))) / Iyy;
 var xd9 = (Ixz * Lb + Ixx * Nb + (Ixz * (Iyy - Ixx - Izz) * x[8] +(Math.pow(Ixz,2) + Ixx * (Ixx - Iyy)) * x[6]) * x[7]) / (Ixx * Izz - Math.pow(Ixz,2));
 var cosPitch = Math.cos(x[10]);

 if( Math.abs(cosPitch) <= 0.00001)
 {
    cosPitch = 0.00001 * sign(cosPitch); // sign muss noch geschrieben werden!!
 }

 var tanPitch = Math.sin(x[10]) / cosPitch;
 var xd10 = x[6] + (Math.sin(x[9]) * x[7] + Math.cos(x[9]) * x[8]) * tanPitch;
 var xd11 = Math.cos(x[9]) * x[7] - Math.sin(x[9]) * x[8];
 var xd12 = (Math.sin(x[9]) * x[7] + Math.cos(x[9]) * x[8]) / cosPitch;
 var xdot = [xd1,xd2,xd3,xd4,xd5,xd6,xd7,xd8,xd9,xd10,xd11,xd12];

 return xdot;
}


// Bis hierhin bin ich beim Durchschauen gekommen!!




EoMQver2 = function(t,x)
{
 /*// FLIGHT Equations of Motion
 // Quaternion Option

 // September 10, 2016
 // ===============================================================
 // Copyright 2006-2016 by ROBERT F. STENGEL. All rights reserved.

 global m Ixx Iyy Izz Ixz S b cBar CONHIS u tuHis deluHis uInc MODEL RUNNING

 // Select Aerodynamic Model

 if MODEL == 0
 AeroModel = @AeroModelAlpha;
 else
 AeroModel = @AeroModelMach;
 end
 */
 var ev = event(t,x);
 var value = ev[1];
 var isterminal = ev[2];
 var direction = ev[3];

 // Earth-to-Body-Axis Transformation Matrix
 var HEB  = RMQ(x[10],x[11],x[12],x[13]);
 // Atmospheric State
 x[6] = min(x[6],0); // Limit x[6] to <= 0 m
 var atmos = Atmos(-x[6]);
 var airDens = atmos[0];
 var airPres = atmos[1];
 var temp = atmos[2];
 var soundSpeed = atmos[3];

 // Body-Axis Wind Field
 var Phi = Math.atan2(2*(x[10]*x[13] + x[11]*x[12]),(1 - 2*(Math.pow(x[10],2) + Math.pow(x[11],2))));
 var Theta = Math.asin(2*(x[11]*x[13] - x[10]*x[12]));
 var Psi = Math.atan2(2*(x[12]*x[13] + x[10]*x[11]),(1 - 2*(Math.pow(x[11],2) + Math.pow(x[12],2))));
 var windb = WindField(x[3],Phi,Theta,Psi);

 // Body-Axis Gravity Components
 var gb = HEB * [0;0;9.80665];

 // Air-Relative Velocity Vector
 x[1] = max(x[1],0); // Limit axial velocity to >= 0 m/s
 var Va = [x[1];x[2];x[3]] + windb;
 var V = Math.sqrt(transpose(Va) * Va); // transpose muss noch geschrieben werden!!
 var alphar = Math.atan(Va[3] / Math.abs(Va[1])); // alphar = min(alphar, (pi/2 - 1e-6)); // Limit angle of attack to <= 90 deg
 var alpha = 57.2957795 * alphar;
 var betar = Math.asin(Va[2] / V);
 var beta = 57.2957795 * betar;
 var Mach = V / soundSpeed;
 var qbar = 0.5 * airDens * Math.pow(V,2);

 // Incremental Flight Control Effects
 var uTotal;

 if(CONHIS >=1 && RUNNING == 1)
 {
    [uInc] = interp1(tuHis,deluHis,t);
    uInc = transpose(uInc); // transpose muss noch geschrieben werden!!
    uTotal = u + uInc;
 }
 else
 {
    uTotal = u;
 }

 // Force and Moment Coefficients; Thrust
 var mod;

 if(MODEL==0)
 {
    mod = AeroModelAlpha(x,uTotal,Mach,alphar,betar,V);
 }
 else
 {
    mod = AeroModelMach(x,uTotal,Mach,alphar,betar,V);
 }

 var CD = mod[1];
 var CL = mod[2];
 var CY = mod[3];
 var Cl = mod[4];
 var Cm = mod[5];
 var Cn = mod[6];
 var Thrust = mod[7];
 var qbarS = qbar * S;
 var CX = -CD * Math.cos(alphar) + CL * Math.sin(alphar); // Body-axis X coefficient
 var CZ =  -CD * Math.sin(alphar) - CL * Math.cos(alphar); // Body-axis Z coefficient

 // State Accelerations
 var Xb = (CX * qbarS + Thrust) / m;
 var Yb = CY * qbarS / m;
 var Zb = CZ * qbarS / m;
 var Lb = Cl * qbarS * b;
 var Mb = Cm * qbarS * cBar;
 var Nb = Cn * qbarS * b;
 var nz = -Zb / 9.80665; // Normal load factor

 // Dynamic Equations
 var xd1 = Xb + gb[1] + x[9] * x[2] - x[8] * x[3];
 var xd2 = Yb + gb[2] - x[9] * x[1] + x[7] * x[3];
 var xd3 = Zb + gb[3] + x[8] * x[1] - x[7] * x[2];
 var y = transpose(HEB) * [x[1];x[2];x[3]]; // transpose muss noch geschrieben werden!!
 var xd4 = y[1];
 var xd5 = y[2];
 var xd6 = y[3];
 var xd7 =  (Izz * Lb + Ixz * Nb - (Ixz * (Iyy - Ixx - Izz) * x[7] + (Math.pow(Ixz,2) + Izz * (Izz - Iyy)) * x[9]) * x[8]) / (Ixx * Izz - Math.pow(Ixz,2));
 var xd8 =  (Mb - (Ixx - Izz) * x[7] * x[9] - Ixz * (Math.pow(x[7],2) - Math.pow(x[9],2))) / Iyy;
 var xd9 = (Ixz * Lb + Ixx * Nb + (Ixz * (Iyy - Ixx - Izz) * x[9] +(Math.pow(Ixz,2) + Ixx * (Ixx - Iyy)) * x[7]) * x[8]) / (Ixx * Izz - Math.pow(Ixz,2));

 // Quaternion Propagation
 var p = x[7];
 var q = x[8];
 var r = x[9];
 var Q = 0.5*[0, r, -q, p,-r, 0, p, q, q, -p, 0, r, -p, -q, -r, 0];
 var qVec = [x[10]; x[11]; x[12]; x[13]];
 var qd = Q*qVec;
 var xdot = [xd1;xd2;xd3;xd4;xd5;xd6;xd7;xd8;xd9;qd[1];qd[2];qd[3];qd[4]];

 return xdot;
}






flight = function()
{
 /*// FLIGHT, version 2.0 -- Generic 6-DOF Trim, Linear Model, and Flight Path Simulation
 // Euler Angle/Quaternion Option for Rotation (Direction Cosine) Matrix

 // November 23, 2016
 // ===============================================================
 // Copyright 2006-2016 by ROBERT F. STENGEL. All rights reserved.

 clear
 global GEAR CONHIS SPOIL u x V uInc tuHis deluHis TrimHist SMI MODEL RUNNING


 // This is the SCRIPT FILE. It contains the Main Program, which:
 //  Defines initial conditions
 //  Contains aerodynamic data tables (if required)
 //  Calculates longitudinal trim condition
 //  Calculates stability-and-control derivatives
 //  Simulates flight path using nonlinear equations of motion

 // Functions used by FLIGHT:
 //  AeroModelAlpha.m High-Alpha, Low-Mach aerodynamic coefficients of the aircraft,
 //  thrust model, and geometric and inertial properties
 //  AeroModelMach.m Low-Alpha, High-Mach aerodynamic coefficients of the aircraft,
 //  thrust model, and geometric and inertial properties
 //  Atmos.m Air density, sound speed
 //  DCM.m  Direction-cosine matrix
 //  EoM.m  Equations of motion for integration
 // 	LinModel.m Equations of motion for linear model definition
 // TrimCost.m Cost function for trim solution
 // WindField.m Wind velocity components

 // DEFINITION OF THE STATE VECTOR
 // With Euler Angle DCM option:
 // x[1] =  Body-axis x inertial velocity, ub, m/s
 // x[2] = Body-axis y inertial velocity, vb, m/s
 // x[3] = Body-axis z inertial velocity, wb, m/s
 // x[4] = North position of center of mass WRT Earth, xe, m
 // x[5] = East position of center of mass WRT Earth, ye, m
 // x[6] = Negative of c.m. altitude WRT Earth, ze = -h, m
 // x[7] = Body-axis roll rate, pr, rad/s
 // x[8] = Body-axis pitch rate, qr, rad/s
 // x[9] = Body-axis yaw rate, rr,rad/s
 // x[10] = Roll angle of body WRT Earth, phir, rad
 // x[11] = Pitch angle of body WRT Earth, thetar, rad
 // x[12] = Yaw angle of body WRT Earth, psir, rad
 // With Quaternion DCM option:
 // x[1] =  Body-axis x inertial velocity, ub, m/s
 // x[2] = Body-axis y inertial velocity, vb, m/s
 // x[3] = Body-axis z inertial velocity, wb, m/s
 // x[4] = North position of center of mass WRT Earth, xe, m
 // x[5] = East position of center of mass WRT Earth, ye, m
 // x[6] = Negative of c.m. altitude WRT Earth, ze = -h, m
 // x[7] = Body-axis roll rate, pr, rad/s
 // x[8] = Body-axis pitch rate, qr, rad/s
 // x[9] = Body-axis yaw rate, rr,rad/s
 // x[10] = q1, x Component of quaternion
 // x[11] = q2, y Component of quaternion
 // x[12] = q3, z Component of quaternion
 // x[13] = q4, cos(Euler) Component of quaternion

 // DEFINITION OF THE CONTROL VECTOR
 // u[1] =  Elevator, dEr, rad, positive: trailing edge down
 // u[2] =  Aileron, dAr, rad, positive: left trailing edge down
 // u[3] =  Rudder, dRr, rad, positive: trailing edge left
 // u[4] =  Throttle, dT,
 // u[5] = Asymmetric Spoiler, dASr, rad
 // u[6] = Flap, dFr, rad
 // u[7] = Stabilator, dSr, rad

 // ======================================================================
 // USER INPUTS
 // ===========
 // FLIGHT Flags (1 = ON, 0 = OFF)
 */

 var MODEL = 1; // Aerodynamic model selection:  0: Incompressible flow, high angle of attack;  1: Compressible flow, low angle of attack
 
 var TRIM =  0; // Trim flag (= 1 to calculate trim @ I.C.)
 var LINEAR =  0; // Linear model flag (= 1 to calculate and store F and G)
 var SIMUL = 1; // Flight path flag (= 1 for nonlinear simulation)
 var GEAR =  0; // Landing gear DOWN (= 1) or UP (= 0)
 var SPOIL = 0; // Symmetric Spoiler DEPLOYED (= 1) or CLOSED (= 0)
 var CONHIS = 1; // Control history ON (= 1) or OFF (= 0)
 var dF =  0; // // Flap setting, deg
 var RUNNING = 0; // internal flag, -

 // Initial Altitude (ft), Indicated Airspeed (kt)
 var hft = 10000; // Altitude above Sea Level, ft
 var VKIAS = 150; // Indicated Airspeed, kt
 var hm = hft * 0.3048; // Altitude above Sea Level, m
 var VmsIAS = VKIAS * 0.5154;// Indicated Airspeed, m/s

 // US Standard Atmosphere, 1976, Table Lookup for I.C.
 var atmos_1 = Atmos(hm);
 var airDens = atmos_1[0];
 var airPres = atmos_1[1];
 var temp = atmos_1[2];
 var soundSpeed = atmos_1[3];

 // Dynamic Pressure (N/m^2), and True Airspeed (m/s)
 var qBarSL = 0.5*1.225*Math.pow(VmsIAS,2); // Dynamic Pressure at sea level, N/m^2
 var V = Math.sqrt(2*qBarSL/airDens); // True Airspeed, TAS, m/s
 var TASms = V;

 // Alphabetical List of Initial Conditions
 var alpha = 0; // Angle of attack, deg (relative to air mass)
 var beta = 0; // Sideslip angle, deg (relative to air mass)
 var dA = 0; // Aileron angle, deg
 var dAS = 0; // Asymmetric spoiler angle, deg
 var dE = 0; // Elevator angle, deg
 var dR = 0; // Rudder angle, deg
 var dS = 0; // Stabilator setting, deg
 var dT = 0; // Throttle setting, // / 100
 var hdot = 0; // Altitude rate, m/s
 var p = 0; // Body-axis roll rate, deg/s
 var phi = 0; // Body roll angle wrt earth, deg
 var psi = 0; // Body yaw angle wrt earth, deg
 var q = 0; // Body-axis pitch rate, deg/sec
 var r = 0; // Body-axis yaw rate, deg/s
 var SMI = 0; // Static margin increment due to center-of-mass variation from reference, ///100
 var tf = 1; // Final time for simulation, sec
 var ti = 0; // Initial time for simulation, sec
 var theta = alpha; // Body pitch angle wrt earth, deg [theta = alpha if hdot = 0]
 var xe = 0; // Initial longitudinal position, m
 var ye =  0; // Initial lateral position, m
 var ze = -hm; // Initial vertical position, m [h: + up, z: + down]

 // Initial Conditions Depending on Prior Initial Conditions
 var gamma = 57.2957795 * Math.atan(hdot / Math.sqrt(Math.pow(V,2) - Math.pow(hdot,2))); // Inertial Vertical Flight Path Angle, deg
 var qbar = 0.5 * airDens * Math.pow(V,2); // Dynamic Pressure, N/m^2
 var IAS = Math.sqrt(2 * qbar / 1.225); // Indicated Air Speed, m/s
 var Mach = V / soundSpeed; // Mach Number
 var uInc = [];

 // Initial Control Perturbation (Test Inputs: rad or 100//)
 var delu = [0,0,0,0,0,0,0];

 // Initial State Perturbation (Test Inputs: m, m/s, rad, or rad/s)
 var delx = [0,0,0,0,0,0,0,0,0,0,0,0];

 // Control Perturbation History (Test Inputs: rad or 100//)
 // =======================================================
 // Each control effector represented by a column
 // Each row contains control increment delta-u(t) at time t:
 var tuHis = [0, 33, 67, 100];
 var deluHis = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

 // State Vector and Control Initialization, rad
 var phir = phi * 0.01745329;
 var thetar = theta * 0.01745329;
 var psir = psi * 0.01745329;
 var windb = WindField(-ze,phir,thetar,psir);
 var alphar = alpha * 0.01745329;
 var betar = beta * 0.01745329;
 var x = [V * Math.cos(alphar) * Math.cos(betar) - windb[0],V * Math.sin(betar) - windb[1],V * Math.sin(alphar) * Math.cos(betar) - windb[2],xe,ye,ze,p * 0.01745329,q * 0.01745329,r * 0.01745329,phir,thetar,psir];
 var u = [dE * 0.01745329,dA * 0.01745329,dR * 0.01745329,dT,dAS * 0.01745329,dF * 0.01745329,dS * 0.01745329];

 // Trim Calculation (for Steady Level Flight at Initial V and h)
 // =============================================================
 // Always use Euler Angles for trim calculation
 // Trim Parameter Vector (OptParam):
 // 1 = Stabilator, rad
 // 2 = Throttle,
 // 3 = Pitch Angle, rad

 if(TRIM >= 1)
 {
    var OptParam = [];
    var TrimHist = [];
    var InitParam = [0.0369,0.1892,0.0986];

  //  [OptParam,J,ExitFlag,Output] = fminsearch('TrimCost',InitParam); ÜBERPRÜFEN!!

    // Optimizing Trim Error Cost with respect to dSr, dT, and Theta
    var TrimHist;
    var Index = [];
    
    for(var i = 0; i < TrimHist.length, i++)
    {
      Index[0] = i;
    }
    
    var TrimStabDeg = 57.2957795*OptParam[0];
    var TrimThrusPer = 100*OptParam[1];
    var TrimPitchDeg = 57.2957795*OptParam[2];
    var TrimAlphaDeg = TrimPitchDeg - gamma;

    // Insert trim values in nominal control and state vectors
    var u = [u[0],u[1],u[2],OptParam[1],u[4],u[5],OptParam[0]];
    var x = [V * Math.cos(OptParam[2]),x[1],V * Math.sin(OptParam[2]),x[3],x[4],x[5],x[6],x[7],x[8],x[9],OptParam[2],x[11]];
 }

 // Stability-and-Control Derivative Calculation
 // ============================================


 // Flight Path Simulation, with Quaternion Option
 // ==============================================

 if(SIMUL >= 1)
 {
    RUNNING = 1;

    var xo;
    
    for(var i = 0; i < x.length; i++)
    {
      xo[i] = x[i] + delx[i];
      u[i] = u[i] + delu[i];
    }

    var kHis;

    x = explizit(0,1,xo,1/20);
    


    var VAirRel = [];
    var vEarth = [];
    var AlphaAR = [];
    var BetaAR = [];
    var windBody = [];
    var airDensHis = [];
    var soundSpeedHis = [];
    var qbarHis = [];
    var GammaHis = [];
    var XiHis  = [];
    var Alphar;

    for(var i = 0; i < kHis; i++)
    {
      var windb  = WindField(-x[i][5],x[i][9],x[i][10],x[i][11]);
      var windBody = [windBody, windb];
      [airDens,airPres,temp,soundSpeed] = Atmos(-x[i][5]);
      airDensHis = [airDensHis, airDens];
      soundSpeedHis = [soundSpeedHis, soundSpeed];
    }
    
    var vbody = [];
    
    for(var i = 0; i < kHis; i++)
    {
      vBody  = [x[i][0], x[i][1], x[i][2]];
    }
    
    vBody  = transpose(vBody); // transpose muss noch geschrieben werden!!
    var vBodyAir = vBody + windBody;

    for(var i = 0; i < kHis; i++)
    {
      var vE = transpose(DCM(x[i][9],x[i][10],x[i][11])) * [vBody[0][i], vBody[1][i], vBody[2][i]]; // ÜBERPRÜFEN!
      var VER = Math.sqrt(Math.pow(vE[0],2) + Math.pow(vE[1],2) + Math.pow(vE[2],2));
      var VAR = Math.sqrt(Math.pow(vBodyAir[0][i],2) + Math.pow(vBodyAir[1][i],2) + Math.pow(vBodyAir[2][i],2));
      var VARB = Math.sqrt(Math.pow(vBodyAir[0][i],2) + Math.pow(vBodyAir[2][i],2));

      if(vBodyAir[0][i] >= 0)
      {
	Alphar = Math.asin(vBodyAir[2][i] / VARB);
      }
      else
      {
	Alphar = Math.PI - Math.asin(vBodyAir[2][i] / VARB);
      }

      AlphaAR = [AlphaAR, Alphar];
      var Betar =  Math.asin(vBodyAir[1][i] / VAR);
      BetaAR = [BetaAR, Betar];
      vEarth = [vEarth, vE];
      var Xir = Math.asin(vEarth[1][i] / Math.sqrt(Math.pow(vEarth[0][i],2) + Math.pow(vEarth[1][i],2)));

      if(vEarth[0][i] <= 0 && vEarth[1][i] <= 0)
      {
	Xir = -Math.PI - Xir;
      }

      if(vEarth[0][i] <= 0 && vEarth[1][i] >= 0)
      {
	Xir = Math.PI - Xir;
      }

      var Gammar = Math.asin(-vEarth[2][i]/ VER);
      var GammaHis = [GammaHis, Gammar];
      XiHis = [XiHis, Xir];
      VAirRel = [VAirRel, VAR];
    }

    var MachHis = elemdiv(VAirRel, soundSpeedHis);
    var AlphaDegHis  = 57.2957795 * AlphaAR;
    var BetaDegHis = 57.2957795 * BetaAR;
    var qbarHis = 0.5 * elemmult(elemmult(airDensHis, VAirRel),VAirRel);   // ÜBERPRÜFEN!!!
    var GammaDegHis = 57.2957795 * GammaHis;
    var XiDegHis = 57.2957795 * XiHis;
 }
}
