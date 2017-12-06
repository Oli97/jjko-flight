/*globale Variablen*/
var Trimhist, x, u, V;
/*%   Supporting Calculations for Geometric, Inertial, and Aerodynamic
%   Properties of BizJet B
%   June 12, 2015
%	Copyright 2006-2015 by ROBERT F. STENGEL.  All rights reserved.

    clear
    disp('==============================================')
    disp('Geometric and Inertial Properties for BizJet B')
    disp('==============================================')
    Date    =   date
    */
//%   Comparable Bizjet Weight Distribution, lb, based on empty weight less engine weight.
//%   (from Stanford AA241 Notes, Ilan Kroo)
var    WingSys         =   1020;
var    TailSys         =   288;
var    BodySys         =   930;
var    GearSys         =   425;
var    NacelleSys      =   241;
var    PropSys         =   340;
var    ControlSys      =   196;
var    InstrSys        =   76;
var    HydrPneuSys     =   94;
var    ElecSys         =   361;
var    AvionSys        =   321;
var    FurnEquipSys    =   794;
var    ACSys           =   188;
var    AntiIceSys      =   101;
var    LoadHandSys     =   2;
var    EmptyStruc      =   (WingSys+TailSys+BodySys+GearSys+PropSys+ControlSys+InstrSys+HydrPneuSys+ElecSys+AvionSys+FurnEquipSys+ACSys+AntiIceSys+LoadHandSys);
                      //  % Less nacelles and engines
var    EngWgt          =   1002;
var    EmptyWgt        =   (EmptyStruc+EngWgt);    //% Less nacelles
var    EmptyStrucWgt   =   EmptyStruc/EmptyWgt;

var    WingRatio       =   (WingSys+GearSys+HydrPneuSys+AntiIceSys)/EmptyStruc;
var    HTRatio         =   0.75*(TailSys)/EmptyStruc;
var    VTRatio         =   0.25*(TailSys)/EmptyStruc;
var    FusRatio        =   (BodySys+PropSys+ControlSys+InstrSys+ElecSys+AvionSys+FurnEquipSys+ACSys+LoadHandSys)/EmptyStruc;
var    Total           =   (WingRatio+HTRatio+VTRatio+FusRatio);
var    NacelleRatio    =   (NacelleSys)/EngWgt;
  //                      % Related to engine weight rather than empty structure

//%	BizJet B Geometric Properties
//%   x measurements from nose along centerline, negative aft
//%   y & z measurements from centerline, positive right and down

var	S		=	19.51			//% Reference Area, m^2
var 	taperw	=	0.5333		//	% Wing Taper Ratio
var	cBar	=	1.56			//% Mean Aerodynamic Chord, m
var    sweep   =   11*0.01745329    //% Wing L.E. sweep angle, rad
var    xcp     =   -5.7473        // % Wing center of pressure, m
var    GamWing =   3*0.01745329 //    % Dihedral angle of the wing, rad

//%	BizJet B Mass and Inertial Properties
var	m		=	3000                    //% Total mass for simulation (USER-specified), kg
var    mEmpty  =   2522                   // % Gross empty mass, kg
var    mEng    =   (1+NacelleRatio)*240    //% Mass of engines + nacelles, kg
var    mStruc  =   mEmpty - mEng           //% Empty structural mass (less engines + nacelles), kg
var    mWing   =   WingRatio*mStruc        //% Wing mass, kg
var    mHT     =   HTRatio*mStruc         // % Horizontal tail mass, kg
var    mVT     =   VTRatio*mStruc         // % Vertical tail mass, kg
var    mFus    =   FusRatio*mStruc      //   % Empty fuselage mass, kg
var    mPay    =   0.5*(m - mEmpty)   //     % Payload mass, kg
var    mFuel   =   0.5*(m - mEmpty) //       % Fuel mass, kg

var    xcm     =   xcp - 0.45*cBar// % Center of mass from nose (USER-specified), m

var    lWing   =   xcm - xcp //  % Horizontal distance between c.m and wing c.p., m
var    zWing   =   -0.557  //    % Vertical distance between c.m and wing c.p., m
var    b		=	13.16     //  % Wing Span, m

var    lenFus  =   10.72          // % Fuselage length, m
var    xcpFus  =   -0.25*lenFus   // % Linear-regime fuselage center of pressure, m
var    xcpFusN =   -0.5*lenFus    // % Newtonian-regime fuselage center of pressure, m
var    lFus    =   xcm - xcpFus;  // % Linear fuselage lift cp offset,m
var    lFusN   =   xcm - xcpFusN; // % Newtonian fuselage lift cp offset, m

var    dFus    =   1.555             // % Fuselage diameter, m
var    Sfus    =   (pi/4)*lenFus*dFus //% Plan or side area of fuselage, m
var    Sbase   =   (pi/4)*dFus^2     // % Fuselage cross-sectional area, m^2

var    bHT     =   5.3;            //% Horizontal tail span, m
var    cHT     =   1.1;            //% Mean horizontal tail chord, m
var    swpHT   =   38*0.0174533;   //% Horizontal tail sweep, rad
var    SHT     =   bHT*cHT         //% Horizontal tail area, m^2
var    xHT     =   -11.3426;       //% Linear xcp of horizontal tail
var    lHT     =   xcm - xHT;     // % Horizontal tail length, m
var    zHT     =   1.5;         //   % zcp of horizontal tail, m

var    xVT     =   -10.044;   //     % Linear xcp of vertical tail, m
var 	lVT		=	xcm - xVT;     // % Vertical tail length, m
var    bVT     =   2.409;         // % Vertical tail span, m
var    cVT     =   1.88;          // % Mean vertical tail chord, m
var    swpVT   =   50*0.0174533;  // % Vertical tail sweep, rad
var    SVT     =   bVT*cVT;       // % Vertical tail area, m^2
var    zVT     =   1.5;           // % zcp of vertical tail, m

var    xEng    =   -7.735;        // % xcm of engine, m
var    lEng    =   xcm - xEng;    // % Engine length, m
var    yEng    =   1.1325;        // % ycm of engine, m
var    zEng    =   0.4038;        // % zcm of engine, m

var    xNac        =   -7.7252;           // % xcp of engine, m
var    lNac        =   xcm - xNac;        // % Nacelle length, m
var    bNac        =   2.5;               // % Nacelle span, m
var    cNac        =   1.82;              // % Nacellechord, m
var    dNac        =   0.73;              // % Nacelle diameter, m
var    SbaseNac    =   0.25*pi*dNac^2;    // % Nacelle base area, m^2
var    Snac        =   bNac*cNac;         // % Nacelle plan area, m^2

var    xVent   =   -9.94;         // % xcp of ventral fin, m
var    lVent   =   xcm - xVent;   // % Ventral fin length, m
var    zVent   =   0;             // % zcp of ventral fin, m
var    bVent   =   1;             // % Ventral fin span, m
var    cVent   =   0.85;          // % Ventral fin chord, m
var    Svent   =   bVent*cVent;   // % Ventral fin area, m^2
var    swpVent =   60*0.0174533;  // % Ventral sweep angle, rad

var    Splan   =   S + Sfus + Snac + SHT + Svent //  % Plan area of airplane
var    Swet    =   2*(Splan + Sfus + SVT)        //  % Wetted area of airplane


var    ARwing  =   (b^2) / S         //  % Wing aspect ratio
var    ARHT    =   (bHT^2) / SHT     //  % Horizontal tail aspect ratio
var    ARnac   =   (bNac^2) / Snac   //  % Engine nacelle aspect ratio
var    ARvent  =   (bVent^2)/ Svent  //  % Ventral fin aspect ratio
var    ARVT    =   (bVT^2) / SVT     //  % Vertical tail aspect ratio

//%   Moments and Product of Inertia
var    Ixx     =   (1/12)*((mWing+mFuel)*b^2 + mHT*bHT^2 + mVT*bVT^2) + (0.25*(mFus+mPay)*dFus^2 + mEng*yEng^2 + mVT*zVT^2)
var    Iyy     =   (1/12)*((mFus+mPay)*lenFus^2 + (mWing+mFuel)*cBar^2 + mVT*cHT^2) + (mEng*lEng^2 + mHT*lHT^2 + mVT*lHT^2)
var    Izz     =   (1/12)*((mFus+mPay)*lenFus^2 + (mWing+mFuel)*b^2 + mHT*bHT^2) + (mEng*lEng^2 + mHT*lHT^2 + mVT*lVT^2)
var    Ixz     =   mHT*lHT*zHT + mVT*lVT*zVT + mEng*lEng*zEng

var    dEmax   =   20 * 0.01745329 //%   Maximum Elevator Deflection is ±20 deg
var	dAmax   =   35 * 0.01745329 //%   Maximum Aileron Deflection is ±35 deg
var	dRmax   =   35 * 0.01745329 //%   Maximum Rudder Deflection is ±35 deg


//%   BizJet B Aero Properties
var    AlphaTable	=	[-10 -8 -6 -4 -2 0 2 4 6 8 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 30 35 40 45 50 55 60 65 70 75 80 85 90];
var    Points      =   length(AlphaTable);
var    AlphaRad    =   0.0174533*AlphaTable;
var    SinAlpha    =   Math.sin(AlphaRad);
var    CosAlpha    =   Math.cos(AlphaRad);

  //  %   Newtonian Coefficients
var    CN          =   2*(Splan/S)*SinAlpha.*SinAlpha;
    CN          =   CN.*sign(AlphaRad);
var    CDNewt      =   2*(Splan/S)*Math.abs(SinAlpha.*SinAlpha.*SinAlpha);
var    CLNewt      =   CN.*CosAlpha;
var    cpNewt      =   (S*lWing + Sfus*lFusN + SHT*lHT + Svent*lVent + Snac*lNac)/ (S + Sfus + SHT + Svent + Snac);

//%   Longitudinal Aerodynamics
//%   =========================
//%   Lift
var    CLaWing =   pi*ARwing / (1 + Math.sqrt(1 + (0.5*ARwing/Math.cos(sweep))^2));
var    deda    =   0; %   T-tail
var    CLaHT   =   (1 - deda)*(pi*ARHT / (1 + Math.sqrt(1 + (0.5*ARHT/Math.cos(swpHT))^2)))*SHT / S;
var    CLaFus  =   2*Sbase / S;
var    CLaNac  =   2*SbaseNac / S;
var    CLaVent =   (pi*ARvent / (1 + Math.sqrt(1 + (0.5*ARvent/Math.cos(swpVent))^2)))*Svent / S;
var    CLaTot  =   CLaWing + CLaHT + CLaFus + CLaNac + CLaVent;
  //  %   Assume CL is linear to Alpha = 10 deg = 0.1745 rad
var    CL10    =   CLaTot*0.174533;
  //  %   Assume CL is symmetrically quadratic about CLmax
var    CLmax   =   1.35;
var    delAlph =   4*0.0174533;   // % Stall occurs at 14 deg

var    CLstatic        =   zeros(1,39);
    CLstatic(1:11)  =   CLaTot*AlphaRad(1:11);
var    kStall          =   (CLmax - CL10)/delAlph^2;
    CLstatic(12:21) =   CLmax - kStall*(AlphaRad[15] - AlphaRad(12:21)).*(AlphaRad[15] - AlphaRad(12:21));
    CLstatic[22]    =   0.73;
    CLstatic[23]    =   0.73;
    CLstatic[24]    =   0.74;
    CLstatic[25]    =   0.76;
    CLstatic[26]    =   0.78;
    CLstatic(27:39) =   CN(27:39).*CosAlpha(27:39);
var    CLTable         =   CLstatic;


//%   Drag
var    Re      =   1.225*100*lenFus / 1.725e-5;    %   Reynolds number at 100 m/s
var    Cf      =   0.46*(Math.log10(Re))^(-2.58);       %   Flat plate friction coefficient
var    CDf     =   Cf*Swet / S;
var    CDbase  =   0.12*Sbase / S;                 %   Base pressure drag
var    CDo     =   CDf + CDbase;
var    OEF     =   1.78*(1 - 0.045*ARwing^0.68) - 0.64;
                //% Oswald Efficiency factor, Raymer (straight wing)
var    CDstatic        =   zeros(1,39);
    CDstatic(1:10)  =   CDo + CLstatic(1:10).*CLstatic(1:10) / (OEF*pi*ARwing);
    CDstatic[11]    =   CDstatic[10]*1.15;
    CDstatic[12]    =   CDstatic[11]*1.15;
    CDstatic[13]    =   CDstatic[12]*1.15;
    CDstatic[14]    =   CDstatic[13]*1.15;
    CDstatic[15]    =   CDstatic[14]*1.1;
    CDstatic[16]    =   CDstatic[15]*1.1;
    CDstatic[17]    =   CDstatic[16]*1.1;
    CDstatic[18]    =   CDstatic[17]*1.1;
    CDstatic[19]    =   CDstatic[18]*1.1;
    CDstatic[20]    =   CDstatic[19]*1.1;
    CDstatic[21]    =   CDstatic[20]*1.1;
    CDstatic[22]    =   CDstatic[21]*1.1;
    CDstatic[23]    =   CDstatic[22]*1.1;
    CDstatic[24]    =   CDstatic[23]*1.1;
    CDstatic[25]    =   CDstatic[24]*1.1;
    CDstatic[26]    =   CDstatic[25]*1.1;
    CDstatic(26:39) =   2*(Splan/S)*Math.abs(SinAlpha(26:39).*SinAlpha(26:39).*SinAlpha(26:39));
var    CDNewt          =   2*(Splan/S)*Math.abs(SinAlpha.*SinAlpha.*SinAlpha);
var    CDTable         =   CDstatic;

//%   Pitching Moment (c.m. @ wing c.p.)
var    CmStatic        =   zeros(1,39);
var    SM              =   (CLaWing*lWing + CLaFus*lFus + CLaHT*lHT + CLaNac*lNac + CLaVent*lVent) /(cBar*(CLaWing + CLaFus + CLaHT + CLaNac + CLaVent));
                        //% Static Margin at Alpha = 0 deg
    CmStatic[1:21]  =   -(CLstatic[1:21].*Math.cos(AlphaRad[1:21]) + CDstatic[1:21].*Math.sin(AlphaRad[1:21]))*SM;
    CmStatic[22]    =   CmStatic[21]*0.9;
    CmStatic[23]    =   CmStatic[22]*0.9;
    CmStatic[24]    =   CmStatic[23]*0.9;
    CmStatic[25]    =   CmStatic[24]*0.9;
    CmStatic[26]    =   0.9*CmStatic[25] - 0.1*CN[26].*Math.sign(AlphaRad[26])*cpNewt/cBar;
    CmStatic[27]    =   0.4*CmStatic[26] - 0.6*CN[27].*Math.sign(AlphaRad[27])*cpNewt/cBar;
    CmStatic[28]    =   0.1*CmStatic[27] - 0.9*CN[28].*Math.sign(AlphaRad[28])*cpNewt/cBar;
var    CN              =   2*(Splan/S)*SinAlpha.*SinAlpha;
    CmStatic[29:39] =   -CN[29:39]*cpNewt/cBar;
var    CmTable         =   CmStatic;



//%   CmdETable & CLdEo, Elevator Effect
var    tauDE           =   0.32;  // % Geometric elevator chord/horizontal tail chord
var    tauCO           =   0.68;  // % Elevator Carryover effect
var    Sigmoid         =   [];
    Sigmoid[1:5]    =   tauCO;
    Sigmoid[6:39]   =   tauDE + (tauCO - tauDE) ./ (1 + Math.exp(-15.*(AlphaRad[26] - AlphaRad[6:39])));
var    CmdETable       =   zeros(1,39);
var    CLdEo           =   tauCO*CLaHT;              //  % CLdE at Alpha = 0
var    CmdEo           =   -tauCO*(lHT/cBar)*CLaHT;   // % CmdE at Alpha = 0
    CmdETable[1:39] =   (CmdEo*Sigmoid).*Math.cos(AlphaRad[1:39]);
                                      // % Elevator effect on moment, per rad
    CLdETable[1:39] =   (CLdEo*Sigmoid).*Math.cos(AlphaRad[1:39]);
                                       //% Elevator effect on lift, per rad

//%   Longitudinal Rotary & Unsteady Derivatives
var    CLqHat      =   2*CLaHT*lHT/cBar;
var    CmqHat      =   -CLqHat*lHT/cBar;

//%   Lateral-Directional Aerodynamics
//%   ================================
//%   CYBetaTable, Side Force Sensitivity to Sideslip Angle
var    EndPlate    =   1.1;   // % End-plate effect of T-tail
var    CYBetaVT    =   -EndPlate*(pi*ARVT / (1 + Math.sqrt(1 + (0.5*ARVT/Math.cos(swpVT))^2)))*SVT / S;
var    CYBetaFus   =   -2*Sbase / S;
var    CDoWing     =   0.005;
var    CYBetaWing  =   -CDoWing - (GamWing^2)*pi*ARwing / (1 + Math.sqrt(1 + ARwing^2));
  var  CYBetaVent  =   -0.4*CLaVent;
var    CYBetao     =   CYBetaVT + CYBetaFus + CYBetaWing + CYBetaVent;
var    CYBetaTable =   CYBetao*Math.cos(AlphaRad);

//%   ClBetaTable, Roll Moment Sensitivity to Sideslip Angle
var    ClBetaWing  =   -((1 + 2*taperw)/(6*(1 + taperw)))*(GamWing*CLaWing + (CLTable.*Math.tan(sweep)));
var    ClBetaWF    =   1.2*Math.sqrt(ARwing)*(2*zWing*dFus/b^2);
var    ClBetaVT    =   -zVT*CYBetaVT/b;
var    ClBetao     =   (ClBetaWing + ClBetaWF + ClBetaVT);
var    ClBetaTable =   ClBetao.*Math.cos(AlphaRad);


//%   CnBetaTable, Yaw Moment Sensitivity to Sideslip Angle
var    CnBetaWing  =   0.075*CLTable*GamWing;
var    CnBetaFus   =   CLaFus*lFusN / b;
var    CnBetaVT    =   -CYBetaVT*lVT / b;
var    CnBetaVent  =   -CYBetaVent*lVent / b;
var    CnBetaTable =   (CnBetaWing + CnBetaFus + CnBetaVT).*Math.cos(AlphaRad);

//%   CldATable, Roll Moment Sensitivity to Aileron Deflection
var    tauDA       =   0.25;
var    kDA         =   0.38;
var    CldAo       =   tauDA*(CLaWing/(1 + taperw))*((1 - kDA^2)/3 - (1 - kDA^3)*(1 - taperw)/3);
var    CldATable   =   CldAo*Math.cos(AlphaRad);

var    CYdAo       =   0;  //%   Side force due to aileron, rad

//%   CndATable, Yaw Moment Sensitivity to Aileron Deflection
//%   Cessna 510 has an Aileron-Rudder Interconnect; assume CndA = 0
var    CndATable   =   zeros(1,39);

//   CldRTable, Roll Moment Sensitivity to Rudder Deflection
var    tauDR       =   0.5;   // % Geometric Rudder chord/horizontal tail chord
var    tauCOR      =   0.8;   // % Rudder Carryover effect
var    CldRo       =   tauCOR*zVT*CYBetaVT / b;
var    CldRTable   =   CldRo*Math.cos(AlphaRad);

//   CndRTable, Yaw Moment Sensitivity to Rudder Deflection
var    CndRo       =   -tauCOR*CnBetaVT;
var    CndRTable   =   CndRo*Math.cos(AlphaRad);

//   Lateral-Directional Rotary & Unsteady Derivatives
var    CYrHat      =   -2*CYBetaVT*lVT/b;
  var  ClpHato     =   -(CLaWing + CLaHT*(SHT/S) - CYBetaVT*(SVT/S))*((1 + 3*taperw)/(1 + taperw))/12;
var    ClpHatTable =   ClpHato*Math.cos(AlphaRad);
var    ClrHato     =   -(CLaWing + CLaHT*(SHT/S) - CYBetaVT*(SVT/S))*(1 + 3 * taperw)/(12 * (1 + taperw));
var    ClrHatTable =   ClrHato*Math.cos(AlphaRad);
var    CnpHatTable =   (- CLTable*((1 + 3*taperw)/(1 + taperw))/12) .* Math.cos(AlphaRad);
var    CnrHatVT    =   -2*CnBetaVT*(lVT/b);
var    CnrHatWing  =   -0.103*CLTable.*CLTable - 0.4*CDoWing; // %   (from Seckel)
var    CnrHato     =   CnrHatVT + CnrHatWing;
var    CnrHatTable =   CnrHato .* Math.cos(AlphaRad);


Windfield=function(height,phir,thetar,psir){
/*%	FLIGHT Wind Field Interpolation for 3-D wind as a Function of Altitude
%	June 12, 2015
%	===============================================================
%	Copyright 2006-2015 by ROBERT F. STENGEL.  All rights reserved.
*/
	var windh	=	[-10 0 100 200 500 1000 2000 4000 8000 16000];	//% Height, m
	var windx	=	[0 0 0 0 0 0 0 0 0 0];	//% Northerly wind, m/s
	var windy	=	[0 0 0 0 0 0 0 0 0 0];	//% Easterly wind, m/s
	var windz	=	[0 0 0 0 0 0 0 0 0 0];	//% Vertical wind. m/s

	var winde	=	[interp1(windh,windx,height);interp1(windh,windy,height);interp1(windh,windz,height)];	//% Earth-relative frame
	var HEB		=	DCM(phir,thetar,psir);
	var windb	=	HEB * winde;					//% Body-axis frame

  return windb;
}

LinModel=function(tj,xj){
/*%	FLIGHT Equations of Motion for Linear Model (Jacobian) Evaluation,
%	with dummy state elements added for controls

%	June 12, 2015
%	===============================================================
%	Copyright 2006 by ROBERT F. STENGEL.  All rights reserved.
*/
	var x		=	xj[1:12];
	u		=	xj[13:19];

	var xdot	=	EoM(tj,x);
	var xdotj	=	[xdot;0;0;0;0;0;0;0];

  return xdotj;
}

Atmos = function(geomAlt){
/*	1976 U.S. Standard Atmosphere Interpolation for FLIGHT

%	June 12, 2015
%	===============================================================
%	Copyright 2006-15 by ROBERT F. STENGEL.  All rights reserved.

%	Note:	Function does not extrapolate outside altitude range
%	Input:	Geometric Altitude, m (positive up)
%	Output:	Air Density, kg/m^3
%			Air Pressure, N/m^2
%			Air Temperature, K
%			Speed of Sound, m/s

%	Values Tabulated by Geometric Altitude*/
	var Z	= [-1000,0,2500,5000,10000,11100,15000,20000,47400,51000];
	var H	= [-1000,0,2499,4996,9984,11081,14965,19937,47049,50594];
	var ppo	= [1,1,0.737,0.533,0.262,0.221,0.12,0.055,0.0011,0.0007];
	var rro	= [1,1,0.781,0.601,0.338,0.293,0.159,0.073,0.0011,0.0007];
	var T	= [288.15,288.15,271.906,255.676,223.252,216.65,216.65,216.65,270.65,270.65];
	var a	= [340.294,340.294,330.563,320.545,299.532,295.069,295.069,295.069,329.799,329.799];
	var R		= 6367435;	// Mean radius of the earth, m
	var Dens	= 1.225;	// Air density at sea level, Kg/m^3
	var Pres	= 101300;	// Air pressure at sea level, N/m^2

//	Geopotential Altitude, m
	var geopAlt	=	R * geomAlt / (R + geomAlt);

//	Linear Interpolation in Geopotential Altitude
//	for Temperature and Speed of Sound
	var temp		=	interp1(Z,T,geopAlt);
	var soundSpeed	=	interp1(Z,a,geopAlt);

//	Exponential Interpolation in Geometric Altitude
//for Air Density and Pressure
var betap, betar, airDens, airPres;
	for(k=2;k<=10;k++){
		if (geomAlt <= Z[k]){
			betap	=	Math.log(ppo[k] / ppo[k-1]) / (Z[k] - Z[k-1]);
			betar	=	Math.log(rro[k] / rro[k-1]) / (Z[k] - Z[k-1]);
			airPres	=	Pres * ppo[k-1] * Math.exp(betap * (geomAlt - Z[k-1]));
			airDens	=	Dens * rro[k-1] * Math.exp(betar * (geomAlt - Z[k-1]));
			break;
		}
  }
  var atmos = [airDens,airPres,temp,soundSpeed];
  return atmos;

}

  TrimCost=function(OptParam){
/*	FLIGHT Cost Function for Longitudinal Trim in Steady Level Flight
%	June 12, 2015
%	===============================================================
%	Copyright 2006-2015 by ROBERT F. STENGEL.  All rights reserved.
*/
var R;
	R[1][1]	=1;
  R[1][2]	=0;
  R[1][3]	=0;
  R[2][1]	=0;
  R[2][2]	=1;
  R[2][3]	=0;
  R[3][1]	=0;
  R[3][2]	=0;
  R[3][3]	=1;

/*
% Optimization Vector:
%	1 = Stabilator, rad
%	2 = Throttle, %
%	3 = Pitch Angle, rad
*/

	u	=	[u[1];u[2];u[3];OptParam[2];u[5];u[6];OptParam[1]];

	x	=	[V * Math.cos(OptParam[3]); x[2]; V * Math.sin(OptParam[3]); x[4]; x[5];x[6];x[7];x[8];x[9];x[10];OptParam[3];x[12]];

	var xdot	=	EoM(1,x);
	var xCost	=	[xdot[1] xdot[3] xdot[8]];
	var TrimCost		=	math.transpose(xCost) * R * xCost;
	var ParamCost	=	[OptParam;J];
	TrimHist	=	[TrimHist ParamCost];

  return TrimCost;
}

DCM = function(Phi,Theta,Psi)
{
  var sinR = Math.sin(Phi);
  var cosR = Math.cos(Phi);
  var sinP = Math.sin(Theta);
  var cosP = Math.cos(Theta);
  var sinY = Math.sin(Psi);
  var cosY = Math.cos(Psi);
  var H;

  H[1][1] = cosP * cosY;
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
  var value       =   x[6];
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
