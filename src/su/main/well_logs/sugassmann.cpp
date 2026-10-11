/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUGASSMAN: $Revision: 1.3 $ ; $Date: 2011/11/12 00:40:42 $	*/

#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									",
" SUGASSMAN - Model reflectivity change with rock/fluid properties	",
"									",
"	sugassman [optional parameters] > data.su			",
"									",
" Optional parameters:							",
" nt=500 	number of time samples					",
" ntr=200	number of traces					",
" dt=0.004 	time sampling interval in seconds			",
" mode=0	model isolated gassmann refl coefficient		",
"		=1 embed gassmann RC in random RC series		",
"		=2 R0 parameter sensitivity output			",
" p=0.15 	parameter sensitivity test range (if mode=2)		",
" .... Environment variables ...					",
" temp=140 	Temperature in degrees C				",
" pres=20 	Pressure in megaPascals					",
" .... Caprock variables ....						",
" v1=37900 	caprock P-wave speed (m/s)				",
" r1=44300 	caprock mass density (g/cc)				",
" .... Reservoir fluid variables ....					",
" g=0.56	Gas specific gravity 0.56 (methane)-1.8 (condensate)	",
" api=50 	Gas specific gravity 10 (heavy)-50 (ultra light)	",
" s=35		Brine salinity in ppm/(1000 000				",
" so=.7 	Oil saturation (0-1)					",
" sg=.2 	Gas saturation (0-1)					",
" .... Reservoir rock frame variables ....				",
" kmin=37900 	Bulk modulus (MPa) of rock frame mineral(s) [default=quartz]",
" mumin=44300 	Shear modulus (MPa) of rock frame mineral(s) [default=quartz]",
" rmin=2.67 	Mass density (g/cc) of rock frame mineral(s) [default=quartz]",
" phi=.24 	Rock frame porosity (0-1)				",
" a=1 		Fitting parameters: Mdry/Mmineral ~ 1/(a + b phi^c)	",
" b=15 		... where M is P-wave modulus and defaults are for	",
" c=1 		... Glenn sandstone [see Liner (2nd ed, table 26.2)]	",
"	h=20 			Reservoir thickness (m)			",
"									",
" Notes:								",
" Creates a reflection coefficient series based on Gassmann		",
" theory of velocity and density for porous elastic media		",
NULL};
#endif


/* 
 * Credits: UHouston: Chris Liner	9/23/2009
 *
 * trace header fields set: 
 */
/**************** end self doc ***********************************/

/* Library version of SUGASSMANN: the rock and fluid physics of the program (Gassmann's fluid substitution of the rock, with the
 * Batzle and Wang properties of brine, oil and gas); its main program, which makes traces from them, is left to the caller.
 *
 * Differences from the program:
 *  - the functions are prefixed gassmann_;
 *  - a mixture with so + sg > 1 gives NaN where the program ended the process.
 * (The program's modes 1 and 2 could not be used: mode 1 does nothing, and mode 2 was changed to 0 by the program.)
 */
/* Function prototypes */
float gassmann_rockT0(float h, float kmin, float mumin, float rmin, 
		float phi, float a, float b, float c, float api, 
		float g, float s, float so, float sg, float t, float p);
float gassmann_rockR0(float v1, float r1, float kmin, float mumin, 
		float rmin, float phi, float a, float b, float c, 
		float api, float g, float s, float so, float sg, 
		float t, float p);
float gassmann_rockIp(float kmin, float mumin, float rmin, float phi, float a,
		float b, float c, float api, float g, float s, 
		float so, float sg, float t, float p);
float gassmann_rockVp(float kmin, float mumin, float rmin, float phi, 
		float a, float b, float c, float api, float g, 
		float s, float so, float sg, float t, float p);
float gassmann_rockDensity(float rmin, float phi, float api, float g, 
			float s, float so, float sg, float t, float p);
float gassmann_fluidVelocity(float api, float g, float s, float so, 
			float sg, float t, float p);
float gassmann_fluidModulus(float api, float g, float s, float so, 
			float sg, float t, float p);
float gassmann_fluidDensity(float api, float g, float s, float so, float sg, 
			float t, float p);
float gassmann_brineVelocity(float s, float t, float p);
float gassmann_brineModulus(float s, float t, float p);
float gassmann_brineDensity(float s, float t, float p);
float gassmann_oilModulus(float api, float g, float t, float p);
float gassmann_oilVelocity(float api, float g, float t, float p);
float gassmann_oilDensity(float api, float g, float t, float p);
float gassmann_gasVelocity(float g, float t, float p);
float gassmann_gasModulus(float g, float t, float p);
float gassmann_gasDensity(float g, float t, float p);

/* C language definitions for use with Mathematica output */
#define Power(x, y)	(pow((double)(x), (double)(y)))
#define Sqrt(x)	 (sqrt((double)(x)))
#define Abs(x)		(fabs((double)(x)))
#define Exp(x)		(exp((double)(x)))
#define Log(x)		(log((double)(x)))
#define Sin(x)		(sin((double)(x)))
#define Cos(x)		(cos((double)(x)))
#define tan(x)		(tan((double)(x)))
#define ArcSin(x)	(asin((double)(x)))
#define ArcCos(x)	(acos((double)(x)))
#define Arctan(x)	(atan((double)(x)))
#define Sinh(x)		(sinh((double)(x)))
#define Cosh(x)		(cosh((double)(x)))
#define tanh(x)		(tanh((double)(x)))
#define E		2.71828182845904523536029
#define pi		3.14159265358979323846264
#define Degree		0.01745329251994329576924

float gassmann_rockT0(float h, float kmin, float mumin, float rmin, 
		float phi, float a, float b, float c, 
		float api, float g, float s, float so, 
		float sg, float t, float p)
/*******************************************************************
rockTO - Calculate P-wave normal incidence reflection coefficient of 
	 a fluid saturated reservoir rock. 
*******************************************************************
Notes:
Caprock velocity is v1 (m/s) and density is r1 (g/cc).
Constituent mineral properties are bulk modulus kmin (MPa),
shear modulus \[Mu]min (MPa), and density rmin (g/cc).	Porosity 
is phi (Range: 0-1) and (a,b,c) are empirical frame stiffness 
constants used as mdry/mmin=(a+b phi^c). The oil component is characterized by oil gravity 
API, gas by specific gravity G, and brine by weight fraction of sodium chloride S (ppm/1,000,000).	Oil 
saturation is So (Range: 0-1) and gas saturation is Sg (Range: 0-1).	Water saturation is calculated 
from 1-So-Sg. T is temperature (deg C) and P is pressure (MPa)
*******************************************************************
Credits: UHouston: Chris Liner	9/23/2009
*******************************************************************/
{
	float v;
	float x;

	v = gassmann_rockVp(kmin, mumin, rmin, phi, a, b, c, api, g, s, so, sg, t, p);
	x = 2 * h / v;

	return x;
}

float gassmann_rockR0(float v1, float r1, float kmin, float mumin,
		float rmin, float phi, float a, float b, float c, 
		float api, float g, float s, float so, 
		float sg, float t, float p)
/*******************************************************************
rockRO - Calculate P-wave normal incidence reflection coefficient of a fluid saturated reservoir rock. 
*******************************************************************
Notes:
Caprock velocity is v1 (m/s) and density is r1 (g/cc).
Constituent mineral properties are bulk modulus kmin (MPa),
shear modulus \[Mu]min (MPa), and density rmin (g/cc).	Porosity 
is phi (Range: 0-1) and (a,b,c) are empirical frame stiffness 
constants used as mdry/mmin=(a+b phi^c). The oil component is characterized by oil gravity 
API, gas by specific gravity G, and brine by weight fraction of sodium chloride S (ppm/1,000,000).	Oil 
saturation is So (Range: 0-1) and gas saturation is Sg (Range: 0-1).	Water saturation is calculated 
from 1-So-Sg. T is temperature (deg C) and P is pressure (MPa)
*******************************************************************
Credits: UHouston: Chris Liner	9/23/2009
*******************************************************************/
{
	float ip1, ip2;
	float x;

	ip1 = v1 * r1;
	ip2 = gassmann_rockIp(kmin, mumin, rmin, phi, a, b, c, api, g, s, so, sg, t, p);
	x = (ip2 - ip1)/(ip2 + ip1);

	return x;
}


float gassmann_rockIp(float kmin, float mumin, float rmin, float phi,
		float a, float b, float c, float api, float g,
		float s, float so, float sg, float t, float p)
/*******************************************************************
gassmann_rockIp - Calculate P-wave acoustic impedance (m/s)*(g/cc) of a fluid 
	 saturated reservoir rock. 
********************************************************************
Notes:
Constituent mineral properties are bulk modulus kmin (MPa),
shear modulus \[Mu]min (MPa), and density rmin (g/cc).	Porosity 
is phi (Range: 0-1) and (a,b,c) are empirical frame stiffness 
constants used as mdry/mmin=(a+b phi^c). The oil component is characterized by oil gravity 
API, gas by specific gravity G, and brine by weight fraction of sodium chloride S (ppm/1,000,000).	Oil 
saturation is So (Range: 0-1) and gas saturation is Sg (Range: 0-1).	Water saturation is calculated 
from 1-So-Sg. T is temperature (deg C) and P is pressure (MPa)
*******************************************************************
Credits: UHouston: Chris Liner	9/23/2009
*******************************************************************/
{
	float rsat, vsat;
	float x;

	rsat = gassmann_rockDensity(rmin, phi, api, g, s, so, sg, t, p);
	vsat = gassmann_rockVp(kmin, mumin, rmin, phi, a, b, c, api, g, s, so, sg, t, p);
	x = rsat * vsat;

	return x;
}

float gassmann_rockVp(float kmin, float mumin, float rmin, float phi,
		float a, float b, float c, float api, float g,
		float s, float so, float sg, float t, float p)
/*******************************************************************
gassmann_rockVp - Calculate P-wave speed (m/s) of a fluid saturated reservoir rock.	
********************************************************************
Notes:
Constituent mineral properties are bulk modulus kmin (MPa), 
shear modulus \[Mu]min (MPa), and density rmin (g/cc).	Porosity 
is phi (Range: 0-1) and (a,b,c) are empirical frame stiffness 
constants used as mdry/mmin=(a+b phi^c). The oil component is characterized by oil gravity 
API, gas by specific gravity G, and brine by weight fraction of sodium chloride S (ppm/1,000,000).	Oil 
saturation is So (Range: 0-1) and gas saturation is Sg (Range: 0-1).	Water saturation is calculated 
from 1-So-Sg. T is temperature (deg C) and P is pressure (MPa)
*******************************************************************
Credits: UHouston: Chris Liner	9/23/2009
*******************************************************************/
{
	float kf, mmin, mdry, rhs, msat, rsat;
	float x;

	kf	= gassmann_fluidModulus(api, g, s, so, sg, t, p);
	mmin = kmin + 4 * mumin / 3;
	mdry = mmin /(a + b * Power(phi,c));
	rhs	= mdry / (mmin - mdry) + kf / ( phi * (mmin - kf) );
	msat = mmin * rhs / (1 + rhs);
	rsat = gassmann_rockDensity(rmin, phi, api, g, s, so, sg, t, p);
	x = Sqrt(1000 * msat / rsat);

	return x;
}

float gassmann_rockDensity(float rmin, float phi, float api, float g, float s, 
			float so, float sg, float t, float p)
/*******************************************************************
gassmann_rockDensity - Calculate mass density (g/cc) of a fluid 
		saturated reservoir rock.	\[Continuation]
********************************************************************
Notes:
The mineral density of the rock is rmin, and the porosity \[Continuation]
is phi (Range: 0-1).	The oil component is characterized by oil gravity 
API, gas by specific gravity G, and brine by weight 
fraction of sodium chloride S (ppm/1,000,000).	Oil 
saturation is So (Range: 0-1) and gas saturation is 
Sg (Range: 0-1).	Water saturation is calculated 
from 1-So-Sg. T is temperature (deg C) and P is pressure (MPa)
*******************************************************************
Credits: UHouston: Chris Liner	9/23/2009
*******************************************************************/
{
	float x1, x;

	x1 = gassmann_fluidDensity(api, g, s, so, sg, t, p);
	x = rmin * (1 - phi) + phi * x1;

	return x;
}

float gassmann_fluidVelocity(float api, float g, float s, float so,
			float sg, float t, float p)
/*******************************************************************
gassmann_fluidVelocity - Calculate sound speed (m/s) for a reservoir fluid mix. 
********************************************************************
The oil component is characterized by oil gravity 
API, gas by specific gravity G, and brine by weight 
fraction of sodium chloride S (ppm/1,000,000).	Oil 
saturation is So (Range: 0-1) and gas saturation is 
Sg (Range: 0-1).	Water saturation is calculated 
from 1-So-Sg. T is temperature (deg C) and P is pressure (MPa)
*******************************************************************
Credits: UHouston: Chris Liner	9/23/2009
*******************************************************************/
{
	float x1, x2, x;

	x1 = gassmann_fluidModulus(api, g, s, so, sg, t, p);
	x2 = gassmann_fluidDensity(api, g, s, so, sg, t, p);
	x = Sqrt( 1000 * x1 / x2 );

	return x;
}

float gassmann_fluidModulus(float api, float g, float s, float so, 
			float sg, float t, float p)
/*******************************************************************
gassmann_fluidModulus - Calculate bulk modulus (MPa) of a reservoir fluid mix. 
********************************************************************
The oil component is characterized by oil gravity 
API, gas by specific gravity G, and brine by weight 
fraction of sodium chloride S (ppm/1,000,000).	Oil 
saturation is So (Range: 0-1) and gas saturation is 
Sg (Range: 0-1).	Water saturation is calculated 
from 1-So-Sg. T is temperature (deg C) and P is pressure (MPa)
*******************************************************************
Credits: UHouston: Chris Liner	9/23/2009
*******************************************************************/
{
	float brine, gas, oil;
	float x;

	if ( so + sg > 1 ) return NAN;

	brine = gassmann_brineModulus(s,t,p);
	gas	= gassmann_gasModulus(g,t,p);
	oil	= gassmann_oilModulus(api,g,t,p);

	x = 1.0 / ( (1 - so - sg)/brine + so/oil + sg/gas );

	return x;
}

float gassmann_fluidDensity(float api, float g, float s, float so,
			float sg, float t, float p)
/*******************************************************************
gassmann_fluidDensity - Calculate mass density (g/cc) of a reservoir fluid mix. 
********************************************************************
Notes:
The oil component is characterized by oil gravity 
API, gas by specific gravity G, and brine by weight 
fraction of sodium chloride S (ppm/1,000,000).	Oil 
saturation is So (Range: 0-1) and gas saturation is 
Sg (Range: 0-1).	Water saturation is calculated 
from 1-So-Sg. T is temperature (deg C) and P is pressure (MPa)
*******************************************************************
Credits: UHouston: Chris Liner	9/23/2009
*******************************************************************/
{
	float brine, gas, oil;
	float x;

	if ( so + sg > 1 ) return NAN;

	brine = gassmann_brineDensity(s,t,p);
	gas	= gassmann_gasDensity(g,t,p);
	oil	= gassmann_oilDensity(api,g,t,p);

	x = (1 - so - sg)*brine + so*oil + sg*gas;

	return x;
}

float gassmann_brineModulus(float s, float t, float p)
/*******************************************************************
gassmann_brineModulus - Calculate bulk modulus (MPa) of saline solution characterized 
	  by weight fraction of sodium chloride S (ppm/1,000,000) as 
	  a function of temperature (deg C) and pressure (MPa)
*******************************************************************
Credits: UHouston: Chris Liner	9/23/2009
*******************************************************************/
{
	float x;

	x = (Power(1402.85 + 1.524*p + 0.003437*Power(p,2) - 
	0.000011970000000000002*Power(p,3) + 
	(780 - 10*p + 0.16*Power(p,2))*Power(s,1.5) - 
	1820*Power(s,2) + 4.871*t - 0.0111*p*t + 
	0.00017390000000000003*Power(p,2)*t - 
	1.628e-6*Power(p,3)*t - 0.04783*Power(t,2) + 
	0.0002747*p*Power(t,2) - 
	2.135e-6*Power(p,2)*Power(t,2) + 
	1.2370000000000001e-8*Power(p,3)*Power(t,2) + 
	0.0001487*Power(t,3) - 
	6.502999999999999e-7*p*Power(t,3) - 
	1.4550000000000001e-8*Power(p,2)*Power(t,3) + 
	1.327e-10*Power(p,3)*Power(t,3) - 
	2.197e-7*Power(t,4) + 7.987e-10*p*Power(t,4) + 
	5.23e-11*Power(p,2)*Power(t,4) - 
	4.614e-13*Power(p,3)*Power(t,4) + 
	s*(1170 + 2.6*p - 0.0476*Power(p,2) - 9.6*t - 
		0.0029*p*t + 0.055*Power(t,2) - 
		0.000085*Power(t,3)),2)*
	(1 + (489*p - 0.333*Power(p,2) - 80*t - 2*p*t - 
		0.002*Power(p,2)*t - 3.3*Power(t,2) + 
		0.016*p*Power(t,2) + 0.00175*Power(t,3) - 
		0.000013000000000000001*p*Power(t,3))/1.e6 + 
	s*(0.668 + 0.44*s + 
		(300*p - 2400*p*s + 
		t*(80 - 13*p - 3300*s + 47*p*s + 3*t))/1.e6
		)))/
	(1000.*(1 + 0.0494*Power(10,
	 -4 - (7.786*s)/Power(17.78 + t,0.306))*
	(3676*Power(p,0.64) + 
		0.712*p*Power(Abs(-76.71 + t),1.5))));

	return x;
}

float gassmann_brineVelocity(float s, float t, float p)
/*******************************************************************
gassmann_brineVelocity - Calculate sound speed (m/s) of saline solution characterized 
	  by weight fraction of sodium chloride S (ppm/1,000,000) as 
	  a function of temperature (deg C) and pressure (MPa)
*******************************************************************
Credits: UHouston: Chris Liner	9/23/2009
*******************************************************************/
{
	float x;

	x = 1402.85 + 1.524*p + 0.003437*Power(p,2) - 
	0.000011970000000000002*Power(p,3) + 
	(780 - 10*p + 0.16*Power(p,2))*Power(s,1.5) - 
	1820*Power(s,2) + 4.871*t - 0.0111*p*t + 
	0.00017390000000000003*Power(p,2)*t - 
	1.628e-6*Power(p,3)*t - 0.04783*Power(t,2) + 
	0.0002747*p*Power(t,2) - 
	2.135e-6*Power(p,2)*Power(t,2) + 
	1.2370000000000001e-8*Power(p,3)*Power(t,2) + 
	0.0001487*Power(t,3) - 
	6.502999999999999e-7*p*Power(t,3) - 
	1.4550000000000001e-8*Power(p,2)*Power(t,3) + 
	1.327e-10*Power(p,3)*Power(t,3) - 
	2.197e-7*Power(t,4) + 7.987e-10*p*Power(t,4) + 
	5.23e-11*Power(p,2)*Power(t,4) - 
	4.614e-13*Power(p,3)*Power(t,4) + 
	s*(1170 + 2.6*p - 0.0476*Power(p,2) - 9.6*t - 
	0.0029*p*t + 0.055*Power(t,2) - 
	0.000085*Power(t,3));

	return x;
}


float gassmann_brineDensity(float s, float t, float p)
/*******************************************************************
gassmann_brineDensity - Calculate mass density (g/cc) of saline solution characterized 
		by weight fraction of sodium chloride S (ppm/1,000,000) as 
		a function of temperature (deg C) and pressure (MPa)
*******************************************************************
Credits: UHouston: Chris Liner	9/23/2009
*******************************************************************/
{
	float x;

	x = 1 + (489*p - 0.333*Power(p,2) - 80*t - 2*p*t - 
	0.002*Power(p,2)*t - 3.3*Power(t,2) + 
	0.016*p*Power(t,2) + 0.00175*Power(t,3) - 
	0.000013000000000000001*p*Power(t,3))/1.e6 + 
	s*(0.668 + 0.44*s + 
	(300*p - 2400*p*s + 
	 t*(80 - 13*p - 3300*s + 47*p*s + 3*t))/1.e6);

	return x;
}


float gassmann_oilModulus(float api, float g, float t, float p)
/*******************************************************************
gassmann_oilModulus - Calculate bulk modulus (MPa) of live oil (including solution gas)
	characterized by oil gravity API, gas specific gravity G, 
	at temperature T (deg C) and pressure P (MPa)
*******************************************************************
Credits: UHouston: Chris Liner	9/23/2009
*******************************************************************/
{
	float x;

	x = gassmann_oilDensity(api,g,t,p) * Power( gassmann_oilVelocity(api,g,t,p), 2) / 1000.;

	return x;
}

float gassmann_oilVelocity(float api, float g, float t, float p)
/*******************************************************************
gassmann_oilVelocity - Calculate sound speed (m/s) of live oil (including solution gas)
	characterized by oil gravity API, gas specific gravity G, 
	at temperature T (deg C) and pressure P (MPa)
*******************************************************************
Credits: UHouston: Chris Liner	9/23/2009
*******************************************************************/
{
	float x;

	x = 4.64*p - 3.7*t + 24932.71072306419*
	Power(1/
	((131.5 + api)*(1 + 
		0.0020299999999999997*g*
		Power(Power(E,0.02878*api - 0.00377*t)*p,
		1.205))*(0.972 + 
		0.00038*Power(17.8 + 
		0.409570869105443*g*Sqrt((131.5 + api)*g)*
		Power(Power(E,0.02878*api - 0.00377*t)*p,
		1.205) + t,1.175))*
	(2.6 - 141.5/
		((131.5 + api)*
		(1 + 0.0020299999999999997*g*
		Power(Power(E,0.02878*api - 0.00377*t)*
			p,1.205))*
		(0.972 + 0.00038*
		Power(17.8 + 
			0.409570869105443*g*
			Sqrt((131.5 + api)*g)*
			Power(Power(E,
			0.02878*api - 0.00377*t)*p,1.205) + 
			t,1.175))))),0.5) + 
	0.0115*p*t*(-1 + 4.12*
	Sqrt(-1 + 0.007632508833922262*(131.5 + api)*
		(1 + 0.0020299999999999997*g*
		Power(Power(E,0.02878*api - 0.00377*t)*p,
		1.205))*(0.972 + 
		0.00038*Power(17.8 + 
		0.409570869105443*g*
		Sqrt((131.5 + api)*g)*
		Power(Power(E,0.02878*api - 0.00377*t)*
			p,1.205) + t,1.175))));

	return x;
}


float gassmann_oilDensity(float api, float g, float t, float p)
/*******************************************************************
gassmann_oilDensity - Calculate mass density (g/cc) of live oil (including solution gas) 
	characterized by oil gravity API, gas specific gravity G, 
	at temperature T (deg C) and pressure P (MPa)
*******************************************************************
Credits: UHouston: Chris Liner	9/23/2009
*******************************************************************/
{
	float x;

	x = 0.000349*p + (141.5/(131.5 + api) + 
	0.0024359999999999994*Power(g,2)*
	Power(Power(E,0.02878*api - 0.00377*t)*p,1.205))/
	(0.972 + 0.00038*Power(17.8 + 
	 0.409570869105443*g*Sqrt((131.5 + api)*g)*
		Power(Power(E,0.02878*api - 0.00377*t)*p,
		1.205) + t,1.175)) + 
	(0.00277*p - 1.7099999999999998e-7*Power(p,3))*
	Power(-1.15 + (141.5/(131.5 + api) + 
	 0.0024359999999999994*Power(g,2)*
		Power(Power(E,0.02878*api - 0.00377*t)*p,
		1.205))/
	(0.972 + 0.00038*
		Power(17.8 + 
		0.409570869105443*g*Sqrt((131.5 + api)*g)*
		Power(Power(E,0.02878*api - 0.00377*t)*p,
		1.205) + t,1.175)),2);

	return x;
}


float gassmann_gasVelocity(float g, float t, float p)
/*******************************************************************
gassmann_gasVelocity - Calculate sound speed (m/s) of natural gas characterized by 
	specific gravity G, at temperature T (deg C) and pressure P (MPa)
*******************************************************************
Credits: UHouston: Chris Liner	9/23/2009
*******************************************************************/
{
	float x;

	x = sqrt( ( 1000. * gassmann_gasModulus(g,t,p) ) / gassmann_gasDensity(g,t,p) );

	return x;
}


float gassmann_gasModulus(float g, float t, float p)
/*******************************************************************
gassmann_gasModulus - Calculate bulk modulus (MPa) of natural gas characterized by 
	specific gravity G, at temperature T (deg C) and pressure P (MPa)
*******************************************************************
Credits: UHouston: Chris Liner	9/23/2009
*******************************************************************/
{
	float x;

	x = (p*(0.85 - 8.7/
	Power(E,0.65*(1 + p/(4.892 - 0.4048*g))) + 
	5.6/(2 + p/(4.892 - 0.4048*g)) + 
	27.1/Power(3.5 + p/(4.892 - 0.4048*g),2)))/
	(1 - (p*(0.03 + 0.00527*
		Power(3.5 - (273.15 + t)/(94.72 + 170.75*g),
		3) + (0.0228*
		Power(E,((94.72 + 170.75*g)*
		 Power(p/(4.892 - 0.4048*g),1.2)*
		 (-0.45 - 
			8*Power(0.56 - 
			(94.72 + 170.75*g)/(273.15 + t),2)))
		/(273.15 + t))*(94.72 + 170.75*g)*
		Power(p/(4.892 - 0.4048*g),0.2)*
		Power(3.85 - 
		(273.15 + t)/(94.72 + 170.75*g),2)*
		(-0.45 - 8*
		Power(0.56 - 
			(94.72 + 170.75*g)/(273.15 + t),2)))/
		(273.15 + t)))/
	((4.892 - 0.4048*g)*
	(-0.52 + (0.642*(273.15 + t))/
		(94.72 + 170.75*g) - 
		(0.007*Power(273.15 + t,4))/
		Power(94.72 + 170.75*g,4) + 
		0.019*Power(3.85 - 
		(273.15 + t)/(94.72 + 170.75*g),2) + 
		(p*(0.03 + 0.00527*
		Power(3.5 - 
			(273.15 + t)/(94.72 + 170.75*g),3)))/
		(4.892 - 0.4048*g))));

	return x;
}


float gassmann_gasDensity(float g, float t, float p)
/*******************************************************************
gassmann_gasDensity - Calculate mass density (g/cc) of natural gas 
	characterized by specific gravity G, at temperature T 
	(deg C) and pressure P (MPa)
*******************************************************************
Credits: UHouston: Chris Liner	9/23/2009
*******************************************************************/
{
	float x;

	x = (3.4638699124410666*g*p)/
	((273.15 + t)*(-0.52 + 
	(0.642*(273.15 + t))/(94.72 + 170.75*g) - 
	(0.007*Power(273.15 + t,4))/
	Power(94.72 + 170.75*g,4) + 
	0.109*Power(E,((94.72 + 170.75*g)*
		Power(p/(4.892 - 0.4048*g),1.2)*
		(-0.45 - 8*
		Power(0.56 - 
		 (94.72 + 170.75*g)/(273.15 + t),2)))/
		(273.15 + t))*
	Power(3.85 - (273.15 + t)/(94.72 + 170.75*g),2)
	+ (p*(0.03 + 0.00527*
		Power(3.5 - 
		(273.15 + t)/(94.72 + 170.75*g),3)))/
	(4.892 - 0.4048*g)));

	return x;
}
