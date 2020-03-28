/************************************************************************
 * ComFluSoM - Simulation kit for Fluid Solid Soil Mechanics            *
 * Copyright (C) 2019 Pei Zhang                                         *
 * Email: peizhang.hhu@gmail.com                                        *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/
class MPM_PARTICLE
{
public:
	MPM_PARTICLE();
	// MPM_PARTICLE(int type, const Vector3d& x, double m, double young, double poisson);
	MPM_PARTICLE(int tag, const Vector3d& x, double m);
	void Elastic(Matrix3d& de);
	void SetElastic(double young, double poisson);
	void SetNewtonian(double miu);
	void SetMohrCoulomb(double young, double poisson, double phi, double psi, double c);
	void SetDruckerPrager(int dptype, double young, double poisson, double phi, double psi, double c);
	void SetTensionCutoff(double pmax);
	void SetHoekBrown(double young, double poisson, double a, double s, double mb, double sigmaci, double mi);	
	void SetModifiedCamClay(double young, double poisson, double mslope, double kappa, double cclamda, double nvalue, double pini, double e0);
	void Newtonian(Matrix3d& de);
	void MohrCoulomb(Matrix3d& de);
	void DruckerPrager(Matrix3d& de);
	void HoekBrown(Matrix3d& de);
	void ModifiedCamClay(Matrix3d& de);
	void EOSMorris(double C);
	void EOSMonaghan(double C);
	

    int 						Type;                       // Type of particle, 0 for elastic 1 for fluid 2 for soil.
	int 						ID; 				    	// Index of particle in the list 
	int 						Tag;				    	// Tag of particle

	double 						M;				            // Mass
	double 						Vol;						// Volume
	double 						Vol0;						// Init Volume
	double 						R;							// Radius for DEMPM
	double 						Arc;						// Arc length for boundary nodes
	double 						Mu;							// Shear modulus (Lame's second parameter) or viscosity for fluid
	double                      Mu1;                        // Cam-Clay parameters
	double						La;							// Lame's first parameter
	double                      G;                          // shear modulus
	double 						K;							// bulk modulus
	double 						H;							// Liniear hardening modulus
	double 						Young;						// Young's modus
	double						Poisson;					// Possion ratio
	double 						C;							// Cohesion coefficient, unit [kg/(m*s^2)] (or Pa)
	double 						Phi;						// Angle of internal friction
	double 						Psi;						// Angle of dilatation
	double 						A_dp;						// Drucker–Prager parameters
	double 						B_dp;						// Drucker–Prager parameters
	double 						Ad_dp;						// Drucker–Prager parameters
	double 						Pmax;						// Drucker–Prager parameters for tension cutoff (max pressure)
	bool 						TensionCut;					// Drucker–Prager parameters for tension cutoff
	double						A;							// Hoek Brown rock mass material constants 
	double						S;							// Hoek Brown rock mass material constants
	double						Mb;							// Hoek Brown rock mass material constants 
	double                      Mi;							// Hoek Brown rock mass material constants
	double 						Sigmaci;					// Hoek Brown unconfined compressive strength
    double 						Mslope; 					// Cam-Clay parameters
    double 						Kappa;   					// Cam-Clay parameters
    double 						CCLamda; 					// Cam-Clay Parameters
    double 						Nvalue;					    // Cam-Clay Parameters 
    double 						E0; 						// Cam-Clay initial void
    double 						Pini; 						// Cam-Clay initial presure or conslidation presure

	double 						P;							// Pressure of fluid

	Vector3d					PSize0;						// Vector of half length of particle domain at init
	Vector3d					PSize;						// Vector of half length of particle domain

	Vector3d 					X;				            // Position
	Vector3d 					X0;				            // Init position
	Vector3d 					DeltaX;				        // Increasement of position
	Vector3d					V;							// Velocity
	Vector3d					Vf;							// Fixed velocity
	Vector3d					B;							// Body force acc
	Vector3d					Fh0;						// Hydro force
	Vector3d					Fh;							// Hydro force
	Vector3d					Fc;							// Contact force
	Vector3d					Nor;						// Normal direction (only non-zero for boundary particles)

	Matrix3d					Strain;						// Strain
	Matrix3d					StrainP;					// Plastic strain tensor
	Matrix3d					Stress;						// Stress
	Matrix3d					StressSmooth;				// Smoothed Stress for visualization
	Matrix3d					L;							// Velocity gradient tensor
	Matrix3d					F;							// Derformation gradient tensor
	Matrix3d					Dp;							// Elastic tensor in principal stress space
	Matrix3d					Dpi;						// Inverse of Dp

	bool						FixV;						// Whether the velocity is fixed
	bool						Removed;					// whether this particle is removed

	vector<int>					Lnei;						// List of neighor nodes indexs, used to calculate arc lengh for FSI problems
	vector<size_t>				Lni;						// List of node indexs
	vector<size_t>				Lgi;						// List of gauss point indexs
	vector<double>				LnN;						// List of shape functions
	vector<Vector3d>			LnGN;						// List of gradient of shape functions
};

inline MPM_PARTICLE::MPM_PARTICLE()
{
    Type	= -1;
	ID		= 0;
	Tag		= 0;
	M 		= 0.;
	X 		= Vector3d::Zero();
	X0 		= Vector3d::Zero();
	V 		= Vector3d::Zero();
	Vf 		= Vector3d::Zero();
	B 		= Vector3d::Zero();
	Fh 		= Vector3d::Zero();
	Fh0 	= Vector3d::Zero();
	Fc 		= Vector3d::Zero();
	Nor 	= Vector3d::Zero();

	Strain 	= Matrix3d::Zero();
	StrainP = Matrix3d::Zero();
	Stress 	= Matrix3d::Zero();
	StressSmooth = Matrix3d::Zero();
	F 		= Matrix3d::Identity();

	Lni.resize(0);
	LnN.resize(0);
	LnGN.resize(0);

	FixV	= false;
	Removed	= false;
}

// inline MPM_PARTICLE::MPM_PARTICLE(int type, const Vector3d& x, double m, double young, double poisson)
// {
//     Type	= type;
// 	ID		= 0;
// 	Tag		= 0;
// 	M 		= m;
// 	X 		= x;
// 	X0 		= x;
// 	V 		= Vector3d::Zero();
// 	Vf 		= Vector3d::Zero();
// 	B 		= Vector3d::Zero();
// 	Fh 		= Vector3d::Zero();
// 	Nor 	= Vector3d::Zero();

// 	Strain 	= Matrix3d::Zero();
// 	StrainP = Matrix3d::Zero();
// 	Stress 	= Matrix3d::Zero();
// 	F 		= Matrix3d::Identity();

// 	Lni.resize(0);
// 	LnN.resize(0);
// 	LnGN.resize(0);

// 	FixV	= false;
// 	Removed	= false;

// 	Young 	= young;
// 	Poisson = poisson;

// 	Mu 		= 0.5*Young/(1.+Poisson);
// 	La 		= Young*Poisson/(1.+Poisson)/(1.-2.*Poisson);
// 	K 		= La+2./3.*Mu;
// 	H 		= 0.;
// 	Dp(0,0) = Dp(1,1) = Dp(2,2) = La+2.*Mu;
// 	Dp(0,1) = Dp(1,0) = Dp(0,2) = Dp(2,0) = Dp(1,2) = Dp(2,1) = La;

// 	Dpi = Dp.inverse();

// 	// Matrix3d dpi;
// 	// dpi(0,0) = dpi(1,1) = dpi(2,2) = 1/Young;
// 	// dpi(0,1) = dpi(1,0) = dpi(0,2) = dpi(2,0) = dpi(1,2) = dpi(2,1) = -Poisson/Young;

// 	// cout << Dpi << endl;
// 	// cout << "========" << endl;
// 	// cout << dpi << endl;
// 	// abort();
// }

inline MPM_PARTICLE::MPM_PARTICLE(int tag, const Vector3d& x, double m)
{
    Type	= -1;
	ID		= 0;
	Tag		= tag;
	M 		= m;
	X 		= x;
	X0 		= x;
	V 		= Vector3d::Zero();
	Vf 		= Vector3d::Zero();
	B 		= Vector3d::Zero();
	Fh 		= Vector3d::Zero();
	Fc 		= Vector3d::Zero();
	Nor 	= Vector3d::Zero();

	Strain 	= Matrix3d::Zero();
	StrainP = Matrix3d::Zero();
	Stress 	= Matrix3d::Zero();
	StressSmooth = Matrix3d::Zero();
	F 		= Matrix3d::Identity();

	Lni.resize(0);
	LnN.resize(0);
	LnGN.resize(0);

	FixV	= false;
	Removed	= false;
}

inline void MPM_PARTICLE::SetElastic(double young, double poisson)
{
	Type 	= 0;
	Young 	= young;
	Poisson = poisson;
	Mu 		= 0.5*Young/(1.+Poisson);
	La 		= Young*Poisson/(1.+Poisson)/(1.-2.*Poisson);
	K 		= La+2./3.*Mu;
	H 		= 0.;
	Dp(0,0) = Dp(1,1) = Dp(2,2) = La+2.*Mu;
	Dp(0,1) = Dp(1,0) = Dp(0,2) = Dp(2,0) = Dp(1,2) = Dp(2,1) = La;
	Dpi = Dp.inverse();
}

inline void MPM_PARTICLE::SetNewtonian(double miu)
{
	Type 	= 1;
	Mu 		= miu;
}

inline void MPM_PARTICLE::SetMohrCoulomb(double young, double poisson, double phi, double psi, double c)
{
	Type 	= 2;
	Young 	= young;
	Poisson = poisson;
	Mu 		= 0.5*Young/(1.+Poisson);
	La 		= Young*Poisson/(1.+Poisson)/(1.-2.*Poisson);
	K 		= La+2./3.*Mu;
	H 		= 0.;
	Dp(0,0) = Dp(1,1) = Dp(2,2) = La+2.*Mu;
	Dp(0,1) = Dp(1,0) = Dp(0,2) = Dp(2,0) = Dp(1,2) = Dp(2,1) = La;
	Dpi 	= Dp.inverse();
	Phi 	= phi;
	Psi 	= psi;
	C 		= c;
}

inline void MPM_PARTICLE::SetDruckerPrager(int dptype, double young, double poisson, double phi, double psi, double c)
{
	Type 	= 3;
	Young 	= young;
	Poisson = poisson;
	Mu 		= 0.5*Young/(1.+Poisson);
	La 		= Young*Poisson/(1.+Poisson)/(1.-2.*Poisson);
	K 		= La+2./3.*Mu;
	H 		= 0.;
	Dp(0,0) = Dp(1,1) = Dp(2,2) = La+2.*Mu;
	Dp(0,1) = Dp(1,0) = Dp(0,2) = Dp(2,0) = Dp(1,2) = Dp(2,1) = La;
	Dpi 	= Dp.inverse();
	Phi 	= phi;
	Psi 	= psi;
	C 		= c;
	TensionCut = false;
	// inner
	if (dptype==0)
	{
		double bot = sqrt(3)*(3.+sin(Phi));
		A_dp = 6.*sin(Phi)/bot;
		B_dp = 6.*cos(Phi)/bot;
		Ad_dp = 6.*sin(Psi)/(sqrt(3)*(3.+sin(Psi)));
	}
	// outer
	else if (dptype==1)
	{
		double bot = sqrt(3)*(3.-sin(Phi));
		A_dp = 6.*sin(Phi)/bot;
		B_dp = 6.*cos(Phi)/bot;
		Ad_dp = 6.*sin(Psi)/(sqrt(3)*(3.-sin(Psi)));
	}
	// plane strain
	else if (dptype==2)
	{
		double bot = sqrt(9.+12*tan(Phi)*tan(Phi));
		A_dp = 3.*tan(Phi)/bot;
		B_dp = 3./bot;
		Ad_dp = 3.*tan(Psi)/sqrt(9.+12*tan(Psi)*tan(Psi));
	}
}
// only works with Drucker Prager (12.3.2019)
inline void MPM_PARTICLE::SetTensionCutoff(double pmax)
{
	TensionCut 	= true;
	Pmax 		= pmax;
}
inline void MPM_PARTICLE::SetHoekBrown(double young, double poisson, double a, double s, double mb, double sigmaci, double mi)
{

	Type 		= 4;
	Young       = young;
	Poisson 	= poisson;
	Mu 			= 0.5*Young/(1.+Poisson);
	La 			= Young*Poisson/(1.+Poisson)/(1.-2.*Poisson);
    A           = a;
	S           = s;
	Mb 			= mb;
	Sigmaci     = sigmaci;
	Mi 			= mi;
	Dp(0,0) 	= Dp(1,1) = Dp(2,2) = La+2.*Mu;
	Dp(0,1) 	= Dp(1,0) = Dp(0,2) = Dp(2,0) = Dp(1,2) = Dp(2,1) = La;
	Dpi 		= Dp.inverse();
}
inline void MPM_PARTICLE::SetModifiedCamClay(double young, double poisson, double mslope, double kappa, double cclamda, double nvalue,  double pini, double e0)
{
	Type 		= 5;
    Young 		= young;
	Poisson 	= poisson;
    Mu1 	    = 0.5*Young/(1.+Poisson);
    Mslope 		= mslope;
	Kappa 		= kappa;
	CCLamda     = cclamda;
	Nvalue      = nvalue;
    Pini 		= pini;   											// initial presure or consolidation presure
    E0          = e0;   											// initial viod
    Mu    		= (3.*(1.-2.*Mu1)*(1.+E0)*Pini)/(2.*(1.+Mu1)*Kappa);
    La 			= ((1.+E0)*Pini)/Kappa;
	Dp(0,0) 	= Dp(1,1) = Dp(2,2) = La+(4./3.)*Mu;
	Dp(0,1) 	= Dp(1,0) = Dp(0,2) = Dp(2,0) = Dp(1,2) = Dp(2,1) = La-(2/3)*Mu;
	Dpi 		= Dp.inverse();
}
// Elastic model
inline void MPM_PARTICLE::Elastic(Matrix3d& de)
{
	Stress += 2.*Mu*de + La*de.trace()*Matrix3d::Identity();
}

// Newtonian fluid model
inline void MPM_PARTICLE::Newtonian(Matrix3d& de)
{
	Stress = 2.*Mu*(de - de.trace()/3.*Matrix3d::Identity()) - P*Matrix3d::Identity();
}

void MPM_PARTICLE::EOSMorris(double Cs)
{
	P = Cs*Cs*M/Vol;
}

void MPM_PARTICLE::EOSMonaghan(double Cs)
{
	P = Cs*Cs*M/Vol0/7.*(pow(Vol0/Vol,7.)-1.);
}

// Mohr-Coulomb model
// Based on "An efficient return algorithm for non-associated plasticity with linear yield criteria in principal stress space"
inline void MPM_PARTICLE::MohrCoulomb(Matrix3d& de)
{
	// Apply elastic model first
	Elastic(de);
	SelfAdjointEigenSolver<Matrix3d> eigensolver(Stress);

	double s1 = eigensolver.eigenvalues()(2);
	double s2 = eigensolver.eigenvalues()(1);
	double s3 = eigensolver.eigenvalues()(0);

	Vector3d sb (s1, s2, s3);

	if (s1<s2 || s1<s3 || s2<s3)
	{
		cout << "wrong order of s" << endl;
		abort();
	}

	double f = (s1-s3) +(s1+s3)*sin(Phi) -2.*C*cos(Phi);		// Eq.28

	if (f>0.)
	{
		Matrix3d v0;

		v0.col(0) = eigensolver.eigenvectors().col(2);
		v0.col(1) = eigensolver.eigenvectors().col(1);
		v0.col(2) = eigensolver.eigenvectors().col(0);

		Matrix3d spp = Matrix3d::Zero();
		spp(0,0) = s1;
		spp(1,1) = s2;
		spp(2,2) = s3;

		// cout << "Young= " << Young << endl;
		// cout << "Poisson= " << Poisson << endl;
		// cout << "La= " << La << endl;
		// cout << "Mu= " << Mu << endl;
		// abort();

		// Matrix3d ss = v0 * spp * v0.inverse();

		// cout << S << endl;
		// cout << "=========" << endl;
		// cout << ss << endl;
		// abort();

		Vector3d sc;

		double k = (1.+sin(Phi)) / (1.-sin(Phi));				// Eq.32
		Vector3d a1 (k, 0., -1.);								// Eq.32
		double m = (1.+sin(Psi)) / (1.-sin(Psi));				// Eq.33
		Vector3d b1 (m, 0., -1.);								// Eq.33

		Vector3d rp = Dp*b1 / (b1.transpose()*Dp*a1);			// Eq.27b

		Vector3d sa (1., 1., 1.);
		// sa *= 2.*C*sqrt(k)/(k-1.);								// Eq.34
		sa *= C/tan(Phi);										// (6.117)

		Vector3d rl1 (1., 1., k);								// Eq.40
		Vector3d rl2 (1., k , k);								// Eq.40

		Vector3d rgl1 (1., 1., m);								// Eq.41
		Vector3d rgl2 (1., m , m);								// Eq.41

		double t1 = rgl1.transpose()*Dpi*(sb-sa); 				// Eq.39
		double t2 = rgl2.transpose()*Dpi*(sb-sa);				// Eq.39

		double t1f = rgl1.transpose()*Dpi*rl1;
		double t2f = rgl2.transpose()*Dpi*rl2;

		t1 /= t1f;
		t2 /= t2f;

		double p12 = rp.cross(rl1).dot(sb-sa);					// Eq.45
		double p13 = rp.cross(rl2).dot(sb-sa);					// Eq.46

		// return to apex
		if (t1>0. && t2>0.)
		{
			sc = sa;											// Eq.42
		}
		// return to plane f=0
		else if (p12>=0. && p13<=0.)
		{
			Vector3d dsp = f*rp;								// Eq.27a
			sc = sb - dsp;										// Eq.6
		}
		// return to l1
		else if (p12<0. && p13<0.)
		{
			sc = t1*rl1 + sa;									// Eq.40
		}
		// return to l2
		else if (p12>0. && p13>0.)
		{
			sc = t2*rl2 + sa;									// Eq.40
		}
		else
		{
			cout << "undefined" << endl;
			abort();
		}

		// cout << "f= " << f << endl;

		// double ff = (sc(0)-sc(2)) +(sc(0)+sc(2))*sin(Phi) -2.*C*cos(Phi);

		// if (ff>1.e-8)
		// {
		// 	cout << "f= " << f << endl;
		// 	cout << "ff= " << ff << endl;
		// 	// abort();
		// }

		Matrix3d sp = Matrix3d::Zero();
		sp(0,0) = sc(0);
		sp(1,1) = sc(1);
		sp(2,2) = sc(2);

		Stress = v0 * sp * v0.inverse();

		// if (ID==5195)
		// {
		// 	cout << "f= " << f << endl;
		// 	cout << "ff= " << ff << endl;
		// 	cout << "===================" << endl;
		// 	cout << S << endl;
		// 	cout << "===================" << endl;
		// 	// abort();
 	// 	}

		// if (s1>1.0e-12)
		// {
		// 	cout << s1 << endl;
		// 	cout << s2 << endl;
		// 	cout << s3 << endl;
		// 	cout << S << endl;
		// 	cout << "-------" << endl;
		// 	cout << v0 * spp * v0.inverse() << endl;
		// 	abort();
		// }

		// S(2,2) = S(0,1) = S(0,2) = S(1,0) = S(2,0) = 0.;
	}
}


//DruckerPrager model 
void MPM_PARTICLE::DruckerPrager(Matrix3d& de)
{
	// auto YieldFunc = [](Matrix3d s)
	// {
	// 	double p = s.trace()/3.;
	// 	Matrix3d ss = s - p*Matrix3d::Identity();
	// 	double j2 = 0.5*(ss.array()*ss.array()).sum();
	// 	return sqrt(j2)+A_dp*p-B_dp*C;
	// };
	// Apply elastic model first
	Elastic(de);
	// pressure
	double p = Stress.trace()/3.;
	Matrix3d ss = Stress - p*Matrix3d::Identity();
	double j2 = 0.5*(ss.array()*ss.array()).sum();
	double f = sqrt(j2)+A_dp*p-B_dp*C;

    //	std::cout <<"matrix： "<< Stress <<std::endl; 
	double tao 		= sqrt(j2);									//  Effective shear stress 	
	double taop 	= B_dp-A_dp*Pmax; 						// define parameter for h function
	double aphaP    = sqrt(1+A_dp*A_dp)-A_dp;					// define parameter for h function
	double h 		= tao-taop-aphaP*(p-Pmax); 				// define h function for the diagonal line between the line fs and ft 
	double detalamdas = f/(Mu+K*A_dp*Ad_dp);											//Eq43
	//|| (TensionCut && p>Pmax)
    // std::cout << Stress <<std::endl;

	if (f>0. || (TensionCut && p>Pmax))
	{
		SelfAdjointEigenSolver<Matrix3d> eigensolver(Stress);
		double s1 = eigensolver.eigenvalues()(2);
		double s2 = eigensolver.eigenvalues()(1);
		double s3 = eigensolver.eigenvalues()(0);

		Matrix3d v0;
		v0.col(0) = eigensolver.eigenvectors().col(2);
		v0.col(1) = eigensolver.eigenvectors().col(1);
		v0.col(2) = eigensolver.eigenvectors().col(0);

		Matrix3d ssp = Matrix3d::Zero();
		double pp = (s1+s2+s3)/3.;
		ssp(0,0) = s1-pp;
		ssp(1,1) = s2-pp;
		ssp(2,2) = s3-pp;

		Vector3d df = (0.5*ssp/sqrt(j2) + A_dp/3.*Matrix3d::Identity()).diagonal();		// 3.2 and (6.156)
		Vector3d dg = (0.5*ssp/sqrt(j2) + Ad_dp/3.*Matrix3d::Identity()).diagonal();	// 3.2 and (6.166)
		double chi = dg.transpose()*Dp*df;	// 3.4
		double lamda = f/chi;				// 3.6
		
		Matrix3d sp = Matrix3d::Zero();
		sp(0,0) = s1-lamda*dg(0);			// 3.5
		sp(1,1) = s2-lamda*dg(1);
		sp(2,2) = s3-lamda*dg(2);

		// moodifed pressure at principal stress space
		double ppm = (sp(0,0)+sp(1,1)+sp(2,2))/3.;
		// pressure at apex
		double p_apex = B_dp*C/A_dp;
		Vector3d apex (p_apex, p_apex, p_apex);
		// if (sp.trace()/3.>p_apex)	 sp.diagonal() << p_apex, p_apex, p_apex;
		if (TensionCut)
		{
			if (ppm>Pmax)
			{
				//if (sqrt(j2)-A_dp*Pmax+B_dp*C>0.)
				if ( h>0.)
				{
					double corrpp     = pp-K*A_dp*detalamdas;											//Eq45
   					double corrtao    = B_dp-A_dp*corrpp;	
					sp = ss * (corrtao/sqrt(j2)) + corrpp * Matrix3d::Identity();

					//Vector3d n = (sp.diagonal()-apex).normalized();
					//sp.diagonal() += 3.*(Pmax-ppm)/(n(0)+n(1)+n(2))*n;
				}
				else 
				{
					
        			Matrix3d corrt = Matrix3d::Zero();   

					corrt = (Pmax-pp)*Matrix3d::Identity();    							                // Eq55
        
       				Matrix3d sp = Matrix3d::Zero();
					sp(0,0) = s1+corrt(0);																// Eq55
					sp(1,1) = s2+corrt(1);
					sp(2,2) = s3+corrt(2);

					//double deltaP = pp-Pmax;
					//sp.diagonal() << s1+deltaP, s2+deltaP, s3+deltaP;
				}
			}
		}
		else
		{
			if (ppm>p_apex)		sp.diagonal() = apex;
		}
		// convert back from principal stress space
		Stress = v0 * sp * v0.inverse();
	}
}

// Matsuoka-Nakai model
//void MPM_PARTICLE::MNModel()
//{
//	double sinPhi2 = sin(Phi)*sin(Phi);
//	double kf = (sinPhi2-9.)/(sinPhi2-1.);
//
	// Matrix3d Sb = S - C*Matrix3d::Identity();

//	double I1 = S.trace();
//	double I2 = 0.5*(I1*I1 - (S*S).trace());
//	double I3 = S.determinant();
//
//	double f = 6.*(I3*kf -I1*I2) + C*(12.*I1*I1 + 18.*I2 - 6.*I2*kf) + C*C*(6.*I1*kf - 54.*I1) + C*C*C*(54. - 6.*kf);
//}
 void MPM_PARTICLE::HoekBrown(Matrix3d& de)
// Based on "an exact implementation of the Hoek Brown criterion for elasto-plastic finite element caculations" Paper 1
//Finite element implementation of the Hoek Brown material model with general sofenting behavior" Paper 2
{
	// Apply elastic model first
	Elastic(de);
	SelfAdjointEigenSolver<Matrix3d> eigensolver(Stress);

	double s1 = eigensolver.eigenvalues()(2);
	double s2 = eigensolver.eigenvalues()(1);
	double s3 = eigensolver.eigenvalues()(0);



   // std::cout <<"matrix-O： "<< Stress <<std::endl; 
   // std::cout <<"matrix-e： "<< Strain <<std::endl;

	Vector3d sb (s1, s2, s3);

	if (s1<s2 || s1<s3 || s2<s3)
	{
		cout << "wrong order of s" << endl;
		abort();
	}

	double f = s1-s3-Sigmaci*pow((S-Mb*(s1/Sigmaci)),A);	//Eq. 7 Paper 1

	if (f>0.)
	{
		Matrix3d v0;

		v0.col(0) = eigensolver.eigenvectors().col(2);
		v0.col(1) = eigensolver.eigenvectors().col(1);
		v0.col(2) = eigensolver.eigenvectors().col(0);
		
		Vector3d sc;
            
        double parf = S-Mb*s1/Sigmaci;			      			//Value of f-parenthesis	Eq.31 Ppaer 1
		double k =  1. + A*Mb * pow(parf,(A-1.));				// First component of yield surface normal
		Vector3d a (k, 0., -1.);								// Eq.30a Ppaer 1

		double Ag = A;
		double Mg = Mb;
		double Sg = S;
		double parg = Sg-Mg*s1/Sigmaci;		      				//Value of g-parenthesis 	Eq.32 Ppaer 1
	    double kg = 1. + Ag*Mg* pow(parg,(Ag-1.));				//First component of plastic potential gradient				
		Vector3d b  (kg, 0., -1.);								// Eq.30b Ppaer 1
		Vector3d b1 (1., 1., -2./kg);                           // Eq.45  Ppaer 1
		Vector3d b2 (2.*kg, 1., 1.);							// Eq.47  Ppaer 1

		Vector3d rp = Dp*b / (b.transpose()*Dp*a);				// Eq.27b another paper the same as MB
		

		Vector3d sa (1., 1., 1.);
		sa *= (S*Sigmaci)/ Mb;								 	//using the new HB failure criterion from the new HB2018 version 

		Vector3d Rpa 			= Dp*b;							// Plastic corrector stress at apex by return from region I    					Eq.48 	Ppaer 1                 
    	Vector3d Ra_II 			= Dp*b1;						// Plastic corrector stress at apex by return along the plane sigma1 = sigma2
    	Vector3d Ra_III 		= Dp*b2;						// Plastic corrector stress at apex by return along the plane sigma2 = sigma3
    	Vector3d N_IV_II 		= Rpa.cross(Ra_II);				// Normal of the region II and IV boundary plane								Eq.35 Ppaer 2 
    	Vector3d N_IV_III 		= Ra_III.cross(Ra_III);			// Normal of the region III and IV boundary plane								Eq.35 Ppaer 2 

		double t1   = N_IV_II.dot(sb-sa);                       // Eq.49 and 50 	Ppaer 1
		double t2   = N_IV_III.dot(sb-sa);                      //

		
		double sigmavalue  = s1-Sigmaci*pow(parg,A);
		Vector3d l1 (s1, s1, sigmavalue);                       // Eq.11 value of line1 surface failure 
		Vector3d l2 (s1, sigmavalue, sigmavalue);               // Eq.12 value of line2 surface failure 
		
		Vector3d sigmael1 (s1, s1, s3);							// Figure 11 Paper 2
		Vector3d sigmael6 (s1, s3, s3);							// Figure 11 Paper 2
		Vector3d ne1 (1.,1.,kg);								// Eq.38 Paper 2
		Vector3d ne6 (1.,kg,kg);								// Eq.38 Paper 2

		double p12 = ne1.dot(sb-sigmael1);						// Eq.39 Paper 2
		double p13 = ne6.dot(sb-sigmael6);					    // Eq.39

		// return to apex
		if (t1>0. && t2>0.)
		{
			sc = sa;											// Table 1 Paper 2
		}
		// return to plane f=0
		else if (p12<0. && p13<0.)
		{
			Vector3d dsp = f*rp;								// From the previous one MC
			sc = sb - dsp;										// 
		}
		// return to l1
		else if (p12>0. && p13<0.)
		{
			sc = l1;											// Eq.31
		}
		// return to l2
		else if (p12<0. && p13>0.)
		{
			sc = l2;											// Eq.32
		}
		else
		{
			cout << "undefined" << endl;
			abort();
		}

		Matrix3d sp = Matrix3d::Zero();
		sp(0,0) = sc(0);
		sp(1,1) = sc(1);
		sp(2,2) = sc(2);

		Stress = v0 * sp * v0.inverse();
	}
}

void MPM_PARTICLE::ModifiedCamClay(Matrix3d& de)
{
 	// Apply elastic model first
    Elastic(de);

   // Volumetric stress and devatoric stress
   // Matrix3d epsilonvolumetric = de.trace()/3.*Matrix3d::Identity();  // volumetirc strain
   // Matrix3d epsilondeviatoric = de-epsilonvolumetric;                // deviatoric strain
 
	SelfAdjointEigenSolver<Matrix3d> eigensolver(Stress);
	double s1 = eigensolver.eigenvalues()(2);
	double s2 = eigensolver.eigenvalues()(1);
	double s3 = eigensolver.eigenvalues()(0);

  //  std::cout <<"matrix-O： "<< Stress <<std::endl; 
  //  std::cout <<"matrix-e： "<< Strain <<std::endl; 

	if (s1<s2 || s1<s3 || s2<s3)
	{
		cout << "wrong order of s" << endl;
		abort();
	}

	double pvol 		= (s1+s2+s3)/3.;   // Mean principal stress 
	double qdev 		= (s1-s3);	       // deviatoric stress

	double En    		= (Nvalue-CCLamda*log10(Pini)+(Kappa*log10(Pini/pvol)))-1.;

	//std::cout <<"matrix： "<< Stress <<std::endl; 
	//double G1   = (3*(1-2*Mu)*(1+E0)*Pvol)/(2*(1+Mu)*kappa);                              // The K and G will change because the void change when it deforms
    //double K1 	= ((1+E0)*Pvol)/kappa;
	
	//Matrix3d sv = K*epsilonvolumetric;              	 									// Volumetric stress
	
 	//Matrix3d sd = 2*G*epsilondeviatoric; 													// Deviatoric stress
   
 	double Pccam		= pvol-((qdev*qdev)/((Mslope*Mslope)*pvol));

 	double f 			= pvol*pvol-pvol*Pccam+((qdev*qdev)/(Mslope*Mslope));      			// Yield function
	

    
    //the following part is for the aniostropic function 

	//Matrix3d alphaij = Matrix3d::Zero();

	//double alpha = sqrt((3/2)*((alphaij).dot(alphaij)));                //双点积应该如何表示
	
	//double sprsults = (sd-pvol*alphaij).dot(sd-pvol*alphaij);

	//double f = pvol*pvol-pvol*pc+(3/(2*(Mslope*Mslope-alpha*alpha))*sprsults;     //yield function
   
    if (f>0.)
	{
		Matrix3d v0;

		v0.col(0) = eigensolver.eigenvectors().col(2);
		v0.col(1) = eigensolver.eigenvectors().col(1);
		v0.col(2) = eigensolver.eigenvectors().col(0);



		Matrix3d dfdsigma   = (((2.*pvol-Pccam)/3.)*Matrix3d::Identity())+(3.*(Stress-pvol*Matrix3d::Identity())/Mslope*Mslope);       

		double dfdp 	  	= 2.*pvol-Pccam;

		

	    double Pcbar  		= ((1.+En)/(CCLamda-Kappa))*(dfdp)*Pccam;

	    Matrix3d Dep1       = (dfdsigma*Dp*dfdsigma.transpose()*Dp);


		Matrix3d Dep2       = (((-dfdp)*Pcbar*Matrix3d::Identity())+dfdsigma*Dp*dfdsigma.transpose());
        
       
   		Matrix3d Dep2i 		= Dep2.inverse();
        
        Matrix3d Dep        = Dep1*Dep2i;

	    //std::cout << Dep<<std::endl; 

		Matrix3d sigmaij    = Stress - Dep*de;
		
   	 	double DeltaR = sigmaij.trace()-f;
   	 	//std::cout << DeltaR <<std::endl; 
        while (DeltaR>= 1.e-2)
	   {
           double Pcbarj=Pcbar+0.1;

	       Matrix3d Dep1       = (dfdsigma*Dp*dfdsigma.transpose()*Dp);

	       Matrix3d Dep2       = (((-dfdp)*Pcbarj*Matrix3d::Identity())+dfdsigma*Dp*dfdsigma.transpose());

           Matrix3d Dep2i 		= Dep2.inverse();
        
	       Matrix3d Dep        = Dep1*Dep2i;

	       Matrix3d sigmaij    = Stress - Dep*de;

           DeltaR = sigmaij.trace()-f;

           //return DeltaR
        }
	    




		 
		
		//Matrix3d  Eep       =((dfdsigma*Dp)/((-dfdp)*Pcbar+dfdsigma*Dp*dfdsigma.transpose()))*((1+En)/(CCLamda-kappa)*dfdp*Pc);

		//double    Pc         = Pc + Eep*de; 

		//double Loadingindex = 

		Matrix3d sp = Matrix3d::Zero();

		sp = sigmaij;

		Stress = v0 * sp * v0.inverse();

       // std::cout <<"matrix： "<< Stress <<std::endl; 

	}
}