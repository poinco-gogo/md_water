#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <ctime>
#include <cmath>
#include "vector3.h"
using namespace std;
double gauss();
class Atom
{	
	public:
	Vector3 position;
	Vector3 velocity;
	Vector3 fold, fnew, vnew, rold, rnew;
	string  atomname;
	double  charge, mass, R;
	Atom(double x, double y, double z)
	{
		position.x = x;
		position.y = y;
		position.z = z;
	}
};
class System
{
	public:
	int nstep, natom, kmax;
	double dt, T, gamma, boxsize, cutoff, volume, ewcoeff, lj, es;
	Vector3 origin, a1, a2, a3, g1, g2, g3, b1, b2, b3;

	System()
	{
		dt = 9999;
		T = 9999;
		gamma = 9999;
		nstep = 0;
		natom = 0;
		boxsize = 1;
		origin.x = 9999; origin.y = 9999; origin.z = 9999;
		cutoff = 1.;
		lj = 0;
		es = 0;
	}

	void setup_box()
	{
		volume = boxsize * boxsize * boxsize;
		a1.x = boxsize; a1.y = 0; a1.z = 0;
		a2.x = 0; a2.y = boxsize; a2.z = 0;
		a3.x = 0; a3.y = 0; a3.z = boxsize;
		double rvol2pi = 2. * M_PI / volume;
		g1 = (a2 % a3) * rvol2pi;
		g2 = (a3 % a1) * rvol2pi;
		g3 = (a1 % a2) * rvol2pi;
	}
};
void print_ene(int istep, System& system, double K, int nfree);
bool shake(vector<Atom>& atomVector, vector<int>& shake_list);
void calc_frc(vector<Atom>& atomVector, vector<int>& lj_pair_list, vector<int>& el_pair_list, System& system, vector<Vector3>& g);
double calc_kin(vector<Atom>& atomVector);
void calc_pot(vector<Atom>& atomVector, vector<int>& lj_pair_list, vector<int>& el_pair_list, System& system, vector<Vector3>& g);
void output(ofstream& fo, vector<Atom>& atomVector, System& system);
void make_lj_pair(vector<Atom>& atomVector, vector<int>& lj_pair_list);
void make_el_pair(vector<Atom>& atomVector, vector<int>& el_pair_list);
void make_shake_pair(vector<Atom>& atomVector, vector<int>& shake_list);
bool load_config(ifstream& fi, System& system);
bool load_psf(ifstream& fs, System& system, vector<Atom>& atomVector);
bool load_pdb(ifstream& fp, vector<Atom>& atomVector);
int main (int argc, char** argv)
{
	if (argc < 4)
	{
		cout << "usage: ./a.out trjout psf pdb config\n";
		return 1;
	}
	
	srand( (unsigned int)time(NULL) );
	
	ofstream fo(argv[1]);
	ifstream fs(argv[2]);
	ifstream fp(argv[3]);
	ifstream fc(argv[4]);
	System system;
	if (!load_config(fc, system))
		return 1;
	vector<Atom> atomVector;
	if (!load_psf(fs, system, atomVector))
		return 1;

	const int natom = system.natom;
	const int nwat  = natom / 3;
	const int nfree = natom * 3 - nwat * 3;

	if (!load_pdb(fp, atomVector))
		return 1;

	int print_energy_step = 1;
	int print_trj_step  = 1;
	int nstep = system.nstep;

	const double ps2asu = 20.455;
        const double dt_ps = system.dt * 0.001; // in ps
        const double gamma_ps = system.gamma; // in ps-1
        const double dt = dt_ps * ps2asu;
        const double gamma = gamma_ps / ps2asu;
	const double dt_div2 = dt * 0.5;
	const double dtdt_div2 = dt * dt * 0.5;
	const double kb = 0.001987191; // kcal/mol/K
	const double  T = system.T; // K
        const double A = 1. - gamma * dt * 0.5;
        const double B = 1. + gamma * dt * 0.5;
        const double inB = 1. / B;
	
	cout << "REMARK Number of atoms " << natom << '\n';
	for (int i = 0; i < natom; i++)
	{	
		Atom& at = atomVector[i];
		at.R = sqrt(2. * kb * T * gamma * at.mass / dt);
	}

	cout << "REMARK Number of water molecules " << nwat << '\n';
	cout << "REMARK Degrees of freedom " << nfree << '\n';
	cout << "REMARK dt[fs] " << dt_ps * 1e3 << '\n';
	cout << "REMARK gamma[ps-1] " << gamma_ps << '\n';
	cout << "REMARK T[K] " << T << '\n';
	cout << "REMARK Cutoff[ang.] " << system.cutoff << '\n';

	const double rOH  = 0.9572;
	const double aHOH = 104.52 / 180. * acos(-1.0);

	// make reciprocal vectors
	vector<Vector3> g;
	int kmax = system.kmax;
	int sqkmax = kmax - 1;
	sqkmax *= sqkmax;
	double dum = 0;
	for (int k = 0; k < kmax; k++)
        for (int i = -kmax + 1; i < kmax; i++)
        for (int j = -kmax + 1; j < kmax; j++)
        {
                if (k * k + j * j + i * i > sqkmax or (i==0 && j==0 && k==0))
                        continue;
                Vector3 vtmp(i * system.g1.x, j * system.g2.y, k * system.g3.z);
                if (vabs(vtmp) > dum)
                        dum = vabs(vtmp);
                g.push_back(vtmp);
        }
        cout << "REMARK gmax = " << dum << '\n';

	double ew_self = 0;
	for (int i = 0; i < atomVector.size(); i++)
	{	
		ew_self += atomVector[i].charge * atomVector[i].charge;
	}
	ew_self *= -system.ewcoeff / sqrt(M_PI) * 332.0636;


	int icnt = 0;
	for (int i = 0; i < natom; i++)
	{
		Atom& at = atomVector[i];
		at.velocity.x = gauss() * sqrt(kb * T / at.mass);
		at.velocity.y = gauss() * sqrt(kb * T / at.mass);
		at.velocity.z = gauss() * sqrt(kb * T / at.mass);
	}
	
	output(fo, atomVector, system);

	vector<int> lj_pair_list, el_pair_list, shake_list;
	make_lj_pair(atomVector, lj_pair_list);
	make_el_pair(atomVector, el_pair_list);
	make_shake_pair(atomVector, shake_list);

	calc_pot(atomVector, lj_pair_list, el_pair_list, system, g);
	double Ktmp = calc_kin(atomVector);
	system.es += ew_self;
	print_ene(0, system, Ktmp, nfree);

	calc_frc(atomVector, lj_pair_list, el_pair_list, system, g);
	for (int i = 0; i < atomVector.size(); i++)
	{
		Atom& at = atomVector[i];
		at.fold = at.fnew;
		at.rold = at.position;
		at.position = at.rold + dt * at.velocity + dtdt_div2 / at.mass * at.fold;
	}

	double accum = 0;

	for (int istep = 1; istep <= nstep; istep++)
	{
		calc_frc(atomVector, lj_pair_list, el_pair_list, system, g);
		for (int i = 0; i < atomVector.size(); i++)
		{
			Vector3 noise( gauss(), gauss(), gauss());
			Atom& at = atomVector[i];
			at.rnew = 2. * at.position - at.rold + gamma * dt_div2 * at.rold + dt * dt / at.mass * (at.fnew + at.R * noise);
			at.rnew = at.rnew * inB;
		}

		for (int ishake = 0; ishake < 100; ishake++)
		{
			if (shake(atomVector, shake_list))
				break;
			if (ishake == 99)
			{
				cerr << "error: shake does not converged.\n";
				return 1;
			}
		}

		for (int i = 0; i < atomVector.size(); i++)
		{
			Atom& at = atomVector[i];
			at.vnew = 0.5 / dt * (at.rnew - at.rold);
		}

		double K = calc_kin(atomVector);
		calc_pot(atomVector, lj_pair_list, el_pair_list, system, g);
		system.es += ew_self;
		if (istep % print_energy_step== 0)
			print_ene(istep, system, K, nfree);

		if (istep % print_trj_step== 0)
		{
			output(fo, atomVector, system);
		}

		for (int i = 0; i < atomVector.size(); i++)
		{
			Atom& at = atomVector[i];
			at.rold = at.position;
			at.position = at.rnew;
			at.velocity = at.vnew;
		}
	}
}
/////////////////////////  end of main program


void print_ene(int istep, System& system, double K, int nfree)
{
	static double kb = 0.001987191; //kcal/mol/K
	double totpot = system.es + system.lj;
	cout 
		<< setprecision(4) << fixed
		<< setw(12) << istep
		<< setw(16) << system.lj
		<< setw(16) << system.es
		<< setw(16) << totpot
		<< setw(16) << K
		<< setw(16) << totpot + K
		<< setw(16) << K * 2. / nfree / kb
		<< '\n';
}

bool shake(vector<Atom>& atomVector, vector<int>& shake_list)
{
	static const double eps = 1e-6;
	static const double eps2 = eps * eps;
	static const double rOH = 0.9572;
	static const double aHOH = 104.52 / 180. * acos(-1.0);
	static const double dOH = rOH * rOH;
	static const double rHH = rOH * sin(aHOH/2.) * 2.;
	static const double dHH = rHH * rHH;
	for (int i = 0; i < shake_list.size() / 2; i++)
	{
		Atom& at1 = atomVector[shake_list[2 * i]];
		Atom& at2 = atomVector[shake_list[2 * i + 1]];
		double gamma;
		if (at1.atomname == "O" && at2.atomname == "H")
		{
			gamma = (dOH - norm(at1.rnew - at2.rnew)) / 
				(2.*(1./at1.mass+1./at2.mass)*((at1.position - at2.position) * (at1.rnew - at2.rnew)));
		}
		else
		{
			gamma = (dHH - norm(at1.rnew - at2.rnew)) / 
				(2.*(1./at1.mass+1./at2.mass)*((at1.position - at2.position) * (at1.rnew - at2.rnew)));

		}
		at1.rnew = at1.rnew + gamma * (at1.position - at2.position) / at1.mass;
		at2.rnew = at2.rnew + gamma * (at2.position - at1.position) / at2.mass;
	}

	for (int i = 0; i < shake_list.size() / 2; i++)
	{
		Atom& at1 = atomVector[shake_list[2 * i]];
		Atom& at2 = atomVector[shake_list[2 * i + 1]];
		double r = vabs(at1.rnew - at2.rnew);
		double error;
		if (at1.atomname  == "O" && at2.atomname == "H")
		{
			error = abs(r - rOH);
		}
		else
		{
			error = abs(r - rHH);
		}
		if (error > eps) return false;
	}

	return true;
}

double calc_kin(vector<Atom>& atomVector)
{
	double k = 0;
	for (int i = 0; i < atomVector.size(); i++)
	{
		Atom& at = atomVector[i];
		k += at.mass * norm(at.vnew);
	}
	return k * 0.5;
}

void calc_frc(vector<Atom>& atomVector, vector<int>& lj_pair_list, vector<int>& el_pair_list, System& system, vector<Vector3>& g)
{
	static double cutoff   = system.cutoff;
	static double cutoff2  = cutoff * cutoff;
	static double ewcoeff  = system.ewcoeff;
	static double ewcoeff2 = ewcoeff * ewcoeff;
	static double boxsize  = system.boxsize;
	static double factor = 1. / (4. * ewcoeff * ewcoeff);
	static double const_intra   = -2. * ewcoeff / sqrt(acos(-1.));
	static double const_recipro = 4. * acos(-1.) / system.volume;
	static double A = 582. * 1e3;
	static double A_2 = A * 2.;
	static double B = 595.0;
	static double C = 332.0636;
	for (int i = 0; i < atomVector.size(); i++)
	{
		atomVector[i].fnew.x = 0.;
		atomVector[i].fnew.y = 0.;
		atomVector[i].fnew.z = 0.;
	}

	for (int i = 0; i < lj_pair_list.size(); i++)
	{
		Atom& at1 = atomVector[lj_pair_list[i]];

		for (int j = i + 1; j < lj_pair_list.size(); j++)
		{
			Atom& at2 = atomVector[lj_pair_list[j]];
			Vector3 del = at1.position - at2.position;
			del.x -= boxsize * floor(del.x / boxsize + 0.5);
			del.y -= boxsize * floor(del.y / boxsize + 0.5);
			del.z -= boxsize * floor(del.z / boxsize + 0.5);
			double r2 = norm(del);
			if (r2 > cutoff2) continue;
			double r6 = r2 * r2 * r2;
			double r8 = r6 * r2;
			Vector3 f12 = 6. / r8 * (B - A_2 / r6) * del;
			at1.fnew -= f12;
			at2.fnew += f12;
		}
	//	if (!i) print(f12);
	}
//	print(atomVector[0].fnew);

	// ewald intra force
	for (int i = 0; i < atomVector.size() / 3; i++)
	{
		Vector3 frc(0., 0., 0.);
		int j = 3 * i;
		Atom& at1 = atomVector[j];
		Atom& at2 = atomVector[j + 1];
		Atom& at3 = atomVector[j + 2];
		Vector3 del1 = at1.position - at2.position;
		del1.x -= boxsize * floor(del1.x / boxsize + 0.5);
		del1.y -= boxsize * floor(del1.y / boxsize + 0.5);
		del1.z -= boxsize * floor(del1.z / boxsize + 0.5);
		double adel1 = vabs(del1);
		double sqadel1 = adel1 * adel1;
		frc = at1.charge * at2.charge * C * 
			(const_intra * exp(-ewcoeff2 * sqadel1) 
			 + erfl(ewcoeff * adel1) / adel1) / sqadel1 * del1;
		at1.fnew -= frc;
		at2.fnew += frc;
		Vector3 del2 = at1.position - at3.position;
		del2.x -= boxsize * floor(del2.x / boxsize + 0.5);
		del2.y -= boxsize * floor(del2.y / boxsize + 0.5);
		del2.z -= boxsize * floor(del2.z / boxsize + 0.5);
		double adel2 = vabs(del2);
		double sqadel2 = adel2 * adel2;
		frc = at1.charge * at3.charge * C * 
			(const_intra * exp(-ewcoeff2 * sqadel2) 
			 + erfl(ewcoeff * adel2) / adel2) / sqadel2 * del2;
		at1.fnew -= frc;
		at3.fnew += frc;
		Vector3 del3 = at3.position - at2.position;
		del3.x -= boxsize * floor(del3.x / boxsize + 0.5);
		del3.y -= boxsize * floor(del3.y / boxsize + 0.5);
		del3.z -= boxsize * floor(del3.z / boxsize + 0.5);
		double adel3 = vabs(del3);
		double sqadel3 = adel3 * adel3;
		frc = at3.charge * at2.charge * C * 
			(const_intra * exp(-ewcoeff2 * sqadel3) 
			 + erfl(ewcoeff * adel3) / adel3) / sqadel3 * del3;
		at3.fnew -= frc;
		at2.fnew += frc;
	}


	// ewald direct space force summation
	for (int i = 0; i < el_pair_list.size() / 2; i++)
	{
		Atom& at1 = atomVector[el_pair_list[2 * i]];
		Atom& at2 = atomVector[el_pair_list[2 * i + 1]];

		Vector3 del = at1.position - at2.position;
		del.x -= boxsize * floor(del.x / boxsize + 0.5);
		del.y -= boxsize * floor(del.y / boxsize + 0.5);
		del.z -= boxsize * floor(del.z / boxsize + 0.5);
		double r2 = norm(del);
		if (r2 > cutoff2) continue;
		double r = sqrt(r2);
		Vector3 frc(0., 0., 0.);
		frc = at1.charge * at2.charge * C *
			(-const_intra * exp(-ewcoeff2 * r2) 
			+ erfc(ewcoeff*r) / r) / r2 *del;
		at1.fnew += frc;
		at2.fnew -= frc;
	}
	
	
	// ewald reciprocal space summation
	for (int ii = 0; ii < atomVector.size(); ii++)
	{
		Atom& iat = atomVector[ii];
		Vector3 frc(0., 0., 0.);
	for (int i = 0; i < g.size(); i++)
	{
		double ag = vabs(g[i]);
		double agag = ag * ag;
		double pre, dtmp;
		pre = 0;
		for (int j = 0; j < atomVector.size(); j++)
		{
			Atom& jat = atomVector[j];
			Vector3 del = iat.position - jat.position;
			del.x -= boxsize * floor(del.x / boxsize + 0.5);
			del.y -= boxsize * floor(del.y / boxsize + 0.5);
			del.z -= boxsize * floor(del.z / boxsize + 0.5);
			double dot = g[i] * del;

			pre += jat.charge * sin(dot);
		}
		g[i].z ? dtmp = 1. : dtmp = 0.5;
		frc = frc + (dtmp * exp(-agag * factor) / agag * pre) * g[i];
	}
		iat.fnew += frc * const_recipro * C * iat.charge;
//		cout << vabs(frc) << '\n';
	}
//	print(atomVector[0].fnew);
}

void calc_pot(vector<Atom>& atomVector, vector<int>& lj_pair_list, vector<int>& el_pair_list, System& system, vector<Vector3>& g)
{
	static double cutoff = system.cutoff;
	static double cutoff2 = cutoff * cutoff;
	static double boxsize = system.boxsize;
	static double ewcoeff = system.ewcoeff;
	static double factor = 1. / (4. * ewcoeff * ewcoeff);
	static double const_recipro = 4. * acos(-1.) / system.volume;
	static double A = 582. * 1e3;
	static double B = 595.0;
	static double qO = -0.834;
	static double qH = 0.417;
	static double C  = 332.0636;

	system.lj = 0;
	system.es = 0;
	double ew_direct  = 0;
	double ew_recipro = 0;
	double ew_intra   = 0;
	double U = 0;

	for (int i = 0; i < lj_pair_list.size(); i++)
	{
		Atom& iat = atomVector[lj_pair_list[i]];
		for (int j = i + 1; j < lj_pair_list.size(); j++)
		{
			Atom& jat = atomVector[lj_pair_list[j]];
			Vector3 del = iat.position - jat.position;
			del.x -= boxsize * floor(del.x / boxsize + 0.5);
			del.y -= boxsize * floor(del.y / boxsize + 0.5);
			del.z -= boxsize * floor(del.z / boxsize + 0.5);
			double roo2 = norm(del);
			if (roo2 > cutoff2) continue;
			double roo6 = roo2 * roo2 * roo2;
			system.lj += A / (roo6 * roo6) - B / roo6;
		}
	}

	// ewald intra energy
	for (int i = 0; i < atomVector.size() / 3; i++)
	{
		int j = 3 * i;
		Atom& at1 = atomVector[j];
		Atom& at2 = atomVector[j + 1];
		Atom& at3 = atomVector[j + 2];
		Vector3 del1 = at1.position - at2.position;
		del1.x -= boxsize * floor(del1.x / boxsize + 0.5);
		del1.y -= boxsize * floor(del1.y / boxsize + 0.5);
		del1.z -= boxsize * floor(del1.z / boxsize + 0.5);
		double adel1 = vabs(del1);
		Vector3 del2 = at1.position - at3.position;
		del2.x -= boxsize * floor(del2.x / boxsize + 0.5);
		del2.y -= boxsize * floor(del2.y / boxsize + 0.5);
		del2.z -= boxsize * floor(del2.z / boxsize + 0.5);
		double adel2 = vabs(del2);
		Vector3 del3 = at3.position - at2.position;
		del3.x -= boxsize * floor(del3.x / boxsize + 0.5);
		del3.y -= boxsize * floor(del3.y / boxsize + 0.5);
		del3.z -= boxsize * floor(del3.z / boxsize + 0.5);
		double adel3 = vabs(del3);
		ew_intra += at1.charge*at2.charge*erfl(ewcoeff*adel1)/adel1;
		ew_intra += at1.charge*at3.charge*erfl(ewcoeff*adel2)/adel2;
		ew_intra += at3.charge*at2.charge*erfl(ewcoeff*adel3)/adel3;
	}

	// ewald direct space summation
	for (int i = 0; i < el_pair_list.size() / 2; i++)
	{
		Atom& at1 = atomVector[el_pair_list[2 * i]];
		Atom& at2 = atomVector[el_pair_list[2 * i + 1]];

		Vector3 del = at1.position - at2.position;
		del.x -= boxsize * floor(del.x / boxsize + 0.5);
		del.y -= boxsize * floor(del.y / boxsize + 0.5);
		del.z -= boxsize * floor(del.z / boxsize + 0.5);
		double r = vabs(del);
		if (r > cutoff) continue;
		ew_direct += at1.charge * at2.charge * erfc(ewcoeff*r)/r;
	}

	// ewald reciprocal space summation
	for (int i = 0; i < g.size(); i++)
	{
		double ag = vabs(g[i]);
		double agag = ag * ag;
		double re, im, dtmp;
		re = 0;
		im = 0;
		for (int j = 0; j < atomVector.size(); j++)
		{
			Vector3 del = atomVector[j].position;
			del.x -= boxsize * floor(del.x / boxsize + 0.5);
			del.y -= boxsize * floor(del.y / boxsize + 0.5);
			del.z -= boxsize * floor(del.z / boxsize + 0.5);
			double dot = g[i] * del;

			re += atomVector[j].charge * cos(dot);
			im += atomVector[j].charge * sin(dot);
		}
		g[i].z ? dtmp = 1. : dtmp = 0.5;
		ew_recipro += dtmp * exp(-agag * factor) / agag * (re * re + im * im);
	}

	ew_recipro *= const_recipro;

//	cout << "el_direct  " << el_direct * 332.0636 << '\n';
//	cout << "el_intra    " << el_intra * 332.0636 << '\n';
//	cout << "el_recipro " << el_recipro * 332.0636<< '\n';

	system.es += ew_direct + ew_recipro - ew_intra;
	system.es *= C;
}

void output(ofstream& fo, vector<Atom>& atomVector, System& system)
{
	static const double boxsize = system.boxsize;
	fo << atomVector.size() << '\n' << '\n';
	fo << setprecision(6) << fixed;
	for (int i = 0; i < atomVector.size() / 3; i++) 
	{
		int j = 3 * i;
		Atom& at1 = atomVector[j];
		Atom& at2 = atomVector[j + 1];
		Atom& at3 = atomVector[j + 2];
		Vector3 p1 = at1.position;
		Vector3 p2 = at2.position;
		Vector3 p3 = at3.position;
		Vector3 vcom = (p1 + p2 + p3) / 3.;

		double d = floor(0.5 + vcom.x / boxsize);
		p1.x -= boxsize * d;
		p2.x -= boxsize * d;
		p3.x -= boxsize * d;
		
		d = floor(0.5 + vcom.y / boxsize);
		p1.y -= boxsize * d;
		p2.y -= boxsize * d;
		p3.y -= boxsize * d;
		
		d = floor(0.5 + vcom.z / boxsize);
		p1.z -= boxsize * d;
		p2.z -= boxsize * d;
		p3.z -= boxsize * d;

		fo << " " << at1.atomname << setw(20) << p1.x << setw(20) << p1.y << setw(20) << p1.z << '\n';
		fo << " " << at2.atomname << setw(20) << p2.x << setw(20) << p2.y << setw(20) << p2.z << '\n';
		fo << " " << at3.atomname << setw(20) << p3.x << setw(20) << p3.y << setw(20) << p3.z << '\n';
	}
}

double gauss()
{
	const double A1 = 3.949846138, A3 = 0.252408784;
	const double A5 = 0.076542912, A7 = 0.008355968;
	const double A9 = 0.029899776;
	double sum = 0.;
	for (int i = 0; i < 12; i++)
	{
		sum += rand() / static_cast<double>(RAND_MAX);
	}

	double R = (sum - 6.) / 4.;
	double R2 = R * R;
	return ((((A9 * R2 + A7) * R2 + A5) * R2 + A3) * R2 + A1) * R;
}

void make_lj_pair(vector<Atom>& atomVector, vector<int>& lj_pair_list)
{
	for (int i = 0; i < atomVector.size(); i++)
	{
		if (atomVector[i].atomname == "O")
			lj_pair_list.push_back(i);
	}
}

void make_el_pair(vector<Atom>& atomVector, vector<int>& el_pair_list)
{
	for (int i = 0; i < atomVector.size(); i++)
	{
		Atom& at1 = atomVector[i];

		for (int j = i + 1; j < atomVector.size(); j++)
		{
			Atom& at2 = atomVector[j];

			if (i / 3 == j / 3) continue;

			el_pair_list.push_back(i);
			el_pair_list.push_back(j);
		}
	}
}

void make_shake_pair(vector<Atom>& atomVector, vector<int>& shake_list)
{
	for (int i = 0; i < atomVector.size() / 3; i++)
	{
		int j = 3 * i;
		shake_list.push_back(j);
		shake_list.push_back(j+1);
		shake_list.push_back(j);
		shake_list.push_back(j+2);
		shake_list.push_back(j+1);
		shake_list.push_back(j+2);
	}
}

bool load_config(ifstream& fc, System& system)
{
	string s;
	while (getline(fc, s))
	{
		istringstream is(s);
		string stmp;
		if (s.find("#", 0) != string::npos)
		{
			// found comment line.
			continue;
		}
		else if (s.find("dt", 0) != string::npos)
		{
			is >> stmp >> system.dt;
		}
		else if (s.find("temperature", 0) != string::npos)
		{
			is >> stmp >> system.T;
		}
		else if (s.find("gamma", 0) != string::npos)
		{
			is >> stmp >> system.gamma;
		}
		else if (s.find("numstep", 0) != string::npos)
		{
			is >> stmp >> system.nstep;
		}
		else if (s.find("boxsize", 0) != string::npos)
		{
			is >> stmp >> system.boxsize;
			system.setup_box();
		}
		else if (s.find("origin", 0) != string::npos)
		{
			is >> stmp >> system.origin.x
				>> system.origin.y
				>> system.origin.z;
		}
		else if (s.find("cutoff", 0) != string::npos)
		{
			is >> stmp >> system.cutoff;
		}
		else if (s.find("ewcoeff", 0) != string::npos)
		{
			is >> stmp >> system.ewcoeff;
		}
		else if (s.find("kmax", 0) != string::npos)
		{
			is >> stmp >> system.kmax;
		}
		else 
		{
			is >> stmp;
			cerr << "unknown config param \"" << stmp
				<< "\" found.\n";
			return false;
		}
	}

	return true;
}

bool load_psf(ifstream& fs, System& system, vector<Atom>& atomVector)
{
	string s;
	while (getline(fs, s))
	{
		if (s.find("NATOM", 10) != string::npos)
		{
			istringstream is(s);
			is >> system.natom;
			for (int i = 0; i < system.natom; i++)
			{
				Atom at(0., 0., 0.);
				getline(fs, s);
				istringstream iss(s);
				int itmp;
				string stmp, atname;
				iss >> itmp >> stmp >> itmp >> stmp >> atname
					>> stmp >> at.charge >> at.mass;
				//cout << atname << '\n';
				//cout << s << '\n';
				at.atomname = atname.substr(0, 1);
				atomVector.push_back(at);
			}
			return true;
		}

	}
	return false;
}

bool load_pdb(ifstream& fp, vector<Atom>& atomVector)
{
	string s;
	int icnt = 0;
	while (getline(fp, s))
	{
		if (s.find("ATOM", 0) != string::npos)
		{
			Atom& at = atomVector[icnt++];
			istringstream is(s.substr(30, 24));
			is >> at.position.x >> at.position.y >> at.position.z;
		}
	}
	return atomVector.size() == icnt ? true : false;
}
