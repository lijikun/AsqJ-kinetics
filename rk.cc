// Compile: g++ -std=c++11 -O3 rk.cc -o rk
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#define NUM_OF_VARS 8
#define NUM_OF_PARAMS 7
#define MAX_TIME 10.0
#define RECORD_EVERY 100

void odeint(const std::vector<double>& k, std::vector<double>& x, const double& dt, std::fstream& results);
void dxdt(const std::vector<double>& k, const std::vector<double>& x, std::vector<double>& d);

int main(int argc, char* argv[]){
	std::fstream parameters;
	std::fstream results;
	std::vector<double> x0(NUM_OF_VARS);
	std::vector<double> k(NUM_OF_PARAMS);
	double dt;

	for (auto i = 1; i != argc; ++i){
		std::string paramFileName(argv[i]);
		std::string resultsFileName;
		resultsFileName = paramFileName + ".dat";
		parameters.open(paramFileName, std::ios::in);
		results.open(resultsFileName, std::ios::out);
		if (parameters.is_open() && results.is_open()) {
			std::cout << "Processing parameters file " << paramFileName 
				<< ", writing results file " << resultsFileName << "....\n";
			results << "# Time  Enzyme  aKG  subst  O2  ES  I  P\n";
			for (auto&& param : k)
				parameters >> param;
			for (auto&& initValue : x0)
				parameters >> initValue;
			parameters >> dt;
			parameters.close();
			results << 0.0 << '\t';
			for (auto&& initValue : x0)
				results << initValue << '\t';
			results << '\n';
			odeint(k, x0, dt, results);
			results.close();
		}
		else 
			std::cout << "Error processing files. \n";
	}

	return 0;
}

void odeint(const std::vector<double>& k, std::vector<double>& x, const double& dt, std::fstream& results){
	static std::vector<double> x1(x);	// Temporary variables for intermediate RK4 derivatives
	static std::vector<double> d(NUM_OF_VARS);	// dx/dt
	static std::vector<double> dx1(NUM_OF_VARS);	
	static std::vector<double> dx2(NUM_OF_VARS);
	static std::vector<double> dx3(NUM_OF_VARS);
	static std::vector<double> dx4(NUM_OF_VARS);
	double time = 0.0;
	int counter = 0;

	while (time < MAX_TIME){
		time += dt;
		dxdt(k, x1, d);
		for (int i = 0; i != NUM_OF_VAR; ++i){
			dx1[i] = dt * d[i];
			x1[i] = x[i] + dx1[i] / 2.0;
		}
		dxdt(k, x1, d);
		for (int i = 0; i != NUM_OF_VAR; ++i){
			dx2[i] = dt * d[i];
			x1[i] = x[i] + dx2[i] / 2.0;
		}
		dxdt(k, x1, d);
		for (int i = 0; i != NUM_OF_VAR; ++i){
			dx3[i] = dt * d[i];
			x1[i] = x[i] + dx3[i];
		}
		dxdt(k, x1, d);
		for (int i = 0; i != NUM_OF_VAR; ++i){
			dx4[i] = dt * d[i];
			x[i] += ( dx1[i] / 6.0 + dx2[i] / 3.0 + dx3[i] / 3.0 + dx4[i] / 6.0 );
		}
		if (++ counter == RECORD_EVERY){
			results << time << '\t';
			for (auto&& val : x)	
				results << val << '\t';
			results << '\n';
			counter = 0;
		}
	}

	return;
}

void dxdt(const std::vector<double>& k, const std::vector<double>& x, std::vector<double>& d){
	static double rate1, rate2, rate3, rate4, rate5;
	rate1 = k[0] * x[0] * x[1] - k[1] * x[4];
	rate2 = k[2] * x[4] * x[2] - k[3] * x[5];
	rate3 = k[4] * x[3] * x[5];
	rate4 = k[5] * x[6];
	rate5 = k[6] * x[7];
	d[0] = rate5 - rate1;
	d[1] = - rate1;
	d[2] = - rate2;
	d[3] = - rate3;
	d[4] = rate1 - rate2;
	d[5] = rate2 - rate3;
	d[6] = rate3 - rate4;
	d[7] = rate4 - rate5;

	return;	
}

