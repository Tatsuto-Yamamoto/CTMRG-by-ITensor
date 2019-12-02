#include <time.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

using namespace std;
const double J = 1.0;// coupling constant(J=1 ferro, J=-1 anti-ferro)

void calc(vector<double> &, double &, const int &,const int &,const double &);

int main(){
	/***********************************************/
	//Call initial condition from file
	/***********************************************/
	FILE* par;
	par = fopen("initial.txt", "r");
	if ( par == NULL ){
		printf("Can't open the input file \"initial.txt\" \n");
		exit(1);
	}
	int N_max_r;
	int dim_r;
	double T_i_r;
	double T_f_r;
	double T_step_r;
	double small_beta_r; 

	fscanf(par, "%d%*[^\n]", &N_max_r);
	fscanf(par, "%d%*[^\n]", &dim_r);
	fscanf(par, "%lf%*[^\n]", &T_i_r);
	fscanf(par, "%lf%*[^\n]", &T_f_r);
	fscanf(par, "%lf%*[^\n]", &T_step_r);
	fscanf(par, "%lf%*[^\n]", &small_beta_r);
	fclose(par);

	const int N = N_max_r;// system size
	const int dim = dim_r;// bond dimension
	const double T_i = T_i_r;// initial temperature
	const double T_f = T_f_r;// final temperature
	const double T_step = T_step_r;// temperature step
	const double small_beta = small_beta_r;// Minute inverse-temperature for differentiation

	if(T_step <= 0){
		cout << "T_step <= 0" << endl;
		exit(1);
	}
	/***********************************************/

	/***********************************************/
	//Output the result to a file
	/***********************************************/
	char filename[100];
	sprintf(filename, "N=%d_dim=%d_Itensor.txt",N,dim);
	ofstream outputfile(filename);
	char value[100];
	sprintf(value, "#temperature\tentropy\t\t\tfree_energy\t\tenergy\t\t\tspecific_heat");
	outputfile << value <<endl;
	/***********************************************/

	clock_t start,end;
	double temperature = T_i;
	double beta;//inverse temperature
	cout << "start"<<endl;

	while(temperature < T_f){
		start = clock();
		cout << "temperature = " << temperature <<endl;
		beta = 1/temperature;

		const int number_of_thermo_quantity = 3;
		vector<double> thermo(number_of_thermo_quantity);//this stores thermodynamical quantities

		calc(thermo,beta,N,dim,J);
		double entropy = thermo[0];
		double log_Z = thermo[1];//partition function
		double free_energy = -temperature*log_Z/(N*N);
		double energy = thermo[2];

		// numerical differentiation
		vector<double> thermo_2(number_of_thermo_quantity);
		double beta_up = beta+small_beta;
		calc(thermo_2,beta_up,N,dim,J);
		double log_Z_up = thermo_2[1];

		double beta_low = beta-small_beta;
		calc(thermo_2,beta_low,N,dim,J);
		double log_Z_low = thermo_2[1];

		// specific heat
		double spe_heat = beta*beta*(log_Z_up - 2*log_Z + log_Z_low) / (small_beta*small_beta)/(N*N);

		// output to a file
		char thermodynamics[1000];
		sprintf(thermodynamics, "%1.6E\t%1.16E\t%1.16E\t%1.16E\t%1.16E",
			temperature,entropy,free_energy,energy,spe_heat);
		outputfile << thermodynamics <<endl;

		// next step
		temperature += T_step;
		end = clock();
		printf("takes %.2f second.\n",(double)(end-start)/CLOCKS_PER_SEC);
	}

	return 0;
}