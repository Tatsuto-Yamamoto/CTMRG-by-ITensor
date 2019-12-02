#include "itensor/all.h"
#include <chrono>
using namespace itensor;
using namespace std;

double Ent(ITensor, Index, Index);
double log_partition_cal(ITensor, const int &, vector<double> &, vector<double> &);
double energy_cal(ITensor, Index, Index, Index, Index, const double &);

const int dim0 = 2;// degree of freedom for Ising model

void calc(vector<double> & thermo,double & beta,const int & N, const int & maxdim, const double & J){
	vector<double> CTM_normalize(N);// normalization factor of CTM (This is needed for calulate a free enrgy)
	vector<double> Hrtm_normalize(N);// normalization factor of HRTM

	Index a(dim0, "a"), b(dim0, "b"), c(dim0, "c"), d(dim0, "d"), i(dim0, "i"), j(dim0, "j"), k(dim0, "k"), l(dim0, "l"), m(dim0, "m");

	//Create W (Boltzmann Weight), (if you want to use this code for another model, you may change the definition of W.)
	auto W = ITensor(a, b, c, d);
	auto Sig = [](int s) { return 1.-2.*(s-1); };
	for(auto sa : range1(dim0))
	for(auto sb : range1(dim0))
	for(auto sc : range1(dim0))
	for(auto sd : range1(dim0)){
	    auto E = J*(Sig(sa)*Sig(sb)+Sig(sb)*Sig(sc)
	            +Sig(sc)*Sig(sd)+Sig(sd)*Sig(sa));
	    auto P = exp(-E*beta);
	    W.set(a(sa),b(sb),c(sc),d(sd),P);
	}

	// initial (CTM and HRTM), (if you change the boundary condition, you may rewrite this.)
	auto C = W*delta(a,i)*delta(b,j)*ITensor(c, d).fill(1.0); // i, j
	auto Hrtm = ITensor(c).fill(1.0)*W*delta(a,k)*delta(b,l)*delta(d,m); // k, l, m
	CTM_normalize[0] = norm(C);
	Hrtm_normalize[0] = norm(Hrtm);
	C /= CTM_normalize[0];
	Hrtm /= Hrtm_normalize[0];

	int N_temp = 1;// Current CTM size 

	// CTMRG loop
	for(auto count : range1(N-1)){
		// expand CTM
		C *= Hrtm*delta(i,m);
		C *= prime(Hrtm)*delta(j,prime(l));
		C *= W*delta(c,k)*delta(d,prime(k));
		C = C*delta(l,i)*delta(prime(m),j);//a,b,i,j
		/***************************/
			// calculate energy
			thermo[2] = energy_cal(C,a,b,i,j,J);
		/***************************/
		auto [Comb, comb] = combiner(b,i);
		C *= Comb;
		auto i_temp = comb;
		auto [Comb2,comb2] = combiner(a,j);
		C *= Comb2;
		auto j_temp = comb2;//i_temp, j_temp

		N_temp += 1;

		// create Density Matrix
		auto rho = C*delta(j_temp,prime(i_temp));
		rho *= prime(C)*delta(prime(j_temp),prime(i_temp,2));
		rho *= prime(C,2)*delta(prime(j_temp,2),prime(i_temp,3));
		rho *= prime(C,3)*delta(prime(j_temp,3),prime(i_temp)); // i_temp, i_temp'
		auto [U,D] = diagPosSemiDef(rho,{"MaxDim=",maxdim,"Truncate=",true, "respectDegenerate=",true});
		auto u = commonIndex(U,D);

		// renormalize CTM
		C *= U;
		C *= prime(U)*delta(j_temp,prime(i_temp));// u, u'
		CTM_normalize[count] = norm(C);
		C /= CTM_normalize[count];

		//calculate entropy
		thermo[0] = Ent(D, u, prime(u));

		// expand Hrtm
		Hrtm = W*delta(c,k)*Hrtm;
		Hrtm *=delta(a,k); // k,b,d,l,m
		auto [Comb3,comb3] = combiner(b,l);
		Hrtm *= Comb3;
		auto l_temp = comb3;
		auto [Comb4,comb4] = combiner(d,m);
		Hrtm *= Comb4;
		auto m_temp = comb4;

		// renormalize Hrtm
		Hrtm *= U*delta(l_temp,i_temp);
		Hrtm *= prime(U)*delta(m_temp,prime(i_temp));//u, u', k
		Hrtm_normalize[count] = norm(Hrtm);
		Hrtm /= Hrtm_normalize[count];

		// calulate a logarithm of partition function
		thermo[1] = log_partition_cal(C,N_temp,CTM_normalize,Hrtm_normalize);//

		// replace indices
		auto new_dim = dim(u);
		auto i_new = Index(new_dim, "i");
		auto j_new = Index(new_dim, "j");
		auto l_new = Index(new_dim, "l");
		auto m_new = Index(new_dim, "m");
		C = delta(prime(u),j_new)*C*delta(u,i_new);
		Hrtm = delta(prime(u),m_new)*Hrtm*delta(u,l_new);
		i = i_new;
		j = j_new;
		l = l_new;
		m = m_new;
	}
}

double Ent(ITensor D, Index a, Index b){
	auto sum = 0.0;
	auto n = dim(a);
	for(auto i : range1(n)){
		sum+=elt(D, a=i, b=i);
	}
	auto norm_D = D/sum;
	auto ent = 0.0;
	for(auto i : range1(n)){
		auto val = elt(norm_D, a=i, b=i);
		if(val>0) ent-=val*log(val);
	}
	return ent;
}

double log_partition_cal(ITensor C,const int & N,vector<double> & CTM_normalize,vector<double> & Hrtm_normalize){//自由エネルギー。繰り込まれた後のCを使う
	auto C2 = C;
	C2 *= prime(C);
	double Z_nol = elt(C2*C2);//Z_nol = Tr(C^4)

	double sum_C = 0;
	double sum_Hrtm = 0;
	for (int i = 0; i < N ; i++){
		sum_C += log(CTM_normalize[i]);
	}
	for (int i = 0; i < N-1 ; i++){
		sum_Hrtm += (N-(i+1))*log(Hrtm_normalize[i]);
	}
	double log_Z = log(Z_nol) + 4*sum_C + 8*sum_Hrtm;

	return log_Z;
}

double energy_cal(ITensor C,Index a,Index b,Index i,Index j, const double & J){
//      a 	  j
//      +    +++
//  	|     |
// b  +-*-----*
//      |     |
//    + |     |
//  i +-*-----+
//    +
//
// C(a,b,i,j)

	auto C3 = C;// i,j,a,b
	C3 *= prime(C)*delta(prime(j),i)*delta(prime(a),b);// i,i',b,b'
	C3 *= C*delta(j,prime(i))*delta(a,prime(b));// i,j,a,b
	double z = elt(C3*C);// partition function

	double energy = 0;
	auto Sig = [](int s) { return 1.-2.*(s-1); };
	for(auto sa : range1(dim(a)))
	for(auto sb : range1(dim(b)))
	for(auto si : range1(dim(i)))
	for(auto sj : range1(dim(j))) {
	    auto E = J*Sig(sa)*Sig(sb) * elt(C3,a=sa,b=sb,i=si,j=sj) * elt(C,a=sa,b=sb,i=si,j=sj);
	    energy += E;
	}
	energy = 2*energy / z;

	return energy;
}