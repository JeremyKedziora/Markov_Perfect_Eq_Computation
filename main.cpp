#include<nr3.h>
#include<iostream>
#include<cmath>
#include<nr3.h>
#include<quadrature.h>
#include<roots.h>
#include<gamma.h>
#include<gauss_wgts.h>
#include<vector>
#include<fstream>
#include<ctime>
#include<mins.h>
#include<Amoeba.h>
#include<Ran.h>
//#include<Scythe_Double_Matrix.h>
//#include<Scythe_Double_Matrix.cc>
//#include<Scythe_Optimize.h>
//#include<Scythe_Simulate.h>

using namespace std;


/*****************************************************************************************************************************
An attempt to compute the equilibrium of my model: note the index outside a function must be different from those used inside!
/*****************************************************************************************************************************
program variables*************************************************************************************************************
*****************************************************************************************************************************/
int i; //a generic index USE THIS ONLY IN THE MAIN PROGRAM
int j; //second generic index USE THIS ONLY IN THE MAIN PROGRAM
int h; //third generic index USE THIS ONLY IN THE SUB-PROGRAMS
int g; //fifth generic index USE THIS ONLY IN THE SUB-PROGRAMS
int t; //fourth generic index USE THIS ONLY IN THE SUB-PROGRAMS
int c; //sixth generic index USE THIS ONLY IN THE SUB-PROGRAMS
int bargain; //determines if we use the bargaining or some arbitrary peace payoffs
int rev; //determines if we allow for revolution after peace or not
int it; //an index for the number of iterations USE ONLY IN THE FIXED POINT ALGORITHM
int nconstraints; //a placeholder for the number of constraints
Doub tol; //a place-holder for the tolerance of the continuation values
VecDoub e(2); //a place-holder for the error between continuation values
VecDoub ind(51); //a simple indicator vector

/******************************
value functions and state space
******************************/
VecDoub K(51); //the state space
MatDoub V(51,2); //the matrix of continuation values; column 1 is the ruler
MatDoub Vn(51,2); //a place-holder matrix for newly computed continuation values
VecDoub VrP(51); //continuation Values for the Ruler for peace
VecDoub VsP(51); //continuation Values for the Subjects for peace
VecDoub VrRev(51); //continuation values for the Ruler for revolution
VecDoub VsRev(51); //continuation values for the subjects for revolution
MatDoub Constraints(51,2); //the actual constraint matrix

/*******************
The model parameters
*******************/
Doub lambda; //upper bound of the feasible set
Doub Step; //stepsize for the grid
Doub eta; //material costs for revolution
int T; //time to rebuild after revolution
Doub p; //political power of the Ruler after peace
Doub dr; //patience of the Ruler
Doub ds; //patience of the subjects
Doub psi; //value of holding office for our Ruler

/*******************
The choice variables
*******************/
VecDoub x(50); //mobilization
VecDoub pi(50); //peace
VecDoub r(50); //revolution

/***********************************************
The payoffs for peace and costs for mobilization
***********************************************/
//the settlement function, assumed exponential in k for now; note that k is the state NOT an index
Doub gam(int k){
	
	Doub gk;
	gk = pow(k,0.5);
	
	return gk;
}

Doub cost(const Doub x){
	return 0.25*(1/(1-x));
}

/**********************************************************************************************************************************
The transition function************************************************************************************************************
**********************************************************************************************************************************/
VecDoub q(51);//vector of transition probabilities on the state space

/********************************************************************
Michael's code for evaluating the binomial distribution: USES INDEX c
********************************************************************/
Doub binomial(int n, int k, Doub p){

	Doub logN = 0;
	int h = 1;
	int c;
	for(c=1;c<=k;c++){
		logN += log(h) - log(c);
		h++;
	}
	for(c=1;c<=n-k;c++){
		logN += log(h) - log(c);
		h++;
	}
	logN += k*log(p) + (n-k)*log(1-p);
	return exp(logN);
}

//the vector version of the transition function: USES INDEX t
VecDoub Q(const Doub x, const Doub k){
	
	Doub p = 1-pow((1/(x+0.03)),0.1)+k/100+0.4;
	//cout << p << endl << endl;
	VecDoub qnew(51);
	
	if(p == 0){
		qnew[0] = 1;
		for(t=1;t<=50;t++){
			q[t] = 0;
		}
	
	}else{
		if(p == 1.0){
			qnew[50] = 1;
			for(t=0;t<50;t++){
				qnew[t] = 0;
			}		
		}else{
			qnew[0] = binomial(50,0,p);
			for(t=1;t<=50;t++){
				qnew[t] = binomial(50,t,p);
			}
		}
	}
	return qnew;
}


//computes new continuation values as a function of the current strategy profile: USES INDEX g
MatDoub Vstar(const MatDoub V, const VecDoub pi, const VecDoub x, const VecDoub r){
	
	VecDoub Vr(51); //continuation Values for the Ruler
	VecDoub Vs(51); //continuation Values for the Subjects
	
	for(g=0;g<=50;g++){
		Vr[g] = V[g][0];
		Vs[g] = V[g][1];
	}
	
	//first we recompute the continuation values
	MatDoub Vnew(51,2);

	//Total defeat
	Vnew[0][0] = 0;
	Vnew[0][1] = 0;
	//Victory
	Vnew[50][0] = Vr[50];
	Vnew[50][1] = Vs[50];
		
	//need to compute everything else
	for(g=1;g<50;g++){
		
		q = Q(x[g],g+1);
		/*seems to work properly when I dont use the dot product program above... annoying...
		for(h=0;h<=20;h++){
			cout << q[h] << ";";
		}*/
		
		Vnew[g][0] = pi[g]*VrP[g]+(1-pi[g])*(r[g]*0+(1-r[g])*((1-dr)*psi+dr*(q[0]*Vr[0]+q[1]*Vr[1]+q[2]*Vr[2]+q[3]*Vr[3]+q[4]*Vr[4]+q[5]*Vr[5]+q[6]*Vr[6]+q[7]*Vr[7]+
		q[8]*Vr[8]+q[9]*Vr[9]+q[10]*Vr[10]+q[11]*Vr[11]+q[12]*Vr[12]+q[13]*Vr[13]+q[14]*Vr[14]+q[15]*Vr[15]+q[16]*Vr[16]+q[17]*Vr[17]+q[18]*Vr[18]+q[19]*Vr[19]+q[20]*Vr[20]+q[21]*Vr[21]
		+q[22]*Vr[22]+q[23]*Vr[23]+q[24]*Vr[24]+q[25]*Vr[25]+q[26]*Vr[26]+q[27]*Vr[27]+q[28]*Vr[28]+q[29]*Vr[29]+q[30]*Vr[30]+q[31]*Vr[31]+q[32]*Vr[32]+q[33]*Vr[33]+q[34]*Vr[34]
		+q[35]*Vr[35]+q[36]*Vr[36]+q[37]*Vr[37]+q[38]*Vr[38]+q[39]*Vr[39]+q[40]*Vr[40]+q[41]*Vr[41]+q[42]*Vr[42]+q[43]*Vr[43]+q[44]*Vr[44]+q[45]*Vr[45]+q[46]*Vr[46]+q[47]*Vr[47]
		+q[48]*Vr[48]+q[49]*Vr[49]+q[50]*Vr[50])));
		Vnew[g][1] = pi[g]*VsP[g]+(1-pi[g])*(r[g]*VsRev[g]+(1-r[g])*((1-ds)*(lambda-cost(x[g]))+ds*(q[0]*Vs[0]+q[1]*Vs[1]+q[2]*Vs[2]+q[3]*Vs[3]+q[4]*Vs[4]+q[5]*Vs[5]+q[6]*Vs[6]+q[7]*Vs[7]+
		q[8]*Vs[8]+q[9]*Vs[9]+q[10]*Vs[10]+q[11]*Vs[11]+q[12]*Vs[12]+q[13]*Vs[13]+q[14]*Vs[14]+q[15]*Vs[15]+q[16]*Vs[16]+q[17]*Vs[17]+q[18]*Vs[18]+q[19]*Vs[19]+q[20]*Vs[20]+q[21]*Vs[21]
		+q[22]*Vs[22]+q[23]*Vs[23]+q[24]*Vs[24]+q[25]*Vs[25]+q[26]*Vs[26]+q[27]*Vs[27]+q[28]*Vs[28]+q[29]*Vs[29]+q[30]*Vs[30]+q[31]*Vs[31]+q[32]*Vs[32]+q[33]*Vs[33]+q[34]*Vs[34]
		+q[35]*Vs[35]+q[36]*Vs[36]+q[37]*Vs[37]+q[38]*Vs[38]+q[39]*Vs[39]+q[40]*Vs[40]+q[41]*Vs[41]+q[42]*Vs[42]+q[43]*Vs[43]+q[44]*Vs[44]+q[45]*Vs[45]+q[46]*Vs[46]+q[47]*Vs[47]
		+q[48]*Vs[48]+q[49]*Vs[49]+q[50]*Vs[50])));
		
	}

	return Vnew;
}


//next I need to compute revolution strategies over the grid, along with mobilization and peace strategies
//first we compute the revolution strategies over the grid, states by mobilization
MatDoub rGrid(const MatDoub V, const int L, const VecDoub G){
	
	MatDoub rgrid(51,L+1);//a grid for all revolution strategies for any value of x
	VecDoub Vs(51); //continuation Values for the Subjects
	
	for(g=0;g<51;g++){
		Vs[g] = V[g][1];		
	}
	
	for(g=0;g<=L;g++){
		for(h=1;h<=49;h++){
			q = Q(G[g],h+1);
			if(max((1-ds)*(lambda-cost(G[g]))+ds*(q[0]*Vs[0]+q[1]*Vs[1]+q[2]*Vs[2]+q[3]*Vs[3]+q[4]*Vs[4]+q[5]*Vs[5]+q[6]*Vs[6]+q[7]*Vs[7]+
		q[8]*Vs[8]+q[9]*Vs[9]+q[10]*Vs[10]+q[11]*Vs[11]+q[12]*Vs[12]+q[13]*Vs[13]+q[14]*Vs[14]+q[15]*Vs[15]+q[16]*Vs[16]+q[17]*Vs[17]+q[18]*Vs[18]+q[19]*Vs[19]+q[20]*Vs[20]+q[21]*Vs[21]
		+q[22]*Vs[22]+q[23]*Vs[23]+q[24]*Vs[24]+q[25]*Vs[25]+q[26]*Vs[26]+q[27]*Vs[27]+q[28]*Vs[28]+q[29]*Vs[29]+q[30]*Vs[30]+q[31]*Vs[31]+q[32]*Vs[32]+q[33]*Vs[33]+q[34]*Vs[34]
		+q[35]*Vs[35]+q[36]*Vs[36]+q[37]*Vs[37]+q[38]*Vs[38]+q[39]*Vs[39]+q[40]*Vs[40]+q[41]*Vs[41]+q[42]*Vs[42]+q[43]*Vs[43]+q[44]*Vs[44]+q[45]*Vs[45]+q[46]*Vs[46]+q[47]*Vs[47]
		+q[48]*Vs[48]+q[49]*Vs[49]+q[50]*Vs[50]),VsRev[h]) == VsRev[h]){rgrid[h][g] = 1;}else{rgrid[h][g] = 0;}
			
			/*appears to be computing the revolution strategies on the grid correctly...
			cout << (1-ds)*(lambda-x[h])+ds*(q[0]*Vs[0]+q[1]*Vs[1]+q[2]*Vs[2]+q[3]*Vs[3]+q[4]*Vs[4]+q[5]*Vs[5]+q[6]*Vs[6]+q[7]*Vs[7]+q[8]*Vs[8]+q[9]*Vs[9]+q[10]*Vs[10]+
			q[11]*Vs[11]+q[12]*Vs[12]+q[13]*Vs[13]+q[14]*Vs[14]+q[15]*Vs[15]+q[16]*Vs[16]+q[17]*Vs[17]+q[18]*Vs[18]+q[19]*Vs[19]+q[20]*Vs[20]) 
			<< setw(12) << VsRev[h] << setw(12) << rgrid[h][g] << endl;*/
		}
	}
	return rgrid;
}

//next we need to compute the best x, so to speak, given the revolution strategies of the subjects
VecDoub xmap(const VecDoub Vr, const MatDoub rg, const VecDoub G, const int L){

	VecDoub xnew(50); //new mobilization offers for the Ruler
	VecDoub Vgrid(L+1); //new continuation values for the Ruler
	Doub ind;
	Doub xmax;

	for(h=1;h<=49;h++){
		for(g=0;g<=L;g++){
			q = Q(G[g],h+1);
			
			/*seems to be spitting out the right probabilities
			for(t=0;t<=20;t++){
				cout << q[t] << ",";
			}
			cout << endl;
			*/
			Vgrid[g] = (1-rg[h][g])*((1-dr)*psi+dr*(q[0]*Vr[0]+q[1]*Vr[1]+q[2]*Vr[2]+q[3]*Vr[3]+q[4]*Vr[4]+q[5]*Vr[5]+q[6]*Vr[6]+q[7]*Vr[7]+
		q[8]*Vr[8]+q[9]*Vr[9]+q[10]*Vr[10]+q[11]*Vr[11]+q[12]*Vr[12]+q[13]*Vr[13]+q[14]*Vr[14]+q[15]*Vr[15]+q[16]*Vr[16]+q[17]*Vr[17]+q[18]*Vr[18]+q[19]*Vr[19]+q[20]*Vr[20]+q[21]*Vr[21]
		+q[22]*Vr[22]+q[23]*Vr[23]+q[24]*Vr[24]+q[25]*Vr[25]+q[26]*Vr[26]+q[27]*Vr[27]+q[28]*Vr[28]+q[29]*Vr[29]+q[30]*Vr[30]+q[31]*Vr[31]+q[32]*Vr[32]+q[33]*Vr[33]+q[34]*Vr[34]
		+q[35]*Vr[35]+q[36]*Vr[36]+q[37]*Vr[37]+q[38]*Vr[38]+q[39]*Vr[39]+q[40]*Vr[40]+q[41]*Vr[41]+q[42]*Vr[42]+q[43]*Vr[43]+q[44]*Vr[44]+q[45]*Vr[45]+q[46]*Vr[46]+q[47]*Vr[47]
		+q[48]*Vr[48]+q[49]*Vr[49]+q[50]*Vr[50]));
			
			/*seems to be computing on the grid correctly
			cout << Vgrid[g] << ",";
			*/
		}
	
		xmax = Vgrid[0];
	
		for(g=1;g<=L;g++){
			if(Vgrid[g]>xmax){
				xmax = Vgrid[g];
				//cout << xmax << endl;
				xnew[h] = G[g];
			}else{xnew[h] = xnew[h];}
		}
		
		
		for(g=1;g<=L;g++){
			if(xnew[g] > lambda){xnew[g] = 0;}else{;}
			
		}
		
		//ind = maxindc(Vgrid);
		//xnew[h] = G[ind];
		//cout << endl;
	}
	return xnew;
}










int main(void){
	
	/****************************************************************************************************
	enter the model parameters and compute peace payoffs*************************************************
	****************************************************************************************************/
	cout << "Let's enter the model parameters..." << endl;
	cout << "Do you want to use bargaining after peace? 2 = yes but skip, 1 = yes, 0 = no.:";
	cin >> bargain;
	if(bargain == 2){
		dr = 0.9;
		ds = 0.9;
		eta = 0.9;
		T = 10;
		p = 0.4;
		psi = 1;
		lambda = 1;
		Step = 0.01;
		tol = 1e-10;
		
		cout << "Do you want to allow the subjects to threaten revolution after peace? 1 = yes, 0 = no.:";
		cin >> rev;
		/******************************************************************************************************
		compute the peace payoffs given this parameterization**************************************************
		******************************************************************************************************/
		if(rev == 1){
			VrP[0] = 0;
			VsP[0] = 0;
			VrRev[0] = 0;
			VsRev[0] = 0;

			VecDoub Damage(T);
			Doub costs = 0;
			for(i=1;i<T;i++){
				Damage[i] = pow(ds,i)*(1-pow(eta,i+1));
				costs += Damage[i];
			}

			for(i=1;i<51;i++){
				VrRev[i] = 0;
				VsRev[i] = (lambda+gam(i))*((1-ds)*costs+pow(ds,T));
				VsP[i] = ((1-dr-p*(1-dr))*gam(i))/(1-dr+p*(dr-ds))+lambda;
				if(max(VsRev[i],VsP[i]) == VsP[i]){
					VrP[i] = (p*(1-ds)*gam(i))/(1-dr+p*(dr-ds))+psi;
					ind[i] = 1;
				}else{
					VrP[i] = p*(gam(i)-VsRev[i]+lambda)/(1-(1-p)*dr)+psi;
					VsP[i] = ((1-p)*(1-dr)*gam(i)+p*(VsRev[i]-lambda))/(1-(1-p)*dr)+lambda;
					ind[i] = 0;
				}
			}
		}else{
			if(rev == 0){
				VrP[0] = 0;
				VsP[0] = 0;
				VrRev[0] = 0;
				VsRev[0] = 0;

				VecDoub Damage(T);
				Doub costs = 0;
				for(i=1;i<T;i++){
					Damage[i] = pow(ds,i)*(1-pow(eta,i+1));
					costs += Damage[i];
				}

				for(i=1;i<51;i++){
					VrRev[i] = 0;
					VsRev[i] = (lambda+gam(i))*((1-ds)*costs+pow(ds,T));
					VsP[i] = ((1-dr-p*(1-dr))*gam(i))/(1-dr+p*(dr-ds))+lambda;
					VrP[i] = (p*(1-ds)*gam(i))/(1-dr+p*(dr-ds))+psi;
					ind[i] = 1;
				}
			}else{;}
		}
	}else{
		if(bargain == 1){
	
			cout << "Please enter the patience of the ruler:";
			cin >> dr;
			cout << "Please enter the patience of the subjects:";
			cin >> ds;
			cout << "Please enter the material damage caused by revolution:";
			cin >> eta;
			cout << "Please enter the time required to rebuild the damage:";
			cin >> T;
			cout << "Please enter the political power of the ruler after peace:";
			cin >> p;
			cout << "Please enter the office priviledges of the Ruler:";
			cin >> psi;
			cout << "Please enter the upper bound of the feasible set, lambda:";
			cin >> lambda;
			cout << "Please enter the stepsize for the Grid:";
			cin >> Step;
			cout << "Please enter the tolerance:";
			cin >> tol;
			
			/******************************************************************************************************
			compute the peace payoffs given this parameterization**************************************************
			******************************************************************************************************/
			VrP[0] = 0;
			VsP[0] = 0;
			VrRev[0] = 0;
			VsRev[0] = 0;
	
			VecDoub Damage(T);
			Doub costs = 0;
			for(i=1;i<T;i++){
				Damage[i] = pow(ds,i)*(1-pow(eta,i+1));
				costs += Damage[i];
			}
	
			for(i=1;i<51;i++){
				VrRev[i] = 0;
				VsRev[i] = (lambda+gam(i))*((1-ds)*costs+pow(ds,T));
				VsP[i] = ((1-dr-p*(1-dr))*gam(i))/(1-dr+p*(dr-ds))+lambda;
				if(max(VsRev[i],VsP[i]) == VsP[i]){
					VrP[i] = (p*(1-ds)*gam(i))/(1-dr+p*(dr-ds))+psi;
					ind[i] = 1;
				}else{
					VrP[i] = p*(gam(i)-VsRev[i]+lambda)/(1-(1-p)*dr)+psi;
					VsP[i] = ((1-p)*(1-dr)*gam(i)+p*(VsRev[i]-lambda))/(1-(1-p)*dr)+lambda;
					ind[i] = 0;
				}
			}
			
			
			}else{
				VrP[0] = 0;
				VsP[0] = 0;
				VrRev[0] = 0;
				VsRev[0] = 0;
	
				VecDoub Damage(T);
				Doub costs = 0;
				for(i=1;i<T;i++){
					Damage[i] = pow(ds,i)*(1-pow(eta,i+1));
					costs += Damage[i];
				}
	
				for(i=1;i<51;i++){
					VrRev[i] = 0;
					VsRev[i] = (lambda+gam(i))*((1-ds)*costs+pow(ds,T));
				}
		
				for(i=1;i<51;i++){
					cout << "Please enter the peace payoff for the Ruler in state "<< i << ":";
					cin >> VrP[i];
					cout << "Please enter the peace payoff for the subjects in state "<< i << ":";
					cin >> VsP[i];
				}
			}
		}
	
	/************************************************************************************************************
	*************************************************************************************************************
	************************************************************************************************************/	
	
	/*******************************************************************************************************
	construct the grid: You were having a problem before because 1 was not 1, it ws just slightly different!
	*******************************************************************************************************/
	int L=lambda/Step;
	
	VecDoub G(L+1); //grid of mobilization values	
	for(i=0;i<=L;i++){
		G[i] = lambda*i/(L);
	}
	
	// check: seems to be correctly outputting the distributions for the first and last entry of the Grid
	/*
	q = Q(G[0]);
	for(j=0;j<=20;j++){cout << q[j] << endl;}
	*/

	/**********************************************
	assign initial policies and continuation values
	**********************************************/
	for(i=0;i<50;i++){
		V[i][0] = 0;
		V[i][1] = 0;
		
		//seems to be putting these correctly
		//cout << V[i][0] << ";";
	}
	V[50][0] = VrP[50];
	V[50][1] = VsP[50];
		

	
	for(i=1;i<=50;i++){
		r[i] = 0;
		pi[i] = 1;
		x[i] = 0;
		
		/*seems to have set initial values for strategies correctly
		cout << r[i] << ";";
		cout << pi[i] << ";";
		cout << x[i] << ";";*/
	}
	MatDoub rg(51,L+1);//a grid for all revolution strategies for any value of x
	VecDoub Vr(51); //continuation Values for the Ruler
	
	//seems to be outputting the probabilities correctly for the initial values
	/*
	for(i=1;i<=50;i++){
		q = Q(x[i],i);
		for(j=0;j<=50;j++){
			cout << q[j] << ",";
		}
		cout << endl << endl;
	}*/





	/************************************************************************************************************
	iterate calling the value function and policy updater until convergence to within the pre-specified tolerance
	************************************************************************************************************/
	it = 1;
	do{
		e[0] = 0;
		e[1] = 0;
		cout << endl;
		cout << it << endl;
		
		/****************************************************************************************
		//compute the new continuation values given the current strategies and assess convergence
		****************************************************************************************/
		Vn = Vstar(V,pi,x,r);
		
		for(i=0;i<=20;i++){
			e[0] += abs(V[i][0]-Vn[i][0]);
			e[1] += abs(V[i][1]-Vn[i][1]); 
		}
		V = Vn;
	
	
		
		/******************************************************************************
		//update the revolution strategies on the grid chosen by the Ruler and subjects
		******************************************************************************/	
		rg = rGrid(V,L,G);
		/*
		for(j=1;j<50;j++){
			for(i=0;i<=L;i++){
				cout << rg[j][i];
			}
			cout << endl;
		}*/
		
		/***************************************************************************************
		//update the mobilization strategies in response to the current revolution over the grid
		***************************************************************************************/
		for(i=0;i<51;i++){
			Vr[i] = V[i][0];
		}
		cout << "x(r,k) = {";		
		x = xmap(Vr,rg,G,L);
		for(i=1;i<=49;i++){
			if(i<49){cout << x[i] << ",";}else{cout << x[i] << "}";}
		}
		cout << endl << endl;

		//update the actual revolution strategies in response to the current x
		cout << endl << "r(x,k) = {";
		for(i=1;i<50;i++){
			for(j=0;j<=L;j++){
				if(x[i] == G[j]){r[i] = rg[i][j];}else{;}
				//if(pi[i] == 1){r[i] = 0;}else{;}
				//if(pi[i] != 0){r[i] = 0;}else{;}
			}
			if(i<49){cout << r[i] << ",";}else{cout << r[i] << "}";}
		}
		cout << endl << endl;



		/**********************************************************************************************
		//update the peace strategies in response to the new revolution strategies and in response to x
		**********************************************************************************************/
		cout << "pi(k) = {";
		for(i=1;i<=49;i++){
			q = Q(x[i],i+1);
			//cout << r[i] << endl;
			if(VrP[i]>(1-r[i])*((1-dr)*psi+dr*(q[0]*Vr[0]+q[1]*Vr[1]+q[2]*Vr[2]+q[3]*Vr[3]+q[4]*Vr[4]+q[5]*Vr[5]+q[6]*Vr[6]+q[7]*Vr[7]+
		q[8]*Vr[8]+q[9]*Vr[9]+q[10]*Vr[10]+q[11]*Vr[11]+q[12]*Vr[12]+q[13]*Vr[13]+q[14]*Vr[14]+q[15]*Vr[15]+q[16]*Vr[16]+q[17]*Vr[17]+q[18]*Vr[18]+q[19]*Vr[19]+q[20]*Vr[20]+q[21]*Vr[21]
		+q[22]*Vr[22]+q[23]*Vr[23]+q[24]*Vr[24]+q[25]*Vr[25]+q[26]*Vr[26]+q[27]*Vr[27]+q[28]*Vr[28]+q[29]*Vr[29]+q[30]*Vr[30]+q[31]*Vr[31]+q[32]*Vr[32]+q[33]*Vr[33]+q[34]*Vr[34]
		+q[35]*Vr[35]+q[36]*Vr[36]+q[37]*Vr[37]+q[38]*Vr[38]+q[39]*Vr[39]+q[40]*Vr[40]+q[41]*Vr[41]+q[42]*Vr[42]+q[43]*Vr[43]+q[44]*Vr[44]+q[45]*Vr[45]+q[46]*Vr[46]+q[47]*Vr[47]
		+q[48]*Vr[48]+q[49]*Vr[49]+q[50]*Vr[50]))){pi[i] = 1;}else{pi[i] = 0;}
			
			if(i<49){cout << pi[i] << ",";}else{cout << pi[i] << "}";}
			if(pi[i] != 0){x[i] = 0;}else{;}
			if(r[i] != 0){pi[i] = 1;}else{;}
		}
		cout << endl << endl;



	//output vectors of continuation values
	cout << "Continuation values for the Ruler are VR(k) = " << endl;
	for(i=0;i<50;i++){
		cout << V[i][0] << ",";
	}
	cout << V[50][0] << endl;
	cout << "Continuation values for the Subjects are VS(k) = " << endl;	
	for(i=0;i<50;i++){
		cout << V[i][1] << ",";
	}
	cout << V[50][1] << endl;

	cout << endl << endl;
	for(j=1;j<50;j++){
		if(j<10){
			cout << j << "  ";}else{cout << j << " ";}
		for(i=0;i<=L;i++){
			cout << rg[j][i];
		}
		cout << " " << pi[j] << endl;
	}
	cout << endl;

	
	it++;
	}while(tol < e[0]|tol< e[1]);

	cout << endl << "Do you wish to view the results now (1 = continue)?";
	cin >> bargain;
	cout << endl;

	/****************************************************************
	First output the current parameter values to keep track of things
	****************************************************************/
	cout << "Ruler patience is (dr): " << dr << endl;
	cout << "Subject patience is (ds): " << ds << endl;
	cout << "Material cost of revolution is: " << eta << endl;
	cout << "Length of time to repair the damage is: " << T << endl;
	cout << "Total societal resources are: " << lambda << endl;
	cout << "Political balance of power after peace: " << p << endl;
	cout << "Political benefits to holding office are: " << psi << endl;
	/************************************************************************************************************
	Results of the fixed point iterations; the continuation values and strategies of a Markov Perfect Equilibrium
	************************************************************************************************************/
	r[50] = 0;	
	for(j=1;j<50;j++){
		if(j<10){
			cout << j << "  ";}else{cout << j << " ";}
		for(i=0;i<=L;i++){
			cout << rg[j][i];
		}
		cout << " " << pi[j] << endl;
	}
	cout << endl;
	for(i=1;i<=49;i++){
			
		if(pi[i] == 1){x[i] = 0;}else{x[i]=x[i];}
			
	}

	for(i=0;i<51;i++){
		K[i] = i;
	}
	
	/**************************************************************
	Need to write out: 1) eq continuation values 2) eq mobilization 
	**************************************************************/
	ofstream out1("MarkovPerfectEq.csv");
	out1 << "VrP" << endl;
	for(i=0;i<51;i++){
		out1 << VrP[i] << ",";
	}
	out1 << endl << endl;
	out1 << "VsP" << endl;
	for(i=0;i<51;i++){
		out1 << VsP[i] << ",";
	}
	out1 << endl << endl;	
	out1 << "pi*" << endl;
	for(i=0;i<51;i++){
		out1 << pi[i] << ",";
	}
	out1 << endl << endl;	
	out1 << "x*" << endl;
	for(i=0;i<51;i++){
		out1 << x[i] << ",";
	}
	out1 << endl << endl;
	out1 << "VR*" << endl;
	for(i=0;i<51;i++){
		out1 << V[i][0] << ",";
	}
	out1 << endl << endl;	
	out1 << "VS*" << endl;
	for(i=0;i<51;i++){
		out1 << V[i][1] << ",";
	}
	out1 << endl << endl;	
	
	cout << "Values for peace for the ruler, and for peace and revolution for the subjects are:" << endl;
	cout << "Ind" << setw(12) << "VrP" << setw(15) << "VsP" << setw(15) << "VsRev" << setw(18) << "pi(k)*" << setw(15) << "x(r,k)*" << setw(15) << "r(x,k)*" << setw(15) << "VR(k)*" 
	<< setw(15) << "VS(k)*" << setw(15) << "State" << endl;
	for(i=1;i<51;i++){
		cout << ind[i] << setw(15) << VrP[i] << setw(15) << VsP[i] << setw(15) << VsRev[i] << setw(15) << pi[i] << setw(15) << x[i] << setw(15) << r[i] << setw(15) << V[i][0] 
		<< setw(15) << V[i][1] << setw(15) << K[i] <<endl;
	}
	
return 0;
}