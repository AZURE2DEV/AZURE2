#include "include/VavilovMoyalFunction.h"
#include "include/ChebyshevPolynomial.h"

// Moyal auxiliary functions
double VavilovMoyalFunction::VMbeta(double beta)
{
 double betamin,betamax,betabar;
 if(beta>=0 && beta<=1)
 {
   betamin = 0.05; betamax = 1.0;
 }
 else
 {
   return 0;
 }
 betabar = 1 + 2 * (beta - betamax) / (betamax - betamin);
 return betabar;
}

double VavilovMoyalFunction::VMk(double k)
{
 double kmin,kmax,kbar;
 if(k>=0.02 && k<0.11)
 {
   kmin = 0.02; kmax = 0.10;
 }
 else if(k>=0.11 && k<0.22)
 {
   kmin = 0.12; kmax = 0.20;
 }
 else if(k>=0.22 && k<0.30)
 {
   kmin = 0.22; kmax = 0.30;
 }
 else
 {
  return 0;
 }
 kbar = 1 + 2 * (k - kmax) / (kmax - kmin) ;
 return kbar;
}

// Moyal Parameters. (An hard-coded table that is loaded by the class constructor)
void VavilovMoyalFunction::MakeTable()
{
 VMT[0][0][0] = 0; VMT[0][0][1] = 0; VMT[0][0][2] = 0.25850874; VMT[0][0][3] = 0.43142617; VMT[0][0][4] = 0.25225964; VMT[0][0][5] = 1.2593243; VMT[0][0][6] = -0.024864376; VMT[0][0][7] = 0.035855696; VMT[0][0][8] = 10.234692; VMT[0][0][9] = 21.487522;
 VMT[0][1][0] = 0; VMT[0][1][1] = 1; VMT[0][1][2] = 0.024880692; VMT[0][1][3] = 0.042127077; VMT[0][1][4] = 0.023834176; VMT[0][1][5] = -0.046865180; VMT[0][1][6] = 0.0018751903; VMT[0][1][7] = -0.00084479939; VMT[0][1][8] = -1.9952390; VMT[0][1][9] = -3.4343169;
 VMT[0][2][0] = 0; VMT[0][2][1] = 2; VMT[0][2][2] = 0.0047404356; VMT[0][2][3] = 0.0073167928; VMT[0][2][4] = 0.0021624675; VMT[0][2][5] = -0.0077222986; VMT[0][2][6] = 0.0012668869; VMT[0][2][7] = 0; VMT[0][2][8] = -0.45679694; VMT[0][2][9] = -1.1063164;
 VMT[0][3][0] = 0; VMT[0][3][1] = 3; VMT[0][3][2] = -0.00074445130; VMT[0][3][3] = -0.0014026047; VMT[0][3][4] = -0.0026865597; VMT[0][3][5] = 0.0032241039; VMT[0][3][6] = -0.00048736023; VMT[0][3][7] = -0.00045675843; VMT[0][3][8] = 0; VMT[0][3][9] = -0.21000819;
 VMT[0][4][0] = 1; VMT[0][4][1] = 0; VMT[0][4][2] = 0.032477982; VMT[0][4][3] = 0.040797543; VMT[0][4][4] = 0.064820468; VMT[0][4][5] = -0.20374501; VMT[0][4][6] = -0.0010368495; VMT[0][4][7] = -0.027542114; VMT[0][4][8] = -3.5619655; VMT[0][4][9] = -11.825253;
 VMT[0][5][0] = 1; VMT[0][5][1] = 1; VMT[0][5][2] = 0.0073225731; VMT[0][5][3] = 0.016195241; VMT[0][5][4] = -0.0054891384; VMT[0][5][5] = 0.0089882920; VMT[0][5][6] = 0.0034850854; VMT[0][5][7] = -0.0069836141; VMT[0][5][8] = 0.50505298; VMT[0][5][9] = 1.7891643;
 VMT[0][6][0] = 1; VMT[0][6][1] = 2; VMT[0][6][2] = 0; VMT[0][6][3] = 0; VMT[0][6][4] = -0.0089439554; VMT[0][6][5] = 0.018786468; VMT[0][6][6] = 0.0019372124; VMT[0][6][7] = 0; VMT[0][6][8] = 0; VMT[0][6][9] = 0.73410606;
 VMT[0][7][0] = 1; VMT[0][7][1] = 3; VMT[0][7][2] = -0.0015727318; VMT[0][7][3] = -0.0025141668; VMT[0][7][4] = -0.0062756944; VMT[0][7][5] = 0.014484097; VMT[0][7][6] = 0.00070761825; VMT[0][7][7] = 0.0015298434; VMT[0][7][8] = 0; VMT[0][7][9] = 0;
 VMT[0][8][0] = 2; VMT[0][8][1] = 0; VMT[0][8][2] = -0.0059020496; VMT[0][8][3] = -0.0091490215; VMT[0][8][4] = -0.023615759; VMT[0][8][5] = 0.095055662; VMT[0][8][6] = 0.0014330117; VMT[0][8][7] = 0.012631023; VMT[0][8][8] = 0.69387764; VMT[0][8][9] = 4.3133087;
 VMT[0][9][0] = 2; VMT[0][9][1] = 1; VMT[0][9][2] = 0; VMT[0][9][3] = 0.0024714789; VMT[0][9][4] = 0.0039800522; VMT[0][9][5] = -0.0067167236; VMT[0][9][6] = 0; VMT[0][9][7] = 0.0039876546; VMT[0][9][8] = 0; VMT[0][9][9] = -0.89601916;
 VMT[0][10][0] = 2; VMT[0][10][1] = 2; VMT[0][10][2] = -0.0011210142; VMT[0][10][3] = -0.0014064022; VMT[0][10][4] = -0.0024655436; VMT[0][10][5] = 0; VMT[0][10][6] = 0.00046898375; VMT[0][10][7] = 0.0019247256; VMT[0][10][8] = 0; VMT[0][10][9] = -0.32454506;
 VMT[0][11][0] = 3; VMT[0][11][1] = 0; VMT[0][11][2] = 0; VMT[0][11][3] = 0; VMT[0][11][4] = 0; VMT[0][11][5] = -0.020771531; VMT[0][11][6] = 0.00020052730; VMT[0][11][7] = -0.0030188807; VMT[0][11][8] = -0.14047599; VMT[0][11][9] = -1.4500543;
 VMT[0][12][0] = 3; VMT[0][12][1] = 1; VMT[0][12][2] = 0.0011668284; VMT[0][12][3] = 0.0020751278; VMT[0][12][4] = 0.0048447456; VMT[0][12][5] = -0.013049241; VMT[0][12][6] = -0.00036597173; VMT[0][12][7] = -0.0036055679; VMT[0][12][8] = 0; VMT[0][12][9] = 0.39120793;
 VMT[1][0][0] = 0; VMT[1][0][1] = 0; VMT[1][0][2] = 0.27827257; VMT[1][0][3] = 0.41421789; VMT[1][0][4] = 0.20191059; VMT[1][0][5] = 1.3206085; VMT[1][0][6] = 0.016435245; VMT[1][0][7] = 0.033432409; VMT[1][0][8] = 5.4529572; VMT[1][0][9] = 9.3841352;
 VMT[1][1][0] = 0; VMT[1][1][1] = 1; VMT[1][1][2] = 0.045091424; VMT[1][1][3] = 0.12693873; VMT[1][1][4] = 0.053921588; VMT[1][1][5] = -0.14986093; VMT[1][1][6] = -0.010775802; VMT[1][1][7] = -0.013346861; VMT[1][1][8] = -1.2218009; VMT[1][1][9] = -1.8160479;
 VMT[1][2][0] = 0; VMT[1][2][1] = 2; VMT[1][2][2] = 0.0080559636; VMT[1][2][3] = 0.022999801; VMT[1][2][4] = 0.0035068740; VMT[1][2][5] = -0.012720568; VMT[1][2][6] = 0.0051476061; VMT[1][2][7] = -0.0017402116; VMT[1][2][8] = -0.32324120; VMT[1][2][9] = -0.50919193;
 VMT[1][3][0] = 0; VMT[1][3][1] = 3; VMT[1][3][2] = -0.0038974523; VMT[1][3][3] = -0.0086792801; VMT[1][3][4] = -0.012621494; VMT[1][3][5] = 0.024972042; VMT[1][3][6] = 0.0056856517; VMT[1][3][7] = 0.0021052496; VMT[1][3][8] = -0.027373591; VMT[1][3][9] = -0.051384654;
 VMT[1][4][0] = 1; VMT[1][4][1] = 0; VMT[1][4][2] = -0.0014227603; VMT[1][4][3] = -0.030061649; VMT[1][4][4] = -0.046831422; VMT[1][4][5] = 0.10036618; VMT[1][4][6] = 0.036051400; VMT[1][4][7] = 0.0060583916; VMT[1][4][8] = -0.90906096; VMT[1][4][9] = -1.6276904;
 VMT[1][5][0] = 1; VMT[1][5][1] = 1; VMT[1][5][2] = 0; VMT[1][5][3] = 0.031875584; VMT[1][5][4] = -0.0054996531; VMT[1][5][5] = -0.0097751962; VMT[1][5][6] = -0.013438433; VMT[1][5][7] = 0.0015528195; VMT[1][5][8] = 0.12173464; VMT[1][5][9] = 0.21413992;
 VMT[1][6][0] = 1; VMT[1][6][1] = 2; VMT[1][6][2] = 0.0054730726; VMT[1][6][3] = 0.019716857; VMT[1][6][4] = 0.018513506; VMT[1][6][5] = -0.048282515; VMT[1][6][6] = -0.0025421507; VMT[1][6][7] = -0.0045124157; VMT[1][6][8] = 0.040917471; VMT[1][6][9] = 0.066596366;
 VMT[1][7][0] = 1; VMT[1][7][1] = 3; VMT[1][7][2] = 0.0019792507; VMT[1][7][3] = 0.0032596742; VMT[1][7][4] = 0.0068332334; VMT[1][7][5] = -0.0098552378; VMT[1][7][6] = 0.0020169108; VMT[1][7][7] = -0.0015629454; VMT[1][7][8] = 0; VMT[1][7][9] = 0;
 VMT[1][8][0] = 2; VMT[1][8][1] = 0; VMT[1][8][2] = 0.0024848327; VMT[1][8][3] = 0.0052249697; VMT[1][8][4] = 0.0096777473; VMT[1][8][5] = -0.022015201; VMT[1][8][6] = 0.0023036520; VMT[1][8][7] = -0.0023381379; VMT[1][8][8] = 0.086122438; VMT[1][8][9] = 0.16571423;
 VMT[1][9][0] = 2; VMT[1][9][1] = 1; VMT[1][9][2] = -0.0030634124; VMT[1][9][3] = -0.0061757928; VMT[1][9][4] = -0.0090029985; VMT[1][9][5] = 0.026087455; VMT[1][9][6] = 0; VMT[1][9][7] = 0.0021900670; VMT[1][9][8] = 0; VMT[1][9][9] = 0;
 VMT[1][10][0] = 2; VMT[1][10][1] = 2; VMT[1][10][2] = 0; VMT[1][10][3] = 0; VMT[1][10][4] = -0.0012940502; VMT[1][10][5] = 0; VMT[1][10][6] = -0.0015144931; VMT[1][10][7] = 0.00022499176; VMT[1][10][8] = 0; VMT[1][10][9] = 0;
 VMT[1][11][0] = 3; VMT[1][11][1] = 0; VMT[1][11][2] = 0; VMT[1][11][3] = 0; VMT[1][11][4] = -0.0017995317; VMT[1][11][5] = 0.0061667091; VMT[1][11][6] = -0.00061666343; VMT[1][11][7] = 0.00083846081; VMT[1][11][8] = 0; VMT[1][11][9] = 0;
 VMT[1][12][0] = 3; VMT[1][12][1] = 1; VMT[1][12][2] = 0.00075633702; VMT[1][12][3] = 0; VMT[1][12][4] = 0.0034958743; VMT[1][12][5] = -0.011399062; VMT[1][12][6] = 0; VMT[1][12][7] = -0.0013202847; VMT[1][12][8] = 0; VMT[1][12][9] = 0;
 VMT[2][0][0] = 0; VMT[2][0][1] = 0; VMT[2][0][2] = 0.29712948; VMT[2][0][3] = 0.40882632; VMT[2][0][4] = 0.16861629; VMT[2][0][5] = 1.3493897; VMT[2][0][6] = 0.10264949; VMT[2][0][7] = 0.029568177; VMT[2][0][8] = 0; VMT[2][0][9] = 6.6184654;
 VMT[2][1][0] = 0; VMT[2][1][1] = 1; VMT[2][1][2] = 0.035707399; VMT[2][1][3] = 0.18719727; VMT[2][1][4] = 0.030144338; VMT[2][1][5] = -0.083447911; VMT[2][1][6] = -0.043097757; VMT[2][1][7] = -0.0048515387; VMT[2][1][8] = 0; VMT[2][1][9] = -1.4540925;
 VMT[2][2][0] = 0; VMT[2][2][1] = 2; VMT[2][2][2] = 0.0096221631; VMT[2][2][3] = 0.056954987; VMT[2][2][4] = 0.013891826; VMT[2][2][5] = -0.048061360; VMT[2][2][6] = -0.0022647176; VMT[2][2][7] = -0.0040797531; VMT[2][2][8] = 0; VMT[2][2][9] = -0.39529833;
 VMT[2][3][0] = 0; VMT[2][3][1] = 3; VMT[2][3][2] = -0.0018402821; VMT[2][3][3] = 0; VMT[2][3][4] = -0.0058030495; VMT[2][3][5] = 0.0076473951; VMT[2][3][6] = 0.0094531290; VMT[2][3][7] = 0.00040403265; VMT[2][3][8] = 0; VMT[2][3][9] = -0.044293243;
 VMT[2][4][0] = 1; VMT[2][4][1] = 0; VMT[2][4][2] = 0.0097572934; VMT[2][4][3] = 0.014474912; VMT[2][4][4] = 0; VMT[2][4][5] = -0.0026863185; VMT[2][4][6] = 0.032738857; VMT[2][4][7] = -0.00163000060; VMT[2][4][8] = 0; VMT[2][4][9] = -0.73866379;
 VMT[2][5][0] = 1; VMT[2][5][1] = 1; VMT[2][5][2] = -0.0049821585; VMT[2][5][3] = 0.023020158; VMT[2][5][4] = -0.0038717547; VMT[2][5][5] = 0.024494430; VMT[2][5][6] = -0.012442571; VMT[2][5][7] = 0.0018200105; VMT[2][5][8] = 0; VMT[2][5][9] = 0.088741049;
 VMT[2][6][0] = 1; VMT[2][6][1] = 2; VMT[2][6][2] = 0.0020301312; VMT[2][6][3] = 0.019300232; VMT[2][6][4] = 0.0082387775; VMT[2][6][5] = -0.047890063; VMT[2][6][6] = -0.0088293329; VMT[2][6][7] = -0.0037432073; VMT[2][6][8] = 0; VMT[2][6][9] = 0;
 VMT[2][7][0] = 1; VMT[2][7][1] = 3; VMT[2][7][2] = -0.0018723311; VMT[2][7][3] = 0; VMT[2][7][4] = -0.010116105; VMT[2][7][5] = 0.017778596; VMT[2][7][6] = 0.0052537299; VMT[2][7][7] = 0.0019950380; VMT[2][7][8] = 0; VMT[2][7][9] = 0;
 VMT[2][8][0] = 2; VMT[2][8][1] = 0; VMT[2][8][2] = 0; VMT[2][8][3] = 0.0025023704; VMT[2][8][4] = 0.0036317285; VMT[2][8][5] = -0.0035216040; VMT[2][8][6] = 0; VMT[2][8][7] = -0.00021119745; VMT[2][8][8] = 0; VMT[2][8][9] = 0.044693973;
 VMT[2][9][0] = 2; VMT[2][9][1] = 1; VMT[2][9][2] = 0.0018831112; VMT[2][9][3] = 0.0050574313; VMT[2][9][4] = 0.0085359607; VMT[2][9][5] = -0.016209200; VMT[2][9][6] = -0.0032283517; VMT[2][9][7] = -0.0014346306; VMT[2][9][8] = 0; VMT[2][9][9] = 0;
 VMT[2][10][0] = 2; VMT[2][10][1] = 2; VMT[2][10][2] = -0.00073403108; VMT[2][10][3] = 0; VMT[2][10][4] = -0.0055135670; VMT[2][10][5] = 0.013179324; VMT[2][10][6] = 0.0013340546; VMT[2][10][7] = 0.0012222675; VMT[2][10][8] = 0; VMT[2][10][9] = 0;
 VMT[2][11][0] = 3; VMT[2][11][1] = 0; VMT[2][11][2] = -0.0015291686; VMT[2][11][3] = -0.0037707379; VMT[2][11][4] = -0.0043657818; VMT[2][11][5] = 0.024434909; VMT[2][11][6] = 0.0043608779; VMT[2][11][7] = 0.0023599053; VMT[2][11][8] = 0; VMT[2][11][9] = 0;
 VMT[2][12][0] = 3; VMT[2][12][1] = 1; VMT[2][12][2] = 0.0043541673; VMT[2][12][3] = 0.0094550140; VMT[2][12][4] = 0.014507659; VMT[2][12][5] = -0.037768479; VMT[2][12][6] = -0.0075640352; VMT[2][12][7] = -0.0039165276; VMT[2][12][8] = 0; VMT[2][12][9] = 0;
 return;
}

// Return the table value from the arguments
double VavilovMoyalFunction::VMTable(unsigned int z, unsigned int i, unsigned int m, unsigned int n)
{
 // Find the correct address
 for (unsigned int x=0; x<13; x++)
 {
   if(VMT[z][x][0] == m && VMT[z][x][1] == n)
    return VMT[z][x][i+1];
 }
 return 0;
}

// Moyal Coefficients.
// i=0: lambda_0, where not need the table.
// i=1,2,3,4,5,6 : a_1, a_2, a_3, a_4, a_5, a_6
// i=7 : lambda_0.95 (2-sigma)
// i=8 : lambda_0.995 (3-sigma)
double VavilovMoyalFunction::VMa(unsigned int i, double k, double beta)
{
 ChebyshevPolynomial Tm,Tn;
 double asum = 0;
 unsigned int z;
 if(k>=0.02 && k<0.11)
 {
   z = 0;
   if(i==0)
    return -3.03;
 }
 else if(k>=0.11 && k<0.22)
 {
   z = 1;
   if(i==0)
    return -3.04;
 }
 else if(k>=0.22 && k<0.30)
 {
   z = 2;
   if(i==0)
    return -3.05;
 }
 else
 {
  return 0;
 }
 for(unsigned int m=0; m<4; m++)
 {
  for(unsigned int n=0; n<4; n++)
  {
    if(m+n<=4)
    {
     double kterm = this->VMk(k);
     double betaterm = this->VMbeta(beta);
     double aterm = this->VMTable(z,i,m,n) * Tm.GetValue(m,kterm) * Tn.GetValue(n,betaterm);
     asum = asum + aterm;
    }
  }
 }
 return asum;
}

// Evaluate the value of the Vavilov distribution by the Moyal approximation
double VavilovMoyalFunction::VMTwoSigma(double k, double beta, double lambda)
{
 // Evaluate the function
 double f = this->VMa(1,k,beta) * std::exp(-1*this->VMa(2,k,beta) * (lambda + this->VMa(5,k,beta) * lambda * lambda) - (this->VMa(3,k,beta) * std::exp(-1*this->VMa(4,k,beta) * (lambda + this->VMa(6,k,beta) * lambda * lambda))));
 return (f>0) ? f : 0.0;
}

double VavilovMoyalFunction::VMThreeSigma(double k, double beta, double lambda)
{
 // Global parameters. Notice that a7 is lambda-TwoSigma, and a8 is Lambda-ThreeSigma
 double f = this->VMTwoSigma(k,beta,this->VMa(7,k,beta));
 double I0 = 0.045;
 double b = (this->VMa(8,k,beta) - this->VMa(7,k,beta)) / (this->VMa(8,k,beta) * this->VMa(7,k,beta));
 double a = 1 / std::log(this->VMa(8,k,beta)/this->VMa(7,k,beta));
 double D = this->VMa(7,k,beta) * this->VMa(7,k,beta) * (f - (a*I0)/this->VMa(7,k,beta)) * (1/(1 - this->VMa(7,k,beta) * a * b)) ;
 double delta = a * (I0/D-b);
 // Evaluate the function
 double x = D * ((1 / (lambda * lambda)) + (delta / lambda));
 return (x>0) ? x : 0.0;
}

// Evaluate the Vavilov-Moyal function
double VavilovMoyalFunction::VMMain(double k, double beta, double lambda)
{
 // Obtain the distribution value
 if(k>=0.02 && k<0.22)
 {
   double lambda0 = this->VMa(0,k,beta); //Absolute minimum value cut-off
   // Function selection
   if(lambda0 <= lambda && lambda < this->VMa(7,k,beta))
   {
     return this->VMTwoSigma(k,beta,lambda);
   }
   else if(this->VMa(7,k,beta) <= lambda && lambda < this->VMa(8,k,beta))
   {
     return this->VMThreeSigma(k,beta,lambda);
   }
   else
   {
    return 0;
   }
 }
 else if(k>=0.22 && k<0.30)
 {
   double lambda0 = this->VMa(0,k,beta); //Absolute minimum value cut-off
   // Function selection
   if(lambda0 <= lambda && lambda < this->VMa(8,k,beta))
   {
     return this->VMTwoSigma(k,beta,lambda);
   }
   else
   {
     return 0;
   }
 }
 else
 {
   return 0; //Out-of-Domain
 }
}


// Define the distribution domain
void VavilovMoyalFunction::SetMoyalStep(double xi, double beta, double k, double DEM, unsigned int numberstep, bool TrimNegative)
{
 double euler = (std::lgamma(0.999999) - std::lgamma(1.000001)) / (0.000002); // Euler's Constant
 double lambdalow = this->VMa(0,k,beta); //lambda_0
 double lambdahigh = this->VMa(8,k,beta); //lambda_0.995
 //Since the distribution domain cannot be bounded analytically, the number of steps are the division by unity.
 double lambdastep = 1.0/numberstep;
 VMCbeta = beta;
 VMCk = k;
 VMDEM = DEM;
 VMxi = xi;
 MoyalStep = xi*lambdastep;
 MoyalMinimum = xi*(lambdalow + 1 - euler + beta*beta + std::log(k)) + DEM;
 MoyalMaximum = xi*(lambdahigh + 1 - euler + beta*beta + std::log(k)) + DEM;
 // Cut energy loss greater than average
 if((MoyalMinimum < 0) && TrimNegative)
 {
  MoyalMinimum = 0;
 }
 return;
}

double VavilovMoyalFunction::GetValue(double AtEnergy)
{
 double euler = (std::lgamma(0.999999) - std::lgamma(1.000001)) / (0.000002); // Euler's Constant
 double MoyalLambda = (AtEnergy - VMDEM)/VMxi - 1 + euler - VMCbeta*VMCbeta - std::log(VMCk);
 return this->VMMain(VMCk,VMCbeta,MoyalLambda);
}