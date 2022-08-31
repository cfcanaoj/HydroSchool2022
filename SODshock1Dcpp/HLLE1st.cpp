#include<iostream>
#include<fstream>

#include<math.h> // pow
using namespace std;

//衝撃波管問題　一次元　最も簡単なやつ//

namespace control {
  int nstep;
  const int Nend=2500;  //grid number of t
  const int Ig=100; //grid number of x
  const int Is= 2;
  const int Ie = Is + Ig - 1;
  const int In = Ie + Is;
  double tsim;
  double dt;
  const double Cc=0.3;
}

namespace grid {
  using namespace control;
  const  double xmin=0.0;
  const  double xmax=1.0;
  double dx;
  double xa[In+1];
  double xb[In];
  void SetupGrid();
  void DetermineTimestep();
}

namespace field{
  using namespace control;
    double   rho[In];
    double    vx[In];
    double     p[In];
    double    ei[In];
    double    cs[In];
    double rhovx[In];
    double    et[In];
  void SetInitialProfile();
  void UpdateFields();
  void EstimateConsvVariable();
  void EstimatePrimvVariable();
}

namespace eos{
const double gam=1.4; //spesific heat ratio
}

namespace flux{
  using namespace control;
  const int ndn=0;
  const int nvx=1;
  const int nei=2;
  const int npr=3;
  const int ncs=4;
  const int nhyd=5;
  double  svc[In][nhyd];
  void EstimateStateVector();

  const int mdn=0;
  const int mvx=1;
  const int met=2;
  const int mflx=3;
  double  numflux[In+1][mflx];
  void EstimateNumFlux();
  void hllflux(double const (&leftst)[2*mflx+2],double const (&rigtst)[2*mflx+2],double (&numfa)[mflx]);
}


namespace dataio{
  const double dtdmp=1.0e-3;
  double tout=0.0;
  int index=0;
  void OutputRadialProfile();
}

int main(void){
  using namespace control;
  using namespace grid;
  using namespace field;
  using namespace flux;
  using namespace dataio;

  SetupGrid();
  SetInitialProfile();
  EstimateConsvVariable();
  tsim=0.0;
  for (nstep=0;nstep<Nend;nstep++){
    DetermineTimestep();
    EstimateStateVector();
    EstimateNumFlux();
    UpdateFields();
    EstimatePrimvVariable();
    tsim=tsim+dt;
    OutputRadialProfile();
  }

  return 0;
}

void grid::SetupGrid(void){
  int i;
  dx=(xmax-xmin)/Ig;
  for (i=0;i<=In;i++){
    xa[i] = (i-Is)*dx;
  }
  for (i=0;i<In;i++){
    xb[i] = (xa[i+1]+xa[i])/2.0;
    //    cout<<i<<" "<<xb[i]<<endl;
  }
}

void grid::DetermineTimestep(void){
  using namespace field;
  int i;
  double dtloc;
  double dtmin;

  dtmin=1.0e99;
  for (i=Is;i<=Ie;i++){
    dtloc = dx/(cs[i]+fabs(vx[i]));
    if(dtloc < dtmin){
      dtmin=dtloc;
    }
  }
  dt=Cc*dtmin;
  dt=min(dt,1.0e-4);
}

void field::SetInitialProfile(void){
  using namespace grid;
  using namespace eos;
  int i;
  //初期値
  for (i =0; i < In; i++) {
    if (xb[i]< 0.5) {
      rho[i] = 1.0;
        p[i] = 1.0;
       vx[i] = 0.0;
       ei[i] = p[i]/(gam-1.0);
       cs[i] = pow(gam*p[i]/rho[i],0.5);
    }else{
      rho[i] = 0.125;
        p[i] = 0.1;
       vx[i] = 0.0;
       ei[i] = p[i]/(gam-1.0);
       cs[i] = pow(gam*p[i]/rho[i],0.5);
    }
    //    cout<<"ei"<<" "<<ei[i]<<endl;
  }
}

void field::EstimateConsvVariable(void){
  using namespace eos;
  int i;
  for (i =0; i < In; i++) {
    rhovx[i]= rho[i]*vx[i];
       et[i]= ei[i] + 0.5*rho[i]*vx[i]*vx[i];
       //       cout<<"et"<<" "<<et[i]<<endl;
  }
}

void field::EstimatePrimvVariable(void){
  using namespace eos;
  int i;
  for (i =0; i < In; i++) {
    vx[i] = rhovx[i]/rho[i];
    ei[i] = et[i] - 0.5*rho[i]*vx[i]*vx[i];
     p[i] = ei[i]*(gam-1.0);
    cs[i] = pow(gam*p[i]/rho[i],0.5);

  }
}

void flux::EstimateStateVector(void){
  using namespace field;
  int i;
  for (i =0; i < In; i++) {
    svc[i][ndn]=rho[i];
    svc[i][nvx]= vx[i];
    svc[i][nei]= ei[i];
    svc[i][npr]=  p[i];
    svc[i][ncs]= cs[i];
  }
}


void flux::EstimateNumFlux(void){
  using namespace grid;
  using namespace field;
  int i,n,m;
  double svcaleft[In+1][nhyd];
  double svcarigt[In+1][nhyd];

  double  ualeft[In+1][mflx];
  double  faleft[In+1][mflx];
  double  csaleft[In+1];
  double  vxaleft[In+1];

  double  uarigt[In+1][mflx];
  double  farigt[In+1][mflx];
  double  csarigt[In+1];
  double  vxarigt[In+1];

  double leftst[2*mflx+2];
  double rigtst[2*mflx+2];
  double numfa[mflx];

  for (i =1; i <= In-1; i++) {
    for (n =0; n < nhyd; n++) {
      svcaleft[i+1][n] = svc[i][n]; //i+1/2,left side
      svcarigt[i  ][n] = svc[i][n]; //i-1/2,right side
    }
  }

  for (i =1; i <= In-1; i++){
    ualeft[i+1][mdn] = svcaleft[i+1][ndn]; // rho
    ualeft[i+1][mvx] = svcaleft[i+1][ndn]*svcaleft[i+1][nvx]; //rho v
    ualeft[i+1][met] = svcaleft[i+1][nei]
    + 0.5*svcaleft[i+1][ndn]*svcaleft[i+1][nvx]*svcaleft[i+1][nvx]; //e_t = e_i + 1/2 rho v^2

    faleft[i+1][mdn] = svcaleft[i+1][ndn]*svcaleft[i+1][nvx]; // rho v
    faleft[i+1][mvx] = svcaleft[i+1][ndn]*svcaleft[i+1][nvx]*svcaleft[i+1][nvx]
                      +svcaleft[i+1][npr]; //rho v v+p
    faleft[i+1][met] = (
                 + 0.5*svcaleft[i+1][ndn]*svcaleft[i+1][nvx]*svcaleft[i+1][nvx]
                      +svcaleft[i+1][nei]
                      +svcaleft[i+1][npr])*svcaleft[i+1][nvx]; // (et+p)*v
    csaleft[i+1] = svcaleft[i+1][ncs];
    vxaleft[i+1] = svcaleft[i+1][nvx];

    uarigt[i  ][mdn] = svcarigt[i  ][ndn]; // rho
    uarigt[i  ][mvx] = svcarigt[i  ][ndn]*svcarigt[i  ][nvx]; //rho v
    uarigt[i  ][met] = svcarigt[i  ][nei]
    + 0.5*svcarigt[i  ][ndn]*svcarigt[i  ][nvx]*svcarigt[i  ][nvx]; //e_t = e_i + 1/2 rho v^2

    farigt[i  ][mdn] = svcarigt[i  ][ndn]*svcarigt[i  ][nvx]; // rho v
    farigt[i  ][mvx] = svcarigt[i  ][ndn]*svcarigt[i  ][nvx]*svcarigt[i  ][nvx]
                      +svcarigt[i  ][npr]; //rho v v+p
    farigt[i  ][met] = (
                  +0.5*svcarigt[i  ][ndn]*svcarigt[i  ][nvx]*svcarigt[i  ][nvx]
                  +    svcarigt[i  ][nei]
                      +svcarigt[i  ][npr])*svcarigt[i  ][nvx]; // (et+p)*v
    csarigt[i ] = svcarigt[i][ncs];
    vxarigt[i ] = svcarigt[i][nvx];

  }

  for (i =Is; i <= Ie+1; i++){
    leftst[0]=ualeft[i][mdn];
    leftst[1]=ualeft[i][mvx];
    leftst[2]=ualeft[i][met];
    leftst[3]=faleft[i][mdn];
    leftst[4]=faleft[i][mvx];
    leftst[5]=faleft[i][met];
    leftst[6]=csaleft[i];
    leftst[7]=vxaleft[i];

    rigtst[0]=uarigt[i][mdn];
    rigtst[1]=uarigt[i][mvx];
    rigtst[2]=uarigt[i][met];
    rigtst[3]=farigt[i][mdn];
    rigtst[4]=farigt[i][mvx];
    rigtst[5]=farigt[i][met];
    rigtst[6]=csarigt[i];
    rigtst[7]=vxarigt[i];

    hllflux(leftst,rigtst,numfa);
    for(m=0;m<mflx;m++){
      numflux[i][m] = numfa[m];
    }
  }

}

void flux::hllflux(double const (&leftst)[2*mflx+2],double const (&rigtst)[2*mflx+2],double (&numfa)[mflx]){
  using namespace std; // for max(A,B), min(A,B)
    double csR,csL;
    double vxR,vxL;
    double sR,sL;
    double uR[mflx];
    double uL[mflx];
    double fR[mflx];
    double fL[mflx];
    int m;

    for (m=0;m<mflx;m++){
      uL[m]=leftst[m];
      fL[m]=leftst[m+mflx];
      uR[m]=rigtst[m];
      fR[m]=rigtst[m+mflx];
    }

    csL=leftst[6];
    vxL=leftst[7];
    csR=rigtst[6];
    vxR=rigtst[7];

    sR = max(vxR + csR ,vxL + csL );
    sL = min(vxR - csR ,vxL - csL );
    
    for (m=0;m<mflx;m++){
      numfa[m]= ( sR * fL[m] - sL * fR[m] + sR * sL* (uR[m] - uL[m]) ) / (sR - sL);
    }
    
}


void field::UpdateFields(){
  using namespace flux;
  using namespace grid;

  int i;
  for(i=Is;i<=Ie;i++){
      rho[i]=   rho[i] - (numflux[i+1][mdn]-numflux[i][mdn])*dt/dx;
    rhovx[i]= rhovx[i] - (numflux[i+1][mvx]-numflux[i][mvx])*dt/dx;
       et[i]=    et[i] - (numflux[i+1][met]-numflux[i][met])*dt/dx;
  }

}



void dataio::OutputRadialProfile(){
  using namespace control;
  using namespace grid;
  using namespace field;

  FILE *ofile;
  int i;
  char outfile[20];
  FILE       *fout;
  int          ret;

  if(tsim > tout+dtdmp){
  cout << "t="<<tsim <<" dt="<<dt<< endl;
  ret=sprintf(outfile,"snap/t%05d.dat",index);
  ofile = fopen(outfile,"w");
  fprintf(ofile,  "# %12.7e\n",tsim);
  fprintf(ofile,  "# %12s %12s %12s %12s\n","x[cm]","rho[g/cm^3]","p[erg/cm^3]","v_x[cm/s]");
  for(i=Is;i<=Ie;i++){
    fprintf(ofile,"  %12.5e %12.5e %12.5e %12.5e \n",xb[i],rho[i],p[i],vx[i]);
  }
  fclose(ofile);
  index=index+1;
  tout = tsim;
  }
}