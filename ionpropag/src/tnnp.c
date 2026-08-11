/*29:*/
#line 829 "tnnp.web"

// tnnp.web $Rev$ 

#define TNNP_C
#define TAB_ROWS 500
#define TAB_LO -100.0
#define TAB_HI 50.0
#define TAB_STEP ((TAB_HI-TAB_LO) /(TAB_ROWS-1) )  \

#define test_param(source_nr,loc) if(source_nr) ifo->param[(loc) = ifo->Nparam++]= (source_nr) ; \

#define NSTATS 33
#define parm(s) (tpar->s[celltype])  \

#define CVcF (ttab_CVcF[celltype]) 
#define Kupsquare (ttab_Kupsquare[celltype]) 
#define exptaufca (ttab_exptaufca[celltype]) 
#define exptaug (ttab_exptaug[celltype]) 
#define Qt (ttab_Qt[celltype])  \

#define sm cell_status[0]
#define sh cell_status[1]
#define sj cell_status[2]
#define sxr1 cell_status[3]
#define sxr2 cell_status[4]
#define sxs cell_status[5]
#define ss cell_status[6]
#define sr cell_status[7]
#define sd cell_status[8] \

#define sf cell_status[9] \

#define sfca cell_status[10] \

#define sg cell_status[11]
#define Cai cell_status[12]
#define CaSR cell_status[13]
#define Nai cell_status[14]
#define Ki cell_status[15]
#define frozen cell_status[16] \
 \

#define os_E_Ka cell_status[TNNP_NSVAR+0]
#define os_E_Na cell_status[TNNP_NSVAR+1]
#define os_E_Ks cell_status[TNNP_NSVAR+2]
#define os_E_Ca cell_status[TNNP_NSVAR+3]
#define os_rec_IK1 cell_status[TNNP_NSVAR+4]
#define os_rec_ipK cell_status[TNNP_NSVAR+5]
#define os_rec_iNaK cell_status[TNNP_NSVAR+6]
#define os_FCa_inf cell_status[TNNP_NSVAR+7]
#define os_G_inf cell_status[TNNP_NSVAR+8] \

#define TABI (1.0/TAB_STEP) 
#define interp_tab(tab) (d1*tab[row][celltype]+d0*tab[row+1][celltype])  \


#line 833 "tnnp.web"

/*42:*/
#line 1184 "tnnp.web"

#include <limits.h> 
#include <math.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>         
#ifdef MATLAB_MEX_FILE
#undef strdup             
#endif
#line 1193 "tnnp.web"
#include <time.h> 

#include "membranes/membrane.h"
#include "membranes/membrane_internal.h"
#include "membranes/tnnp_p.h"    
#include "membranes/tnnp.h"



/*:42*/
#line 834 "tnnp.web"

/*35:*/
#line 1049 "tnnp.web"


static char*statnames[NSTATS];


/*:35*//*44:*/
#line 1226 "tnnp.web"


static TNNP_param*tpar;

static int Ko_param_loc,GKs_param_loc,GKr_param_loc,GNa_param_loc,
Gto_param_loc;


static int freeze_param_loc;


/*:44*//*47:*/
#line 1287 "tnnp.web"


static double ttab_CVcF[TNNP_NTYPES];
static double ttab_Kupsquare[TNNP_NTYPES];
static double ttab_exptaufca[TNNP_NTYPES];
static double ttab_exptaug[TNNP_NTYPES];
static double ttab_Qt[TNNP_NTYPES];


/*:47*//*48:*/
#line 1300 "tnnp.web"


static double mtab_exp_M[TAB_ROWS][TNNP_NTYPES];
static double mtab_exp_H[TAB_ROWS][TNNP_NTYPES];
static double mtab_exp_J[TAB_ROWS][TNNP_NTYPES];
static double mtab_exp_F[TAB_ROWS][TNNP_NTYPES];
static double mtab_exp_Xr1[TAB_ROWS][TNNP_NTYPES];
static double mtab_exp_Xr2[TAB_ROWS][TNNP_NTYPES];
static double mtab_exp_Xs[TAB_ROWS][TNNP_NTYPES];
static double mtab_exp_R[TAB_ROWS][TNNP_NTYPES];
static double mtab_exp_S[TAB_ROWS][TNNP_NTYPES];
static double mtab_exp_D[TAB_ROWS][TNNP_NTYPES];
static double mtab_M_inf[TAB_ROWS][TNNP_NTYPES];
static double mtab_H_inf[TAB_ROWS][TNNP_NTYPES];
static double mtab_J_inf[TAB_ROWS][TNNP_NTYPES];
static double mtab_F_inf[TAB_ROWS][TNNP_NTYPES];
static double mtab_Xr1_inf[TAB_ROWS][TNNP_NTYPES];
static double mtab_Xr2_inf[TAB_ROWS][TNNP_NTYPES];
static double mtab_Xs_inf[TAB_ROWS][TNNP_NTYPES];
static double mtab_R_inf[TAB_ROWS][TNNP_NTYPES];
static double mtab_S_inf[TAB_ROWS][TNNP_NTYPES];
static double mtab_D_inf[TAB_ROWS][TNNP_NTYPES];




/*:48*/
#line 835 "tnnp.web"

/*54:*/
#line 1498 "tnnp.web"


static inline double pow3(double x)
{
return(x*x*x);
}




/*:54*/
#line 836 "tnnp.web"


double tnnp_step(double vm,
yyy_t*__restrict__ cell_status,
Membrane_cell_info*ctype,
double Isd,
double dt,
float*dtime,
float simtime,
long elm)
{

/*45:*/
#line 1240 "tnnp.web"


int celltype;

double INa;
double ICaL;
double Ito;
double IKs;
double IKr;
double IK1;
double INaCa;
double INaK;
double IbNa,IbCa,IbK;
double IpCa;
double IpK;
double Frel= 0.0;
double Fleak= 0.0;
double Fup= 0.0;
double frozen_until;


/*:45*//*46:*/
#line 1266 "tnnp.web"


double E_Ka,E_Na,E_Ks,E_Ca;

double rec_IK1,rec_ipK,rec_iNaK;


double FCa_inf,G_inf;



/*:46*/
#line 848 "tnnp.web"

double Itotal,mh0,mh1;

celltype= ctype->ccode;
/*68:*/
#line 2015 "tnnp.web"

{
if(vm<TAB_LO){
Warning(2,"clipped vm=%f to %f at %s",vm,TAB_LO,node_to_string(elm));
vm= TAB_LO;
}
else if(vm> TAB_HI-TAB_STEP){
Warning(2,"clipped vm=%f to %f at %s",vm,TAB_HI,node_to_string(elm));
vm= TAB_HI-TAB_STEP;
}
}


/*:68*/
#line 852 "tnnp.web"

/*51:*/
#line 1379 "tnnp.web"

{
double R= 8.314472;
double F= 96485.3415;
double T= parm(temperature)+273.15;
double RToverF= 1000.*R*T/F;
double VmFoverRT= vm/RToverF;
double cKo;
double cGKs;
double cGKr;
double cGNa;
double cGto;

if(tpar->Ko_param)cKo= ctype->param[Ko_param_loc];
else cKo= parm(Ko);

if(tpar->GKs_param)cGKs= ctype->param[GKs_param_loc];
else cGKs= parm(GKs);

if(tpar->GKr_param)cGKr= ctype->param[GKr_param_loc];
else cGKr= parm(GKr);

if(tpar->GNa_param)cGNa= ctype->param[GNa_param_loc];
else cGNa= parm(GNa);

if(tpar->Gto_param)cGto= ctype->param[Gto_param_loc];
else cGto= parm(Gto);

if(tpar->freeze_param)frozen_until= ctype->param[freeze_param_loc];
else frozen_until= 0;
if(simtime<frozen_until)dt= 0.0;

/*52:*/
#line 1422 "tnnp.web"

{
double Ak1,Bk1;
long ts,te;

if(dt==0.0){
ts= te= 1;
}else{
ts= 0.5+simtime/dt;
te= 0.5+tpar->eint/dt;
}

if(tpar->eint==0||ts%te==0){
E_Ka= RToverF*log(cKo/Ki);
E_Na= RToverF*log(parm(Nao)/Nai);
E_Ks= RToverF*log((cKo+parm(pKNa)*parm(Nao))/(Ki+parm(pKNa)*Nai));
E_Ca= 0.5*RToverF*log(parm(Cao)/Cai);
Ak1= 0.1/(1.+exp(0.06*(vm-E_Ka-200)));
Bk1= (3.*exp(0.0002*(vm-E_Ka+100))+
exp(0.1*(vm-E_Ka-10)))/(1.+exp(-0.5*(vm-E_Ka)));
rec_IK1= Ak1/(Ak1+Bk1);
rec_iNaK= (1./(1.+0.1245*exp(-0.1*VmFoverRT)+0.0353*exp(-VmFoverRT)));
rec_ipK= 1./(1.+exp((25-vm)/5.98));
}else{
E_Ka= os_E_Ka;
E_Na= os_E_Na;
E_Ks= os_E_Ks;
E_Ca= os_E_Ca;
rec_IK1= os_rec_IK1;
rec_ipK= os_rec_ipK;
rec_iNaK= os_rec_iNaK;
}
if(tpar->eint!=0&&ts%te==0){
os_E_Ka= E_Ka;
os_E_Na= E_Na;
os_E_Ks= E_Ks;
os_E_Ca= E_Ca;
os_rec_IK1= rec_IK1;
os_rec_ipK= rec_ipK;
os_rec_iNaK= rec_iNaK;
}
}


/*:52*/
#line 1411 "tnnp.web"

/*53:*/
#line 1471 "tnnp.web"

{
double exv= exp(2*VmFoverRT);
double env= exp((parm(n)-1)*VmFoverRT);

INa= cGNa*pow3(sm)*sh*sj*(vm-E_Na);
ICaL= parm(GCaL)*sd*sf*sfca*4*(F*VmFoverRT)*
(exv*Cai-0.341*parm(Cao))/(exv-1.);
Ito= cGto*sr*ss*(vm-E_Ka);
IKr= cGKr*sqrt(cKo/5.4)*sxr1*sxr2*(vm-E_Ka);
IKs= cGKs*sxs*sxs*(vm-E_Ks);
IK1= parm(GK1)*rec_IK1*(vm-E_Ka);
INaCa= parm(knaca)*(1./(pow3(parm(KmNai))+pow3(parm(Nao))))*
(1./(parm(KmCa)+parm(Cao)))*
(1./(1+parm(ksat)*env))*
(exp(parm(n)*VmFoverRT)*pow3(Nai)*parm(Cao)-
env*pow3(parm(Nao))*Cai*2.5);
INaK= parm(knak)*(cKo/(cKo+parm(KmK)))*
(Nai/(Nai+parm(KmNa)))*rec_iNaK;
IpCa= parm(GpCa)*Cai/(parm(KpCa)+Cai);
IpK= parm(GpK)*rec_ipK*(vm-E_Ka);
IbNa= parm(GbNa)*(vm-E_Na);
IbCa= parm(GbCa)*(vm-E_Ca);
IbK= parm(GbK)*(vm-E_Ka);
}


/*:53*/
#line 1412 "tnnp.web"

}


/*:51*/
#line 853 "tnnp.web"

Itotal= IKr+IKs+IK1+Ito+INa+IbNa+ICaL+IbCa+
INaK+INaCa+IpCa+IpK+IbK;
if(!finite(Itotal)){
fprintf(stderr,"IKr   = %f\n",IKr);
fprintf(stderr,"IKs   = %f\n",IKs);
fprintf(stderr,"IK1   = %f\n",IK1);
fprintf(stderr,"Ito   = %f\n",Ito);
fprintf(stderr,"IbK   = %f\n",IbK);
fprintf(stderr,"INa   = %f\n",INa);
fprintf(stderr,"IbNa  = %f\n",IbNa);
fprintf(stderr,"ICaL  = %f\n",ICaL);
fprintf(stderr,"IbCa  = %f\n",IbCa);
fprintf(stderr,"INaK  = %f\n",INaK);
fprintf(stderr,"INaCa = %f\n",INaCa);
fprintf(stderr,"IpCa  = %f\n",IpCa);
fprintf(stderr,"IpK   = %f\n",IpK);
fprintf(stderr,"Istim = %f\n",Isd);
Error(1,"tnnp_step: Itotal=%f at %s\n",Itotal,node_to_string(elm));
}
mh0= sm*sh;
if(dt> 0.0){
/*55:*/
#line 1514 "tnnp.web"

{
double Xr1_inf,Xr2_inf,Xs_inf,R_inf,S_inf,D_inf;
double F_inf,M_inf,H_inf,J_inf;

double exp_M,exp_H,exp_J,exp_Xr1,exp_Xr2,exp_D,exp_F;
double exp_Xs,exp_S,exp_R;

if(frozen&&vm> -50.0)frozen= 0;
if(!frozen)/*57:*/
#line 1592 "tnnp.web"

{
double sumK;
double dNai,dKi;
double FCaM;
double FCaSR;

FCaM= -(ICaL+IbCa+IpCa-2*INaCa)*0.5*CVcF;
/*58:*/
#line 1625 "tnnp.web"

{
double CaSRsquare,A;

if(CaSR<=0.0){
Frel= 0;
CaSR= 0.0;
}else{
CaSRsquare= CaSR*CaSR;
A= 0.016464*CaSRsquare/(0.0625+CaSRsquare)+0.008232;
Frel= A*sd*sg;
}
}



/*:58*/
#line 1600 "tnnp.web"

Fleak= 0.00008*(CaSR-Cai);
Fup= parm(Fmaxup)/(1.+(Kupsquare/(Cai*Cai)));
FCaSR= Fup-Frel-Fleak;
/*60:*/
#line 1700 "tnnp.web"

{
double b,c;
double total;
double delta_CaSR;
double CaCSQN;

CaCSQN= parm(Bufsr)*CaSR/(CaSR+parm(Kbufsr));
delta_CaSR= dt*(parm(Vc)/parm(Vsr))*FCaSR;
total= CaCSQN+delta_CaSR+CaSR;

b= parm(Bufsr)-total+parm(Kbufsr);
c= parm(Kbufsr)*total;
CaSR= 0.5*(sqrt(b*b+4*c)-b);
}


/*:60*/
#line 1604 "tnnp.web"

/*59:*/
#line 1681 "tnnp.web"

{
double b,c;
double Cai_total;
double delta_Cai;
double CaBuf;

CaBuf= parm(Bufc)*Cai/(Cai+parm(Kbufc));
delta_Cai= dt*(FCaM-FCaSR);
Cai_total= CaBuf+Cai+delta_Cai;
b= parm(Bufc)-Cai_total+parm(Kbufc);
c= parm(Kbufc)*Cai_total;
Cai= 0.5*(sqrt(b*b+4*c)-b);
}


/*:59*/
#line 1605 "tnnp.web"



dNai= -(INa+IbNa+3*INaK+3*INaCa+parm(pKNa)*IKs)*CVcF;
Nai+= dt*dNai;


sumK= IK1+Ito+IKr+(1.0-parm(pKNa))*IKs-2*INaK+IpK+IbK;
if(tpar->use_Isd)sumK+= Isd;
dKi= -sumK*CVcF;
Ki+= dt*dKi;
}



/*:57*/
#line 1523 "tnnp.web"

/*67:*/
#line 1976 "tnnp.web"

{
int row;
double m,d0,d1;

m= (vm-TAB_LO)*TABI;
row= m;
d0= m-row;
d1= 1.0-d0;
if(row<0||row>=TAB_ROWS-1)Error(1,"this can't happen");

exp_M= interp_tab(mtab_exp_M);
exp_H= interp_tab(mtab_exp_H);
exp_J= interp_tab(mtab_exp_J);
exp_Xr1= interp_tab(mtab_exp_Xr1);
exp_Xr2= interp_tab(mtab_exp_Xr2);
exp_Xs= interp_tab(mtab_exp_Xs);
exp_R= interp_tab(mtab_exp_R);
exp_S= interp_tab(mtab_exp_S);
exp_D= interp_tab(mtab_exp_D);
exp_F= interp_tab(mtab_exp_F);

M_inf= interp_tab(mtab_M_inf);
H_inf= interp_tab(mtab_H_inf);
J_inf= H_inf;
Xr1_inf= interp_tab(mtab_Xr1_inf);
Xr2_inf= interp_tab(mtab_Xr2_inf);
Xs_inf= interp_tab(mtab_Xs_inf);
R_inf= interp_tab(mtab_R_inf);
S_inf= interp_tab(mtab_S_inf);
D_inf= interp_tab(mtab_D_inf);
F_inf= interp_tab(mtab_F_inf);
}


/*:67*/
#line 1524 "tnnp.web"

/*66:*/
#line 1943 "tnnp.web"

{
long ts,te;

ts= 0.5+simtime/dt;
te= 0.5+tpar->eint/dt;
if(tpar->eint==0||ts%te==0){
FCa_inf= (1./(1.+pow((Cai/0.000325),8))+
0.1/(1.+exp((Cai-0.0005)/0.0001))+
0.20/(1.+exp((Cai-0.00075)/0.0008))+
0.23)/1.46;
if(Cai<0.00035)G_inf= 1./(1.+pow((Cai/0.00035),6));
else G_inf= 1./(1.+pow((Cai/0.00035),16));
}else{
FCa_inf= os_FCa_inf;
G_inf= os_G_inf;
}
if(tpar->eint!=0&&ts%te==0){
os_FCa_inf= FCa_inf;
os_G_inf= G_inf;
}
}


/*:66*/
#line 1525 "tnnp.web"

/*69:*/
#line 2040 "tnnp.web"

{
double fcaold,gold;

sm= M_inf-(M_inf-sm)*exp_M;

if(!frozen){
sh= H_inf-(H_inf-sh)*exp_H;
sj= J_inf-(J_inf-sj)*exp_J;
sxr1= Xr1_inf-(Xr1_inf-sxr1)*exp_Xr1;
sxr2= Xr2_inf-(Xr2_inf-sxr2)*exp_Xr2;
sxs= Xs_inf-(Xs_inf-sxs)*exp_Xs;
ss= S_inf-(S_inf-ss)*exp_S;
sr= R_inf-(R_inf-sr)*exp_R;
sd= D_inf-(D_inf-sd)*exp_D;

sf= F_inf-(F_inf-sf)*exp_F;


fcaold= sfca;
sfca= FCa_inf-(FCa_inf-sfca)*exptaufca;

#if 1
if(sfca> fcaold&&(vm)> -60)sfca= fcaold;

#endif
#line 2066 "tnnp.web"

gold= sg;
sg= G_inf-(G_inf-sg)*exptaug;

#if 1
if(sg> gold&&(vm)> -60)sg= gold;
#endif
#line 2073 "tnnp.web"
}
}


/*:69*/
#line 1526 "tnnp.web"

}


/*:55*/
#line 875 "tnnp.web"

/*70:*/
#line 2107 "tnnp.web"

{
double threshold= 0.05;

mh1= sm*sh;
if(mh0<=threshold&&mh1> threshold){
*dtime= simtime+dt*(threshold-mh0)/(mh1-mh0);
}
}



/*:70*/
#line 876 "tnnp.web"

}
return Itotal;
}


/*:29*//*31:*/
#line 914 "tnnp.web"


void tnnp_init(double dt)
{
int celltype;
int row;
double vm;

/*36:*/
#line 1061 "tnnp.web"

{
#define X(nr, name, value) statnames[nr] =  strdup(name)
/*38:*/
#line 1079 "tnnp.web"

{
X(0,"m",sm);
X(1,"h",sh);
X(2,"j",sj);
X(3,"xr1",sxr1);
X(4,"xr2",sxr2);
X(5,"xs",sxs);
X(6,"s",ss);
X(7,"r",sr);
X(8,"d",sd);
X(9,"f",sf);
X(10,"fca",sfca);
X(11,"g",sg);
X(12,"Cai",Cai);
X(13,"CaSR",CaSR);
X(14,"Nai",Nai);
X(15,"Ki",Ki);
X(16,"frozen",frozen);
X(17+0,"INa",INa);
X(17+1,"IKr",IKr);
X(17+2,"IKs",IKs);
X(17+3,"IK1",IK1);
X(17+4,"Ito",Ito);
X(17+5,"IbNa",IbNa);
X(17+6,"ICaL",ICaL);
X(17+7,"IbCa",IbCa);
X(17+8,"INaCa",INaCa);
X(17+9,"IpCa",IpCa);
X(17+10,"IpK",IpK);
X(17+11,"INaK",INaK);
X(17+12,"IbK",IbK);
X(17+13,"Frel",Frel);
X(17+14,"Fleak",Fleak);
X(17+15,"Fup",Fup);
}


/*:38*/
#line 1064 "tnnp.web"

#undef X
}


/*:36*/
#line 922 "tnnp.web"

for(celltype= 0;celltype<TNNP_NTYPES;celltype++){
/*56:*/
#line 1537 "tnnp.web"

{
double F= 96485.3415;
double Q10= 3.0;

Qt= pow(Q10,(parm(temperature)-37)*0.1);

CVcF= parm(capacitance)/(parm(Vc)*F);


Kupsquare= parm(Kup)*parm(Kup);

exptaufca= exp(-dt*Qt/(parm(taufca)));
exptaug= exp(-dt*Qt/(parm(taug)));
}


/*:56*/
#line 924 "tnnp.web"

for(row= 0;row<TAB_ROWS;row++){
vm= TAB_LO+row*TAB_STEP;
/*62:*/
#line 1756 "tnnp.web"

{
double am,bm,ah,bh,aj,bj,tau_m,tau_h,tau_j,tmpa;

am= 1./(1.+exp((-60.-vm)/5.));
bm= 0.1/(1.+exp((vm+35.)/5.))+0.10/(1.+exp((vm-50.)/200.));
tau_m= am*bm/Qt;
mtab_exp_M[row][celltype]= exp(-dt/tau_m);
tmpa= 1.+exp((parm(Vhalf_m)-vm)/9.03);
mtab_M_inf[row][celltype]= 1./(tmpa*tmpa);


if(vm>=-40.){
ah= 0.;
bh= (0.77/(0.13*(1.+exp(-(vm+10.66)/11.1))));
tau_h= 1.0/((ah+bh)*Qt);
}
else{
ah= (0.057*exp(-(vm+80.)/6.8));
bh= (2.7*exp(0.079*vm)+(3.1e5)*exp(0.3485*vm));
tau_h= 1.0/((ah+bh)*Qt);
}
mtab_exp_H[row][celltype]= exp(-dt/tau_h);
tmpa= 1.+exp((vm-parm(Vhalf_h))/7.43);
mtab_H_inf[row][celltype]= 1./(tmpa*tmpa);


if(vm>=-40.){
aj= 0.;
bj= (0.6*exp((0.057)*vm)/(1.+exp(-0.1*(vm+32.))));
tau_j= 1.0/((aj+bj)*Qt);
}
else{
aj= (((-2.5428e4)*exp(0.2444*vm)-(6.948e-6)*
exp(-0.04391*vm))*(vm+37.78)/
(1.+exp(0.311*(vm+79.23))));
bj= (0.02424*exp(-0.01052*vm)/(1.+exp(-0.1378*(vm+40.14))));
tau_j= 1.0/((aj+bj)*Qt);
}
mtab_exp_J[row][celltype]= exp(-dt/tau_j);
}


/*:62*//*63:*/
#line 1814 "tnnp.web"

{
double axr1,axr2,bxr1,bxr2,Axs,Bxs,tau_Xr1,tau_Xr2,tau_Xs;

mtab_Xr1_inf[row][celltype]= 1./(1.+exp((-26.-vm)/7.));
axr1= 450./(1.+exp((-45.-vm)/10.));
bxr1= 6./(1.+exp((vm-(-30.))/11.5));
tau_Xr1= axr1*bxr1/Qt;

mtab_Xr2_inf[row][celltype]= 1./(1.+exp((vm-(-88.))/24.));
axr2= 3./(1.+exp((-60.-vm)/20.));
bxr2= 1.12/(1.+exp((vm-60.)/20.));
tau_Xr2= axr2*bxr2/Qt;

mtab_Xs_inf[row][celltype]= 1./(1.+exp((-5.0-vm+parm(VshiftXs))/14.0));
Axs= 1100./(sqrt(1.+exp((-10.0-vm+parm(VshiftXs))/6)));
Bxs= 1./(1.+exp((vm-parm(VshiftXs)-60.0)/20.));
tau_Xs= parm(FtXs)*Axs*Bxs/Qt;

mtab_exp_Xr1[row][celltype]= exp(-dt/tau_Xr1);
mtab_exp_Xr2[row][celltype]= exp(-dt/tau_Xr2);
mtab_exp_Xs[row][celltype]= exp(-dt/tau_Xs);
}


/*:63*//*64:*/
#line 1849 "tnnp.web"

{
double r_inf= 0,s_inf= 0,tau_r= 0,tau_s= 0;


switch(parm(basetype)){
case TNNP_EPIC:case TNNP_RVEPI:
r_inf= 1./(1.+exp((20-vm)/6.));
s_inf= 1./(1.+exp((vm+20)/5.));
tau_r= (9.5*exp(-(vm+40.)*(vm+40.)/1800.)+0.8)/Qt;
tau_s= (85.*exp(-(vm+45.)*(vm+45.)/320.)+5./(1.+exp((vm-20.)/5.))+3.)/Qt;
break;
case TNNP_ENDO:
r_inf= 1./(1.+exp((20-vm)/6.));
s_inf= 1./(1.+exp((vm+28)/5.));
tau_r= (9.5*exp(-(vm+40.)*(vm+40.)/1800.)+0.8)/Qt;
tau_s= (1000.*exp(-(vm+67)*(vm+67)/1000.)+8.)/Qt;
break;
case TNNP_MCEL:case TNNP_RVMC:
r_inf= 1./(1.+exp((20-vm)/6.));
s_inf= 1./(1.+exp((vm+20)/5.));
tau_r= (9.5*exp(-(vm+40.)*(vm+40.)/1800.)+0.8)/Qt;
tau_s= (85.*exp(-(vm+45.)*(vm+45.)/320.)+5./(1.+exp((vm-20.)/5.))+3.)/Qt;
break;
case TNNP_T5:case TNNP_T6:
Error(1,"this can't happen");
break;
default:Error(1,"unknown cell type %d",celltype);
}
mtab_R_inf[row][celltype]= r_inf;
mtab_S_inf[row][celltype]= s_inf;
mtab_exp_R[row][celltype]= exp(-dt/tau_r);
mtab_exp_S[row][celltype]= exp(-dt/tau_s);
}


/*:64*//*65:*/
#line 1885 "tnnp.web"

{
double Ad,Bd,Cd,tau_d,tau_f;

mtab_D_inf[row][celltype]= 1./(1.+exp((-5-vm)/7.5));
Ad= 1.4/(1.+exp((-35-vm)/13))+0.25;
Bd= 1.4/(1.+exp((vm+5)/5));
Cd= 1./(1.+exp((50-vm)/20));
tau_d= (Ad*Bd+Cd)/Qt;
mtab_exp_D[row][celltype]= exp(-dt/tau_d);

mtab_F_inf[row][celltype]= 1./(1.+exp((vm+20)/7));
tau_f= (1125*exp(-(vm+27)*(vm+27)/240)+80+165/(1.+exp((25-vm)/10)))/Qt;
mtab_exp_F[row][celltype]= exp(-dt/tau_f);
}


/*:65*/
#line 927 "tnnp.web"

}
}
}



/*:31*//*32:*/
#line 945 "tnnp.web"


void tnnp_info(Membrane_cell_info*ct,
Membrane_info*ifo,
void*prm)
{
int i;
tpar= (TNNP_param*)prm;

ct->mcode= 0;
ct->ccode= TNNP_ENDO;
ct->param= NULL;

if(tpar->eint)ifo->Nsvar= TNNP_NSVAR+TNNP_NSVAR_EINT;
else ifo->Nsvar= TNNP_NSVAR;

ifo->Ntypes= TNNP_NTYPES;

ifo->Tname= malloc(ifo->Ntypes*sizeof(char*));
for(i= 0;i<ifo->Ntypes;i++)ifo->Tname[i]= strdup(tnnp_cname[i]);

/*33:*/
#line 989 "tnnp.web"

{
int max_param= 6;

ifo->Nparam= 0;
ifo->param= calloc(max_param,sizeof(int));
test_param(tpar->Ko_param,Ko_param_loc);
test_param(tpar->GKs_param,GKs_param_loc);
test_param(tpar->GKr_param,GKr_param_loc);
test_param(tpar->GNa_param,GNa_param_loc);
test_param(tpar->Gto_param,Gto_param_loc);
test_param(tpar->freeze_param,freeze_param_loc);
if(ifo->Nparam> max_param)Error(1,"coding error");
}


/*:33*/
#line 966 "tnnp.web"


ifo->info= &tnnp_info;
ifo->init= &tnnp_init;
ifo->infinite= &tnnp_infinite;
ifo->step= &tnnp_step;
ifo->get_status= &tnnp_get_status;
ifo->ccode= &tnnp_ccode;

}



/*:32*//*34:*/
#line 1010 "tnnp.web"


void tnnp_infinite(Membrane_cell_info*ctype,vm_t*Vm,
yyy_t*cell_status,int inittype)
{
int celltype= ctype->ccode;

frozen= tpar->freeze;
if(inittype==MEMBRANE_INIT_FROM_STATE)return;

*Vm= parm(V_init);
sm= 0.0;
sh= 0.775;
sj= 0.765;
sxr1= 0.0;
sxr2= 0.483;
sxs= 0.0;
ss= parm(ss_init);
sr= 0.0;
sd= 0.0;
sf= 0.988;
sfca= 1.0;
sg= 1.0;
Cai= parm(Cai_init);
CaSR= parm(CaSR_init);
Nai= parm(Nai_init);
Ki= parm(Ki_init);
}



/*:34*//*39:*/
#line 1130 "tnnp.web"


void tnnp_get_status(float vm,yyy_t*cell_status,Membrane_cell_info*ctype,
float*Stats,int*Nstats,char***Names)
{
/*45:*/
#line 1240 "tnnp.web"


int celltype;

double INa;
double ICaL;
double Ito;
double IKs;
double IKr;
double IK1;
double INaCa;
double INaK;
double IbNa,IbCa,IbK;
double IpCa;
double IpK;
double Frel= 0.0;
double Fleak= 0.0;
double Fup= 0.0;
double frozen_until;


/*:45*//*46:*/
#line 1266 "tnnp.web"


double E_Ka,E_Na,E_Ks,E_Ca;

double rec_IK1,rec_ipK,rec_iNaK;


double FCa_inf,G_inf;



/*:46*/
#line 1135 "tnnp.web"

double dt= 0.0,simtime= 1.0;
long elm= 0;

celltype= ctype->ccode;
/*68:*/
#line 2015 "tnnp.web"

{
if(vm<TAB_LO){
Warning(2,"clipped vm=%f to %f at %s",vm,TAB_LO,node_to_string(elm));
vm= TAB_LO;
}
else if(vm> TAB_HI-TAB_STEP){
Warning(2,"clipped vm=%f to %f at %s",vm,TAB_HI,node_to_string(elm));
vm= TAB_HI-TAB_STEP;
}
}


/*:68*/
#line 1140 "tnnp.web"

/*51:*/
#line 1379 "tnnp.web"

{
double R= 8.314472;
double F= 96485.3415;
double T= parm(temperature)+273.15;
double RToverF= 1000.*R*T/F;
double VmFoverRT= vm/RToverF;
double cKo;
double cGKs;
double cGKr;
double cGNa;
double cGto;

if(tpar->Ko_param)cKo= ctype->param[Ko_param_loc];
else cKo= parm(Ko);

if(tpar->GKs_param)cGKs= ctype->param[GKs_param_loc];
else cGKs= parm(GKs);

if(tpar->GKr_param)cGKr= ctype->param[GKr_param_loc];
else cGKr= parm(GKr);

if(tpar->GNa_param)cGNa= ctype->param[GNa_param_loc];
else cGNa= parm(GNa);

if(tpar->Gto_param)cGto= ctype->param[Gto_param_loc];
else cGto= parm(Gto);

if(tpar->freeze_param)frozen_until= ctype->param[freeze_param_loc];
else frozen_until= 0;
if(simtime<frozen_until)dt= 0.0;

/*52:*/
#line 1422 "tnnp.web"

{
double Ak1,Bk1;
long ts,te;

if(dt==0.0){
ts= te= 1;
}else{
ts= 0.5+simtime/dt;
te= 0.5+tpar->eint/dt;
}

if(tpar->eint==0||ts%te==0){
E_Ka= RToverF*log(cKo/Ki);
E_Na= RToverF*log(parm(Nao)/Nai);
E_Ks= RToverF*log((cKo+parm(pKNa)*parm(Nao))/(Ki+parm(pKNa)*Nai));
E_Ca= 0.5*RToverF*log(parm(Cao)/Cai);
Ak1= 0.1/(1.+exp(0.06*(vm-E_Ka-200)));
Bk1= (3.*exp(0.0002*(vm-E_Ka+100))+
exp(0.1*(vm-E_Ka-10)))/(1.+exp(-0.5*(vm-E_Ka)));
rec_IK1= Ak1/(Ak1+Bk1);
rec_iNaK= (1./(1.+0.1245*exp(-0.1*VmFoverRT)+0.0353*exp(-VmFoverRT)));
rec_ipK= 1./(1.+exp((25-vm)/5.98));
}else{
E_Ka= os_E_Ka;
E_Na= os_E_Na;
E_Ks= os_E_Ks;
E_Ca= os_E_Ca;
rec_IK1= os_rec_IK1;
rec_ipK= os_rec_ipK;
rec_iNaK= os_rec_iNaK;
}
if(tpar->eint!=0&&ts%te==0){
os_E_Ka= E_Ka;
os_E_Na= E_Na;
os_E_Ks= E_Ks;
os_E_Ca= E_Ca;
os_rec_IK1= rec_IK1;
os_rec_ipK= rec_ipK;
os_rec_iNaK= rec_iNaK;
}
}


/*:52*/
#line 1411 "tnnp.web"

/*53:*/
#line 1471 "tnnp.web"

{
double exv= exp(2*VmFoverRT);
double env= exp((parm(n)-1)*VmFoverRT);

INa= cGNa*pow3(sm)*sh*sj*(vm-E_Na);
ICaL= parm(GCaL)*sd*sf*sfca*4*(F*VmFoverRT)*
(exv*Cai-0.341*parm(Cao))/(exv-1.);
Ito= cGto*sr*ss*(vm-E_Ka);
IKr= cGKr*sqrt(cKo/5.4)*sxr1*sxr2*(vm-E_Ka);
IKs= cGKs*sxs*sxs*(vm-E_Ks);
IK1= parm(GK1)*rec_IK1*(vm-E_Ka);
INaCa= parm(knaca)*(1./(pow3(parm(KmNai))+pow3(parm(Nao))))*
(1./(parm(KmCa)+parm(Cao)))*
(1./(1+parm(ksat)*env))*
(exp(parm(n)*VmFoverRT)*pow3(Nai)*parm(Cao)-
env*pow3(parm(Nao))*Cai*2.5);
INaK= parm(knak)*(cKo/(cKo+parm(KmK)))*
(Nai/(Nai+parm(KmNa)))*rec_iNaK;
IpCa= parm(GpCa)*Cai/(parm(KpCa)+Cai);
IpK= parm(GpK)*rec_ipK*(vm-E_Ka);
IbNa= parm(GbNa)*(vm-E_Na);
IbCa= parm(GbCa)*(vm-E_Ca);
IbK= parm(GbK)*(vm-E_Ka);
}


/*:53*/
#line 1412 "tnnp.web"

}


/*:51*/
#line 1141 "tnnp.web"

/*37:*/
#line 1072 "tnnp.web"

{
#define X(nr, name, value)  Stats[nr] =  value
/*38:*/
#line 1079 "tnnp.web"

{
X(0,"m",sm);
X(1,"h",sh);
X(2,"j",sj);
X(3,"xr1",sxr1);
X(4,"xr2",sxr2);
X(5,"xs",sxs);
X(6,"s",ss);
X(7,"r",sr);
X(8,"d",sd);
X(9,"f",sf);
X(10,"fca",sfca);
X(11,"g",sg);
X(12,"Cai",Cai);
X(13,"CaSR",CaSR);
X(14,"Nai",Nai);
X(15,"Ki",Ki);
X(16,"frozen",frozen);
X(17+0,"INa",INa);
X(17+1,"IKr",IKr);
X(17+2,"IKs",IKs);
X(17+3,"IK1",IK1);
X(17+4,"Ito",Ito);
X(17+5,"IbNa",IbNa);
X(17+6,"ICaL",ICaL);
X(17+7,"IbCa",IbCa);
X(17+8,"INaCa",INaCa);
X(17+9,"IpCa",IpCa);
X(17+10,"IpK",IpK);
X(17+11,"INaK",INaK);
X(17+12,"IbK",IbK);
X(17+13,"Frel",Frel);
X(17+14,"Fleak",Fleak);
X(17+15,"Fup",Fup);
}


/*:38*/
#line 1075 "tnnp.web"

#undef X
}

/*:37*/
#line 1142 "tnnp.web"

/*40:*/
#line 1153 "tnnp.web"

{
int i;
for(i= 0;i<NSTATS;i++){
if(!(finite(Stats[i]))){
Error(1,"tnnp_get_status: %s=%f",statnames[i],Stats[i]);
}
}
}


/*:40*/
#line 1143 "tnnp.web"

*Nstats= NSTATS;
*Names= statnames;
}


/*:39*//*41:*/
#line 1168 "tnnp.web"


int tnnp_ccode(char*name)
{
int i;

for(i= 0;i<TNNP_NTYPES;i++){
if(strcmp(name,tnnp_cname[i])==0)return i;
}
Error(1,"I don't know about cell type \"%s\"",name);
return-1;
}


/*:41*/
