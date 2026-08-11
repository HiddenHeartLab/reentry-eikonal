/*33:*/
#line 984 "tp06.web"

// tp06.web $Rev$ 

#define TP06_C
#define TAB_ROWS 500
#define TAB_LO -100.0
#define TAB_HI 50.0
#define TAB_STEP ((TAB_HI-TAB_LO) /(TAB_ROWS-1) )  \

#define NVARS 43
#define parm(s) (tpar->s[celltype])  \

#define C_over_Vc_F (ttab_C_over_Vc_F[celltype]) 
#define C_over_Vss_F (ttab_C_over_Vss_F[celltype]) 
#define Kupsquare (ttab_Kupsquare[celltype]) 
#define Qt (ttab_Qt[celltype]) 
#define Nao3 (ttab_pow3_Nao[celltype])  \

#define sm cell_status[0]
#define sh cell_status[1]
#define sj cell_status[2]
#define sxr1 cell_status[3]
#define sxr2 cell_status[4]
#define sxs cell_status[5]
#define ss cell_status[6]
#define sr cell_status[7] \

#define sd cell_status[8]
#define sf cell_status[9]
#define f2 cell_status[10]
#define f_CaSS cell_status[11] \

#define k_CaSR cell_status[12]
#define k_1 cell_status[13]
#define k_2 cell_status[14]
#define R_rel cell_status[15]
#define O_rel cell_status[16] \

#define Cai cell_status[17]
#define CaSR cell_status[18]
#define CaSS cell_status[19]
#define Nai cell_status[20]
#define Ki cell_status[21]
#define frozen cell_status[22] \
 \

#define os_E_Ka cell_status[TP06_NSVAR+0]
#define os_E_Na cell_status[TP06_NSVAR+1]
#define os_E_Ks cell_status[TP06_NSVAR+2]
#define os_E_Ca cell_status[TP06_NSVAR+3]
#define os_rec_IK1 cell_status[TP06_NSVAR+4]
#define os_rec_ipK cell_status[TP06_NSVAR+5]
#define os_rec_iNaK cell_status[TP06_NSVAR+6] \

#define TABI (1.0/TAB_STEP) 
#define interp_tab(tab) (d1*tab[row][celltype]+d0*tab[row+1][celltype])  \


#line 988 "tp06.web"

/*45:*/
#line 1358 "tp06.web"

#include <limits.h> 
#include <float.h> 
#include <math.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>         
#ifdef MATLAB_MEX_FILE
#undef strdup             
#endif
#line 1368 "tp06.web"
#include <time.h> 

#include "membranes/membrane.h"
#include "membranes/membrane_internal.h"
#include "membranes/tp06_p.h"    
#include "membranes/tp06.h"



/*:45*/
#line 989 "tp06.web"

/*38:*/
#line 1218 "tp06.web"


static char*statnames[NVARS];


/*:38*//*47:*/
#line 1401 "tp06.web"


static TP06_param*tpar;

static int Ko_param_loc,GKs_param_loc,GKr_param_loc,GNa_param_loc,
GCaL_param_loc;




/*:47*//*50:*/
#line 1459 "tp06.web"


static double Faraday= 96485.3415;
static double Q10= 3.0;

static double ttab_C_over_Vc_F[TP06_NTYPES];
static double ttab_C_over_Vss_F[TP06_NTYPES];
static double ttab_Kupsquare[TP06_NTYPES];
static double ttab_Qt[TP06_NTYPES];
static double ttab_pow3_Nao[TP06_NTYPES];


/*:50*//*51:*/
#line 1475 "tp06.web"


static double mtab_exp_M[TAB_ROWS][TP06_NTYPES];
static double mtab_exp_H[TAB_ROWS][TP06_NTYPES];
static double mtab_exp_J[TAB_ROWS][TP06_NTYPES];
static double mtab_exp_F[TAB_ROWS][TP06_NTYPES];
static double mtab_exp_F2[TAB_ROWS][TP06_NTYPES];
static double mtab_exp_Xr1[TAB_ROWS][TP06_NTYPES];
static double mtab_exp_Xr2[TAB_ROWS][TP06_NTYPES];
static double mtab_exp_Xs[TAB_ROWS][TP06_NTYPES];
static double mtab_exp_R[TAB_ROWS][TP06_NTYPES];
static double mtab_exp_S[TAB_ROWS][TP06_NTYPES];
static double mtab_exp_D[TAB_ROWS][TP06_NTYPES];
static double mtab_M_inf[TAB_ROWS][TP06_NTYPES];
static double mtab_HJ_inf[TAB_ROWS][TP06_NTYPES];
static double mtab_F_inf[TAB_ROWS][TP06_NTYPES];
static double mtab_F2_inf[TAB_ROWS][TP06_NTYPES];
static double mtab_Xr1_inf[TAB_ROWS][TP06_NTYPES];
static double mtab_Xr2_inf[TAB_ROWS][TP06_NTYPES];
static double mtab_Xs_inf[TAB_ROWS][TP06_NTYPES];
static double mtab_R_inf[TAB_ROWS][TP06_NTYPES];
static double mtab_S_inf[TAB_ROWS][TP06_NTYPES];
static double mtab_D_inf[TAB_ROWS][TP06_NTYPES];




/*:51*/
#line 990 "tp06.web"

/*58:*/
#line 1695 "tp06.web"


static inline double pow3(double x)
{
return(x*x*x);
}




/*:58*/
#line 991 "tp06.web"


double tp06_step(double vm,
yyy_t*cell_status,
Membrane_cell_info*ctype,
double Isd,
double dt,
float*dtime,
float simtime,
long elm)
{
/*48:*/
#line 1414 "tp06.web"


int celltype;

double INa;
double ICaL;
double Ito;
double IKs;
double IKr;
double IK1;
double INaCa;
double INaK;
double IbNa,IbCa;
double IpCa;
double IpK;
double Frel= 0.0;
double Fleak= 0.0;
double Fup= 0.0;
double Fxfer= 0.0;


/*:48*//*49:*/
#line 1439 "tp06.web"


double E_Ka,E_Na,E_Ks,E_Ca;

double rec_IK1,rec_ipK,rec_iNaK;





/*:49*/
#line 1002 "tp06.web"

double Itotal;
double mh0,mh1;

celltype= ctype->ccode;
/*77:*/
#line 2253 "tp06.web"

{
if(vm<TAB_LO){
Warning(2,"clipped vm=%f to %f",vm,TAB_LO);
vm= TAB_LO;
}
else if(vm> TAB_HI-TAB_STEP){
Warning(2,"clipped vm=%f to %f",vm,TAB_HI);
vm= TAB_HI-TAB_STEP;
}
}


/*:77*/
#line 1007 "tp06.web"

/*54:*/
#line 1558 "tp06.web"

{
double RToverF,FoverRT;
double VmFoverRT;
double cKo;
double cGKs;
double cGKr;
double cGNa;
double cGCaL;

if(tpar->Ko_param)cKo= ctype->param[Ko_param_loc];
else cKo= parm(Ko);

if(tpar->GKs_param)cGKs= ctype->param[GKs_param_loc];
else cGKs= parm(GKs);

if(tpar->GKr_param)cGKr= ctype->param[GKr_param_loc];
else cGKr= parm(GKr);

if(tpar->GNa_param)cGNa= ctype->param[GNa_param_loc];
else cGNa= parm(GNa);

if(tpar->GCaL_param)cGCaL= ctype->param[GCaL_param_loc];
else cGCaL= parm(GCaL);

/*66:*/
#line 1906 "tp06.web"

{
Frel= parm(Fmax_rel)*O_rel*(CaSR-CaSS);
}


/*:66*/
#line 1583 "tp06.web"

/*67:*/
#line 1917 "tp06.web"

{
Fup= parm(Fmaxup)/(1.+(Kupsquare/(Cai*Cai)));
}


/*:67*/
#line 1584 "tp06.web"

Fleak= parm(Fmax_leak)*(CaSR-Cai);
Fxfer= parm(Fmax_xfer)*(CaSS-Cai);

/*55:*/
#line 1595 "tp06.web"

{
double R= 8.314472;
double T= parm(temperature)+273.15;

RToverF= 1000.*R*T/Faraday;
FoverRT= Faraday/(R*T*1000.);
VmFoverRT= vm*FoverRT;
}


/*:55*/
#line 1588 "tp06.web"

/*56:*/
#line 1612 "tp06.web"

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
E_Na= RToverF*log(parm(Nao)/Nai);
E_Ka= RToverF*log(cKo/Ki);
E_Ks= RToverF*log((cKo+parm(pKNa)*parm(Nao))/(Ki+parm(pKNa)*Nai));
E_Ca= 0.5*RToverF*log(parm(Cao)/Cai);
Ak1= 0.1/(1.+exp(0.06*(vm-E_Ka-200)));
Bk1= (3.*exp(0.0002*(vm-E_Ka+100))+
exp(0.1*(vm-E_Ka-10)))/(1.+exp(-0.5*(vm-E_Ka)));
rec_IK1= Ak1/(Ak1+Bk1);
rec_iNaK= (1./(1.+0.1245*exp(-0.1*VmFoverRT)+0.0353*exp(-VmFoverRT)));
rec_ipK= 1./(1.+exp((25-vm)/5.98));
}else{
E_Na= os_E_Na;
E_Ka= os_E_Ka;
E_Ks= os_E_Ks;
E_Ca= os_E_Ca;
rec_IK1= os_rec_IK1;
rec_ipK= os_rec_ipK;
rec_iNaK= os_rec_iNaK;
}
if(tpar->eint!=0&&ts%te==0){
os_E_Na= E_Na;
os_E_Ka= E_Ka;
os_E_Ks= E_Ks;
os_E_Ca= E_Ca;
os_rec_IK1= rec_IK1;
os_rec_ipK= rec_ipK;
os_rec_iNaK= rec_iNaK;
}
}


/*:56*/
#line 1589 "tp06.web"

/*57:*/
#line 1661 "tp06.web"

{
double env= exp((parm(n)-1)*VmFoverRT);
double a;

INa= cGNa*pow3(sm)*sh*sj*(vm-E_Na);

a= exp(2.0*(vm-15.0)*FoverRT);
if(fabs(vm-15.0)<=DBL_EPSILON*15.0){
ICaL= cGCaL*sd*sf*f2*f_CaSS*4.0*Faraday*FoverRT*
((0.25*CaSS-parm(Cao))/(2*FoverRT));
}else{
ICaL= cGCaL*sd*sf*f2*f_CaSS*4.0*(vm-15.0)*Faraday*FoverRT*
((0.25*CaSS*a-parm(Cao))/(a-1.0));
}

Ito= parm(Gto)*sr*ss*(vm-E_Ka);
IKr= cGKr*sqrt(cKo/5.4)*sxr1*sxr2*(vm-E_Ka);
IKs= cGKs*sxs*sxs*(vm-E_Ks);
IK1= parm(GK1)*rec_IK1*(vm-E_Ka);
INaCa= parm(knaca)*
(1./(pow3(parm(KmNai))+pow3(parm(Nao))))*
(1./(parm(KmCa)+parm(Cao)))*
(1./(1+parm(ksat)*env))*
(exp(parm(n)*VmFoverRT)*pow3(Nai)*parm(Cao)-env*Nao3*Cai*2.5);
INaK= parm(knak)*(cKo/(cKo+parm(KmK)))*
(Nai/(Nai+parm(KmNa)))*rec_iNaK;
IpCa= parm(GpCa)*Cai/(parm(KpCa)+Cai);
IpK= parm(GpK)*rec_ipK*(vm-E_Ka);
IbNa= parm(GbNa)*(vm-E_Na);
IbCa= parm(GbCa)*(vm-E_Ca);
}


/*:57*/
#line 1590 "tp06.web"

}


/*:54*/
#line 1008 "tp06.web"

Itotal= IKr+IKs+IK1+Ito+INa+IbNa+ICaL+IbCa+
INaK+INaCa+IpCa+IpK;
if(!finite(Itotal)){
printf("IKr   = %f\n",IKr);
printf("IKs   = %f\n",IKs);
printf("IK1   = %f\n",IK1);
printf("Ito   = %f\n",Ito);
printf("INa   = %f\n",INa);
printf("IbNa  = %f\n",IbNa);
printf("ICaL  = %f\n",ICaL);
printf("IbCa  = %f\n",IbCa);
printf("INaK  = %f\n",INaK);
printf("INaCa = %f\n",INaCa);
printf("IpCa  = %f\n",IpCa);
printf("IpK   = %f\n",IpK);
printf("Istim = %f\n",Isd);
Error(1,"tp06_step: Itotal=%f at %s\n",Itotal,node_to_string(elm));
}
mh0= sm*sh;
/*59:*/
#line 1711 "tp06.web"

{
double Xr1_inf,Xr2_inf,Xs_inf,R_inf,S_inf,D_inf;
double F_inf,M_inf,H_inf,J_inf,F2_inf;

double exp_M,exp_H,exp_J,exp_Xr1,exp_Xr2,exp_D,exp_F,exp_F2;
double exp_Xs,exp_S,exp_R;

if(frozen&&vm> -50.0)frozen= 0;
/*61:*/
#line 1793 "tp06.web"

{
double FNai;
double FKi;
double FCai;
double FCaSR;
double FCaSS;

/*64:*/
#line 1859 "tp06.web"

{
FCai= -(IbCa+IpCa-2*INaCa)*0.5*C_over_Vc_F+
parm(Vsr)/parm(Vc)*(Fleak-Fup)+
Fxfer;
FCaSR= Fup-Fleak-Frel;
FCaSS= -ICaL*0.5*C_over_Vss_F
+parm(Vsr)/parm(Vss)*Frel
-parm(Vc)/parm(Vss)*Fxfer;
}


/*:64*/
#line 1801 "tp06.web"


/*63:*/
#line 1844 "tp06.web"

{
double s= 0.0;

if(!frozen)s= IK1+Ito+IKr+(1.0-parm(pKNa))*IKs-2*INaK+IpK;
if(tpar->use_Isd)s+= Isd;
FKi= -s*C_over_Vc_F;
Ki+= dt*FKi;
}


/*:63*/
#line 1803 "tp06.web"

if(!frozen){
/*62:*/
#line 1833 "tp06.web"

{
FNai= -(INa+IbNa+3*INaK+3*INaCa+parm(pKNa)*IKs)*C_over_Vc_F;
Nai+= dt*FNai;
}


/*:62*/
#line 1805 "tp06.web"

/*68:*/
#line 1965 "tp06.web"

{
double b,c;
double total;
double buf;

buf= parm(Bufc)*Cai/(Cai+parm(Kbufc));
total= Cai+buf+dt*FCai;

b= parm(Bufc)-total+parm(Kbufc);
c= parm(Kbufc)*total;
Cai= 0.5*(sqrt(b*b+4*c)-b);
}


/*:68*/
#line 1806 "tp06.web"

/*69:*/
#line 1983 "tp06.web"

{
double b,c;
double total;
double buf;

buf= parm(Bufsr)*CaSR/(CaSR+parm(Kbufsr));
total= CaSR+buf+dt*FCaSR;

b= parm(Bufsr)-total+parm(Kbufsr);
c= parm(Kbufsr)*total;
CaSR= 0.5*(sqrt(b*b+4*c)-b);
}


/*:69*/
#line 1807 "tp06.web"

/*70:*/
#line 1998 "tp06.web"

{
double b,c;
double total;
double buf;

buf= parm(Bufss)*CaSS/(CaSS+parm(Kbufss));
total= CaSS+buf+dt*FCaSS;

b= parm(Bufss)-total+parm(Kbufss);
c= parm(Kbufss)*total;
CaSS= 0.5*(sqrt(b*b+4*c)-b);
}



/*:70*/
#line 1808 "tp06.web"

}
}



/*:61*/
#line 1720 "tp06.web"

/*76:*/
#line 2217 "tp06.web"

{
int row;
double d0,d1;

/*75:*/
#line 2202 "tp06.web"

{
double m;

m= (vm-TAB_LO)*TABI;
row= m;
d0= m-row;
d1= 1.0-d0;
if(row<0||row>=TAB_ROWS-1)Error(1,"this can't happen");
}


/*:75*/
#line 2222 "tp06.web"

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
exp_F2= interp_tab(mtab_exp_F2);

M_inf= interp_tab(mtab_M_inf);
H_inf= interp_tab(mtab_HJ_inf);
J_inf= H_inf;
Xr1_inf= interp_tab(mtab_Xr1_inf);
Xr2_inf= interp_tab(mtab_Xr2_inf);
Xs_inf= interp_tab(mtab_Xs_inf);
R_inf= interp_tab(mtab_R_inf);
S_inf= interp_tab(mtab_S_inf);
D_inf= interp_tab(mtab_D_inf);
F_inf= interp_tab(mtab_F_inf);
F2_inf= interp_tab(mtab_F2_inf);
}


/*:76*/
#line 1721 "tp06.web"

/*65:*/
#line 1889 "tp06.web"

{
if(!frozen){
double power= parm(KkCaSR)/CaSR;
double temp;
double ma= parm(max_sr),mi= parm(min_sr);

k_CaSR= ma-(ma-mi)/(1.0+power*power);
k_1= parm(k_1f)/k_CaSR;
k_2= parm(k_2f)*k_CaSR;
R_rel+= dt*(parm(k_4)*(1.0-R_rel)-k_2*CaSS*R_rel);
temp= k_1*CaSS*CaSS;
O_rel= (temp*R_rel)/(parm(k_3)+temp);
}
}


/*:65*//*78:*/
#line 2278 "tp06.web"

{
double a,b,tau,inf;

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

f2= F2_inf-(F2_inf-f2)*exp_F2;

a= CaSS*CaSS*400.0;
b= 1.0/(1.0+a);
if(tpar->variant_severi_2009){
tau= 80.0*b+1.0;
inf= 0.9/(1.0+exp((CaSS-1.95)/0.15))+0.1;
if(inf> f_CaSS&&vm> -60.){
f_CaSS= f_CaSS;
}else{
f_CaSS= inf-(inf-f_CaSS)*exp(-(dt/tau));
}
}else{
tau= 80.0*b+2.0;
inf= 0.6*b+0.4;
f_CaSS= inf+(f_CaSS-inf)*exp(-(dt/tau));
}
}
}


/*:78*/
#line 1722 "tp06.web"

}


/*:59*/
#line 1028 "tp06.web"

/*79:*/
#line 2338 "tp06.web"

{
double threshold= 0.05;

mh1= sm*sh;
if(mh0<=threshold&&mh1> threshold){
*dtime= simtime+dt*(threshold-mh0)/(mh1-mh0);
}
}



/*:79*/
#line 1029 "tp06.web"

return Itotal;
}


/*:33*//*35:*/
#line 1065 "tp06.web"


void tp06_init(double dt)
{
int celltype;
int row;
double vm;

/*40:*/
#line 1253 "tp06.web"

{
#define X(nr, name, value) statnames[nr] =  strdup(name)
/*42:*/
#line 1273 "tp06.web"

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
X(10,"f2",f2);
X(11,"f_CaSS",f_CaSS);
X(12,"k_CaSR",k_CaSR);
X(13,"k_1",k_1);
X(14,"k_2",k_2);
X(15,"R_Rel",R_rel);
X(16,"O_Rel",O_rel);
X(17,"Cai",Cai);
X(18,"CaSR",CaSR);
X(19,"CaSS",CaSS);
X(20,"Nai",Nai);
X(21,"Ki",Ki);
X(22,"frozen",frozen);
X(23+0,"INa",INa);
X(23+1,"IKr",IKr);
X(23+2,"IKs",IKs);
X(23+3,"IK1",IK1);
X(23+4,"Ito",Ito);
X(23+5,"IbNa",IbNa);
X(23+6,"ICaL",ICaL);
X(23+7,"IbCa",IbCa);
X(23+8,"INaCa",INaCa);
X(23+9,"IpCa",IpCa);
X(23+10,"IpK",IpK);
X(23+11,"INaK",INaK);
X(23+12,"Frel",Frel);
X(23+13,"Fleak",Fleak);
X(23+14,"Fup",Fup);
X(23+15,"Fxfer",Fxfer);
X(23+16,"E_Na",E_Na);
X(23+17,"E_Ca",E_Ca);
X(23+18,"E_Ka",E_Ka);
X(23+19,"E_Ks",E_Ks);
}



/*:42*/
#line 1256 "tp06.web"

#undef X
}


/*:40*/
#line 1073 "tp06.web"

for(celltype= 0;celltype<TP06_NTYPES;celltype++){
/*60:*/
#line 1735 "tp06.web"

{

Qt= pow(Q10,(parm(temperature)-parm(base_temperature))*0.1);

C_over_Vc_F= parm(capacitance)/(parm(Vc)*Faraday);


C_over_Vss_F= parm(capacitance)/(parm(Vss)*Faraday);

Kupsquare= parm(Kup)*parm(Kup);

Nao3= pow3(parm(Nao));
}


/*:60*/
#line 1075 "tp06.web"

for(row= 0;row<TAB_ROWS;row++){
vm= TAB_LO+row*TAB_STEP;
/*71:*/
#line 2023 "tp06.web"

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
}
else{
ah= (0.057*exp(-(vm+80.)/6.8));
bh= (2.7*exp(0.079*vm)+(3.1e5)*exp(0.3485*vm));
}
tau_h= 1.0/((ah+bh)*Qt);
mtab_exp_H[row][celltype]= exp(-dt/tau_h);
tmpa= 1.+exp((vm-parm(Vhalf_h))/7.43);
mtab_HJ_inf[row][celltype]= 1./(tmpa*tmpa);


if(vm>=-40.){
aj= 0.;
bj= (0.6*exp((0.057)*vm)/(1.+exp(-0.1*(vm+32.))));
}
else{
aj= (((-2.5428e4)*exp(0.2444*vm)-(6.948e-6)*
exp(-0.04391*vm))*(vm+37.78)/
(1.+exp(0.311*(vm+79.23))));
bj= (0.02424*exp(-0.01052*vm)/(1.+exp(-0.1378*(vm+40.14))));
}
tau_j= parm(Ftj)/((aj+bj)*Qt);
mtab_exp_J[row][celltype]= exp(-dt/tau_j);
}


/*:71*//*72:*/
#line 2081 "tp06.web"

{
double axr1,axr2,bxr1,bxr2,Axs,Bxs,tau_Xr1,tau_Xr2,tau_Xs;

mtab_Xr1_inf[row][celltype]= 1./(1.+exp((-26.-vm)/7.));
axr1= 450./(1.+exp((-45.-vm)/10.));
bxr1= 6./(1.+exp((vm-(-30.))/11.5));
tau_Xr1= axr1*bxr1/Qt;


mtab_Xr2_inf[row][celltype]= 1./(1.+exp((vm-(-88.+parm(xr2_off)))/24.));
axr2= 3./(1.+exp((-60.-vm)/20.));
bxr2= 1.12/(1.+exp((vm-60.)/20.));
tau_Xr2= axr2*bxr2/Qt;


mtab_Xs_inf[row][celltype]= 1./(1.+exp((-5.0-vm+parm(VshiftXs))/14.0));
Axs= 1400./(sqrt(1.+exp((5.0-vm+parm(VshiftXs))/6.)));
Bxs= 1./(1.+exp((vm-parm(VshiftXs)-35.0)/15.));
tau_Xs= parm(FtXs)*(Axs*Bxs+80.0)/Qt;


mtab_exp_Xr1[row][celltype]= exp(-dt/tau_Xr1);
mtab_exp_Xr2[row][celltype]= exp(-dt/tau_Xr2);
mtab_exp_Xs[row][celltype]= exp(-dt/tau_Xs);
}


/*:72*//*73:*/
#line 2119 "tp06.web"

{
double a,r_inf,s_inf= 0,tau_r,tau_s= 0;


switch(parm(basetype)){
case TP06_EPIC:case TP06_RVEPI:
s_inf= 1./(1.+exp((vm+20)/5.));
tau_s= (85.*exp(-(vm+45.)*(vm+45.)/320.)+5./(1.+exp((vm-20.)/5.))+3.)/Qt;
break;
case TP06_ENDO:
s_inf= 1./(1.+exp((vm+28)/5.));
tau_s= (1000.*exp(-(vm+67)*(vm+67)/1000.)+8.)/Qt;
break;
case TP06_MCEL:case TP06_RVMC:
s_inf= 1./(1.+exp((vm+20)/5.));
tau_s= (85.*exp(-(vm+45.)*(vm+45.)/320.)+5./(1.+exp((vm-20.)/5.))+3.)/Qt;
break;
case TP06_T5:case TP06_T6:
Error(1,"this can't happen");
break;
default:Error(1,"unknown cell type %d",celltype);
}
r_inf= 1./(1.+exp((20-vm)/6.));
a= vm+40.;
tau_r= (9.5*exp(-a*a/1800.)+0.8)/Qt;
mtab_R_inf[row][celltype]= r_inf;
mtab_S_inf[row][celltype]= s_inf;
mtab_exp_R[row][celltype]= exp(-dt/tau_r);
mtab_exp_S[row][celltype]= exp(-dt/tau_s);
}


/*:73*//*74:*/
#line 2153 "tp06.web"

{
double power,alpha,beta,gamma,tau;

mtab_D_inf[row][celltype]= 1./(1.+exp((-8+parm(d_CaL_off)-vm)/7.5));
alpha= 1.4/(1.+exp((-35-vm)/13))+0.25;
beta= 1.4/(1.+exp((vm+5)/5));
gamma= 1.0/(1.+exp((50-vm)/20));
tau= (alpha*beta+gamma)/Qt;
mtab_exp_D[row][celltype]= exp(-dt/tau);


mtab_F_inf[row][celltype]= 1./(1.+exp((vm+20)/7));
power= vm+27;
alpha= 1102.5*exp(-(power*power)/225);
beta= 200/(1+exp((13-vm)/10));
gamma= 180/(1+exp((vm+30)/10))+20;
tau= parm(FtF)*(alpha+beta+gamma);
mtab_exp_F[row][celltype]= exp(-dt/tau);


if(tpar->variant_severi_2009){
mtab_F2_inf[row][celltype]= 0.3/(1.0+exp((vm+35.0)/7.0))+0.7;

}else{
mtab_F2_inf[row][celltype]= 0.67/(1.0+exp((vm+35.0)/7.0))+0.33;
}
power= vm+25.;
alpha= 600.*exp(-(power*power)/170.);
beta= 31./(1.+exp((25.-vm)/10.));
gamma= 16./(1.+exp((vm+30)/10.));
tau= alpha+beta+gamma;
mtab_exp_F2[row][celltype]= exp(-dt/tau);



}



/*:74*/
#line 1078 "tp06.web"

}
}
}



/*:35*//*36:*/
#line 1096 "tp06.web"


void tp06_info(Membrane_cell_info*ct,
Membrane_info*ifo,
void*prm)
{
int i;
tpar= prm;

ct->mcode= 0;
ct->ccode= TP06_ENDO;
ct->param= NULL;

if(tpar->eint)ifo->Nsvar= TP06_NSVAR+TP06_NSVAR_EINT;
else ifo->Nsvar= TP06_NSVAR;
ifo->Nvar= NVARS;
ifo->Ntypes= TP06_NTYPES;

ifo->Tname= malloc(ifo->Ntypes*sizeof(char*));
for(i= 0;i<ifo->Ntypes;i++)ifo->Tname[i]= strdup(tp06_cname[i]);

ifo->Nparam= 0;
ifo->param= calloc(6,sizeof(int));
if(tpar->Ko_param){
Ko_param_loc= ifo->Nparam++;
ifo->param[Ko_param_loc]= tpar->Ko_param;
}
if(tpar->GKs_param){
GKs_param_loc= ifo->Nparam++;
ifo->param[GKs_param_loc]= tpar->GKs_param;
}
if(tpar->GKr_param){
GKr_param_loc= ifo->Nparam++;
ifo->param[GKr_param_loc]= tpar->GKr_param;
}
if(tpar->GNa_param){
GNa_param_loc= ifo->Nparam++;
ifo->param[GNa_param_loc]= tpar->GNa_param;
}
if(tpar->GCaL_param){
GCaL_param_loc= ifo->Nparam++;
ifo->param[GCaL_param_loc]= tpar->GCaL_param;
}

ifo->info= &tp06_info;
ifo->init= &tp06_init;
ifo->infinite= &tp06_infinite;
ifo->step= &tp06_step;
ifo->get_status= &tp06_get_status;
ifo->ccode= &tp06_ccode;

}


/*:36*//*37:*/
#line 1157 "tp06.web"


void tp06_infinite(Membrane_cell_info*ctype,vm_t*Vm,
yyy_t*cell_status,int inittype)
{
int celltype= ctype->ccode;
double ma= parm(max_sr),mi= parm(min_sr);
double power;
int row;
double d0,d1;
double vm;

frozen= tpar->freeze;
if(inittype==MEMBRANE_INIT_FROM_STATE)return;

*Vm= parm(V_init);
vm= parm(V_init);
/*75:*/
#line 2202 "tp06.web"

{
double m;

m= (vm-TAB_LO)*TABI;
row= m;
d0= m-row;
d1= 1.0-d0;
if(row<0||row>=TAB_ROWS-1)Error(1,"this can't happen");
}


/*:75*/
#line 1174 "tp06.web"


Cai= parm(Cai_init);
CaSR= parm(CaSR_init);
CaSS= parm(CaSS_init);
Nai= parm(Nai_init);
Ki= parm(Ki_init);

sm= interp_tab(mtab_M_inf);
sh= interp_tab(mtab_HJ_inf);
sj= interp_tab(mtab_HJ_inf);
sxr1= interp_tab(mtab_Xr1_inf);
sxr2= interp_tab(mtab_Xr2_inf);
sxs= interp_tab(mtab_Xs_inf);
ss= interp_tab(mtab_S_inf);
sr= interp_tab(mtab_R_inf);

sd= interp_tab(mtab_D_inf);
sf= interp_tab(mtab_F_inf);
f2= interp_tab(mtab_F2_inf);
if(tpar->variant_severi_2009){
f_CaSS= 0.9/(1.0+exp((CaSS-1.95)/0.15))+0.1;
}else{
f_CaSS= 0.6/(1.0+(CaSS/0.05)*(CaSS/0.05))+0.4;
}

power= parm(KkCaSR)/CaSR;
k_CaSR= ma-(ma-mi)/(1+power*power);
k_1= parm(k_1f)/k_CaSR;
k_2= parm(k_2f)*k_CaSR;
if(parm(k_4)> 0&&k_2> 0){
R_rel= parm(k_4)/(k_2*CaSS+parm(k_4));
}else{
R_rel= 0.0;
}
O_rel= (k_1*CaSS*CaSS*R_rel)/(parm(k_3)+k_1*CaSS*CaSS);
}


/*:37*//*39:*/
#line 1229 "tp06.web"

void tp06_get_status(float vm,yyy_t*cell_status,Membrane_cell_info*ctype,
float*Stats,int*Nstats,char***Names)
{
/*48:*/
#line 1414 "tp06.web"


int celltype;

double INa;
double ICaL;
double Ito;
double IKs;
double IKr;
double IK1;
double INaCa;
double INaK;
double IbNa,IbCa;
double IpCa;
double IpK;
double Frel= 0.0;
double Fleak= 0.0;
double Fup= 0.0;
double Fxfer= 0.0;


/*:48*//*49:*/
#line 1439 "tp06.web"


double E_Ka,E_Na,E_Ks,E_Ca;

double rec_IK1,rec_ipK,rec_iNaK;





/*:49*/
#line 1233 "tp06.web"

double dt= 0.0,simtime= 1.0;

celltype= ctype->ccode;
/*77:*/
#line 2253 "tp06.web"

{
if(vm<TAB_LO){
Warning(2,"clipped vm=%f to %f",vm,TAB_LO);
vm= TAB_LO;
}
else if(vm> TAB_HI-TAB_STEP){
Warning(2,"clipped vm=%f to %f",vm,TAB_HI);
vm= TAB_HI-TAB_STEP;
}
}


/*:77*/
#line 1237 "tp06.web"

/*54:*/
#line 1558 "tp06.web"

{
double RToverF,FoverRT;
double VmFoverRT;
double cKo;
double cGKs;
double cGKr;
double cGNa;
double cGCaL;

if(tpar->Ko_param)cKo= ctype->param[Ko_param_loc];
else cKo= parm(Ko);

if(tpar->GKs_param)cGKs= ctype->param[GKs_param_loc];
else cGKs= parm(GKs);

if(tpar->GKr_param)cGKr= ctype->param[GKr_param_loc];
else cGKr= parm(GKr);

if(tpar->GNa_param)cGNa= ctype->param[GNa_param_loc];
else cGNa= parm(GNa);

if(tpar->GCaL_param)cGCaL= ctype->param[GCaL_param_loc];
else cGCaL= parm(GCaL);

/*66:*/
#line 1906 "tp06.web"

{
Frel= parm(Fmax_rel)*O_rel*(CaSR-CaSS);
}


/*:66*/
#line 1583 "tp06.web"

/*67:*/
#line 1917 "tp06.web"

{
Fup= parm(Fmaxup)/(1.+(Kupsquare/(Cai*Cai)));
}


/*:67*/
#line 1584 "tp06.web"

Fleak= parm(Fmax_leak)*(CaSR-Cai);
Fxfer= parm(Fmax_xfer)*(CaSS-Cai);

/*55:*/
#line 1595 "tp06.web"

{
double R= 8.314472;
double T= parm(temperature)+273.15;

RToverF= 1000.*R*T/Faraday;
FoverRT= Faraday/(R*T*1000.);
VmFoverRT= vm*FoverRT;
}


/*:55*/
#line 1588 "tp06.web"

/*56:*/
#line 1612 "tp06.web"

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
E_Na= RToverF*log(parm(Nao)/Nai);
E_Ka= RToverF*log(cKo/Ki);
E_Ks= RToverF*log((cKo+parm(pKNa)*parm(Nao))/(Ki+parm(pKNa)*Nai));
E_Ca= 0.5*RToverF*log(parm(Cao)/Cai);
Ak1= 0.1/(1.+exp(0.06*(vm-E_Ka-200)));
Bk1= (3.*exp(0.0002*(vm-E_Ka+100))+
exp(0.1*(vm-E_Ka-10)))/(1.+exp(-0.5*(vm-E_Ka)));
rec_IK1= Ak1/(Ak1+Bk1);
rec_iNaK= (1./(1.+0.1245*exp(-0.1*VmFoverRT)+0.0353*exp(-VmFoverRT)));
rec_ipK= 1./(1.+exp((25-vm)/5.98));
}else{
E_Na= os_E_Na;
E_Ka= os_E_Ka;
E_Ks= os_E_Ks;
E_Ca= os_E_Ca;
rec_IK1= os_rec_IK1;
rec_ipK= os_rec_ipK;
rec_iNaK= os_rec_iNaK;
}
if(tpar->eint!=0&&ts%te==0){
os_E_Na= E_Na;
os_E_Ka= E_Ka;
os_E_Ks= E_Ks;
os_E_Ca= E_Ca;
os_rec_IK1= rec_IK1;
os_rec_ipK= rec_ipK;
os_rec_iNaK= rec_iNaK;
}
}


/*:56*/
#line 1589 "tp06.web"

/*57:*/
#line 1661 "tp06.web"

{
double env= exp((parm(n)-1)*VmFoverRT);
double a;

INa= cGNa*pow3(sm)*sh*sj*(vm-E_Na);

a= exp(2.0*(vm-15.0)*FoverRT);
if(fabs(vm-15.0)<=DBL_EPSILON*15.0){
ICaL= cGCaL*sd*sf*f2*f_CaSS*4.0*Faraday*FoverRT*
((0.25*CaSS-parm(Cao))/(2*FoverRT));
}else{
ICaL= cGCaL*sd*sf*f2*f_CaSS*4.0*(vm-15.0)*Faraday*FoverRT*
((0.25*CaSS*a-parm(Cao))/(a-1.0));
}

Ito= parm(Gto)*sr*ss*(vm-E_Ka);
IKr= cGKr*sqrt(cKo/5.4)*sxr1*sxr2*(vm-E_Ka);
IKs= cGKs*sxs*sxs*(vm-E_Ks);
IK1= parm(GK1)*rec_IK1*(vm-E_Ka);
INaCa= parm(knaca)*
(1./(pow3(parm(KmNai))+pow3(parm(Nao))))*
(1./(parm(KmCa)+parm(Cao)))*
(1./(1+parm(ksat)*env))*
(exp(parm(n)*VmFoverRT)*pow3(Nai)*parm(Cao)-env*Nao3*Cai*2.5);
INaK= parm(knak)*(cKo/(cKo+parm(KmK)))*
(Nai/(Nai+parm(KmNa)))*rec_iNaK;
IpCa= parm(GpCa)*Cai/(parm(KpCa)+Cai);
IpK= parm(GpK)*rec_ipK*(vm-E_Ka);
IbNa= parm(GbNa)*(vm-E_Na);
IbCa= parm(GbCa)*(vm-E_Ca);
}


/*:57*/
#line 1590 "tp06.web"

}


/*:54*/
#line 1238 "tp06.web"

/*41:*/
#line 1265 "tp06.web"

{
#define X(nr, name, value)  Stats[nr] =  value
/*42:*/
#line 1273 "tp06.web"

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
X(10,"f2",f2);
X(11,"f_CaSS",f_CaSS);
X(12,"k_CaSR",k_CaSR);
X(13,"k_1",k_1);
X(14,"k_2",k_2);
X(15,"R_Rel",R_rel);
X(16,"O_Rel",O_rel);
X(17,"Cai",Cai);
X(18,"CaSR",CaSR);
X(19,"CaSS",CaSS);
X(20,"Nai",Nai);
X(21,"Ki",Ki);
X(22,"frozen",frozen);
X(23+0,"INa",INa);
X(23+1,"IKr",IKr);
X(23+2,"IKs",IKs);
X(23+3,"IK1",IK1);
X(23+4,"Ito",Ito);
X(23+5,"IbNa",IbNa);
X(23+6,"ICaL",ICaL);
X(23+7,"IbCa",IbCa);
X(23+8,"INaCa",INaCa);
X(23+9,"IpCa",IpCa);
X(23+10,"IpK",IpK);
X(23+11,"INaK",INaK);
X(23+12,"Frel",Frel);
X(23+13,"Fleak",Fleak);
X(23+14,"Fup",Fup);
X(23+15,"Fxfer",Fxfer);
X(23+16,"E_Na",E_Na);
X(23+17,"E_Ca",E_Ca);
X(23+18,"E_Ka",E_Ka);
X(23+19,"E_Ks",E_Ks);
}



/*:42*/
#line 1268 "tp06.web"

#undef X
}


/*:41*/
#line 1239 "tp06.web"

/*43:*/
#line 1326 "tp06.web"

{
int i;
for(i= 0;i<NVARS;i++){
if(!(finite(Stats[i]))){
Error(1,"tnnp_get_status: %s=%f",statnames[i],Stats[i]);
}
}
}



/*:43*/
#line 1240 "tp06.web"

*Nstats= NVARS;
*Names= statnames;
}


/*:39*//*44:*/
#line 1342 "tp06.web"


int tp06_ccode(char*name)
{
int i;

for(i= 0;i<TP06_NTYPES;i++){
if(strcmp(name,tp06_cname[i])==0)return i;
}
Error(1,"I don't know about cell type \"%s\"",name);
return-1;
}


/*:44*/
