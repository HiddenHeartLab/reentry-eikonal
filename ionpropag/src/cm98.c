/*22:*/
#line 554 "cm98.web"

// cm98.web $Rev$ 

#define CM98_C
/*66:*/
#line 1575 "cm98.web"

#include <limits.h> 
#include <float.h> 
#include <math.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>         
#ifdef MATLAB_MEX_FILE
#undef strdup             
#endif
#line 1585 "cm98.web"
#include <time.h> 

#ifdef __cplusplus
extern"C"{
#endif
#line 1590 "cm98.web"

#include "membranes/membrane.h"
#include "membranes/membrane_internal.h"
#include "membranes/cm98_p.h"    
#include "membranes/cm98.h"

#ifdef __cplusplus
}
#endif
#line 1599 "cm98.web"

/*:66*/
#line 558 "cm98.web"

#define Fn cell_status[0]
#define h cell_status[1]
#define d cell_status[2]
#define x_r cell_status[3]
#define Na_i cell_status[4]
#define K_i cell_status[5]
#define Ca_rel cell_status[6]
#define o_i cell_status[7]
#define u_i cell_status[8]
#define Cmdn_Ca_i cell_status[9]
#define Csqn_Ca_i cell_status[10]
#define v cell_status[11]
#define m cell_status[12]
#define j cell_status[13]
#define f cell_status[14]
#define x_s cell_status[15]
#define Ca_i cell_status[16]
#define Ca_up cell_status[17]
#define o_a cell_status[18]
#define u_a cell_status[19]
#define f_Ca cell_status[20]
#define Trpn_Ca_i cell_status[21]
#define u cell_status[22]
#define w cell_status[23]
#define frozen cell_status[24]
#define parm(s) (DEV_OR_HOST(tpar)->s[celltype])
#define smosi (DEV_OR_HOST(tpar) ->smooth_singularities)
#define linear_interpolate(x,x0,x1,f0,f1) ((((x1)-(x))*f0(x0)+((x)-(x0))*f1(x1))/((x1)-(x0)))
#define SAFE_EPSILON 1e-8
#define float_equals_constant(x,a) (fabs((x)-(a))<=fabs(a)*SAFE_EPSILON)
#define float_equals_zero(x) (fabs(x)<=DBL_EPSILON)
#define alpha_m(v) (0.32*(((v)+47.13)/(1-exp(-0.1*((v)+47.13)))))
#define alpha_h_low(v) (0.135*exp(-((v)+80)/6.8))
#define alpha_h_high(v) 0.0
#define beta_h_low(v) (3.56*exp(0.079*(v))+3.1e05*exp(0.35*(v)))
#define beta_h_high(v) (1/(0.13*(1+exp(-((v)+10.66)/11.1))))
#define alpha_j_low(v) ((-127140*exp(0.2444*(v))-3.474e-05*exp(-0.04391*(v)))* \
(((v)+37.78)/(1+exp(0.311*((v)+79.23)))))
#define alpha_j_high(v) 0.0
#define beta_j_low(v) (0.1212*(exp(-0.01052*(v))/(1+exp(-0.1378*((v)+40.14)))))
#define beta_j_high(v) (0.3*(exp(-2.535e-07*(v))/(1+exp(-0.1*((v)+32)))))
#define alpha_xr(v) (0.0003*(((v)+14.1)/(1.0-exp(-0.2*((v)+14.1)))))
#define beta_xr(v) (7.3898e-05*(((v)-3.3328)/ \
(-1.0+exp(((v)-3.3328)/5.1237))))
#define alpha_xs(v) (4.0e-05*(((v)-19.9)/(1.0-exp(-((v)-19.9)/17))))
#define beta_xs(v) (3.5e-05*(((v)-19.9)/(-1.0+exp(((v)-19.9)/9))))
#define NVARS 45

#line 559 "cm98.web"


/*23:*/
#line 580 "cm98.web"

static CM98_param*tpar;
static int Ko_param_loc;
static int freeze_param_loc;

CUDA_DEVICE CUDA_CONSTANT static CM98_param*DEV_PREFIX(tpar);
CUDA_DEVICE CUDA_CONSTANT static int DEV_PREFIX(Ko_param_loc);
CUDA_DEVICE CUDA_CONSTANT static int DEV_PREFIX(freeze_param_loc);

/*:23*//*61:*/
#line 1464 "cm98.web"


static char*statnames[NVARS];


/*:61*/
#line 561 "cm98.web"

/*39:*/
#line 1096 "cm98.web"


CUDA_HOSTDEVICE
static inline double pow2(double x)
{
return(x*x);
}

CUDA_HOSTDEVICE
static inline double pow3(double x)
{
return(x*x*x);
}


/*:39*//*50:*/
#line 1243 "cm98.web"


CUDA_HOSTDEVICE
static inline double tau_f(double Vmem)
{
double evm= exp(-(Vmem+10)/6.24);
double tau= (1.0-evm)/(0.035*(Vmem+10)*(1+evm));
return tau;
}

/*:50*/
#line 562 "cm98.web"

/*24:*/
#line 600 "cm98.web"


void cm98_info(Membrane_cell_info*ct,Membrane_info*ifo,
void*prm)
{
int i;

tpar= (CM98_param*)prm;

ct->mcode= 0;
ct->ccode= 1;
ct->param= NULL;

ifo->Nsvar= CM98_NSVAR;
ifo->Nvar= NVARS;

ifo->Ntypes= CM98_NTYPES;

#if HAVE_CUDA
ifo->can_use_cuda= 1;
#else
#line 621 "cm98.web"
 ifo->can_use_cuda= 0;
#endif
#line 623 "cm98.web"

ifo->Tname= (char**)malloc(ifo->Ntypes*sizeof(char*));
for(i= 0;i<ifo->Ntypes;i++)ifo->Tname[i]= strdup(cm98_cname[i]);

ifo->Nparam= 0;
ifo->param= (int*)calloc(3,sizeof(int));
if(tpar->Ko_param){
Ko_param_loc= ifo->Nparam;
ifo->Nparam++;
ifo->param[Ko_param_loc]= tpar->Ko_param;
}
if(tpar->freeze_param){
freeze_param_loc= ifo->Nparam;
ifo->Nparam++;
ifo->param[freeze_param_loc]= tpar->freeze_param;
}

ifo->info= &cm98_info;
ifo->init= &cm98_init;
ifo->infinite= &cm98_infinite;
ifo->step= &cm98_step;
ifo->get_status= &cm98_get_status;
ifo->ccode= &cm98_ccode;
ifo->init_cuda= &cm98_init_cuda;

}

CUDA_GLOBAL
void cm98_info_cuda_dev(Membrane_info_cuda*ifo)
{



ifo->Nsvar= CM98_NSVAR;
ifo->Ntypes= CM98_NTYPES;

ifo->Nparam= 0;

ifo->dtime_inout= 2;


ifo->step= &cm98_step_cuda;
}

void cm98_info_cuda(Membrane_info_cuda*ifo)
{
#ifdef HAVE_CUDA


/*25:*/
#line 707 "cm98.web"


{



CM98_param*tmp;
CUDA_CHECK_RETURN(cudaMalloc((void**)&tmp,sizeof(CM98_param)));
CUDA_CHECK_RETURN(cudaMemcpy(tmp,tpar,sizeof(CM98_param),
cudaMemcpyHostToDevice));
CUDA_CHECK_RETURN(cudaMemcpyToSymbol(DEV_PREFIX(tpar),&tmp,sizeof(CM98_param*),
size_t(0),cudaMemcpyHostToDevice));

CUDA_CHECK_RETURN(cudaMemcpyToSymbol(DEV_PREFIX(Ko_param_loc),
&Ko_param_loc,sizeof(int)));
CUDA_CHECK_RETURN(cudaMemcpyToSymbol(DEV_PREFIX(freeze_param_loc),
&freeze_param_loc,sizeof(int)));

}


/*:25*/
#line 672 "cm98.web"

















cm98_info_cuda_dev<<<1,1>>> (ifo);

CUDA_CHECK_RETURN(cudaDeviceSynchronize());
CUDA_CHECK_RETURN(cudaGetLastError());
#endif
}


/*:24*/
#line 563 "cm98.web"

/*26:*/
#line 731 "cm98.web"


void cm98_init(double dt)
{
/*62:*/
#line 1477 "cm98.web"

{
#define X(nr, name, value) statnames[nr] =  strdup(name)
/*64:*/
#line 1500 "cm98.web"

{
X(0,"Fn",Fn);
X(1,"h",h);
X(2,"d",d);
X(3,"x_r",x_r);
X(4,"Na_i",Na_i);
X(5,"K_i",K_i);
X(6,"Ca_rel",Ca_rel);
X(7,"o_i",o_i);
X(8,"u_i",u_i);
X(9,"Cmdn_Ca_i",Cmdn_Ca_i);
X(10,"Csqn_Ca_i",Csqn_Ca_i);
X(11,"v",v);
X(12,"m",m);
X(13,"j",j);
X(14,"f",f);
X(15,"x_s",x_s);
X(16,"Ca_i",Ca_i);
X(17,"Ca_up",Ca_up);
X(18,"o_a",o_a);
X(19,"u_a",u_a);
X(20,"f_Ca",f_Ca);
X(21,"Trpn_Ca_i",Trpn_Ca_i);
X(22,"u",u);
X(23,"w",w);
X(24,"frozen",frozen);
X(25,"INa",INa);
X(26,"IK1",Ik1);
X(27,"Ito",Ito);
X(28,"IKur",IKur);
X(29,"IKr",IKr);
X(30,"IKs",IKs);
X(31,"ICaL",ICaL);
X(32,"I_pCa",I_pCa);
X(33,"INaK",INaK);
X(34,"INaCa",INaCa);
X(35,"IbNa",IbNa);
X(36,"IbCa",IbCa);
X(37,"Iion",Iion);
X(38,"Irel",Irel);
X(39,"Itr",Itr);
X(40,"Iup",Iup);
X(41,"Iup_leak",Iup_leak);
X(42,"E_Ca",E_Ca);
X(43,"E_K",E_Ka);
X(44,"E_Na",E_Na);
}



/*:64*/
#line 1480 "cm98.web"

#undef X
}

/*:62*/
#line 735 "cm98.web"

}


/*:26*/
#line 564 "cm98.web"

/*28:*/
#line 752 "cm98.web"


void cm98_infinite(Membrane_cell_info*ctype,vm_t*Vm,
yyy_t*cell_status,int inittype)
{
frozen= tpar->freeze;
if(inittype==MEMBRANE_INIT_FROM_STATE)return;

*Vm= -81.2;
h= 9.65e-01;
d= 1.37e-04;
x_r= 3.29e-05;
Na_i= 1.12e01;
K_i= 1.39e02;
Ca_rel= 1.49;
o_i= 9.99e-01;
u_i= 9.99e-01;
Cmdn_Ca_i= 2.05e-03;
Csqn_Ca_i= 6.51;
v= 1.00;

m= 2.91e-03;
j= 9.78e-01;
f= 9.99e-01;

x_s= 1.87e-02;
Ca_i= 1.02e-04;
Ca_up= 1.49;
o_a= 3.04e-02;
u_a= 4.96e-03;
f_Ca= 7.75e-01;

Trpn_Ca_i= 1.18e-02;
u= 0.00;
w= 9.99e-01;

Fn= 0.0;

}


/*:28*/
#line 565 "cm98.web"

/*29:*/
#line 807 "cm98.web"


double cm98_step(double Vmem,
yyy_t*cell_status,
Membrane_cell_info*ctype,
double Isd,
double deltaT,
float*dtime,
float simtime,
long elm)
{
double Iion;
Iion= courtemanche(Vmem,cell_status,ctype,Isd,deltaT,
dtime,simtime,NULL);

if(tpar->test_Itotal_finite){
if(!finite(Iion)){
Error(1,"cm98_step: Iion=%f Vmem=%f at %s simtime=%f at %s",
Iion,Vmem,simtime,node_to_string(elm));
}
}

return Iion;
}

CUDA_HOSTDEVICE
void cm98_step_cuda(double Vmem,
yyy_t*cell_status,
Membrane_cell_info*ctype,
double Isd,
double deltaT,
float*dtime,
float simtime,
long elm,
double*Iion,
int*status)
{
*Iion= courtemanche(Vmem,cell_status,ctype,Isd,deltaT,
dtime,simtime,NULL);
}

void cm98_get_status(float vm,yyy_t*cell_status,Membrane_cell_info*ctype,
float*Stats,int*Nstats,char***Names)
{
double Iion;
Iion= courtemanche(vm,cell_status,ctype,0.0,0.0,
NULL,0,Stats);
if(tpar->test_Itotal_finite){
if(!finite(Iion)){
print_stats(Stats);
Error(1,"cm98_get_status: Iion=%f vm=%f",Iion,vm);
}
}
*Nstats= NVARS;
*Names= statnames;
}

void print_stats(float*Stats)
{
int i;
for(i= 0;i<NVARS;i++){
fprintf(stderr,"cm98 stat %02d = %s = %f\n",i,statnames[i],Stats[i]);
}
}



/*:29*/
#line 566 "cm98.web"

/*30:*/
#line 894 "cm98.web"


CUDA_HOSTDEVICE
double courtemanche(double Vmem,
yyy_t*cell_status,
Membrane_cell_info*ctype,
double Isd,
double deltaT,
float*dtime,float simtime,
float*Stats)
{

/*6:*/
#line 215 "cm98.web"


double INa= NAN;
double Ik1= NAN;
double Ito= NAN;
double IKur= NAN;
double IKr= NAN;
double IKs= NAN;
double ICaL= NAN;
double I_pCa= NAN;
double INaK= NAN;
double INaCa= NAN;
double IbNa= NAN;
double IbCa= NAN;
double Iion= NAN;
double Irel= NAN;
double Itr= NAN;
double Iup= NAN;
double Iup_leak= NAN;


/*:6*//*7:*/
#line 239 "cm98.web"


double E_Na= NAN;
double E_Ka= NAN;
double E_Ca= NAN;


/*:7*//*8:*/
#line 249 "cm98.web"


double alpha= NAN;
double beta= NAN;
double tau= NAN;
double inf= NAN;
double VFoverRT= NAN;


/*:8*//*59:*/
#line 1412 "cm98.web"


int celltype= -1;
double cKo= NAN;
double frozen_until= NAN;
double mh0= NAN;


/*:59*/
#line 906 "cm98.web"

/*9:*/
#line 261 "cm98.web"


double R= 8.3143;
double F= 96.4867;
double z_Na= 1.0;
double z_K= 1.0;
double z_Ca= 2.0;


/*:9*/
#line 907 "cm98.web"


celltype= ctype->ccode;
VFoverRT= Vmem*F/(R*parm(T));
if(DEV_OR_HOST(tpar)->Ko_param)
cKo= ctype->param[DEV_OR_HOST(Ko_param_loc)];
else
cKo= parm(K_o);
if(DEV_OR_HOST(tpar)->freeze_param)
frozen_until= ctype->param[DEV_OR_HOST(freeze_param_loc)];
else
frozen_until= 0;

if((frozen_until> 0&&simtime<frozen_until)||
(frozen_until<0&&simtime> -frozen_until))
deltaT= 0.0;

/*31:*/
#line 951 "cm98.web"

{
double RToverF= R*parm(T)/F;

E_Na= RToverF/z_Na*log(parm(Na_o)/Na_i);
E_Ka= RToverF/z_K*log(cKo/K_i);
E_Ca= RToverF/z_Ca*log(parm(Ca_o)/Ca_i);
}


/*:31*/
#line 924 "cm98.web"

mh0= m*h;
/*35:*/
#line 1022 "cm98.web"

{
if(!smosi&&float_equals_constant(Vmem,-47.13)){
alpha= 3.2;
}else if(smosi&&Vmem> -47.2&&Vmem<-47.1){
alpha= linear_interpolate(Vmem,-47.2,-47.1,alpha_m,alpha_m);
}else{
alpha= alpha_m(Vmem);
}
beta= 0.08*exp(-Vmem/11);
tau= 1/(alpha+beta);
inf= alpha*tau;
m= inf+(m-inf)*exp(-(deltaT/tau));
}


/*:35*/
#line 926 "cm98.web"

if(!frozen)/*36:*/
#line 1047 "cm98.web"

{
if(Vmem<-40){
alpha= alpha_h_low(Vmem);
beta= beta_h_low(Vmem);
}else if(Vmem<-39&&DEV_OR_HOST(tpar)->smooth_steps){
alpha= linear_interpolate(Vmem,-40,-39,alpha_h_low,alpha_h_high);
beta= linear_interpolate(Vmem,-40,-39,beta_h_low,beta_h_high);
}else{
alpha= alpha_h_high(Vmem);
beta= beta_h_high(Vmem);
}
tau= 1/(alpha+beta);
inf= alpha*tau;
h= inf+(h-inf)*exp(-(deltaT/tau));
}


/*:36*//*37:*/
#line 1073 "cm98.web"

{
if(Vmem<-40){
alpha= alpha_j_low(Vmem);
beta= beta_j_low(Vmem);
}else if(Vmem<-39&&DEV_OR_HOST(tpar)->smooth_steps){
alpha= linear_interpolate(Vmem,-40,-39,alpha_j_low,alpha_j_high);
beta= linear_interpolate(Vmem,-40,-39,beta_j_low,beta_j_high);
}else{
alpha= alpha_j_high(Vmem);
beta= beta_j_high(Vmem);
}
tau= 1/(alpha+beta);
inf= alpha*tau;
j= inf+(j-inf)*exp(-(deltaT/tau));
}

/*:37*//*42:*/
#line 1121 "cm98.web"

{
alpha= 0.65/(exp(-(Vmem+10)/8.5)+exp(-(Vmem-30)/59.0));
beta= 0.65/(2.5+exp((Vmem+82)/17.0));
tau= 1/((alpha+beta)*parm(K_Q10));
inf= 1/(1+exp(-(Vmem+20.47)/17.54));
o_a= inf+(o_a-inf)*exp(-(deltaT/tau));


alpha= 1/(18.53+exp((Vmem+113.7)/10.95));
beta= 1/(35.56+exp(-(Vmem+1.26)/7.44));
tau= 1/((alpha+beta)*parm(K_Q10));
inf= 1/(1+exp((Vmem+43.1)/5.3));
o_i= inf+(o_i-inf)*exp(-(deltaT/tau));
}

/*:42*//*44:*/
#line 1146 "cm98.web"

{
alpha= 0.65/(exp(-(Vmem+10)/8.5)+exp(-(Vmem-30)/59.0));
beta= 0.65/(2.5+exp((Vmem+82)/17.0));
tau= 1/((alpha+beta)*parm(K_Q10));
inf= 1/(1+exp(-(Vmem+30.3)/9.6));
u_a= inf+(u_a-inf)*exp(-(deltaT/tau));


alpha= 1/(21+exp(-(Vmem-185)/28));
beta= exp((Vmem-158)/16);
tau= 1/((alpha+beta)*parm(K_Q10));
inf= 1/(1+exp((Vmem-99.45)/27.48));
u_i= inf+(u_i-inf)*exp(-(deltaT/tau));
}


/*:44*//*46:*/
#line 1180 "cm98.web"

{
if(!smosi&&float_equals_constant(Vmem,-14.1)){
alpha= 0.0015;
}else if(smosi&&Vmem> -14.2&&Vmem<-14.0){
alpha= linear_interpolate(Vmem,-14.2,-14.0,alpha_xr,alpha_xr);
}else{
alpha= alpha_xr(Vmem);
}

if(!smosi&&float_equals_constant(Vmem,3.3328)){
beta= 0.00037836118;
}else if(smosi&&Vmem> 3.3&&Vmem<3.4){
beta= linear_interpolate(Vmem,3.3,3.4,beta_xr,beta_xr);
}else{
beta= beta_xr(Vmem);
}

tau= 1/(alpha+beta);
inf= 1/(1+exp(-(Vmem+14.1)/6.5));
x_r= inf+(x_r-inf)*exp(-(deltaT/tau));
}


/*:46*//*48:*/
#line 1216 "cm98.web"

{
if(!smosi&&float_equals_constant(Vmem,19.9)){
alpha= 0.00068;
beta= 0.000315;
}else if(smosi&&Vmem> 19.8&&Vmem<20.0){
alpha= linear_interpolate(Vmem,19.8,20.0,alpha_xs,alpha_xs);
beta= linear_interpolate(Vmem,19.8,20.0,beta_xs,beta_xs);
}else{
alpha= alpha_xs(Vmem);
beta= beta_xs(Vmem);
}
tau= 0.5/(alpha+beta);
inf= 1/sqrt(1+exp(-(Vmem-19.9)/12.7));
x_s= inf+(x_s-inf)*exp(-(deltaT/tau));
}


/*:48*//*51:*/
#line 1253 "cm98.web"

{
if(Ca_i<=0){
Error(1,"courtemanche: Ca_i=%f Vmem=%f",Ca_i,Vmem);
}

if(!smosi&&float_equals_constant(Vmem,-10.0)){

tau= 1.0/(2.0*6.24*0.035);
}else if(smosi&&Vmem> -10.1&&Vmem<-9.9){
tau= linear_interpolate(Vmem,-10.1,-9.9,tau_f,tau_f);
}else{
tau= tau_f(Vmem);
}
inf= 1.0/(1.0+exp(-(Vmem+10)/8));
d= inf+(d-inf)*exp(-(deltaT/tau));


tau= 9.0/(0.0197*exp(-pow2(0.0337*(Vmem+10)))+0.02);
inf= 1.0/(1+exp((Vmem+28)/6.9));
f= inf+(f-inf)*exp(-(deltaT/tau));


tau= 2.0;
inf= 1.0/(1+(Ca_i/0.00035));
f_Ca= inf+(f_Ca-inf)*exp(-(deltaT/tau));

}


/*:51*//*55:*/
#line 1333 "cm98.web"

{
tau= 8.0;
inf= 1/(1+exp(-(Fn-3.4175e-13)/13.67e-16));
u= inf+(u-inf)*exp(-(deltaT/tau));


tau= 1.91+2.09/(1+exp(-(Fn-3.4175e-13)/13.67e-16));
inf= 1-1/(1+exp(-(Fn-6.835e-14)/13.67e-16));
v= inf+(v-inf)*exp(-(deltaT/tau));



tau= 6.0*(1-exp(-(Vmem-7.9)/5))/
((1+0.3*exp(-(Vmem-7.9)/5))*(Vmem-7.9));
inf= 1-1/(1+exp(-(Vmem-40)/17));
w= inf+(w-inf)*exp(-(deltaT/tau));

}


/*:55*/
#line 927 "cm98.web"

/*38:*/
#line 1090 "cm98.web"

{
INa= parm(g_Na)*pow3(m)*h*j*(Vmem-E_Na);
}


/*:38*//*41:*/
#line 1112 "cm98.web"

{
Ik1= parm(g_K1)*(Vmem-E_Ka)/(1+exp(0.07*(Vmem+80)));
}


/*:41*//*43:*/
#line 1137 "cm98.web"

{
Ito= parm(g_to)*pow3(o_a)*o_i*(Vmem-E_Ka);
}


/*:43*//*45:*/
#line 1163 "cm98.web"

{
double g;

g= 0.005+(0.05/(1+exp(-(Vmem-15)/13)));
IKur= parm(f_Kur)*g*pow3(u_a)*u_i*(Vmem-E_Ka);
}


/*:45*//*47:*/
#line 1204 "cm98.web"

{
IKr= (parm(g_Kr)*x_r*(Vmem-E_Ka))/(1.0+exp((Vmem+15)/22.4));
}


/*:47*//*49:*/
#line 1234 "cm98.web"

{
IKs= parm(g_Ks)*pow2(x_s)*(Vmem-E_Ka);
}


/*:49*//*52:*/
#line 1283 "cm98.web"

{
ICaL= parm(g_CaL)*d*f*f_Ca*(Vmem-65);
I_pCa= parm(I_pCa_max)*(Ca_i/(0.0005+Ca_i));
}


/*:52*//*53:*/
#line 1294 "cm98.web"

{
double f_NaK;
double sigma= (exp(parm(Na_o)/63.7)-1.0)/7.0;

if(Na_i<=0)Error(1,"courtemanche: Na_i=%f Vmem=%f",Na_i,Vmem);

f_NaK= 1.0/(1+0.1245*exp(-0.1*VFoverRT)+
0.0365*sigma*exp(-VFoverRT));

INaK= parm(INaK_max)*f_NaK*(1/(1+pow((parm(K_mNa_i)/Na_i),1.5)))*
(cKo/(cKo+parm(K_mK_o)));
}


/*:53*//*54:*/
#line 1313 "cm98.web"

{
INaCa= (parm(INaCa_max)*
(exp(parm(Gamma)*VFoverRT)*pow3(Na_i)*parm(Ca_o)-
exp((parm(Gamma)-1.0)*VFoverRT)*pow3(parm(Na_o))*Ca_i
)
)/
((pow3(parm(K_mNa))+pow3(parm(Na_o)))*(parm(K_mCa)+parm(Ca_o))*
(1+parm(k_sat)*exp((parm(Gamma)-1.0)*VFoverRT))
);

IbNa= parm(g_bNa)*(Vmem-E_Na);

IbCa= parm(g_bCa)*(Vmem-E_Ca);
}


/*:54*//*56:*/
#line 1355 "cm98.web"

{
Irel= parm(k_rel)*u*u*v*w*(Ca_rel-Ca_i);
}


/*:56*/
#line 928 "cm98.web"

if(!frozen)/*57:*/
#line 1365 "cm98.web"

{
double sumK;

Itr= (Ca_up-Ca_rel)/180;

Iup= parm(Iup_max)/(1+(parm(K_up)/Ca_i));

Iup_leak= (Ca_up/parm(Ca_up_max))*parm(Iup_max);

Fn= 1e-12*parm(V_rel)*Irel-
(5e-13/F)*parm(C_m)*(0.5*ICaL-0.2*INaCa);


Na_i+= deltaT*parm(C_m)*((-3*INaK-3*INaCa-IbNa-INa)/(F*parm(V_i)));

sumK= -2*INaK+Ik1+Ito+IKur+IKr+IKs;
if(DEV_OR_HOST(tpar)->use_Isd)sumK+= Isd;

K_i-= deltaT*parm(C_m)*sumK/(F*parm(V_i));
}

/*:57*//*58:*/
#line 1390 "cm98.web"

{
double b1,b2;

b1= parm(C_m)*((2*INaCa-I_pCa-ICaL-IbCa)/(2*F*parm(V_i)))+
((parm(V_up)*(Iup_leak-Iup)+Irel*parm(V_rel))/parm(V_i));
b2= 1+((parm(Trpn_max)*parm(K_mTrpn))/pow2(Ca_i+parm(K_mTrpn)))+
((parm(Cmdn_max)*parm(K_mCmdn))/pow2(Ca_i+parm(K_mCmdn)));
Ca_i+= deltaT*(b1/b2);

Ca_up+= deltaT*(Iup-Iup_leak-
Itr*(parm(V_rel)/parm(V_up)));

Ca_rel+= deltaT*(Itr-Irel)/
(1+((parm(Csqn_max)*parm(K_mCsqn))/pow2(Ca_rel+parm(K_mCsqn))));


Cmdn_Ca_i= parm(Cmdn_max)*(Ca_i/(Ca_i+parm(K_mCmdn)));
Trpn_Ca_i= parm(Trpn_max)*(Ca_i/(Ca_i+parm(K_mTrpn)));
Csqn_Ca_i= parm(Csqn_max)*(Ca_i/(Ca_i+parm(K_mCsqn)));
}

/*:58*/
#line 929 "cm98.web"

if(deltaT==0.0)
Iion= 0.0;
else
Iion= INa+Ik1+Ito+IKur+IKr+IKs+ICaL+
I_pCa+INaK+INaCa+IbNa+IbCa;

#if 0
if(!finite(Iion)){
Error(1,"courtemanche: INa=%f Ik1=%f Ito=%f IKur=%f IKr=%f IKs=%f ICaL=%f"
"I_pCa=%f INaK=%f INaCa=%f IbNa=%f IbCa=%f",
INa,Ik1,Ito,IKur,IKr,IKs,ICaL,
I_pCa,INaK,INaCa,IbNa,IbCa);
}
#endif
#line 944 "cm98.web"

if(Stats)/*63:*/
#line 1488 "cm98.web"

{
#define X(nr, name, value)  Stats[nr] =  value
/*64:*/
#line 1500 "cm98.web"

{
X(0,"Fn",Fn);
X(1,"h",h);
X(2,"d",d);
X(3,"x_r",x_r);
X(4,"Na_i",Na_i);
X(5,"K_i",K_i);
X(6,"Ca_rel",Ca_rel);
X(7,"o_i",o_i);
X(8,"u_i",u_i);
X(9,"Cmdn_Ca_i",Cmdn_Ca_i);
X(10,"Csqn_Ca_i",Csqn_Ca_i);
X(11,"v",v);
X(12,"m",m);
X(13,"j",j);
X(14,"f",f);
X(15,"x_s",x_s);
X(16,"Ca_i",Ca_i);
X(17,"Ca_up",Ca_up);
X(18,"o_a",o_a);
X(19,"u_a",u_a);
X(20,"f_Ca",f_Ca);
X(21,"Trpn_Ca_i",Trpn_Ca_i);
X(22,"u",u);
X(23,"w",w);
X(24,"frozen",frozen);
X(25,"INa",INa);
X(26,"IK1",Ik1);
X(27,"Ito",Ito);
X(28,"IKur",IKur);
X(29,"IKr",IKr);
X(30,"IKs",IKs);
X(31,"ICaL",ICaL);
X(32,"I_pCa",I_pCa);
X(33,"INaK",INaK);
X(34,"INaCa",INaCa);
X(35,"IbNa",IbNa);
X(36,"IbCa",IbCa);
X(37,"Iion",Iion);
X(38,"Irel",Irel);
X(39,"Itr",Itr);
X(40,"Iup",Iup);
X(41,"Iup_leak",Iup_leak);
X(42,"E_Ca",E_Ca);
X(43,"E_K",E_Ka);
X(44,"E_Na",E_Na);
}



/*:64*/
#line 1491 "cm98.web"

#undef X
}

/*:63*/
#line 945 "cm98.web"

if(dtime)/*60:*/
#line 1448 "cm98.web"

{
double threshold= parm(mh_threshold);

if(m*h> threshold){
if(frozen)frozen= 0;
if(mh0<=threshold)*dtime= simtime+deltaT*(threshold-mh0)/(m*h-mh0);
}
}



/*:60*/
#line 946 "cm98.web"

return(Iion);
}


/*:30*/
#line 567 "cm98.web"

/*65:*/
#line 1555 "cm98.web"


int cm98_ccode(char*name)
{
int i;

for(i= 0;i<CM98_NTYPES;i++){
if(strcmp(name,cm98_cname[i])==0)return i;
}
Error(1,"I don't know about cell type \"%s\"",name);
return-1;
}




/*:65*/
#line 568 "cm98.web"

/*27:*/
#line 741 "cm98.web"


void cm98_init_cuda()
{

}

/*:27*/
#line 569 "cm98.web"



/*:22*/
