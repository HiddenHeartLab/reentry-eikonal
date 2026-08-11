/*8:*/
#line 217 "ms.web"

// ms.web $Rev$ 

#define MS_C
/*24:*/
#line 547 "ms.web"

#include <limits.h> 
#include <float.h> 
#include <math.h> 
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>         
#ifdef MATLAB_MEX_FILE
#undef strdup             
#endif
#line 557 "ms.web"
#include <time.h> 

#include "membranes/membrane.h"
#include "membranes/membrane_internal.h"
#include "membranes/ms_p.h"    
#include "membranes/ms.h"



/*:24*/
#line 221 "ms.web"

#define h cell_status[0] \
 \

#define parm(s) (tpar->s[celltype])  \

#define NVARS 1

#line 222 "ms.web"


/*9:*/
#line 241 "ms.web"


static MS_param*tpar;


/*:9*//*19:*/
#line 483 "ms.web"


static char*statnames[NVARS];


/*:19*/
#line 224 "ms.web"

/*10:*/
#line 257 "ms.web"


void ms_info(Membrane_cell_info*ct,Membrane_info*ifo,
void*prm)
{
int i;
tpar= prm;

ct->mcode= 0;
ct->ccode= 1;
ct->param= NULL;

ifo->Nsvar= MS_NSVAR;
ifo->Nvar= NVARS;

ifo->Ntypes= MS_NTYPES;

ifo->Tname= malloc(ifo->Ntypes*sizeof(char*));
for(i= 0;i<ifo->Ntypes;i++)ifo->Tname[i]= strdup(ms_cname[i]);

ifo->Nparam= 0;
ifo->param= calloc(3,sizeof(int));

ifo->info= &ms_info;
ifo->init= &ms_init;
ifo->infinite= &ms_infinite;
ifo->step= &ms_step;
ifo->get_status= &ms_get_status;
ifo->ccode= &ms_ccode;
}


/*:10*/
#line 225 "ms.web"

/*11:*/
#line 292 "ms.web"


void ms_init(double dt)
{
/*20:*/
#line 496 "ms.web"

{
#define X(nr, name, value) statnames[nr] =  strdup(name)
/*22:*/
#line 519 "ms.web"

{
X(0,"h",h);
}


/*:22*/
#line 499 "ms.web"

#undef X
}

/*:20*/
#line 296 "ms.web"

}


/*:11*/
#line 226 "ms.web"

/*12:*/
#line 303 "ms.web"


void ms_infinite(Membrane_cell_info*ctype,vm_t*Vm,
yyy_t*cell_status,int inittype)
{
*Vm= tpar->V_init;

h= tpar->h_init;
}


/*:12*/
#line 227 "ms.web"

/*14:*/
#line 331 "ms.web"


double ms_step(double Vmem,
yyy_t*cell_status,
Membrane_cell_info*ctype,
double Isd,
double deltaT,
float*dtime,
float simtime,
long elm)
{
double Jtotal;
Jtotal= mitchellschaeffer(Vmem,cell_status,ctype,Isd,deltaT,
dtime,simtime,NULL);
if(tpar->test_Itotal_finite){
if(!finite(Jtotal)){
Error(1,"ms_step: Jtotal=%f Vmem=%f at %s",
Jtotal,Vmem,node_to_string(elm));
}
}
return Jtotal;
}


void ms_get_status(float vm,yyy_t*cell_status,Membrane_cell_info*ctype,
float*Stats,int*Nstats,char***Names)
{
double Jtotal;
Jtotal= mitchellschaeffer(vm,cell_status,ctype,0.0,0.0,
NULL,0,Stats);
if(tpar->test_Itotal_finite){
if(!finite(Jtotal)){
print_stats_ms(Stats);
Error(1,"ms_get_status: Jtotal=%f vm=%f",Jtotal,vm);
}
}
*Nstats= NVARS;
*Names= statnames;
}


void print_stats_ms(float*Stats)
{
int i;
for(i= 0;i<NVARS;i++){
fprintf(stderr,"ms stat %02d = %s = %f\n",i,statnames[i],Stats[i]);
}
}


/*:14*/
#line 228 "ms.web"

/*15:*/
#line 392 "ms.web"


double mitchellschaeffer(double Vmem,
yyy_t*cell_status,
Membrane_cell_info*ctype,
double Isd,
double deltaT,
float*dtime,float simtime,
float*Stats)
{
/*17:*/
#line 433 "ms.web"


int celltype;
double Jtotal;


/*:17*/
#line 402 "ms.web"

/*13:*/
#line 317 "ms.web"


double Vmin= tpar->V_min;
double Vmax= tpar->V_max;


/*:13*/
#line 403 "ms.web"

double V= (Vmem-Vmin)/(Vmax-Vmin);

celltype= ctype->ccode;

Jtotal= (h*V*V*(1-V))/(tpar->tau_in)-V/(tpar->tau_out);
Jtotal*= (Vmax-Vmin);



/*16:*/
#line 423 "ms.web"

{
if(V<tpar->v_gate){
h= h+deltaT*((1-h)/tpar->tau_open);
}else if(V> tpar->v_gate){
h= h-deltaT*(h/tpar->tau_close);
}
}


/*:16*/
#line 413 "ms.web"


if(Stats)/*21:*/
#line 507 "ms.web"

{
#define X(nr, name, value)  Stats[nr] =  value
/*22:*/
#line 519 "ms.web"

{
X(0,"h",h);
}


/*:22*/
#line 510 "ms.web"

#undef X
}

/*:21*/
#line 415 "ms.web"

if(dtime)/*18:*/
#line 462 "ms.web"

{
double threshold= 0.5;
double Vmem_threshold= (1-threshold)*Vmin+threshold*Vmax;
double Vmem_new= Vmem+deltaT*(Jtotal-Isd);

if(Vmem_new> Vmem_threshold){
if(Vmem<=Vmem_threshold){
*dtime= simtime+deltaT*(Vmem_threshold-Vmem)/(Vmem_new-Vmem);
}else if(*dtime<0){
*dtime= simtime;
}
}
}



/*:18*/
#line 416 "ms.web"

return-Jtotal;

}

/*:15*/
#line 229 "ms.web"

/*23:*/
#line 529 "ms.web"


int ms_ccode(char*name)
{
int i;

for(i= 0;i<MS_NTYPES;i++){
if(strcmp(name,ms_cname[i])==0)return i;
}
Error(1,"I don't know about cell type \"%s\"",name);
return-1;
}


/*:23*/
#line 230 "ms.web"



/*:8*/
