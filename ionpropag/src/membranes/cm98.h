/*1:*/
#line 90 "cm98.web"


// cm98.web $Rev$ 

#ifdef __cplusplus
extern"C"{
#endif
#line 97 "cm98.web"

/*3:*/
#line 157 "cm98.web"


#define CM98_DEFAULT 0
#define CM98_T1 1
#define CM98_T2 2
#define CM98_T3 3
#define CM98_T4 5
#define CM98_T5 6
#define CM98_NTYPES 6


/*:3*/
#line 98 "cm98.web"

/*5:*/
#line 208 "cm98.web"

#define CM98_NSVAR 25


/*:5*/
#line 99 "cm98.web"

/*4:*/
#line 168 "cm98.web"


#ifdef CM98_C
char cm98_cname[CM98_NTYPES][20]= 
{"default","type1","type2","type3","type4","type5"};
#else
#line 174 "cm98.web"
 extern char cm98_cname[CM98_NTYPES][20];
#endif
#line 176 "cm98.web"


/*:4*/
#line 100 "cm98.web"


void cm98_info(Membrane_cell_info*ct,
Membrane_info*ifo,
void*prm);

void cm98_info_cuda(Membrane_info_cuda*ifo);

void cm98_init(double dt);

void cm98_init_cuda();

void cm98_infinite(Membrane_cell_info*ctype,vm_t*Vm,
yyy_t*cell_status,int inittype);

double cm98_step(double vm,yyy_t*cell_status,
Membrane_cell_info*ctype,double Isd,
double dt,float*dtime,
float simtime,long elm);

CUDA_HOSTDEVICE
void cm98_step_cuda(double vm,yyy_t*cell_status,
Membrane_cell_info*ctype,double Isd,
double dt,float*dtime,
float simtime,long elm,double*Iia,int*status);

void cm98_get_status(float vm,yyy_t*cell_status,Membrane_cell_info*ctype,
float*Stats,int*Nstats,char***Names);

int cm98_ccode(char*name);

/*:1*//*2:*/
#line 136 "cm98.web"


CUDA_HOSTDEVICE
double courtemanche(double Vmem,
yyy_t*cell_status,
Membrane_cell_info*ctype,
double Isd,
double deltaT,
float*dtime,float simtime,
float*Stats);

void print_stats(float*Stats);

#ifdef __cplusplus
}
#endif
#line 152 "cm98.web"

/*:2*/
