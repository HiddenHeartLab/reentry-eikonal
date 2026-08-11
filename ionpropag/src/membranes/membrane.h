/*2:*/
#line 37 "membrane.web"


#ifndef MEMBRANE_H_INCLUDED
#define MEMBRANE_H_INCLUDED 1

#include "membranes/membrane.d"
#include <math.h> 

#ifdef HAVE_CUDA
#include "cuda_runtime.h"
#include <stdlib.h> 
#include <stddef.h> 
#endif
#line 50 "membrane.web"

/*:2*//*4:*/
#line 63 "membrane.web"


typedef struct{
unsigned char mcode;
unsigned char ccode;
float*param;
}Membrane_cell_info;


/*:4*//*5:*/
#line 93 "membrane.web"


typedef double yyy_t;
typedef double vm_t;

#define VM_T_IS_DOUBLE (1)
#define VM_T_IS_FLOAT (0)


/*:5*//*6:*/
#line 110 "membrane.web"


#define MEMBRANE_INIT_FROM_SCRATCH 1
#define MEMBRANE_INIT_FROM_STATE   2


/*:6*//*7:*/
#line 126 "membrane.web"


/*9:*/
#line 203 "membrane.web"


struct Membrane_info;


typedef void(*Membrane_info_fp)(Membrane_cell_info*ct,
struct Membrane_info*ifo,
void*prm);





typedef void(*Membrane_init_fp)(double dt);


typedef void(*Membrane_init_cuda_fp)();

typedef void(*Membrane_infinite_fp)(Membrane_cell_info*ctype,
vm_t*Vm,
yyy_t*cell_status,
int inittype);




typedef double(*Membrane_step_fp)(double vm,yyy_t*cell_status,
Membrane_cell_info*ctype,
double Isd,
double dt,float*dtime,
float simtime,
long elm);



typedef void(*Membrane_step_cuda_fp)(double vm,yyy_t*cell_status,
Membrane_cell_info*ctype,
double Isd,
double dt,float*dtime,
float simtime,
long elm,double*Iion,int*status);

typedef void(*Membrane_get_status_fp)(float vm,
yyy_t*cell_status,
Membrane_cell_info*ctype,
float*stats,
int*Nstats,
char***names);




typedef int(*Membrane_ccode_fp)(char*name);





/*:9*/
#line 128 "membrane.web"


typedef struct Membrane_info{

int Nsvar;
int Nvar;
int Ntypes;
char**Tname;

int Nparam;
int*param;


Membrane_info_fp info;
Membrane_init_fp init;
Membrane_infinite_fp infinite;
Membrane_step_fp step;
Membrane_get_status_fp get_status;
Membrane_ccode_fp ccode;


int can_use_cuda;
Membrane_init_cuda_fp init_cuda;


}Membrane_info;

/*:7*//*8:*/
#line 159 "membrane.web"


typedef struct Membrane_info_cuda{

int Nsvar;
int Ntypes;

int Nparam;


int dtime_inout;


Membrane_step_cuda_fp step;


}Membrane_info_cuda;

#endif
#line 178 "membrane.web"

/*:8*/
