/*1:*/
#line 70 "ms.web"


// ms.web $Rev$ 

/*2:*/
#line 111 "ms.web"


#define MS_DEFAULT 0
#define MS_NTYPES 1


/*:2*/
#line 74 "ms.web"

/*4:*/
#line 134 "ms.web"

#define MS_NSVAR 1


/*:4*/
#line 75 "ms.web"

/*3:*/
#line 117 "ms.web"


#ifdef MS_C
char ms_cname[MS_NTYPES][20]= 
{"default"};
#else
#line 123 "ms.web"
 extern char ms_cname[MS_NTYPES][20];
#endif
#line 125 "ms.web"


/*:3*/
#line 76 "ms.web"


void ms_info(Membrane_cell_info*ct,
Membrane_info*ifo,
void*prm);

void ms_init(double dt);

void ms_infinite(Membrane_cell_info*ctype,vm_t*Vm,
yyy_t*cell_status,int inittype);

double ms_step(double vm,yyy_t*cell_status,
Membrane_cell_info*ctype,double Isd,
double dt,float*dtime,
float simtime,long elm);

void ms_get_status(float vm,yyy_t*cell_status,Membrane_cell_info*ctype,
float*Stats,int*Nstats,char***Names);

int ms_ccode(char*name);

double mitchellschaeffer(double Vmem,
yyy_t*cell_status,
Membrane_cell_info*ctype,
double Isd,
double deltaT,
float*dtime,float simtime,
float*Stats);

void print_stats_ms(float*Stats);


/*:1*/
