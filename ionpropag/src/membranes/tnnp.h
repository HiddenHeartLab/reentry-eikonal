/*1:*/
#line 87 "tnnp.web"


/*4:*/
#line 223 "tnnp.web"


#define TNNP_EPIC 0     
#define TNNP_MCEL 1     
#define TNNP_ENDO 2     
#define TNNP_RVMC 3     
#define TNNP_RVEPI 4     
#define TNNP_T5    5     
#define TNNP_T6    6     
#define TNNP_NTYPES 7


/*:4*/
#line 89 "tnnp.web"

/*49:*/
#line 1351 "tnnp.web"


#define TNNP_NSVAR 17


/*:49*//*50:*/
#line 1368 "tnnp.web"


#define TNNP_NSVAR_EINT 9



/*:50*/
#line 90 "tnnp.web"

/*5:*/
#line 242 "tnnp.web"


#ifdef TNNP_C
char tnnp_cname[TNNP_NTYPES][20]= 
{"epic","mcel","endo","rvmc","rvepi","type5","type6"};
#else
#line 248 "tnnp.web"
 extern char tnnp_cname[TNNP_NTYPES][20];
#endif
#line 250 "tnnp.web"




/*:5*/
#line 91 "tnnp.web"


void tnnp_info(Membrane_cell_info*ct,
Membrane_info*ifo,
void*prm);

void tnnp_init(double dt);

void tnnp_infinite(Membrane_cell_info*ctype,vm_t*Vm,
yyy_t*cell_status,int inittype);

double tnnp_step(double vm,yyy_t*cell_status,
Membrane_cell_info*ctype,double Isd,
double dt,
float*dtime,
float simtime,
long elm);

void tnnp_get_status(float vm,yyy_t*cell_status,Membrane_cell_info*ctype,
float*Stats,int*Nstats,char***Names);

int tnnp_ccode(char*name);




/*:1*/
