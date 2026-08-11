/*1:*/
#line 91 "tp06.web"


// tp06.web $Rev$ 

/*4:*/
#line 240 "tp06.web"


#define TP06_EPIC 0     
#define TP06_MCEL 1     
#define TP06_ENDO 2     
#define TP06_RVMC 3     
#define TP06_RVEPI 4     
#define TP06_T5    5     
#define TP06_T6    6     
#define TP06_NTYPES 7


/*:4*/
#line 95 "tp06.web"

/*52:*/
#line 1533 "tp06.web"


#define TP06_NSVAR 23


/*:52*//*53:*/
#line 1548 "tp06.web"


#define TP06_NSVAR_EINT 7



/*:53*/
#line 96 "tp06.web"

/*5:*/
#line 259 "tp06.web"


#ifdef TP06_C
char tp06_cname[TP06_NTYPES][20]= 
{"epic","mcel","endo","rvmc","rvepi","type5","type6"};
#else
#line 265 "tp06.web"
 extern char tp06_cname[TP06_NTYPES][20];
#endif
#line 267 "tp06.web"




/*:5*/
#line 97 "tp06.web"


void tp06_info(Membrane_cell_info*ct,
Membrane_info*ifo,
void*prm);

void tp06_init(double dt);

void tp06_infinite(Membrane_cell_info*ctype,vm_t*Vm,
yyy_t*cell_status,int inittype);

double tp06_step(double vm,yyy_t*cell_status,
Membrane_cell_info*ctype,double Isd,
double dt,float*dtime,
float simtime,long elm);

void tp06_get_status(float vm,yyy_t*cell_status,Membrane_cell_info*ctype,
float*Stats,int*Nstats,char***Names);

int tp06_ccode(char*name);




/*:1*/
