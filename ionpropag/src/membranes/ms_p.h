/*
    ms_p.h

    Declares the functions
      ms_prm()
      PrMsavefile()
    adapted to read/write parameters for program ms.

    This file was generated automatically by the PRM utility.

    Laboratoire de Modelisation Biomedicale
    Institut de Genie Biomedical
    Ecole Polytechnique / Faculte de Medecine
    Universite de Montreal

*/

#define PARAMETERS \
	&dummy

#define PrMNOPARAM
#include <PrM/Defines.h>
#include <PrM/Typedefs.h>
#include <PrM/Functions.h>

typedef struct MS_param MS_param;
typedef struct PrMBMS_param PrMBMS_param;

struct MS_param {
    BooleaN	 test_Itotal_finite;
    Double	 h_init;
    Double	 V_max;
    Double	 V_min;
    Double	 V_init;
    Double	 v_gate;
    Double	 tau_close;
    Double	 tau_open;
    Double	 tau_out;
    Double	 tau_in;
};
struct PrMBMS_param {
    BooleaN	 test_Itotal_finite;
    BooleaN	 h_init;
    BooleaN	 V_max;
    BooleaN	 V_min;
    BooleaN	 V_init;
    BooleaN	 v_gate;
    BooleaN	 tau_close;
    BooleaN	 tau_open;
    BooleaN	 tau_out;
    BooleaN	 tau_in;
};


#ifdef PrMGLOBAL
/* Global parameter variables */
#endif


int	PrMsavefile(char *PrMfilename, char *PrMprogname);

int	ms_prm(
	/* Variables */
	MS_param *PrMp000, 
	int *PrMpargc, char *PrMargv[]);

