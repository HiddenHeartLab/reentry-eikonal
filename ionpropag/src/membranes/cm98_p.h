/*
    cm98_p.h

    Declares the functions
      cm98_prm()
      PrMsavefile()
    adapted to read/write parameters for program cm98.

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

typedef struct CM98_param CM98_param;
typedef struct PrMBCM98_param PrMBCM98_param;

struct CM98_param {
    BooleaN	 use_Isd;
    BooleaN	 freeze;
    Float	 mh_threshold[6];
    BooleaN	 smooth_steps;
    BooleaN	 smooth_singularities;
    BooleaN	 test_Itotal_finite;
    Float	 Gamma[6];
    Float	 K_mCsqn[6];
    Float	 K_mTrpn[6];
    Float	 K_mCmdn[6];
    Float	 Csqn_max[6];
    Float	 Trpn_max[6];
    Float	 Cmdn_max[6];
    Float	 Ca_up_max[6];
    Float	 K_up[6];
    Float	 k_sat[6];
    Float	 K_mCa[6];
    Float	 K_mNa[6];
    Float	 K_mK_o[6];
    Float	 K_mNa_i[6];
    Float	 k_rel[6];
    Float	 Iup_max[6];
    Float	 I_pCa_max[6];
    Float	 INaCa_max[6];
    Float	 INaK_max[6];
    Float	 g_bNa[6];
    Float	 g_bCa[6];
    Float	 g_CaL[6];
    Float	 g_Ks[6];
    Float	 g_Kr[6];
    Float	 f_Kur[6];
    Float	 g_to[6];
    Float	 g_K1[6];
    Float	 g_Na[6];
    Float	 Ca_o[6];
    Float	 Na_o[6];
    Float	 K_o[6];
    Float	 V_rel[6];
    Float	 V_up[6];
    Float	 V_i[6];
    Float	 V_cell[6];
    Float	 C_m[6];
    Float	 K_Q10[6];
    Float	 T[6];
    Int	 freeze_param;
    Int	 Ko_param;
};
struct PrMBCM98_param {
    BooleaN	 use_Isd;
    BooleaN	 freeze;
    BooleaN	 mh_threshold[6];
    BooleaN	 smooth_steps;
    BooleaN	 smooth_singularities;
    BooleaN	 test_Itotal_finite;
    BooleaN	 Gamma[6];
    BooleaN	 K_mCsqn[6];
    BooleaN	 K_mTrpn[6];
    BooleaN	 K_mCmdn[6];
    BooleaN	 Csqn_max[6];
    BooleaN	 Trpn_max[6];
    BooleaN	 Cmdn_max[6];
    BooleaN	 Ca_up_max[6];
    BooleaN	 K_up[6];
    BooleaN	 k_sat[6];
    BooleaN	 K_mCa[6];
    BooleaN	 K_mNa[6];
    BooleaN	 K_mK_o[6];
    BooleaN	 K_mNa_i[6];
    BooleaN	 k_rel[6];
    BooleaN	 Iup_max[6];
    BooleaN	 I_pCa_max[6];
    BooleaN	 INaCa_max[6];
    BooleaN	 INaK_max[6];
    BooleaN	 g_bNa[6];
    BooleaN	 g_bCa[6];
    BooleaN	 g_CaL[6];
    BooleaN	 g_Ks[6];
    BooleaN	 g_Kr[6];
    BooleaN	 f_Kur[6];
    BooleaN	 g_to[6];
    BooleaN	 g_K1[6];
    BooleaN	 g_Na[6];
    BooleaN	 Ca_o[6];
    BooleaN	 Na_o[6];
    BooleaN	 K_o[6];
    BooleaN	 V_rel[6];
    BooleaN	 V_up[6];
    BooleaN	 V_i[6];
    BooleaN	 V_cell[6];
    BooleaN	 C_m[6];
    BooleaN	 K_Q10[6];
    BooleaN	 T[6];
    BooleaN	 freeze_param;
    BooleaN	 Ko_param;
};


#ifdef PrMGLOBAL
/* Global parameter variables */
#endif


int	PrMsavefile(char *PrMfilename, char *PrMprogname);

int	cm98_prm(
	/* Variables */
	CM98_param *PrMp000, 
	int *PrMpargc, char *PrMargv[]);

