/*
    tp06_p.h

    Declares the functions
      tp06_prm()
      PrMsavefile()
    adapted to read/write parameters for program tp06.

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

typedef struct TP06_param TP06_param;
typedef struct PrMBTP06_param PrMBTP06_param;

struct TP06_param {
    BooleaN	 variant_severi_2009;
    Int	 basetype[7];
    BooleaN	 use_Isd;
    Float	 eint;
    Float	 base_temperature[7];
    Float	 temperature[7];
    BooleaN	 freeze;
    Float	 Ki_init[7];
    Float	 Nai_init[7];
    Float	 CaSS_init[7];
    Float	 CaSR_init[7];
    Float	 Cai_init[7];
    Float	 V_init[7];
    Float	 GpK[7];
    Float	 KpCa[7];
    Float	 GpCa[7];
    Float	 n[7];
    Float	 ksat[7];
    Float	 knaca[7];
    Float	 KmCa[7];
    Float	 KmNai[7];
    Float	 d_CaL_off[7];
    Float	 GbCa[7];
    Float	 FtF[7];
    Int	 GCaL_param;
    Float	 GCaL[7];
    Float	 knak[7];
    Float	 KmNa[7];
    Float	 KmK[7];
    Float	 Ftj[7];
    Float	 GbNa[7];
    Float	 Vhalf_h[7];
    Float	 Vhalf_m[7];
    Int	 GNa_param;
    Float	 GNa[7];
    Float	 GK1[7];
    Float	 pKNa[7];
    Float	 xr2_off[7];
    Int	 GKr_param;
    Float	 GKr[7];
    Float	 VshiftXs[7];
    Float	 FtXs[7];
    Int	 GKs_param;
    Float	 GKs[7];
    Float	 Gto[7];
    Float	 min_sr[7];
    Float	 max_sr[7];
    Float	 KkCaSR[7];
    Float	 k_4[7];
    Float	 k_3[7];
    Float	 k_2f[7];
    Float	 k_1f[7];
    Float	 Fmax_xfer[7];
    Float	 Fmax_leak[7];
    Float	 Fmax_rel[7];
    Float	 Kup[7];
    Float	 Fmaxup[7];
    Float	 Kbufss[7];
    Float	 Bufss[7];
    Float	 Kbufsr[7];
    Float	 Bufsr[7];
    Float	 Kbufc[7];
    Float	 Bufc[7];
    Float	 capacitance[7];
    Float	 Vss[7];
    Float	 Vsr[7];
    Float	 Vc[7];
    Float	 Nao[7];
    Float	 Cao[7];
    Int	 Ko_param;
    Float	 Ko[7];
};
struct PrMBTP06_param {
    BooleaN	 variant_severi_2009;
    BooleaN	 basetype[7];
    BooleaN	 use_Isd;
    BooleaN	 eint;
    BooleaN	 base_temperature[7];
    BooleaN	 temperature[7];
    BooleaN	 freeze;
    BooleaN	 Ki_init[7];
    BooleaN	 Nai_init[7];
    BooleaN	 CaSS_init[7];
    BooleaN	 CaSR_init[7];
    BooleaN	 Cai_init[7];
    BooleaN	 V_init[7];
    BooleaN	 GpK[7];
    BooleaN	 KpCa[7];
    BooleaN	 GpCa[7];
    BooleaN	 n[7];
    BooleaN	 ksat[7];
    BooleaN	 knaca[7];
    BooleaN	 KmCa[7];
    BooleaN	 KmNai[7];
    BooleaN	 d_CaL_off[7];
    BooleaN	 GbCa[7];
    BooleaN	 FtF[7];
    BooleaN	 GCaL_param;
    BooleaN	 GCaL[7];
    BooleaN	 knak[7];
    BooleaN	 KmNa[7];
    BooleaN	 KmK[7];
    BooleaN	 Ftj[7];
    BooleaN	 GbNa[7];
    BooleaN	 Vhalf_h[7];
    BooleaN	 Vhalf_m[7];
    BooleaN	 GNa_param;
    BooleaN	 GNa[7];
    BooleaN	 GK1[7];
    BooleaN	 pKNa[7];
    BooleaN	 xr2_off[7];
    BooleaN	 GKr_param;
    BooleaN	 GKr[7];
    BooleaN	 VshiftXs[7];
    BooleaN	 FtXs[7];
    BooleaN	 GKs_param;
    BooleaN	 GKs[7];
    BooleaN	 Gto[7];
    BooleaN	 min_sr[7];
    BooleaN	 max_sr[7];
    BooleaN	 KkCaSR[7];
    BooleaN	 k_4[7];
    BooleaN	 k_3[7];
    BooleaN	 k_2f[7];
    BooleaN	 k_1f[7];
    BooleaN	 Fmax_xfer[7];
    BooleaN	 Fmax_leak[7];
    BooleaN	 Fmax_rel[7];
    BooleaN	 Kup[7];
    BooleaN	 Fmaxup[7];
    BooleaN	 Kbufss[7];
    BooleaN	 Bufss[7];
    BooleaN	 Kbufsr[7];
    BooleaN	 Bufsr[7];
    BooleaN	 Kbufc[7];
    BooleaN	 Bufc[7];
    BooleaN	 capacitance[7];
    BooleaN	 Vss[7];
    BooleaN	 Vsr[7];
    BooleaN	 Vc[7];
    BooleaN	 Nao[7];
    BooleaN	 Cao[7];
    BooleaN	 Ko_param;
    BooleaN	 Ko[7];
};


#ifdef PrMGLOBAL
/* Global parameter variables */
#endif


int	PrMsavefile(char *PrMfilename, char *PrMprogname);

int	tp06_prm(
	/* Variables */
	TP06_param *PrMp000, 
	int *PrMpargc, char *PrMargv[]);

