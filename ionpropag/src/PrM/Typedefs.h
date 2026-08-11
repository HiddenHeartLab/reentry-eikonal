#ifndef PrMTYPES

#define PrMTYPES

typedef unsigned char		Byte;
#if SYSTEME==SUNOS
typedef char			Char;
#else
typedef signed char		Char;
#endif
typedef struct S_Complex	S_Complex;
typedef double			Double;
typedef struct D_Complex	D_Complex;
typedef float			Float;
typedef int			Int;
typedef long			Long;
typedef short			Short;
typedef int			BooleaN;
typedef int			Flag;
typedef char			*RDir;
typedef char			*RFile;
typedef char			*RWDir;
typedef char			*RWFile;
#ifndef _XtIntrinsic_h
typedef char			*String;
#endif
typedef char			*WDir;
typedef char			*WFile;

typedef char			**Text;

struct S_Complex {
    Float	real, imag;
};

struct D_Complex {
    Double	real, imag;
};

#endif
