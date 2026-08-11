/*10:*/
#line 266 "membrane.web"


/*11:*/
#line 299 "membrane.web"

#if defined(HAVE_CUDA) && defined(__CUDACC__)
#define CUDA_DEVICE __device__
#define CUDA_HOSTDEVICE __device__ __host__
#define CUDA_GLOBAL __global__
#define CUDA_CONSTANT __constant__
#define DEV_PREFIX(ss) (ddCUDA ## ss)
#else
#line 307 "membrane.web"
#define CUDA_DEVICE
#define CUDA_HOSTDEVICE
#define CUDA_GLOBAL
#define CUDA_CONSTANT
#define DEV_PREFIX(ss) (ss)
#endif
#line 313 "membrane.web"

#define CUDA_CHECK_RETURN(value) {         \
  cudaError_t _m_cudaStat =  value;         \
  if (_m_cudaStat != cudaSuccess) {        \
    Error(EXIT_FAILURE,"CUDA Error = %s", \
          cudaGetErrorString(_m_cudaStat)); \
  } }

#ifdef __CUDA_ARCH__
#define DEV_OR_HOST(ss) DEV_PREFIX(ss)
#else
#line 324 "membrane.web"
#define DEV_OR_HOST(ss) (ss)
#endif
#line 326 "membrane.web"

/*:11*/
#line 268 "membrane.web"


#if defined(HAVE_CUDA) && defined(__CUDA_ARCH__)
#define Error err_locate_cuda( __LINE__, __FILE__,  __DATE__,\
                            __TIME__, __func__ ), err_exit_cuda

extern CUDA_DEVICE void err_exit_cuda(int code,const char*fmt,...);
extern CUDA_DEVICE void err_locate_cuda(int l,const char*f,const char*d,const char*t,const char*u);


#else
#line 279 "membrane.web"
#define Error err_locate( __LINE__, __FILE__,  __DATE__,\
                            __TIME__, __func__ ), err_exit
#endif
#line 282 "membrane.web"

extern void err_exit(int code,const char*fmt,...);
extern void err_locate(int l,const char*f,const char*d,const char*t,const char*u);
extern void Warning(int code,const char*fmt,...);

char*node_to_string(long n);


/*:10*/
