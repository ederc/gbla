#include <fflas-ffpack/fflas-ffpack.h>

#ifdef __cplusplus
extern "C" {
#endif

uint32_t RowReduce_int32_t ( int32_t p, int32_t * A, uint32_t m, uint32_t n, uint32_t lda) ;
uint32_t RowReduce_double  ( double p, double * A, uint32_t m, uint32_t n, uint32_t lda) ;

#ifdef __cplusplus
}
#endif


uint32_t RowReduce_int32_t( uint32_t p, int32_t * A, uint32_t m, uint32_t n, uint32_t lda)
{
	size_t * P, * Qt ;
	FFPACK::Modular<int32_t> F(p);
	return FFPACK::LUdivine(F,FFLAS::FflasUnit,FFLAS::FflasNoTrans,m,n,A,lda,P,Qt);
}

uint32_t RowReduce_double( double p, double * A, uint32_t m, uint32_t n, uint32_t lda)
{
	size_t * P, * Qt ;
	FFPACK::Modular<double> F(p);
	return FFPACK::LUdivine(F,FFLAS::FflasUnit,FFLAS::FflasNoTrans,m,n,A,lda,P,Qt);
}
