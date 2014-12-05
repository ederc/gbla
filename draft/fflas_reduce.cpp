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
	size_t * P  = FFLAS::fflas_new<size_t>(max(m,n)) ;
	size_t * Qt = FFLAS::fflas_new<size_t>(max(m,n)) ;

	FFPACK::Modular<int32_t> F(p);
	uint32_t r = FFPACK::LUdivine(F,FFLAS::FflasUnit,FFLAS::FflasNoTrans,m,n,A,lda,P,Qt);

	FFLAS::fflas_delete(P);
	FFLAS::fflas_delete(Qt);
	return r ;
}

uint32_t RowReduce_double( double p, double * A, uint32_t m, uint32_t n, uint32_t lda)
{
	size_t * P  = FFLAS::fflas_new<size_t>(max(m,n)) ;
	size_t * Qt = FFLAS::fflas_new<size_t>(max(m,n)) ;

	FFPACK::Modular<double> F(p);
	uint32_t r = FFPACK::LUdivine(F,FFLAS::FflasUnit,FFLAS::FflasNoTrans,m,n,A,lda,P,Qt);

	FFLAS::fflas_delete(P);
	FFLAS::fflas_delete(Qt);
	return r ;

}
