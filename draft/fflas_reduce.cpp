#include <fflas-ffpack/fflas-ffpack.h>

int RowReduce_int32( uint16_t p, int32_t * A, uint32_t m, uint32_t n, uint32_t lda)
{
	size_t * P, * Qt ;
	FFPACK::Modular<int32_t> F(p);
	return FFPACK::LUdivine(F,FFLAS::FflasUnit,FFLAS::FflasNoTrans,m,n,A,lda,P,Qt);
}
