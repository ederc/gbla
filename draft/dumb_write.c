#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

/*
 * a b 0 0 0 0 c
 * a 0 0 0 b 0 c
 * 0 c 0 0 d 0 0
 * 0 0 c 0 0 d 0
 *
 * p_abc = 2 3 4
 * p_cd = 9 8
 *
 * becomes
 * a b 0 0 0 0 c
 * 0 c 0 0 d 0 0
 * 0 0 c 0 0 d 0
 * -------------
 * a 0 0 0 b 0 c
 *
 * becomes
 *
 * a 0 0 | 0 0 b c
 * 0 d 0 | 0 d c 0
 * 0 0 c | 0 0 0 0
 * ---------------
 * a 0 0 | 0 b 0 c
 *
 *
 *
 */


int main()
{
	FILE * toto =fopen("test.gb","wb");
	assert(toto);

	const uint32_t mask = (1U<<31);

	const uint32_t un = 1 ;
	const uint32_t row = 4 ;
	const uint32_t col = 7 ;
	const uint64_t nnz = 10 ;
	const uint32_t mod = 7 ;

	const uint64_t start_zo[] = { 0 , 3, 6 , 8 ,10 };
	/* uint32_t colid_zo[nnz] = { 0, 1, 7, 0, 4, 6, 1, 4, 2, 5 }; */
	const uint32_t colid_zo[] = { 0 , 2 , 7 | mask, 0 | mask, 4 | mask, 6 | mask, 1 |mask, 4 | mask, 2 | mask, 5| mask};
	const uint64_t colid_size = 10 ;
	const uint32_t map_zo_pol[] = { 0, 0, 1, 1 };
	const uint32_t pol_nb = 2 ;
	const uint32_t pol_start [] = { 0 , 3 , 5 } ;
	const uint32_t pol_nnz = pol_start[pol_nb] ;
	const uint32_t pol_data [] = { 2, 3, 4, 9, 8 };

	fwrite(&un,sizeof(uint32_t),1,toto);
	fwrite(&row,sizeof(uint32_t),1,toto);
	fwrite(&col,sizeof(uint32_t),1,toto);
	fwrite(&mod,sizeof(uint32_t),1,toto);
	fwrite(&nnz,sizeof(uint64_t),1,toto);
	fwrite(start_zo,sizeof(uint64_t),row+1,toto);
	fwrite(map_zo_pol,sizeof(uint32_t),row,toto);
	fwrite(&colid_size,sizeof(uint64_t),1,toto);
	fwrite(colid_zo,sizeof(uint32_t),colid_size,toto);
	fwrite(&pol_nb,sizeof(uint32_t),1,toto);
	fwrite(pol_start,sizeof(uint32_t),pol_nb+1,toto);
	fwrite(pol_data,sizeof(uint32_t),pol_nnz,toto);


	fclose(toto);
	return 0;
}

