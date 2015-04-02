#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

typedef uint32_t header_t;
typedef uint32_t dimen_t;
typedef uint64_t index_t;
typedef uint16_t elt_t;

int main (int argc, char *argv[])
{
  if (argc != 2)
    return 1;
  else {
    char *fni = argv[1];
    char *out = ".gbm";
    char fno[100];
    strcpy(fno, fni);
    strcat(fno, out);

    char buffer[10000]; 
    int position[20000];
    int position2[20000];
    int ctr, i, last;
    char *pos;
    char *pos2;
    char tmp[1000];
    FILE* fhi = fopen(fni,"r"); 
    FILE* fho = fopen(fno,"wb"); 
    /*
     * only needed for new format
     *
    header_t header = 0xFFFF0000;
    fwrite(&header, sizeof(header_t), 1, fho);
    */
    fgets(buffer, 10000, fhi);
    printf("%s\n",buffer);
    dimen_t nrows = (dimen_t) atoi(buffer);
    fgets(buffer, 10000, fhi);
    printf("%s\n",buffer);
    dimen_t ncols = (dimen_t) atoi(buffer);
#if 0
    pos = strchr(buffer,',');
    ctr = 0;
    while (pos != NULL) {
      position[ctr] = pos-buffer+1;
      ctr++;
      pos=strchr(pos+1,',');
    }

    // get nrows
    strncpy(tmp, buffer+position[0]+1,  position[1]-position[0]-2);
    dimen_t nrows = (dimen_t) atoi(tmp);
    printf("%u\n", nrows);
    memset(tmp,0,strlen(tmp));
    // get ncols
    strncpy(tmp, buffer+position[1]+1,  position[2]-position[1]-2);
    dimen_t ncols = (dimen_t) atoi(tmp);
    printf("%u\n", ncols);
    memset(tmp,0,strlen(tmp));
#endif
    // set modulus
    dimen_t modulus = 31;

    fwrite(&nrows, sizeof(dimen_t), 1, fho);
    fwrite(&ncols, sizeof(dimen_t), 1, fho);
    fwrite(&modulus, sizeof(dimen_t), 1, fho);

    // get data, assume density of <=10%
    index_t nelts = (index_t)nrows * ncols / 2;
    index_t idx   = 0;
    index_t nnz   = 0;
    elt_t data[nelts];    // elements in matrix
    dimen_t colid[nelts]; // column id of elements
    dimen_t rlen[nrows];  // length of rows
    for (ctr=0; ctr<nrows; ++ctr) {
      rlen[ctr] = 0;
    }
    dimen_t ridx  = 0;
    while (fgets(buffer, sizeof(buffer), fhi) != NULL) {
      pos2  = strchr(buffer,'<');
      ctr = 0;
      while (pos2 != NULL) {
        position2[ctr] = pos2-buffer+1;
        ctr++;
        pos2=strchr(pos2+1,'<');
      }
      printf("d1\n");
      pos = strchr(buffer,',');
      ctr = 0;
      while (pos != NULL) {
        position[ctr] = pos-buffer+1;
        ctr++;
        pos=strchr(pos+1,',');
      }
      if (ctr%3 != 0) { // we are in the last row
        // we add the last position mark by searching for the last ">" in the
        // line and set the position mark one behind the index of that ">"
        pos = strrchr(buffer,'>');
        position[ctr] = pos-buffer+2;
        ctr++;
      }
      if (ctr != 0) {
        int ctr2  = 0;
        for (i=0; i<ctr; i=i+3) {
          printf("idx %u\n",idx);
          printf("%s || %d\n",buffer,position2[ctr2]);
          strncpy(tmp, buffer+position2[ctr2], 1);
          //strncpy(tmp, buffer+position2[ctr2], position[i]-position2[ctr2]-2);
          printf("%s\n",tmp);
          ridx = (index_t) atoi(tmp) -1;
          rlen[ridx]++;
          nnz++;
          printf("%d\n",ridx);
          memset(tmp,0,strlen(tmp));
          strncpy(tmp, buffer+position[i]+1, position[i+1]-position[i]-2);
          printf("%s\n",tmp);
          colid[idx]  = (dimen_t) atoi(tmp) - 1; // magma starts at column 1
          printf("colid[%d] = %u\n",idx, colid[idx]);
          memset(tmp,0,strlen(tmp));
          strncpy(tmp, buffer+position[i+1]+1, position[i+2]-position[i+1]-2);
          printf("%s\n",tmp);
          data[idx]  = (elt_t) atoi(tmp);
          printf("data[%d] = %u\n",idx, data[idx]);
          memset(tmp,0,strlen(tmp));
          idx++;
          ctr2++;
        }
      }
    printf("data[0] = %d\n",data[0]);
    }
    printf("nnz %d\n",nnz);
    for (i=0; i<nnz; ++i)
      printf("%u | %u\n",data[i],colid[i]);
    for (i=0; i<nrows; ++i)
      printf("%u -- %u\n",i,rlen[i]);

    fwrite(&nnz, sizeof(index_t), 1, fho);
    fwrite(data, sizeof(elt_t), nnz, fho);
    fwrite(colid, sizeof(dimen_t), nnz, fho);
    fwrite(rlen, sizeof(dimen_t), nrows, fho);

    fclose(fhi);
    fclose(fho);

    return 0;
  }
}
