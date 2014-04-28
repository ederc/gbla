#include "gb.h"
void print_help() {
  printf("\n");
  printf("DESCRIPTION\n");
  printf("       Computes Gaussian Elimination of a structured hybrid sparse-\n");
  printf("       dense matrix coming from Groebner basis computations.\n");

  printf("OPTIONS\n");
  printf("       -f FILE    file the input matrix is stored in\n");
  printf("       -h         print help\n");
  printf("       -m MEM     free memory on the go\n");
  printf("                  default: 0\n");
  printf("       -t THRDS   number of threads\n");
  printf("                  default: 1\n");
  printf("       -v VALID   validate results with structured Gaussian Elimination\n");
  printf("                  default: 0\n");
  printf("\n");

  return; 
}


int main(int argc, char *argv[])
{  
  const char *file_name = NULL;
  int free_mem          = 0;
  int validate_results  = 0;
  int method            = 0;
  int nthrds            = 1;
  
  int opt;

  opterr  = 0;

  while ((opt = getopt(argc, argv, "f:hmtv")) != -1)
  {
    switch (opt)
    {
      case 'h':
        print_help();
        return 0;
      case 'f':
        file_name = optarg;
        break;
      case 'm': 
        free_mem = atoi(optarg);
        break;
      case 't': 
        nthrds  = atoi(optarg);
        break;
      case 'v': 
        validate_results  = atoi(optarg);
        break;
      case '?':
        if (optopt == 'f')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
                    "Unknown option character `\\x%x'.\n",
                    optopt);
        return 1;
      default:
        abort ();
    }
  }

  if (file_name == NULL)
  {
    fprintf(stderr, "File name is required.\nSee help using '-h' option.\n");
    return 1;
  }

  printf("File %c\nNbr Threads %d\n",file_name,nthrds);

  return 0;
}
