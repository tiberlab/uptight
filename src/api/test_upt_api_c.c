#include <stdio.h>
#include "uptight.h"

int main(int *argc, char **argv)
{
  int handler[UPT_HSIZE];
  int i;


  printf("Handler size in bytes= %d\n",sizeof(handler) );
  
  printf("Initialising UPTIGHT\n");
  f77_upt_initsession(handler);

  printf("Handler received:");
  for(i=0;i<UPT_HSIZE;i++)
  {  
     printf(" %d",handler[i]);
  }
  printf("\n");

  printf("Destructing UPTIGHT\n");

  f77_upt_destructsession(handler);

  printf("Done\n");

}
