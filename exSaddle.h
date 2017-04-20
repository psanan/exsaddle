#ifndef EXSADDLE_H_
#define EXSADDLE_H_


/* Note that NSD (number of spatial dimensions) is expected via -D 
   and is assumed to be either 2 or 3 */
#if (NSD!=2 && NSD!=3)
#error "NSD must be equal to 2 or 3"
#endif

#define U_DOFS     NSD       /* degrees of freedom per velocity node. This can only be NSD with this code.*/
#define P_DOFS     1         /* degrees of freedom per pressure node. This can only be 1 with this code. */
#if NSD==2
#define MAX_QUADP  9
#define Q2_BASIS   9         /* nodes per element */
#define Q1_BASIS   4         /* nodes per element */
#elif NSD==3
#define MAX_QUADP  27 
#define Q2_BASIS   27
#define Q1_BASIS   8
#endif
#define U_BASIS    Q2_BASIS
#define P_BASIS    Q1_BASIS
#define MG_DEPTH   10        /* Maximum MG levels */

#endif
