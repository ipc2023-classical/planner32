//If DEBUG IS SET, debug messages will be print
#ifndef DEBUG_MSG

#include "globals.h"


#ifdef DEBUG
   #define DEBUG_PHPDBS(str) do { str } while ( false )
   #define DEBUG_MSG(str) do { str } while( false )
   #define DEBUG_MAS(str) do { str } while( false )
#else
//#define DEBUG_MSG(str) do {cout << "  " << g_timer() << ": "; str } while ( false )
#define DEBUG_MSG(str) do {if(false) { str }} while(false)
#define DEBUG_PHPDBS(str) do {if(false) { str }} while(false)
//#define DEBUG_PHPDBS(str) do { str } while ( false )

//   #define DEBUG_MAS(str) do { str } while( false )
#define DEBUG_MAS(str) do { if(false) { str } } while( false )

#endif
#endif
//#define DEBUG_PHPDBS(str) do {cout << "  " << g_timer() << ": "; str } while ( false )
//#define DEBUG_TIMES(str) do {if(false) { str }}while(false)
