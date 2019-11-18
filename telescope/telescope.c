#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "telescope.h"


int main (int argc, char **argv) {
	
	shockrad *item = NULL;
	double i;

	for (i = 0; i < 10; i++) {

		//item = (shockrad*) malloc( sizeof(shockrad) );
		//item->dObsTime = i;
		//HASH_ADD(hh, map, dObsTime, sizeof(double), item);

		sendToTelescope(i, 1, 1);
		sendToTelescope(i, 2, 2);
		sendToTelescope(i, 3, 3);
		sendToTelescope(i, 3, 3);
	}

	displayTelescopeData();


	return(0);
}
