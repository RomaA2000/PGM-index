#include <iostream>
#include <stdlib.h>

void StartTests ();

int main ()
{
	StartTests();

	return 0;
}


void StartTests ()
{	
	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s1000000000 1000 10000000000 >> 10in9ld.out");
	// 10^9    numbers
	// e^434  exp
	
	// 10     tests
	// ~10000 sec		(~2.5 hour)

	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s1000000000 1000 10000000000 >> 10in9ld.out");
	// 10^9    numbers
	// e^4340  exp
	
	// 10     tests
	// ~10000 sec		(~2.5 hour)

	// About 5 hours in total
	return;
}