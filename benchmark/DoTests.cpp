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
	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s100000 1000 100000000 >> 10in5ld.out");
	// 10^5    numbers
	// e^434   exp

	// 1000    tests
	// ~100    sec	 	(~1.5 min)
	
	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s1000000 1000 1000000000 >> 10in6ld.out");
	// 10^6    numbers
	// e^434   exp
	
	// 1000    tests
	// ~1000   sec	 	(~20 min)

	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s10000000 1000 1000000000 >> 10in7ld.out");
	// 10^7    numbers
	// e^434   exp

	// 100     tests
	// ~1000   sec	 	(~20 min)
	
	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s100000000 1000 10000000000 >> 10in8ld.out");
	// 10^8    numbers
	// e^434   exp
	
	// 100     tests
	// ~10000  sec		(~2.5 hour)





	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s100000 10000 100000000 >> 10in5ld.out");
	// 10^5    numbers
	// e^4340  exp

	// 1000    tests
	// ~100    sec	 	(~1.5 min)
	
	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s1000000 10000 1000000000 >> 10in6ld.out");
	// 10^6    numbers
	// e^4340  exp
	
	// 1000    tests
	// ~1000   sec	 	(~20 min)

	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s10000000 10000 1000000000 >> 10in7ld.out");
	// 10^7    numbers
	// e^4340  exp

	// 100     tests
	// ~1000   sec	 	(~20 min)
	
	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s100000000 10000 10000000000 >> 10in8ld.out");
	// 10^8    numbers
	// e^4340  exp
	
	// 100     tests
	// ~10000  sec		(~2.5 hour)




	// About 6 hours in total
	return;
}