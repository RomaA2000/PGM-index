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
	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s100000 100 >> 10in5ld.out");
	// 10^5    numbers
	// e^43    exp

	// 1e8     querries
	// ~1000   sec	 	(~1 hour)
	
	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s1000000 100 >> 10in6ld.out");
	// 10^6    numbers
	// e^43    exp
	
	// 1e8     querries
	// ~1000   sec	 	(~1 hour)

	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s10000000 100 >> 10in7ld.out");
	// 10^7    numbers
	// e^43    exp

	// 1e8     querries
	// ~1000   sec	 	(~1 hour)
	
	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s100000000 100 >> 10in8ld.out");
	// 10^8    numbers
	// e^43    exp
	
	// 1e8     querries
	// ~1000   sec		(~1 hour)




/*
	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s100000 1000 >> 10in5ld.out");
	// 10^5    numbers
	// e^434   exp

	// 1e8     querries
	// ~1000   sec	 	(~1 hour)
	
	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s1000000 1000 >> 10in6ld.out");
	// 10^6    numbers
	// e^434   exp
	
	// 1e8     querries
	// ~1000   sec	 	(~1 hour)

	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s10000000 1000 >> 10in7ld.out");
	// 10^7    numbers
	// e^434   exp

	// 1e8     querries
	// ~1000   sec	 	(~1 hour)
	
	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s100000000 1000 >> 10in8ld.out");
	// 10^8    numbers
	// e^434   exp
	
	// 1e8     querries
	// ~1000   sec		(~1 hour)





	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s100000 1000 >> 10in5ld.out");
	// 10^5    numbers
	// e^4340  exp

	// 1e8     querries
	// ~1000   sec	 	(~1 hour)
	
	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s1000000 1000 >> 10in6ld.out");
	// 10^6    numbers
	// e^4340  exp
	
	// 1e8     querries
	// ~1000   sec	 	(~1 hour)

	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s10000000 1000 >> 10in7ld.out");
	// 10^7    numbers
	// e^4340  exp

	// 1e8     querries
	// ~1000   sec	 	(~1 hour)
	
	system ("cmake . -DCMAKE_BUILD_TYPE=Release ; make benchmark -j8 ; benchmark/benchmark -s100000000 1000 >> 10in8ld.out");
	// 10^8    numbers
	// e^4340  exp
	
	// 1e8     querries
	// ~1000   sec		(~1 hour)

*/


	// About 6 hours in total
	return;
}