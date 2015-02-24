#include "3dbpp.cpp"
// Not entirely working. Take a look at: 
// https://bitbucket.org/krris/3d-bin-packing

int main(int argc, char const *argv[])
{


	// int Nrbc=80, Nplt=10, Nwbc=5, Nsickle=30;
	int Nrbc=103, Nplt=0, Nwbc=0, Nsickle=0;
	int N[4] = {Nrbc, Nplt, Nwbc, Nsickle};
	int number = Nrbc + Nplt + Nwbc + Nsickle;
	int fac=1;
	double Widths[4]  = {7.88, 2, 10, 6};
	double Heights[4] = {2.7, 2, 10, 6};
	double Depths[4]  = {7.88, 1, 10, 4};
	int w[number], h[number], d[number];
	int x[number], y[number], z[number], bno[number];
	int i,j,s=0;
	for (i = 0; i < 4; ++i)
	{
		for (j = 0; j < N[i]; ++j)
		{
			w[s] = int(ceil(Widths[i] * fac));
			h[s] = int(ceil(Heights[i] * fac));
			d[s] = int(ceil(Depths[i] * fac));
			s++;
		}
	}
	int W=32, H=32, D=32;
	int nodelimit=0, iterlimit=0, timelimit=3;
	int packingtype=0;

	int ub, lb;
	int nodeused, iterused, timeused;
	try {
		binpack3d(number, W * fac, H * fac, D * fac,
	               w, h, d, 
	               x, y, z, bno,
	               &lb, &ub, 
	              nodelimit, iterlimit, timelimit, 
	              &nodeused, &iterused, &timeused, 
	              packingtype);
	} catch (int e) {
		printf("An exception occurred. Exception Nr. %d \n", e);
	}

	for (i = 0; i < number; ++i)
	{
		printf("pack (%f,%f,%f) \n", x[i] * 1.0/fac, y[i] * 1.0/fac, z[i] * 1.0/fac);
	}

	printf("n           = %d\n", number);
	printf("nodelimit   = %d\n", nodelimit);
	printf("iterlimit   = %d\n", iterlimit);
	printf("timelimit   = %d\n", timelimit);
	printf("nodeused = %d\n", nodeused);
	printf("iterused = %d\n", iterused);
	printf("timeused = %d\n", timeused);
	printf("packingtype = %d\n", packingtype);

	return 0;
}