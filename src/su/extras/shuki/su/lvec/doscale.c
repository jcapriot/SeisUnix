doscale(scale,n,p)
float scale,*p;
{
	nn;
	float *pp,t;
	pp = p;
	t = scale;
	nn = n;
	do {
		*(pp++) *= t;
	} while(--nn);
}
