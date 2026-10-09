float powerof(p,n)
float *p;
{
	nn;
	float *pp,ppower;

	nn = n;
	pp = p;
	ppower = 0.0;
	while(nn--) {
		ppower += (*pp) * (*pp);
		pp++;
	}
	return(ppower);
}
