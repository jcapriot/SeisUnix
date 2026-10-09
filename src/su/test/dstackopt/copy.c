copy(p,q,n)
char *p,*q;
{
	char *pp,*qq;
	pp = p; qq = q;
	while(n--) *(qq++) = *(pp++);
}
