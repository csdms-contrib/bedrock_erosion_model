#include<malloc.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#define FREE_ARG char*
#define NR_END 1

#define PI 3.141592653589793
#define sqrt2 1.414213562373
#define oneoversqrt2 0.707106781186
#define fillincrement 0.01

float delrho,alpha,*load,**erode,**slope,S,averosionrate,avreboundrate,avelevation;
float oneoverdeltax,oneoverdeltax2,**U,K,D,X,duration,timestep;
float *topovec,thresh,**deltah,**channel,**area,**sedflux,**sedfluxold,**flow1;
float **flow2,**flow3,**flow4,**flow5,**flow6,**flow7,**flow8,**flow,**topo,**topoold;
float **topo2,deltax,*ax,*ay,*bx,*by,*cx,*cy,*ux,*uy,*rx,*ry;
int *nn,*topovecind,lattice_size_x,lattice_size_y,*iup,*idown,*jup,*jdown;

float *vector(nl,nh)
long nh,nl;
/* allocate a float vector with subscript range v[nl..nh] */
{
        float *v;

        v=(float *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(float)));
        return v-nl+NR_END;
}

int *ivector(nl,nh)
long nh,nl;
/* allocate an int vector with subscript range v[nl..nh] */
{
        int *v;

        v=(int *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(int)));
        return v-nl+NR_END;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

void free_vector(float *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

int **imatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
      int  i,**m;

       /*allocate pointers to rows */
        m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
      m -= nrl;

       /*allocate rows and set pointers to them */
        for(i=nrl;i<=nrh;i++) {
                      m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
      }
       /* return pointer to array of pointers to rows */
        return m;
}

float **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    int i;
    float **m;

        /*allocate pointers to rows */
        m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
    m -= nrl;

   /*allocate rows and set pointers to them */
      for(i=nrl;i<=nrh;i++) {
                      m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float)
);
            m[i] -= ncl;
    }
      /* return pointer to array of pointers to rows */
      return m;
}

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

float ran3(idum)
int *idum;
{
        static int inext,inextp;
        static long ma[56];
        static int iff=0;
        long mj,mk;
        int i,ii,k;

        if (*idum < 0 || iff == 0) {
                iff=1;
                mj=MSEED-(*idum < 0 ? -*idum : *idum);
                mj %= MBIG;
                ma[55]=mj;
                mk=1;
                for (i=1;i<=54;i++) {
                        ii=(21*i) % 55;
                        ma[ii]=mk;
                        mk=mj-mk;
                        if (mk < MZ) mk += MBIG;
                        mj=ma[ii];
                }
                for (k=1;k<=4;k++)
                        for (i=1;i<=55;i++) {
                                ma[i] -= ma[1+(i+30) % 55];
                                if (ma[i] < MZ) ma[i] += MBIG;
                        }
                inext=0;
                inextp=31;
                *idum=1;
        }
        if (++inext == 56) inext=1;
        if (++inextp == 56) inextp=1;
        mj=ma[inext]-ma[inextp];
        if (mj < MZ) mj += MBIG;
        ma[inext]=mj;
        return mj*FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

float gasdev(idum)
int *idum;
{
        static int iset=0;
        static float gset;
        float fac,r,v1,v2;
        float ran3();

        if  (iset == 0) {
                do {
                        v1=2.0*ran3(idum)-1.0;
                        v2=2.0*ran3(idum)-1.0;
                        r=v1*v1+v2*v2;
                } while (r >= 1.0);
                fac=sqrt(-2.0*log(r)/r);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        } else {
                iset=0;
                return gset;
        }
}

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 100000

void indexx(n,arr,indx)
float arr[];
int indx[],n;
{
        unsigned long i,indxt,ir=n,itemp,j,k,l=1;
        int jstack=0,*istack;
        float a;

        istack=ivector(1,NSTACK);
        for (j=1;j<=n;j++) indx[j]=j;
        for (;;) {
                if (ir-l < M) {
                        for (j=l+1;j<=ir;j++) {
                                indxt=indx[j];
                                a=arr[indxt];
                                for (i=j-1;i>=1;i--) {
                                        if (arr[indx[i]] <= a) break;
                                        indx[i+1]=indx[i];
                                }
                                indx[i+1]=indxt;
                        }
                        if (jstack == 0) break;
                        ir=istack[jstack--];
                        l=istack[jstack--];
                } else {
                        k=(l+ir) >> 1;
                        SWAP(indx[k],indx[l+1]);
                        if (arr[indx[l+1]] > arr[indx[ir]]) {
                                SWAP(indx[l+1],indx[ir])
                        }
                        if (arr[indx[l]] > arr[indx[ir]]) {
                                SWAP(indx[l],indx[ir])
                        }
                        if (arr[indx[l+1]] > arr[indx[l]]) {
                                SWAP(indx[l+1],indx[l])
                        }
                        i=l+1;
                        j=ir;
                        indxt=indx[l];
                        a=arr[indxt];
                        for (;;) {
                                do i++; while (arr[indx[i]] < a);
                                do j--; while (arr[indx[j]] > a);
                                if (j < i) break;
                                SWAP(indx[i],indx[j])
                        }
                        indx[l]=indx[j];
                        indx[j]=indxt;
                        jstack += 2;
                        if (ir-i+1 >= j-l) {
                                istack[jstack]=ir;
                                istack[jstack-1]=i;
                                ir=j-1;
                        } else {
                                istack[jstack]=j-1;
                                istack[jstack-1]=l;
                                l=i;
                        }
                }
        }
        free_ivector(istack,1,NSTACK);
}
#undef M
#undef NSTACK
#undef SWAP

void tridag(float a[], float b[], float c[], float r[], float u[],
	unsigned long n)
{
	unsigned long j;
	float bet,*gam;

	gam=vector(1,n);
	u[1]=r[1]/(bet=b[1]);
	for (j=2;j<=n;j++) {
		gam[j]=c[j-1]/bet;
		bet=b[j]-a[j]*gam[j];
		u[j]=(r[j]-a[j]*u[j-1])/bet;
	}
	for (j=(n-1);j>=1;j--)
		u[j] -= gam[j+1]*u[j+1];
	free_vector(gam,1,n);
}

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void fourn(float data[], unsigned long nn[], int ndim, int isign)
{
	int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	float tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;

	for (ntot=1,idim=1;idim<=ndim;idim++)
		ntot *= nn[idim];
	nprev=1;
	for (idim=ndim;idim>=1;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
		ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for (i2=1;i2<=ip2;i2+=ip1) {
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					for (i2=i1;i2<=ip3;i2+=ifp2) {
						k1=i2;
						k2=k1+ifp1;
						tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
						tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}
#undef SWAP

void setupgridneighbors()
{    int i,j;

     idown=ivector(1,lattice_size_x);
     iup=ivector(1,lattice_size_x);
     jup=ivector(1,lattice_size_y);
     jdown=ivector(1,lattice_size_y);
     for (i=1;i<=lattice_size_x;i++)
      {idown[i]=i-1;
       iup[i]=i+1;}
     idown[1]=1;
     iup[lattice_size_x]=lattice_size_x;
     for (j=1;j<=lattice_size_y;j++)
      {jdown[j]=j-1;
       jup[j]=j+1;}
     jdown[1]=1;
     jup[lattice_size_y]=lattice_size_y;
}

void fillinpitsandflats(i,j)
int i,j;
{    float min;

     min=topo[i][j];
     if (topo[iup[i]][j]<min) min=topo[iup[i]][j];
     if (topo[idown[i]][j]<min) min=topo[idown[i]][j];
     if (topo[i][jup[j]]<min) min=topo[i][jup[j]];
     if (topo[i][jdown[j]]<min) min=topo[i][jdown[j]];
     if (topo[iup[i]][jup[j]]<min) min=topo[iup[i]][jup[j]];
     if (topo[idown[i]][jup[j]]<min) min=topo[idown[i]][jup[j]];
     if (topo[idown[i]][jdown[j]]<min) min=topo[idown[i]][jdown[j]];
     if (topo[iup[i]][jdown[j]]<min) min=topo[iup[i]][jdown[j]];
     if ((topo[i][j]<=min)&&(i>1)&&(j>1)&&(i<lattice_size_x)&&(j<lattice_size_y))
      {topo[i][j]=min+fillincrement;
       fillinpitsandflats(i,j);
       fillinpitsandflats(iup[i],j);
       fillinpitsandflats(idown[i],j);
       fillinpitsandflats(i,jup[j]);
       fillinpitsandflats(i,jdown[j]);
       fillinpitsandflats(iup[i],jup[j]);
       fillinpitsandflats(idown[i],jup[j]);
       fillinpitsandflats(idown[i],jdown[j]);
       fillinpitsandflats(iup[i],jdown[j]);}
}

void mfdflowroute(i,j)
int i,j;
{    float tot;

     tot=0;
     if (topo[i][j]>topo[iup[i]][j])
      tot+=pow(topo[i][j]-topo[iup[i]][j],1.1);
     if (topo[i][j]>topo[idown[i]][j])
      tot+=pow(topo[i][j]-topo[idown[i]][j],1.1);
     if (topo[i][j]>topo[i][jup[j]])
      tot+=pow(topo[i][j]-topo[i][jup[j]],1.1);
     if (topo[i][j]>topo[i][jdown[j]])
      tot+=pow(topo[i][j]-topo[i][jdown[j]],1.1);
     if (topo[i][j]>topo[iup[i]][jup[j]])
      tot+=pow((topo[i][j]-topo[iup[i]][jup[j]])*oneoversqrt2,1.1);
     if (topo[i][j]>topo[iup[i]][jdown[j]])
      tot+=pow((topo[i][j]-topo[iup[i]][jdown[j]])*oneoversqrt2,1.1);
     if (topo[i][j]>topo[idown[i]][jup[j]])
      tot+=pow((topo[i][j]-topo[idown[i]][jup[j]])*oneoversqrt2,1.1);
     if (topo[i][j]>topo[idown[i]][jdown[j]])
      tot+=pow((topo[i][j]-topo[idown[i]][jdown[j]])*oneoversqrt2,1.1);
     if (topo[i][j]>topo[iup[i]][j])
      flow1[i][j]=pow(topo[i][j]-topo[iup[i]][j],1.1)/tot;
       else flow1[i][j]=0;
     if (topo[i][j]>topo[idown[i]][j])
      flow2[i][j]=pow(topo[i][j]-topo[idown[i]][j],1.1)/tot;
       else flow2[i][j]=0;
     if (topo[i][j]>topo[i][jup[j]])
      flow3[i][j]=pow(topo[i][j]-topo[i][jup[j]],1.1)/tot;
       else flow3[i][j]=0;
     if (topo[i][j]>topo[i][jdown[j]])
      flow4[i][j]=pow(topo[i][j]-topo[i][jdown[j]],1.1)/tot;
       else flow4[i][j]=0;
     if (topo[i][j]>topo[iup[i]][jup[j]])
      flow5[i][j]=pow((topo[i][j]-topo[iup[i]][jup[j]])*oneoversqrt2,1.1)/tot;
       else flow5[i][j]=0;
     if (topo[i][j]>topo[iup[i]][jdown[j]])
      flow6[i][j]=pow((topo[i][j]-topo[iup[i]][jdown[j]])*oneoversqrt2,1.1)/tot;
       else flow6[i][j]=0;
     if (topo[i][j]>topo[idown[i]][jup[j]])
      flow7[i][j]=pow((topo[i][j]-topo[idown[i]][jup[j]])*oneoversqrt2,1.1)/tot;
       else flow7[i][j]=0;
     if (topo[i][j]>topo[idown[i]][jdown[j]])
      flow8[i][j]=pow((topo[i][j]-topo[idown[i]][jdown[j]])*oneoversqrt2,1.1)/tot;
       else flow8[i][j]=0;
     flow[iup[i]][j]+=flow[i][j]*flow1[i][j];
     flow[idown[i]][j]+=flow[i][j]*flow2[i][j];
     flow[i][jup[j]]+=flow[i][j]*flow3[i][j];
     flow[i][jdown[j]]+=flow[i][j]*flow4[i][j];
     flow[iup[i]][jup[j]]+=flow[i][j]*flow5[i][j];
     flow[iup[i]][jdown[j]]+=flow[i][j]*flow6[i][j];
     flow[idown[i]][jup[j]]+=flow[i][j]*flow7[i][j];
     flow[idown[i]][jdown[j]]+=flow[i][j]*flow8[i][j];
}

void calculatealongchannelslope(i,j)
int i,j;
{    float down;

     down=0;
     if (topo[iup[i]][j]-topo[i][j]<down) down=topo[iup[i]][j]-topo[i][j];
     if (topo[idown[i]][j]-topo[i][j]<down) down=topo[idown[i]][j]-topo[i][j];
     if (topo[i][jup[j]]-topo[i][j]<down) down=topo[i][jup[j]]-topo[i][j];
     if (topo[i][jdown[j]]-topo[i][j]<down) down=topo[i][jdown[j]]-topo[i][j];
     if ((topo[iup[i]][jup[j]]-topo[i][j])*oneoversqrt2<down)
      down=(topo[iup[i]][jup[j]]-topo[i][j])*oneoversqrt2;
     if ((topo[idown[i]][jup[j]]-topo[i][j])*oneoversqrt2<down)
      down=(topo[idown[i]][jup[j]]-topo[i][j])*oneoversqrt2;
     if ((topo[iup[i]][jdown[j]]-topo[i][j])*oneoversqrt2<down)
      down=(topo[iup[i]][jdown[j]]-topo[i][j])*oneoversqrt2;
     if ((topo[idown[i]][jdown[j]]-topo[i][j])*oneoversqrt2<down)
      down=(topo[idown[i]][jdown[j]]-topo[i][j])*oneoversqrt2;
     slope[i][j]=fabs(down)/deltax;
}

void computeflexure()
{    int i,j,index;
     float fact;

     for (j=1;j<=lattice_size_y;j++)
       for (i=1;i<=lattice_size_x;i++)
        {load[2*(i-1)*lattice_size_y+2*j-1]=erode[i][j];
         load[2*(i-1)*lattice_size_y+2*j]=0.0;}
     fourn(load,nn,2,1);
     load[1]*=1/delrho;
     load[2]*=1/delrho;
     for (j=1;j<=lattice_size_y/2;j++)
       {fact=1/(delrho+pow(alpha*j*PI/lattice_size_y,4.0));
        load[2*j+1]*=fact;
        load[2*j+2]*=fact;
        load[2*lattice_size_y-2*j+1]*=fact;
        load[2*lattice_size_y-2*j+2]*=fact;}
     for (i=1;i<=lattice_size_x/2;i++)
       {fact=1/(delrho+pow(alpha*i*PI/lattice_size_x,4.0));
        load[2*i*lattice_size_y+1]*=fact;
        load[2*i*lattice_size_y+2]*=fact;
        load[2*lattice_size_x*lattice_size_y-2*i*lattice_size_y+1]*=fact;
        load[2*lattice_size_x*lattice_size_y-2*i*lattice_size_y+2]*=fact;}
     for (i=1;i<=lattice_size_x/2;i++)
       for (j=1;j<=lattice_size_y/2;j++)
         {fact=1/(delrho+pow(alpha*sqrt((PI*PI*i*i)/(lattice_size_x*lattice_size_x)
        +(PI*PI*j*j)/(lattice_size_y*lattice_size_y)),4.0));
          index=2*i*lattice_size_y+2*j+1;
          load[index]*=fact;
          load[index+1]*=fact;
          index=2*i*lattice_size_y+2*lattice_size_y-2*j+1;
          load[index]*=fact;
          load[index+1]*=fact;
          index=2*lattice_size_x*lattice_size_y-2*(i-1)*lattice_size_y-2*(lattice_size_y-j)+1;
          load[index]*=fact;
          load[index+1]*=fact;
          index=2*lattice_size_x*lattice_size_y-2*lattice_size_y-
           2*(i-1)*lattice_size_y+2*(lattice_size_y-j)+1;
          load[index]*=fact;
          load[index+1]*=fact;}
     fourn(load,nn,2,-1);
     for (j=2;j<=lattice_size_y-1;j++)
       for (i=2;i<=lattice_size_x-1;i++)
        {U[i][j]=delrho*(1-delrho)*(load[2*(i-1)*lattice_size_y+2*j-1]
          /(lattice_size_x*lattice_size_y);}
}

void setupmatrices()
{    int i,j;

     ax=vector(1,lattice_size_x);
     ay=vector(1,lattice_size_y);
     bx=vector(1,lattice_size_x);
     by=vector(1,lattice_size_y);
     cx=vector(1,lattice_size_x);
     cy=vector(1,lattice_size_y);
     ux=vector(1,lattice_size_x);
     uy=vector(1,lattice_size_y);
     rx=vector(1,lattice_size_x);
     ry=vector(1,lattice_size_y);
     topo=matrix(1,lattice_size_x,1,lattice_size_y);
     topo2=matrix(1,lattice_size_x,1,lattice_size_y);
     topoold=matrix(1,lattice_size_x,1,lattice_size_y);
     sedflux=matrix(1,lattice_size_x,1,lattice_size_y);
     sedfluxold=matrix(1,lattice_size_x,1,lattice_size_y);
     area=matrix(1,lattice_size_x,1,lattice_size_y);
     slope=matrix(1,lattice_size_x,1,lattice_size_y);
     deltah=matrix(1,lattice_size_x,1,lattice_size_y);
     load=vector(1,2*lattice_size_x*lattice_size_y);
     erode=matrix(1,lattice_size_x,1,lattice_size_y);
     U=matrix(1,lattice_size_x,1,lattice_size_y);
     channel=matrix(1,lattice_size_x,1,lattice_size_y);
     flow=matrix(1,lattice_size_x,1,lattice_size_y);
     flow1=matrix(1,lattice_size_x,1,lattice_size_y);
     flow2=matrix(1,lattice_size_x,1,lattice_size_y);
     flow3=matrix(1,lattice_size_x,1,lattice_size_y);
     flow4=matrix(1,lattice_size_x,1,lattice_size_y);
     flow5=matrix(1,lattice_size_x,1,lattice_size_y);
     flow6=matrix(1,lattice_size_x,1,lattice_size_y);
     flow7=matrix(1,lattice_size_x,1,lattice_size_y);
     flow8=matrix(1,lattice_size_x,1,lattice_size_y);
     topovec=vector(1,lattice_size_x*lattice_size_y);
     topovecind=ivector(1,lattice_size_x*lattice_size_y);
}

main()
{    FILE *fp0,*fp1,*fp2,*fp3,*fp4,*fp5;
     float change,maxelevation,capacity,time,max,hillslopeerosionrate;
     int printinterval,i,j,t;

     fp0=fopen("sedfluxdriventopoinitial","r");
     fp1=fopen("sedfluxdriventopo","w");
     fp2=fopen("sedfluxdrivenerosion","w");
     fp3=fopen("sedfluxdrivensedflux","w");
     fp4=fopen("sedfluxdrivenrebound","w");
     fp5=fopen("modeltimehistory","w");
     lattice_size_x=256;
     lattice_size_y=256;
     deltax=500;   /* m */
     delrho=0.28;  /* (rho_m-rho_c)/rho_c */
     oneoverdeltax=1.0/deltax;
     oneoverdeltax2=1.0/(deltax*deltax);
     timestep=1.0; /* kyr */
     alpha=200;    /* multiples of deltax */
     hillslopeerosionrate=0.01; /* m/kyr */
     K=0.0005;     /* m^1/2/kyr */
     D=1.0;        /* m^2/kyr */
     X=0.005;      /* m^-1 */
     alpha/=deltax;
     setupmatrices();
     setupgridneighbors();
     for (j=1;j<=lattice_size_y;j++)
      for (i=1;i<=lattice_size_x;i++)
       {fscanf(fp0,"%f",&topo[i][j]);topo[i][j]/=5;
        if ((i==1)||(j==1)||(i==lattice_size_x)||(j==lattice_size_y)) topo[i][j]=0;
        topoold[i][j]=topo[i][j];
        flow[i][j]=deltax*deltax;
        area[i][j]=0;
        U[i][j]=0;
        erode[i][j]=0;
        channel[i][j]=0;
        sedflux[i][j]=0;
        sedfluxold[i][j]=0;}
     time=0;
     nn=ivector(1,2);
     nn[1]=lattice_size_x;
     nn[2]=lattice_size_y;
     duration=40000;         /* 40 Myr */
     printinterval=40000;
     while (time<duration)
      {for (j=1;j<=lattice_size_y;j++)
        for (i=1;i<=lattice_size_x;i++)
         topoold[i][j]=topo[i][j];
       /*compute uplift rate*/
       if (time>1000) computeflexure();  /* compute isostatic uplift */
       else
        {/* or prescribe a tectonic uplift rate */
         for (j=2;j<=lattice_size_y-1;j++)
          for (i=2;i<=lattice_size_x-1;i++)
           U[i][j]=1.0;  /* m/kyr */}
       avelevation=0;avreboundrate=0;maxelevation=0;
       for (j=2;j<=lattice_size_y-1;j++)
        for (i=2;i<=lattice_size_x-1;i++)
         {avelevation+=topo[i][j];
          if (topo[i][j]>maxelevation) maxelevation=topo[i][j];
          avreboundrate+=U[i][j];
          topoold[i][j]+=U[i][j]*timestep;
          topo[i][j]+=U[i][j]*timestep;
          deltah[i][j]=0;
          channel[i][j]=0;}
       avelevation/=(lattice_size_x-2)*(lattice_size_y-2);
       avreboundrate/=(lattice_size_x-2)*(lattice_size_y-2);
       printf("%f %f %f %f %f\n",time,avelevation,averosionrate,avreboundrate,maxelevation);
       fprintf(fp5,"%f %f %f %f %f\n",time,avelevation,averosionrate,avreboundrate,maxelevation);
       /* hillslope erosion can occur by prescribing a uniform rate, as done here
	      or using the ADI technique to solve the diffusion equation on hillslopes */
       for (j=2;j<=lattice_size_y-1;j++)
        for (i=2;i<=lattice_size_x-1;i++)
         if (channel[i][j]==0) topo[i][j]-=hillslopeerosionrate*timestep;
       for (j=2;j<=lattice_size_y-1;j++)
        for (i=2;i<=lattice_size_x-1;i++)
         erode[i][j]=(topoold[i][j]-topo[i][j])/timestep;
       for (j=1;j<=lattice_size_y;j++)
        for (i=1;i<=lattice_size_x;i++)
         fillinpitsandflats(i,j);
       /*route water from highest gridpoint to lowest*/
       for (j=1;j<=lattice_size_y;j++)
        for (i=1;i<=lattice_size_x;i++)
         flow[i][j]=deltax*deltax;
       for (j=1;j<=lattice_size_y;j++)
        for (i=1;i<=lattice_size_x;i++)
         topovec[(j-1)*lattice_size_x+i]=topo[i][j];
       indexx(lattice_size_x*lattice_size_y,topovec,topovecind);
       t=lattice_size_x*lattice_size_y+1;
       while (t>1)
        {t--;
         i=(topovecind[t])%lattice_size_x;
         if (i==0) i=lattice_size_x;
         j=(topovecind[t])/lattice_size_x+1;
         if (i==lattice_size_x) j--;
         mfdflowroute(i,j);}
       for (i=2;i<=lattice_size_x-1;i++)
        for (j=2;j<=lattice_size_y-1;j++)
         area[i][j]=flow[i][j];
       /*perform upwind differencing */
       max=0;
       for (i=2;i<=lattice_size_x-1;i++)
        for (j=2;j<=lattice_size_y-1;j++)
         {calculatealongchannelslope(i,j);
          capacity=slope[i][j]*sqrt(area[i][j]);
          if (capacity>1/X)
           {change=timestep*K*sqrt(fabs(sedflux[i][j]))*deltax*slope[i][j];
            deltah[i][j]+=change;
            erode[i][j]+=change/timestep;
            if (deltah[i][j]<0) deltah[i][j]=0;
            channel[i][j]=1;}
          topo[i][j]-=deltah[i][j];
          if (topo[i][j]<0) topo[i][j]=0;
          if (K*sqrt(fabs(sedflux[i][j]))*deltax>max)
           max=K*sqrt(fabs(sedflux[i][j]))*deltax;}
       for (j=1;j<=lattice_size_y;j++)
        for (i=1;i<=lattice_size_x;i++)
         fillinpitsandflats(i,j);
       averosionrate=0;
       for (i=2;i<=lattice_size_x-1;i++)
        for (j=2;j<=lattice_size_y-1;j++)
         {flow[i][j]=erode[i][j];
          averosionrate+=erode[i][j];}
       averosionrate/=(lattice_size_x-2)*(lattice_size_y-2);
       for (j=1;j<=lattice_size_y;j++)
        for (i=1;i<=lattice_size_x;i++)
         topovec[(j-1)*lattice_size_x+i]=topo[i][j];
       indexx(lattice_size_x*lattice_size_y,topovec,topovecind);
       t=lattice_size_x*lattice_size_y+1;
       while (t>1)
        {t--;
         i=(topovecind[t])%lattice_size_x;
         if (i==0) i=lattice_size_x;
         j=(topovecind[t])/lattice_size_x+1;
         if (i==lattice_size_x) j--;
         mfdflowroute(i,j);}
       for (i=2;i<=lattice_size_x-1;i++)
        for (j=2;j<=lattice_size_y-1;j++)
         {sedfluxold[i][j]=sedflux[i][j];
          sedflux[i][j]=flow[i][j];}
       time+=timestep;
       if (max>5.0*deltax/timestep)
        {time-=timestep;
         timestep/=2.0;
         for (i=2;i<=lattice_size_x-1;i++)
          for (j=2;j<=lattice_size_y-1;j++)
           {sedflux[i][j]=sedfluxold[i][j];
            topo[i][j]=topoold[i][j]-U[i][j]*timestep;}}
        else
         {if (max<0.5*deltax/timestep) timestep*=1.2;}
       if (time>=printinterval)
        {printinterval+=40000;
         for (i=1;i<=lattice_size_x;i++)
          for (j=1;j<=lattice_size_y;j++)
           {fprintf(fp1,"%f\n",topo[i][j]);
            fprintf(fp2,"%f\n",erode[i][j]);
            fprintf(fp3,"%f\n",sedflux[i][j]);
            fprintf(fp4,"%f\n",U[i][j]);}}}
}
