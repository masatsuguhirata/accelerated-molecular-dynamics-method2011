//実行方法
//mpicc -O0 -o filename filename.c -lm
//mpirun -np (procs) ./filename -a (N) (TMAX) (blocksize) (cutoff distance)
//exp)mpirun -np 8 ./Ar -a 128 0.5 16 20 

//???cutoffNの削除

//mpicc -DN=粒子数 -fast -o filename filename.c -lm

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define	LAYER_NUMBER	14		//[layer number of cluster]
#define	EPS	16.7			//[x10^-15 erg]
#define	SIGM	3.4				//[A]
#define	MASS	6.636			//[x10^-23 g]
#define	d_t 0.005				// [x10^-12 s] Δt
//#define TMAX 0.5			//[x10^-12 s]
#define TAVE 2000

//#define FILEOUT
//#ifndef N
//#define	N 128
//#endif
//#define bs 64
//#define drand//通常は隠す,確認用

void ArConf (double x0[], double y0[], double z0[],int N);// 初期配置を生成する関数

int main(int argc,char *argv[]){

    int N;//argv
    int bs;//argv
    double TMAX,cutoffL;//argv
    if (*argv[1] == '-')
	if (*(argv[1]+1) == 'a'){
	    N = atoi(argv[2]);
	    TMAX = atof(argv[3]);
	    bs = atoi(argv[4]);
	    cutoffL=atof(argv[5]);
	}
    cutoffL=cutoffL*cutoffL;
    double rri,ri6,x,y,z,tmp,coeff;
    int i,j;
    int ib,jb,nb,jstart,iend,jend;
    //int ij;
    int step;
    int t_ave=1;
    double cutoffR;		//cutoff rate
    long int cutoffN=0;//cutoff number
    long int loopN;		//total number of loop 
    double *x0,*y0,*z0,*x1,*y1,*z1;// t,t+Δtでの位置
    double *xAve,*yAve,*zAve;// t,t+Δtでの位置
    double *vx0,*vy0,*vz0,*vx1,*vy1,*vz1;// t,t+Δtでの速度
    double *ax0,*ay0,*az0;
    double *master;//double *work;//不要
    double *ax1,*ay1,*az1;//double  ax1[N],ay1[N],az1[N];// t,t+Δtでの加速度
    double t=0.0;
    double a,b,c,d,ddt,dt2,m2,sigm3,sigm6,rndx,rndy,rndz;
    double E,U=0.0,T;
    double t1,t2,t3,t4,t5,t6,t33;
    double t16,t26,t34,t45;
    int MAXSTEP;
    double maxbox,minbox;
    FILE * fp;//ich
    FILE * fp2;//v
    FILE * fp3;//EK
    FILE * fp4;//U
    FILE * fp5;//Uk+U
    FILE * fp6;//T
    int fflush(FILE *steram);
    master=(double*)malloc(sizeof(double)*(21*N)+1);
    ax1=master+1;ay1=master+N+1;az1=master+2*N+1;
    x0=master+3*N+1;y0=master+4*N+1;z0=master+5*N+1;
    x1=master+6*N+1;y1=master+7*N+1;z1=master+8*N+1;
    xAve=master+9*N+1;yAve=master+10*N+1;zAve=master+11*N+1;
    vx0=master+12*N+1;vy0=master+13*N+1;vz0=master+14*N+1;
    vx1=master+15*N+1;vy1=master+16*N+1;vz1=master+17*N+1;
    ax0=master+18*N+1;ay0=master+19*N+1;az0=master+20*N+1;


    sigm3 = SIGM*SIGM*SIGM;
    sigm6 = sigm3*sigm3;
    a = 2.0*sigm6;// 2σ^66
    b = 12.0*a*EPS/MASS;// 24εσ^6/m
    dt2	= d_t*0.5;
    ddt = dt2*d_t;
    m2 = MASS*0.5;
    c = 4.0*EPS*sigm6;
    d = 20/(4.14*(N-2));

//前準備　クエンチ・アニール
    for(i=0;i<N;i++){// 初期化
	ax0[i] = vx0[i] = x0[i] = xAve[i] = 0.0;
	ay0[i] = vy0[i] = y0[i] = yAve[i] = 0.0;
	az0[i] = vz0[i] = z0[i] = zAve[i] = 0.0;
    }

    ArConf(x0,y0,z0,N);// 初期位置の計算

/*
#ifdef drand//確認 
    srand48(123);
#else
    srand(time(NULL));// 乱数の初期化
#endif	
    for(i = 1; i < N; i++){
#ifdef drand
	rndx=(drand48()-0.5)*0.2;//0.1未満の乱数計算
	rndy=(drand48()-0.5)*0.2;
	rndz=(drand48()-0.5)*0.2;
#else
	rndx = (double)rand() / 2147483647;// 1未満の乱数計算
	rndy = (double)rand() / 2147483647;
	rndz = (double)rand() / 2147483647;
#endif
	x0[i] += rndx;
	y0[i] += rndy;
	z0[i] += rndz;	
    }
*/
    
    for(i=0;i<N;i++){// a(t)=F(t)/m
	for(j=i+1;j<N;j++){
	    x = x0[j]-x0[i];
	    y = y0[j]-y0[i];
	    z = z0[j]-z0[i];
	    rri = 1.0/(x*x + y*y + z*z);
	    ri6 = rri*rri*rri;
	    coeff = b*ri6*rri*(a*ri6-1);

	    ax0[i] -= tmp = coeff*x;
	    ax0[j] += tmp;

	    ay0[i] -= tmp = coeff*y;
	    ay0[j] += tmp;

	    az0[i] -= tmp = coeff*z;
	    az0[j] += tmp;

	    U += c*ri6*(sigm6*ri6-1);// U(0)の計算    
	}
    }
/*
#ifdef FILEOUT
    fp = fopen("ich-cutoff.dat","w");
    fp2 = fopen("v-cutoff.dat","w");
    fp3 = fopen("Ek-cutoff.dat","w");
    fp4 = fopen("U-cutoff.dat","w");
    fp5 = fopen("Ek_U-cutoff.dat","w");
    fp6 = fopen("T-cutoff.dat","w");

    fprintf(fp,"%d\n",N);	// r(0)の出力
    for (i=0; i<N; ++i) {
	fprintf(fp,"Ar %14.6f %14.6f %14.6f\n",x0[i],y0[i],z0[i]);
    }
    fprintf(fp,"t=%14.6f\n",t);
    fprintf(fp,"\n");
    fprintf(fp2,"t=%14.6f\nAr	x	y	z\n",t);// v(0)の出力
    for (i=0; i<N; ++i) {
	fprintf(fp2,"%d	%14.6f	%14.6f	%14.6f\n",i+1,vx0[i],vy0[i],vz0[i]);
    }
    fprintf(fp2,"\n");
    fprintf(fp3,"t	Ek\n0	0\n");	// Ek(0)の出力
    fprintf(fp4,"t	U\n%14.6f	%14.6f\n",t,U);// U(0)の出力
    fprintf(fp5,"t	Ek+U\n%14.6f	%14.6f\n",t,U);// Ek(0)+U(0)の出力
    fprintf(fp6,"t	T\n%14.6f	0\n",t);// Tの出力
#endif
*/
//前準備おわり

//時間計測開始
  double ttime;
  struct timeval tv1,tv2,tv3,tv4,tv5;
  gettimeofday(&tv1, NULL);


//二次元分割準備
    nb=(N-1)/bs+1;//二次元領域分割
    int counter,Ctmp,b1,b2,nblks1,nblks2;
    counter = ((nb+1)*nb)*0.5;

    int *COUNT;
    int *ibtbl,*jbtbl,*bs1,*be1,*bs2i,*be2i,*bs2j,*be2j;
    COUNT=(int*)malloc(sizeof(int)*(counter*6-nb*2));
    ibtbl=COUNT;
    jbtbl=COUNT+counter;
    bs1=COUNT+2*counter;
    be1=COUNT+2*counter+nb;
    bs2i=COUNT+2*counter+2*nb;
    be2i=COUNT+3*counter+nb;
    bs2j=COUNT+4*counter;
    be2j=COUNT+5*counter-nb;

    counter = 0;
    for (ib =0;ib < nb; ib++){
	for (jb=ib;jb <nb; jb++){
	    ibtbl[counter]=ib;
	    jbtbl[counter]=jb;
	    counter++;
	}
    }
    b1=0;b2=0;					
    //for(Ctmp =my_rank; Ctmp <(nb+1)*nb/2; Ctmp+=p){   //Ctmpはn+(n-1)...+1.counterまでする。
    for(Ctmp =0; Ctmp <(nb+1)*nb/2; Ctmp++){   //Ctmpはn+(n-1)...+1.counterまでする。
	if(ibtbl[Ctmp]==jbtbl[Ctmp]){
	    bs1[b1] = ibtbl[Ctmp]*bs;//b1の配列の数はi行のブロック数
	    be1[b1] = bs1[b1]+bs;
	    if(be1[b1] > N)
		be1[b1]=N;
	    b1++;
	}else{
	    bs2i[b2]=ibtbl[Ctmp]*bs;//b2の配列の数はcounter-b1
	    be2i[b2]=ibtbl[Ctmp]*bs+bs;
	    if (be2i[b2] > N){
		be2i[b2] = N;
	    }

	    bs2j[b2] = jbtbl[Ctmp]*bs;
	    be2j[b2] = jbtbl[Ctmp]*bs+bs;
	    if (be2j[b2] > N)
		be2j[b2] = N;
	    b2++;
	}
    }
    nblks1=b1;
    nblks2=b2;
//二次元分割準備終わり	


//3重ループ開始部分
    maxbox=100;
    minbox=maxbox*(-1);
    MAXSTEP=TMAX/d_t;

    gettimeofday(&tv2, NULL);
    for(step=1; step<=MAXSTEP;step++){
	t=(step-1)*d_t;
	for(i=0;i<N;i++){// 初期化
	    ax1[i] = 0.0;
	    ay1[i] = 0.0;
	    az1[i] = 0.0;
	}

	E = U = 0.0;

//速度ヴェルレ:位置
	for(i=0;i<N;i++){// r(t+Δt)
	    x1[i] = x0[i] + d_t*vx0[i] + ax0[i]*ddt;
	    y1[i] = y0[i] + d_t*vy0[i] + ay0[i]*ddt;
	    z1[i] = z0[i] + d_t*vz0[i] + az0[i]*ddt;
	}

//境界条件位置(t+dt)
	for(i=0;i<N;i++){
	    if(x1[i]>maxbox)
		x1[i]=maxbox*2-x1[i];
	    if(x1[i]<minbox)
		x1[i]=minbox*2-x1[i];
	    if(y1[i]>maxbox)
		y1[i]=maxbox*2-y1[i];
	    if(y1[i]<minbox)
		y1[i]=minbox*2-y1[i];
	    if(z1[i]>maxbox)
		z1[i]=maxbox*2-z1[i];
	    if(z1[i]<minbox)
		z1[i]=minbox*2-z1[i];
	}

//Lennard-Jones,二次元分割
	gettimeofday(&tv3, NULL);
	for (b1=0;b1 < nblks1; b1++){
	    iend=be1[b1];
	    for (i=bs1[b1];i <iend; i++)
		for(j= i+1;j <iend;j++){
		    x = x1[j]-x1[i];
		    y = y1[j]-y1[i];
		    z = z1[j]-z1[i];

		    rri = (x*x + y*y + z*z);
		    rri = 1.0/rri;
		    ri6 = rri*rri*rri;
		    coeff = b*ri6*rri*(a*ri6-1);

		    ax1[i] -= tmp = coeff*x;
		    ax1[j] += tmp;


		    ay1[i] -= tmp = coeff*y;
		    ay1[j] += tmp;

		    az1[i] -= tmp = coeff*z;
		    az1[j] += tmp;

		    U += c*ri6*(sigm6*ri6-1);// ポテンシャルU
		}
	}

	for (b2=0;b2 < nblks2; b2++){
	    jstart=bs2j[b2];
	    iend =be2i[b2];jend=be2j[b2];
	    for (i=bs2i[b2];i <iend; i++)
		for(j= jstart;j < jend;j++){	
		    x = x1[j]-x1[i];
		    y = y1[j]-y1[i];
		    z = z1[j]-z1[i];

		    rri = (x*x + y*y + z*z);
		    rri = 1.0/rri;
		    ri6 = rri*rri*rri;
		    coeff = b*ri6*rri*(a*ri6-1);

		    ax1[i] -= tmp = coeff*x;
		    ax1[j] += tmp;

		    ay1[i] -= tmp = coeff*y;
		    ay1[j] += tmp;

		    az1[i] -= tmp = coeff*z;
		    az1[j] += tmp;

		    U += c*ri6*(sigm6*ri6-1);// ポテンシャルU
		}
	}
	gettimeofday(&tv4, NULL);
	
//加速度,ポテンシャルの共有		
	t34 += tv4.tv_sec-tv3.tv_sec+(tv4.tv_usec - tv3.tv_usec)/1000000.0;

//境界条件速度・加速度
	for(i=0;i<N;i++){
	    if(x1[i]>maxbox||x1[i]<minbox){
		vx0[i]=-vx0[i];
		ax0[i]=-ax0[i];
	    }
	    if(y1[i]>maxbox||y1[i]<minbox){
		vy0[i]=-vy0[i];
		ay0[i]=-ay0[i];
	    }
	    if(z1[i]>maxbox||z1[i]<minbox){
		vz0[i]=-vz0[i];
		az0[i]=-az0[i];
	    }
	}

//速度ヴェルレ：速度
	for(i=0;i<N;i++){// v(t+Δt)
	    vx1[i] = vx0[i] + (ax0[i]+ax1[i])*dt2;
	    vy1[i] = vy0[i] + (ay0[i]+ay1[i])*dt2;
	    vz1[i] = vz0[i] + (az0[i]+az1[i])*dt2;
	}

	for(i=0;i<N;i++)// 運動エネルギーEkの計算
	    E += m2*(vx1[i]*vx1[i] + vy1[i]*vy1[i] + vz1[i]*vz1[i]);

	T = d*E;// 温度の計算

//ファイル出力
/*
#ifdef FILEOUT
	    if ( t_ave < TAVE ){
		for (i=0; i<N; ++i) {
		    xAve[i] += x1[i];
		    yAve[i] += y1[i];
		    zAve[i] += z1[i];
		}
		++t_ave;
	    } else {
		fprintf(fp,"%d\n",N);	//TAVEstep毎の平均を出力
		for (i=0; i<N; ++i) {
		    xAve[i] /= (double)TAVE;
		    yAve[i] /= (double)TAVE;
		    zAve[i] /= (double)TAVE;
		    fprintf(fp,"Ar_Ave %14.6f %14.6f %14.6f\n",xAve[i],yAve[i],zAve[i]);
		}
		fprintf(fp,"t=%14.6f\n",t+d_t);
		fprintf(fp,"\n");
		t_ave = 1;
	    }    
	    fflush(fp);

//各step毎の位置
	    for (i=0; i<N; ++i){
		fprintf(fp,"Ar %14.6f %14.6f %14.6f\n",x1[i],y1[i],z1[i]);
	    }
	    fprintf(fp,"t=%14.6f\n",t+d_t);
	    fprintf(fp,"\n");

//その他の出力		
	    fprintf(fp2,"t=%14.6f\nAr	x	y	z\n",t+d_t);//v(t+Δt)の出力
	    for (i=0; i<N; ++i) {		
		fprintf(fp2,"%d	%14.6f	%14.6f	%14.6f\n",i+1,vx1[i],vy1[i],vz1[i]);	
	    }
	    fprintf(fp2,"\n");
	    fflush(fp2);
	    fprintf(fp3,"%14.6f	%14.6f\n",t+d_t,E);// Ek(t+Δt)の出力
	    fflush(fp3);
	    fprintf(fp4,"%14.6f	%14.6f\n",t+d_t,U);// U(t+Δt)の出力
	    fflush(fp4);
	    fprintf(fp5,"%14.6f	%14.6f\n",t+d_t,E+U);// Ek(t+Δt)+U(t+Δt)の出力
	    fflush(fp5);
	    fprintf(fp6,"%14.6f	%14.6f\n",t+d_t,T);// T(t+Δt)の出力
	    fflush(fp6);
#endif
*/
//更新
	for(i=0;i<N;i++){// 更新
	    x0[i] = x1[i];
	    y0[i] = y1[i];
	    z0[i] = z1[i];
	    ax0[i] = ax1[i];
	    ay0[i] = ay1[i];
	    az0[i] = az1[i];
	    vx0[i] = vx1[i];
	    vy0[i] = vy1[i];
	    vz0[i] = vz1[i];
	}
    }	
    gettimeofday(&tv5, NULL);
    t16 = tv5.tv_sec-tv1.tv_sec+(tv5.tv_usec - tv1.tv_usec)/1000000.0;
    t26 = tv5.tv_sec-tv2.tv_sec+(tv5.tv_usec - tv2.tv_usec)/1000000.0;
    
//各プロセスのランク番号と実行時間等の出力
	printf("N=%d,TMAX=%f,bs=%d,cutoffL=%f,プロセス数=ⅰ\n全体の実行時間 3重loop実行時間 2重loop実行時間 集団通信時間 my_rank 計算量\n",N,TMAX,bs,sqrt(cutoffL));    

	    printf("%f	%f	%f	0	0	%Ld\n",t16,t26,t34,cutoffN);

//ループの回した回数とカットオフ率の表示
    loopN=(long int)MAXSTEP*N*(N-1)*0.5;
    cutoffR=(double)cutoffN/loopN;
    cutoffR=1.0-cutoffR;
	printf("Total number of loop=%Ld,cutoff number=%Ld,cutoff length=%14.6f,cutoff rate=%f \n\n",loopN,cutoffN,sqrt(cutoffL),cutoffR);

    return 0;
}


void ArConf (double x0[], double y0[], double z0[],int N){

        int i,k,jl,ng;
	double r;
	r = SIGM*pow(2.0,1.0/6.0);

  /****  Calcuration of whole particle number  ****/

	ng = 1;
	for(k=1;k<LAYER_NUMBER;k++){
		for(i=1,jl=0;i<=k+1;i++){
			jl+=i;
		}
		ng += 20*jl-30*k-18;
	}

  /****  Calcuration of Base Vector  ****/
	{
		double cos1[6],sin1[6],x1,y1,z1,x2,y2,z2;
		double  pi,x3,y3,z3,a,b,c,d,l,alpha,alphac,alphas,betac,betas;
		int m, mInit, mEnd, lnum;

		pi = 4.0*atan(1.0);
		d = 1.0/(sin(pi/5.0)*sin(pi/5.0));
		a = (d-2.0)*r/2.0;
		l = sqrt(4.0-d)*r;
		b = l/(2.0*sin(pi/5.0));
		c = b*cos(pi/5.0);

		x1 = 0.0;   y1 = 0.0;   z1 = r;
		x2 = b;     y2 = 0.0;   z2 = a;
		x3 = c;     y3 = l/2.0; z3 = -a;

    /**** MAIN LOOP ****/

    //xgBefor[0]=ygBefor[0]=zgBefor[0]=0.0;  //Core

		for(i=0;i<=5;i++){
			cos1[i] = cos(0.4*pi*i);
			sin1[i] = sin(0.4*pi*i);
		}

		for(m=mInit=lnum=1;(lnum<LAYER_NUMBER)&&(m<N);lnum++,mInit=m){
    	/**** configurate Northern Band ****/
			x0[m] = lnum*x1;
			y0[m] = lnum*y1;
			z0[m] = lnum*z1;
			m++;
			for(i=0;i<=4;i++){
				int n;
				for(n=1;n<=lnum;n++){
					alpha = lnum-n;
					for(k=1;(k<=n)&&(m<N);k++,m++){
						betac = k*cos1[i]+(n-k)*cos1[i+1];
						betas = k*sin1[i]+(n-k)*sin1[i+1];
						x0[m] = alpha*x1 + betac*x2 - betas*y2;
						y0[m] = alpha*y1 + betas*x2 + betac*y2;
						z0[m] = alpha*z1 + n*z2;
					}
				}
			}

      /**** configurate Equator Band ****/
			for(i=0;i<=4;i++){
				int n;
				for(n=1;n<=lnum-1;n++){
					for(k=1;(k<=n)&&(m<N);k++,m++){
						alphac = (lnum-n)*cos1[i]+(n-k)*cos1[i+1];
						alphas = (lnum-n)*sin1[i]+(n-k)*sin1[i+1];
						betac = k*cos1[i];
						betas = k*sin1[i];
						x0[m] = alphac*x2 - alphas*y2 + betac*x3 - betas*y3;
						y0[m] = alphas*x2 + alphac*y2 + betas*x3 + betac*y3;
						z0[m] = (lnum-k)*z2 + k*z3;
					}
				}
			}

      /**** configurate Southern Hemisphere ****/
			alphac = cos(pi/5.0);
			alphas = sin(pi/5.0);
			for(mEnd=m-1;(mEnd>=mInit)&&(m<N);mEnd--,m++){
				x0[m] = alphac*x0[mEnd] - alphas*y0[mEnd];
				y0[m] = alphas*x0[mEnd] + alphac*y0[mEnd];
				z0[m] = -z0[mEnd];
			}
		}
	}
}

