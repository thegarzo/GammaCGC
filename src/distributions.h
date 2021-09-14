#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <iostream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

double MV_DIP_J0(double x,void * p)
{
	double k = *(double *) p;
	double K1= gsl_sf_bessel_K1(x*LamQCD);
	double B2 = beta*(1-x*LamQCD*K1);
	return  (exp(-B2))*gsl_sf_bessel_J0(k*x)*x;
}

double MV_DIP_J0_0(double x,void * p)
{
	double k = *(double *) p;
	double polfactor = pow(Q2P*x*x/4, gammap);
	double logfactor = log( pow(x*LamQCD,-2) + exp(1) );
	double B2 = polfactor*logfactor;
	return  (exp(-B2))*gsl_sf_bessel_J0(k*x)*x;
}

class Distributions{
	
	public:
    Distributions(){
        
        A1=1;
        A2=1;
        
        int rep1 = 0;
        int rep2 = 0;
        
        char wf1[256];
        char wf2[256];
        char lgx[256];
        
        if (wfmode == 0){
            sprintf(wf1,"../wfC/ft_g1119_qs02_0168_af1_399_N%d.dat",A1);
            sprintf(wf2,"../wfC/ft_g1119_qs02_0168_af1_399_N%d.dat",A2);
            sprintf(lgx,"../wfC/LargeXG1119.dat");
        }
        if (wfmode == 1){
            sprintf(wf1,"../wfC/ft_mv_qs02_02_af1_N%d.dat",A1);
            sprintf(wf2,"../wfC/ft_mv_qs02_02_af1_N%d.dat",A2);
            sprintf(lgx,"../wfC/LargeX.dat");
        }
        if (wfmode == 2){
            QS21 = (double)(A1)*QS02;
            QS22 = (double)(A2)*QS02;
        }
        
        double temp;
        
        FILE *table1 = fopen(wf1,"r");
        FILE *table2 = fopen(wf2,"r");
        FILE *table3 = fopen(lgx,"r");
        
        double uGD1f[NY][Nkt+1];
        double uGD2f[NY][Nkt+1];
        double uGD1a[NY][Nkt+1];
        double uGD2a[NY][Nkt+1];
        double Y[NY];
        double kT[NY][Nkt+1];
        
        for (int iY=0; iY < NY; iY++)
        {
            for (int ikT=1; ikT < Nkt+1; ikT++)
            {
                
                if (fscanf(table1,"%lf %lf %lf %lf", &Y[iY], &kT[iY][ikT], &uGD1f[iY][ikT], &uGD1a[iY][ikT]) != 4)
                {
                    printf("Error reading table 1!\n");
                    exit(1);
                }
                

                if (fscanf(table2,"%lf %lf %lf %lf", &Y[iY], &kT[iY][ikT], &uGD2f[iY][ikT], &uGD2a[iY][ikT]) != 4)
                {
                    printf("Error reading table 2!\n");
                    exit(1);
                }
            }
            
            kT[iY][0] = 0;
            uGD1f[iY][0] = get_N_lin_interp(kT[iY][1], kT[iY][2], uGD1f[iY][1], uGD1f[iY][2], 0);
            uGD2f[iY][0] = get_N_lin_interp(kT[iY][1], kT[iY][2], uGD2f[iY][1], uGD2f[iY][2], 0);
            uGD1a[iY][0] = get_N_lin_interp(kT[iY][1], kT[iY][2], uGD1a[iY][1], uGD1a[iY][2], 0);
            uGD2a[iY][0] = get_N_lin_interp(kT[iY][1], kT[iY][2], uGD2a[iY][1], uGD2a[iY][2], 0);
            //initialize spline at each rapidity
            // rep = 0 => fundamental
            // rep = 1 => Adjoint
            
            if(rep1 == 0 )
            {
                accWF1[iY] = gsl_interp_accel_alloc();
                WF1[iY] = gsl_spline_alloc(gsl_interp_cspline, Nkt+1);
                gsl_spline_init (WF1[iY], &kT[iY][0], &uGD1f[iY][0], Nkt+1);
                exponent_rep[0] = 0.5;
            }
            else
            {
                accWF1[iY] = gsl_interp_accel_alloc();
                WF1[iY] = gsl_spline_alloc(gsl_interp_cspline, Nkt+1);
                gsl_spline_init (WF1[iY], &kT[iY][0], &uGD1a[iY][0], Nkt+1);
                exponent_rep[0] = 1;
            }
            
            if(rep2 == 0 )
            {
                accWF2[iY] = gsl_interp_accel_alloc();
                WF2[iY] = gsl_spline_alloc(gsl_interp_cspline, Nkt+1);
                gsl_spline_init (WF2[iY], &kT[iY][0], &uGD2f[iY][0], Nkt+1);
                exponent_rep[1] = 0.5;
            }
            else
            {
                accWF2[iY] = gsl_interp_accel_alloc();
                WF2[iY] = gsl_spline_alloc(gsl_interp_cspline, Nkt+1);
                gsl_spline_init (WF2[iY], &kT[iY][0], &uGD2a[iY][0], Nkt+1);
                exponent_rep[1] = 1;
            }
//            cout<< "Here, at iY = "<< iY<< "\n" ;

        }

        double xvals[100];
        double LxTable[100];
        
        for (int iLx = 0; iLx < 100; iLx++) {
            if (fscanf(table3,"%lf %lf", &xvals[iLx], &LxTable[iLx]) != 2)
            {
                printf("Error reading table Large X!\n");
                exit(1);
            }
        }
        
        accLargeX = gsl_interp_accel_alloc();
        LargeX = gsl_spline_alloc(gsl_interp_cspline, 100);
        gsl_spline_init (LargeX, &xvals[0], &LxTable[0], 100);
        
        fclose(table1);
        fclose(table2);
        fclose(table3);
        
        //Initialisation of internal variables.
        
        dY = Y[1]-Y[0];  // width of rapidity bin
    
        int index_t = Nkt+1-1;
        while (kT[0][index_t] > kT_max) index_t--;
        
//        DK = kT[0][1]-kT[0][0];
        
        i_max = index_t;
        kT1 = kT[0][index_t];
        kT2 = kT[0][index_t-1];
//        cout << DK <<"\n";
        
    }
    
	Distributions(int AA1, int rep1, int AA2, int rep2): A1(AA1), A2(AA2)
	{
		
		char wf1[256];
		char wf2[256];
		char lgx[256];
		
		if (wfmode == 0){
			sprintf(wf1,"../wfC/ft_g1119_qs02_0168_af1_399_N%d.dat",A1);
			sprintf(wf2,"../wfC/ft_g1119_qs02_0168_af1_399_N%d.dat",A2);
			sprintf(lgx,"../wfC/LargeX.dat");
		}
		if (wfmode == 1){
			sprintf(wf1,"../wfC/ft_mv_qs02_02_af1_N%d.dat",A1);
			sprintf(wf2,"../wfC/ft_mv_qs02_02_af1_N%d.dat",A2);
			sprintf(lgx,"../wfC/LargeX.dat");
		}
		if (wfmode == 2){
			QS21 = (double)(A1)*QS02;
			QS22 = (double)(A2)*QS02;
		}
		
		double temp;
		
		FILE *table1 = fopen(wf1,"r");
		FILE *table2 = fopen(wf2,"r");
		FILE *table3 = fopen(lgx,"r");
		
		double uGD1f[NY][Nkt+1];
		double uGD2f[NY][Nkt+1];
		double uGD1a[NY][Nkt+1];
		double uGD2a[NY][Nkt+1];
		double Y[NY];
		double kT[NY][Nkt+1];
		
		for (int iY=0; iY < NY; iY++)
		{
			for (int ikT=1; ikT < Nkt+1; ikT++)
			{
				
				if (fscanf(table1,"%lf %lf %lf %lf", &Y[iY], &kT[iY][ikT], &uGD1f[iY][ikT], &uGD1a[iY][ikT]) != 4)
				{
					printf("Error reading table 1!\n");
					exit(1);
				}
				

				if (fscanf(table2,"%lf %lf %lf %lf", &Y[iY], &kT[iY][ikT], &uGD2f[iY][ikT], &uGD2a[iY][ikT]) != 4)
				{
					printf("Error reading table 2!\n");
					exit(1);
				}
			}
			
			kT[iY][0] = 0;
			uGD1f[iY][0] = get_N_lin_interp(kT[iY][1], kT[iY][2], uGD1f[iY][1], uGD1f[iY][2], 0);
			uGD2f[iY][0] = get_N_lin_interp(kT[iY][1], kT[iY][2], uGD2f[iY][1], uGD2f[iY][2], 0);
			uGD1a[iY][0] = get_N_lin_interp(kT[iY][1], kT[iY][2], uGD1a[iY][1], uGD1a[iY][2], 0);
			uGD2a[iY][0] = get_N_lin_interp(kT[iY][1], kT[iY][2], uGD2a[iY][1], uGD2a[iY][2], 0);
			//initialize spline at each rapidity
			// rep = 0 => fundamental
			// rep = 1 => Adjoint
			
			if(rep1 == 0 )
			{
				accWF1[iY] = gsl_interp_accel_alloc();
				WF1[iY] = gsl_spline_alloc(gsl_interp_cspline, Nkt+1);
				gsl_spline_init (WF1[iY], &kT[iY][0], &uGD1f[iY][0], Nkt+1);
				exponent_rep[0] = 0.5;
			}
			else
			{
				accWF1[iY] = gsl_interp_accel_alloc();
				WF1[iY] = gsl_spline_alloc(gsl_interp_cspline, Nkt+1);
				gsl_spline_init (WF1[iY], &kT[iY][0], &uGD1a[iY][0], Nkt+1);
				exponent_rep[0] = 1;
			}
			
			if(rep2 == 0 )
			{
				accWF2[iY] = gsl_interp_accel_alloc();
				WF2[iY] = gsl_spline_alloc(gsl_interp_cspline, Nkt+1);
				gsl_spline_init (WF2[iY], &kT[iY][0], &uGD2f[iY][0], Nkt+1);
				exponent_rep[1] = 0.5;
			}
			else
			{
				accWF2[iY] = gsl_interp_accel_alloc();
				WF2[iY] = gsl_spline_alloc(gsl_interp_cspline, Nkt+1);
				gsl_spline_init (WF2[iY], &kT[iY][0], &uGD2a[iY][0], Nkt+1);
				exponent_rep[1] = 1;
			}
//			cout<< "Here, at iY = "<< iY<< "\n" ;

		}

		double xvals[100];
		double LxTable[100];
		
		for (int iLx = 0; iLx < 100; iLx++) {
			if (fscanf(table3,"%lf %lf", &xvals[iLx], &LxTable[iLx]) != 2)
			{
				printf("Error reading table Large X!\n");
				exit(1);
			}
		}
		
		accLargeX = gsl_interp_accel_alloc();
		LargeX = gsl_spline_alloc(gsl_interp_cspline, 100);
		gsl_spline_init (LargeX, &xvals[0], &LxTable[0], 100);
		
		fclose(table1);
		fclose(table2);
		fclose(table3);
		
		//Initialisation of internal variables.
		
		dY = Y[1]-Y[0];  // width of rapidity bin
	
		int index_t = Nkt+1-1;
		while (kT[0][index_t] > kT_max) index_t--;
		
//		DK = kT[0][1]-kT[0][0];
		
		i_max = index_t;
		kT1 = kT[0][index_t];
		kT2 = kT[0][index_t-1];
//		cout << DK <<"\n";

		
	}
	
	virtual ~Distributions() 
	{
		for(int ij =0; ij< NY; ij++ )
		{
			gsl_spline_free (WF1[ij]);
			gsl_spline_free (WF2[ij]);
			gsl_interp_accel_free (accWF1[ij]);
			gsl_interp_accel_free (accWF2[ij]);
		}

		gsl_spline_free (LargeX);
		gsl_interp_accel_free (accLargeX);

	}
    
    void ReDist(int AA1, int rep1, int AA2, int rep2)
    {
        
        A1=AA1; A2=AA2;
        
        char wf1[256];
        char wf2[256];
        char lgx[256];
        
        if (wfmode == 0){
            sprintf(wf1,"../wfC/ft_g1119_qs02_0168_af1_399_N%d.dat",A1);
            sprintf(wf2,"../wfC/ft_g1119_qs02_0168_af1_399_N%d.dat",A2);
            sprintf(lgx,"../wfC/LargeX.dat");
        }
        if (wfmode == 1){
            sprintf(wf1,"../wfC/ft_mv_qs02_02_af1_N%d.dat",A1);
            sprintf(wf2,"../wfC/ft_mv_qs02_02_af1_N%d.dat",A2);
            sprintf(lgx,"../wfC/LargeX.dat");
        }
        if (wfmode == 2){
            QS21 = (double)(A1)*QS02;
            QS22 = (double)(A2)*QS02;
        }
        
        double temp;
        
        FILE *table1 = fopen(wf1,"r");
        FILE *table2 = fopen(wf2,"r");
        FILE *table3 = fopen(lgx,"r");
        
        double uGD1f[NY][Nkt+1];
        double uGD2f[NY][Nkt+1];
        double uGD1a[NY][Nkt+1];
        double uGD2a[NY][Nkt+1];
        double Y[NY];
        double kT[NY][Nkt+1];
        
        for (int iY=0; iY < NY; iY++)
        {
            for (int ikT=1; ikT < Nkt+1; ikT++)
            {
                
                if (fscanf(table1,"%lf %lf %lf %lf", &Y[iY], &kT[iY][ikT], &uGD1f[iY][ikT], &uGD1a[iY][ikT]) != 4)
                {
                    printf("Error reading table 1!\n");
                    exit(1);
                }
                

                if (fscanf(table2,"%lf %lf %lf %lf", &Y[iY], &kT[iY][ikT], &uGD2f[iY][ikT], &uGD2a[iY][ikT]) != 4)
                {
                    printf("Error reading table 2!\n");
                    exit(1);
                }
            }
            
            kT[iY][0] = 0;
            uGD1f[iY][0] = get_N_lin_interp(kT[iY][1], kT[iY][2], uGD1f[iY][1], uGD1f[iY][2], 0);
            uGD2f[iY][0] = get_N_lin_interp(kT[iY][1], kT[iY][2], uGD2f[iY][1], uGD2f[iY][2], 0);
            uGD1a[iY][0] = get_N_lin_interp(kT[iY][1], kT[iY][2], uGD1a[iY][1], uGD1a[iY][2], 0);
            uGD2a[iY][0] = get_N_lin_interp(kT[iY][1], kT[iY][2], uGD2a[iY][1], uGD2a[iY][2], 0);
            //initialize spline at each rapidity
            // rep = 0 => fundamental
            // rep = 1 => Adjoint
            
            if(rep1 == 0 )
            {
                accWF1[iY] = gsl_interp_accel_alloc();
                WF1[iY] = gsl_spline_alloc(gsl_interp_cspline, Nkt+1);
                gsl_spline_init (WF1[iY], &kT[iY][0], &uGD1f[iY][0], Nkt+1);
                exponent_rep[0] = 0.5;
            }
            else
            {
                accWF1[iY] = gsl_interp_accel_alloc();
                WF1[iY] = gsl_spline_alloc(gsl_interp_cspline, Nkt+1);
                gsl_spline_init (WF1[iY], &kT[iY][0], &uGD1a[iY][0], Nkt+1);
                exponent_rep[0] = 1;
            }
            
            if(rep2 == 0 )
            {
                accWF2[iY] = gsl_interp_accel_alloc();
                WF2[iY] = gsl_spline_alloc(gsl_interp_cspline, Nkt+1);
                gsl_spline_init (WF2[iY], &kT[iY][0], &uGD2f[iY][0], Nkt+1);
                exponent_rep[1] = 0.5;
            }
            else
            {
                accWF2[iY] = gsl_interp_accel_alloc();
                WF2[iY] = gsl_spline_alloc(gsl_interp_cspline, Nkt+1);
                gsl_spline_init (WF2[iY], &kT[iY][0], &uGD2a[iY][0], Nkt+1);
                exponent_rep[1] = 1;
            }
//            cout<< "Here, at iY = "<< iY<< "\n" ;

        }

        double xvals[100];
        double LxTable[100];
        
        for (int iLx = 0; iLx < 100; iLx++) {
            if (fscanf(table3,"%lf %lf", &xvals[iLx], &LxTable[iLx]) != 2)
            {
                printf("Error reading table Large X!\n");
                exit(1);
            }
        }
        
        accLargeX = gsl_interp_accel_alloc();
        LargeX = gsl_spline_alloc(gsl_interp_cspline, 100);
        gsl_spline_init (LargeX, &xvals[0], &LxTable[0], 100);
        
        fclose(table1);
        fclose(table2);
        fclose(table3);
        
        //Initialisation of internal variables.
        
        dY = Y[1]-Y[0];  // width of rapidity bin
    
        int index_t = Nkt+1-1;
        while (kT[0][index_t] > kT_max) index_t--;
        
//        DK = kT[0][1]-kT[0][0];
        
        i_max = index_t;
        kT1 = kT[0][index_t];
        kT2 = kT[0][index_t-1];
//        cout << DK <<"\n";

        
    }
    
    int getNT(){return A2;}

	////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////
	
	
	double N_K_X(int nucleus, double xx, double kk)
	{
		#if(FIXDX==0)
		double x = xx;
		#endif
		#if(FIXDX==1)
		double x = fixedx;
		#endif
		
		double k = kk;
	
		if (wfmode == 2)
		{
			return wfsimple(nucleus, x, k);
		}

		if ( x <= XMIN || x >= 1.0 )
		{
			return 0.0;
		}
		
		if ( k < 0 )
		{
			return 0;
		}
		if (k > kT_max)
		{
			//get index corresponding to max kT
//			double b = log( N_K_X(nucleus, x, kT1 )/N_K_X(nucleus, x, kT2 ) ) / log( kT1/kT2 );
//			double a = N_K_X(nucleus, x, kT1)/pow(kT1,b);
//			return a*pow(k,b);
			return 0;

		}

//		double norm = NC/4.*LamQCD*LamQCD/alpha(k) ;
//		double norm = k*k;
		double norm = 1;
		double myY = gsl_max(0.0,log(x0/x));
		int iY = (int) (myY/dY);    // rapidity bin

		if (iY+1 < NY)
		{
			double tmp,tmp1,tmp2;
			
			if (nucleus == 1)
			{
				tmp1 = gsl_spline_eval(WF1[iY], k, accWF1[iY]);
				tmp2 = gsl_spline_eval(WF1[iY+1], k, accWF1[iY+1]);
			}
			else if (nucleus == 2)
			{
				tmp1 = gsl_spline_eval(WF2[iY], k, accWF2[iY]);
				tmp2 = gsl_spline_eval(WF2[iY+1], k, accWF2[iY+1]);
			} else
			{
				tmp=0; tmp1=0; tmp2=0;
				printf("Error: invalid Nucleus\n");
				return 0.;
			}
	
			tmp = tmp1 + (myY-dY*iY)/dY * (tmp2-tmp1);  // lin interpolation in Y
			/* large-mak correction */
			return norm*tmp*(x > x0 ? pow(abs(gsl_spline_eval(LargeX, x, accLargeX)), exponent_rep[nucleus-1]) : 1.0);
			//return norm*tmp*(x > x0 ? pow( (1. - x)/(1. - x0), 4.)  : 1.0);
		} else
		{
			//printf("x too small in Phi(x)!\n");
			return 0.;
		}

		return 0.;
	}
	
	int PrintWF(double x, const char *out)
	{
		char fname[256];
		sprintf(fname,"%s_wf_x_%0.5f.dat",out,x);
		printf("Writing wf to file %s\n",fname);
		FILE *wfout = fopen(fname,"w");
		
		for(double k = 0; k <= 4.*kT_max ; k += 0.01)
		{
			fprintf(wfout,"%10.3e\t%10.5e\t%10.5e\n",k, N_K_X(1,x,k), N_K_X(2,x,k) );
		}
		
		
		fclose(wfout);
		return 0;
	}
    double IMatch(int nucleus, double xx, double Qi,int Nk){
        double dk = Qi/(Nk-1);
        double sum = 0;
        
        for(int ik=1;ik<Nk-1; ik++){
            sum += pow(ik*dk,3) * N_K_X(nucleus,xx,ik*dk);
        }
        
        sum += pow(Qi,3) * N_K_X(nucleus,xx,Qi)/2.;
        
        return 2 * dk * sum;
    }
    
    double dIMatchdx(int nucleus, double xx,double dx, double Qi,int Nk){

        return ( IMatch(nucleus,xx,Qi,Nk)-IMatch(nucleus,xx-dx,Qi,Nk) )/dx  ;
    }
    
	
private:

	gsl_interp_accel *accWF1[NY];
	gsl_spline *WF1[NY];
	gsl_interp_accel *accWF2[NY];
	gsl_spline *WF2[NY];
	gsl_interp_accel *accLargeX;
	gsl_spline *LargeX;

	
	
	int NM;
	int NP;

	int A1;
	int A2;
	
	double dY; // width of rapidity bin
	double DK; // width of rapidity bin
	int i_max;
	double kT1;
	double kT2;

	
	double QS21;
	double QS22;
	
	double irR;
	double uvR;
	double xi, yi;
	
	double exponent_rep[2];
	
	
	double get_N_lin_interp(double x0, double x1, double y0, double y1, double x)
	{
		return y0 + (x-x0)*(y1-y0)/(x1-x0);
	}
	
    
	
	double wfsimple(int nucleus, double x, double k)
	{
		double Q2;
		if (nucleus == 1)	Q2 = QS21;
		if (nucleus == 2)	Q2 = QS22;
		
		return 8.*M_PI*Q2/(k*k)*( 1.0 - exp(-k*k*k*k/Q2/Q2/4.) );
	}
	
};

#endif
