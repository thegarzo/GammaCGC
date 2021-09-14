//
//  pdfs.h
//  
//  Created by Oscar Garcia-Montero on 18.10.17.
//

#ifndef pdfs_h
#define pdfs_h

#include <iostream>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

class PDFs{
	
	public:
	
    PDFs(){
        iP=0;
        char pdfname[256];
        sprintf(pdfname,"../pdf_files/cteq6m_gluon.dat");
            
        double temp;
        
        FILE *table = fopen(pdfname,"r");
        
        double pdfA[NQQ][NXX];
        double XA[NQQ][NXX];
        double QA[NQQ];
        //CHANGE SO THAT X GETS CUBIC SPLINE AND Q GETS THE SHITTY LINEAR ONE. THIS IS NECESSARY FOR DIFFERENTIATION
        //CHANGE FORTRAN CODE TO GET TABS!
    
        for (int iQ=0; iQ < NQQ; iQ++)
        {
            for (int iX=0; iX < NXX; iX++)
            {
                if (fscanf(table,"%lf %lf %lf", &QA[iQ], &XA[iQ][iX], &pdfA[iQ][iX]) != 3)
                {
                    printf("Error reading distribution table!\n");
                    exit(1);
                }
//                cout << QA[iQ] << "\t"<< XA[iQ][iX] << "\n";
            }
            //initialize spline at each scale Q
            
            accPDF[iQ] = gsl_interp_accel_alloc();
            PDF[iQ] = gsl_spline_alloc(gsl_interp_cspline, NXX);
            gsl_spline_init (PDF[iQ], &XA[iQ][0], &pdfA[iQ][0], NXX);
        }
        fclose(table);
//        cout<< "gluon distribution loaded!\n";
    }
        
	PDFs(int iparton)
	{
        iP=iparton;
		char pdfname[256];
		
		switch (iparton) {
				case -5:
					sprintf(pdfname,"../pdf_files/cteq6m_bbar.dat");
					break;
				case -4:
					sprintf(pdfname,"../pdf_files/cteq6m_cbar.dat");
					break;
				case -3:
					sprintf(pdfname,"../pdf_files/cteq6m_sbar.dat");
					break;
				case -2:
					sprintf(pdfname,"../pdf_files/cteq6m_dbar.dat");
					break;
				case -1:
					sprintf(pdfname,"../pdf_files/cteq6m_ubar.dat");
					break;
				case 0:
					sprintf(pdfname,"../pdf_files/cteq6m_gluon.dat");
					break;
				case 1:
					sprintf(pdfname,"../pdf_files/cteq6m_u.dat");
					break;
				case 2:
					sprintf(pdfname,"../pdf_files/cteq6m_d.dat");
					break;
				case 3:
					sprintf(pdfname,"../pdf_files/cteq6m_s.dat");
					break;
				case 4:
					sprintf(pdfname,"../pdf_files/cteq6m_c.dat");
					break;
				case 5:
					sprintf(pdfname,"../pdf_files/cteq6m_b.dat");
					break;
				default:
					cout << "Parton identifier, "<< iparton <<" is not valid\n";
					exit(1);
					break;
		}
		
		double temp;
		
		FILE *table = fopen(pdfname,"r");
		
		double pdfA[NQQ][NXX];
		double XA[NQQ][NXX];
		double QA[NQQ];
		//CHANGE SO THAT X GETS CUBIC SPLINE AND Q GETS THE SHITTY LINEAR ONE. THIS IS NECESSARY FOR DIFFERENTIATION
		//CHANGE FORTRAN CODE TO GET TABS!
	
		for (int iQ=0; iQ < NQQ; iQ++)
		{
			for (int iX=0; iX < NXX; iX++)
			{
				if (fscanf(table,"%lf %lf %lf", &QA[iQ], &XA[iQ][iX], &pdfA[iQ][iX]) != 3)
				{
					printf("Error reading distribution table!\n");
					exit(1);
				}
//				cout << QA[iQ] << "\t"<< XA[iQ][iX] << "\n";
			}
			//initialize spline at each scale Q
			
			accPDF[iQ] = gsl_interp_accel_alloc();
			PDF[iQ] = gsl_spline_alloc(gsl_interp_cspline, NXX);
			gsl_spline_init (PDF[iQ], &XA[iQ][0], &pdfA[iQ][0], NXX);
		}
		fclose(table);
		
		
		switch (iparton) {
				case -5:
				cout<< "anti b quark distribution loaded!\n";
				break;
				case -4:
				cout<< "anti c quark distribution loaded!\n";
				break;
				case -3:
				cout<< "anti s quark distribution loaded!\n";
				break;
				case -2:
				cout<< "anti d quark distribution loaded!\n";
				break;
				case -1:
				cout<< "anti u quark distribution loaded!\n";
				break;
				case 0:
				cout<< "gluon distribution loaded!\n";
				break;
				case 1:
				cout<< "u quark distribution loaded!\n";
				break;
				case 2:
				cout<< "d quark distribution loaded!\n";
				break;
				case 3:
				cout<< "s quark distribution loaded!\n";
				break;
				case 4:
				cout<< "c quark distribution loaded!\n";
				break;
				case 5:
				cout<< "b quark distribution loaded!\n";
				break;
		}
		
		
	}
    
    void ReDist(int iparton)
    {
        
        iP=iparton;
        char pdfname[256];
        
        switch (iparton) {
                case -5:
                    sprintf(pdfname,"../pdf_files/cteq6m_bbar.dat");
                    break;
                case -4:
                    sprintf(pdfname,"../pdf_files/cteq6m_cbar.dat");
                    break;
                case -3:
                    sprintf(pdfname,"../pdf_files/cteq6m_sbar.dat");
                    break;
                case -2:
                    sprintf(pdfname,"../pdf_files/cteq6m_dbar.dat");
                    break;
                case -1:
                    sprintf(pdfname,"../pdf_files/cteq6m_ubar.dat");
                    break;
                case 0:
                    sprintf(pdfname,"../pdf_files/cteq6m_gluon.dat");
                    break;
                case 1:
                    sprintf(pdfname,"../pdf_files/cteq6m_u.dat");
                    break;
                case 2:
                    sprintf(pdfname,"../pdf_files/cteq6m_d.dat");
                    break;
                case 3:
                    sprintf(pdfname,"../pdf_files/cteq6m_s.dat");
                    break;
                case 4:
                    sprintf(pdfname,"../pdf_files/cteq6m_c.dat");
                    break;
                case 5:
                    sprintf(pdfname,"../pdf_files/cteq6m_b.dat");
                    break;
                default:
                    cout << "Parton identifier, "<< iparton <<" is not valid\n";
                    exit(1);
                    break;
        }
        
        double temp;
        
        FILE *table = fopen(pdfname,"r");
        
        double pdfA[NQQ][NXX];
        double XA[NQQ][NXX];
        double QA[NQQ];
        //CHANGE SO THAT X GETS CUBIC SPLINE AND Q GETS THE SHITTY LINEAR ONE. THIS IS NECESSARY FOR DIFFERENTIATION
        //CHANGE FORTRAN CODE TO GET TABS!
    
        for (int iQ=0; iQ < NQQ; iQ++)
        {
            for (int iX=0; iX < NXX; iX++)
            {
                if (fscanf(table,"%lf %lf %lf", &QA[iQ], &XA[iQ][iX], &pdfA[iQ][iX]) != 3)
                {
                    printf("Error reading distribution table!\n");
                    exit(1);
                }
//                cout << QA[iQ] << "\t"<< XA[iQ][iX] << "\n";
            }
            //initialize spline at each scale Q
            
            accPDF[iQ] = gsl_interp_accel_alloc();
            PDF[iQ] = gsl_spline_alloc(gsl_interp_cspline, NXX);
            gsl_spline_init (PDF[iQ], &XA[iQ][0], &pdfA[iQ][0], NXX);
        }
        fclose(table);
        
        
        switch (iparton) {
                case -5:
                cout<< "anti b quark distribution loaded!\n";
                break;
                case -4:
                cout<< "anti c quark distribution loaded!\n";
                break;
                case -3:
                cout<< "anti s quark distribution loaded!\n";
                break;
                case -2:
                cout<< "anti d quark distribution loaded!\n";
                break;
                case -1:
                cout<< "anti u quark distribution loaded!\n";
                break;
                case 0:
                cout<< "gluon distribution loaded!\n";
                break;
                case 1:
                cout<< "u quark distribution loaded!\n";
                break;
                case 2:
                cout<< "d quark distribution loaded!\n";
                break;
                case 3:
                cout<< "s quark distribution loaded!\n";
                break;
                case 4:
                cout<< "c quark distribution loaded!\n";
                break;
                case 5:
                cout<< "b quark distribution loaded!\n";
                break;
        }
        
        
    }
    
    
	virtual ~PDFs()
	{
		for(int ij =0; ij< NY; ij++ )
		{
			gsl_spline_free (PDF[ij]);
			gsl_interp_accel_free (accPDF[ij]);
		}
	}
	
	// --------------------------------------- EVALUATION OF PDFS --------------------------------------- //
	
	double fpart(double xx, double qq)
	{
		double fxQ;
		
		if ( xx <= XpdfMIN || xx >= 1.0 )
		{
			return 0.0;
		}
		
		if ( qq < QpdfMIN ||  qq > QpdfMAX )
		{
			return 0.0;
		}
		
		//		double norm = NC/4.*LamQCD*LamQCD/alpha(k) ;
		//		double norm = k*k;
		
		int jQ;
		double QQ1,QQ2;
		double tmp,tmp1,tmp2;
		
		jQ = (int) ( log10(qq/QpdfMIN)/dQ );
		QQ1= QpdfMIN * pow(10, jQ * dQ);
		QQ2= QpdfMIN * pow(10, (jQ+1) * dQ);
		
		tmp1 = gsl_spline_eval(PDF[jQ], xx, accPDF[jQ]);
		tmp2 = gsl_spline_eval(PDF[jQ+1], xx, accPDF[jQ+1]);
		
		tmp = tmp1 + (tmp2-tmp1) * (qq-QQ1)/(QQ2-QQ1) ;  // lin interpolation in Q

		return tmp;
		
	}
	
    int getiP(){return iP;}
	

//	
	double xfpart(double xx, double QQ)
	{
		return xx * fpart(xx,QQ);
	}
	
    double dxfpartdx(double xx, double QQ, double dx)
    {
        return  ( (xx+dx) * fpart(xx+dx,QQ)  - (xx) * fpart(xx,QQ) )/dx ;
    }
    
	private:
    
    int iP;
	gsl_interp_accel *accPDF[NQQ];
	gsl_spline *PDF[NQQ];
	
};

#endif /* pdfs_h */
