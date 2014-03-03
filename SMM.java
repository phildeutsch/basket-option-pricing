/*
 *	Calculates the price of a basket option
 * 	Philipp Deutsch, philipp.g.deutsch@gmail.com
 * 	April 2012
 */


import java.util.Random;

import Jama.CholeskyDecomposition;
import Jama.Matrix;

public class SMM {

	private static double[] S;      // fwd swap rates
	private static double[] S0;     // initial fwd swap rates
	private static double[] s;  	// implied swaption volas
	private static double[] B;    	// zero bonds
	private static double[] B0;    	// initial zero bonds
	private static double[] C;		// value of the annuities
	private static double[] D;		// discount factor
	private static double[] m;		// drift of the fwd swap rates (0...N-1)
	private static double[] notio;	// notional of swap
	private static double rho;		// correlation between fwd swap rates
	private static double[][] COR;	// Correlation Matrix
	private static double K;		// Strike of the swaption
	private static int M;       	// time to maturity
	private static int X;           // time to expiry
	
	public static double[] price;	// prices of european swaptions
	
	public SMM(int time_to_mat, int time_to_exp, double cor) {
		X = time_to_exp;
		M = time_to_mat - X;
		
		S		= new double[M];
		S0		= new double[M];
		s   	= new double[M];
		m		= new double[M];
		B		= new double[M+1];
		B0		= new double[M+1];
		C		= new double[M];
		D		= new double[M+1];
		notio 	= new double[M];
		price	= new double[M];
		COR	 	= new double[M][M];
		
		
		for (int i=0; i<M; i++) s[i] = 0.4;
		for (int i=0; i<M+1; i++) B[i] = 0.99 - 0.01*i;
		for (int i=0; i<M+1; i++) B0[i] = B[i];
		
		double annu = B[M];
		for (int i=M-1; i>=0; i--) {
			S[i] = (B[i] - B[M])/annu;
			C[i] = annu;
			price[i] = 0;
			S0[i] = S[i];
			annu += B[i];
		}
		
		for (int i=0; i<M; i++) notio[i] = 10000.;
		
		K = 0.01;
		rho = cor;
		
		for(int i=0; i<M; i++){
			for (int j=0; j<M; j++){
				if (i == j) COR[i][j] = 1;
				else COR[i][j] = rho;
			}
		}
		
		/*
		Matrix mat = new Matrix(COR);
		mat.print(3, 3);
		*/
	}
	
	public static void calcBonds(int t) {
		int n = M-t;
		double[] sol = new double[n];
		double[][] A = new double[n][n];
		for (int i=0; i<n; i++){
			sol[i] = 0;
			A[i][0] = S[M-i-1] + 1;
			for (int j=i+2; j<n; j++) A[i][j] = 0;
			if (i+1<n) A[i][i+1] = -1;
			for (int j=1; j<i+1; j++) A[i][j] = S[M-1-1];
		}
		sol[n-1] = 1.;
		Matrix matA = new Matrix(A);
		Matrix matsol = new Matrix(sol,n);
		Matrix b = matA.solve(matsol);
		
		/*
		matA.print(3, 3);
		matsol.print(3, 3);
		b.print(3, 3);
		*/
		
		C[M-1] = B[M];
		for (int i=0; i<M+1; i++){
			if (i<=t) B[i] = 1;
			else B[i] = b.get(n-(i-t), 0);
		}
		for (int i=M-2; i>t; i--) C[i] = C[i+1] + B[i+1];
		D[t] = B0[M]/B[M];
	}
	public static void calcDrifts(int t){
		double prod = 1.;
		for (int n=t; n<M; n++){
			m[n] = 0;
			prod = 1.;
			for (int i=n; i<M-1; i++){
				if (i>n+1) {
					for (int j=n+1; j<=i; j++) prod *= 1 + S[j];
				}
				if (i+1 == n) m[n] += prod * C[i+1]/C[M-1] * S[n] * s[n] * S[i+1] * s[i+1];
				else m[n] += prod * C[i+1]/C[M-1] * rho * S[n] * s[n] * S[i+1] * s[i+1];
			}
		}
	}
	
	public static void simSwap(int t) {
		Random gen = new Random();
		double[] W = new double[M];
		int runs = 1;
		double SR = 0;
		
		for (int i=0; i<M; i++){
			W[i] = gen.nextGaussian();
		}
		Matrix matW = new Matrix(W,1);
		Matrix matC = new Matrix(COR);
	    CholeskyDecomposition matChol = new CholeskyDecomposition(matC);
	    Matrix matL = matChol.getL();
	    Matrix matWcor = matL.times(matW.transpose());
	    
	    double expo = 0;
	    for (int i=t; i<M; i++){
	    	SR = 0;
	    	for (int k=0; k<runs; k++){
	    		expo  = s[i] * matWcor.get(i,0) + (m[i] - 0.5 * s[i] * s[i]);
	    		SR   += S[i] * Math.exp(expo);
	    	}
	    	S[i] = SR / runs;
	    }
	}
	
	public static void calcEuropeanPVs() {
		int paths = 10000;
		double PV = 0;
		for( int n=0; n<M; n++){
			PV = 0;
			for (int p=0; p<paths; p++){
				for (int i=0; i<M; i++) {
					S[i] = S0[i];
					B[i] = B0[i];
				}
				B[M] = B0[M];
				for (int i=0; i<n+1; i++) {
					calcDrifts(i);
					simSwap(i);
					calcBonds(i);
//					for (int j=0; j<M; j++) System.out.format("%.10f ",m[j]);
//					System.out.println();
				}
				PV  += Math.max(S[n]-K,0) * C[n] * D[n] * notio[n];
			}
			price[n] = PV/paths;
		}
	}

	public static void main(String[] argv) {
		int timeToMat = 3;
		int timeToExp = 1;
		
		for (int cor=1; cor<100; cor++){
			new SMM(timeToMat, timeToExp, cor/100.);
			calcEuropeanPVs();
			System.out.format("%.2f ",cor/100.);
			for (int i=0; i<M; i++) {
				System.out.format("%.10f ",price[i]);
			}
			System.out.println();
		}
	}
	
}
