import Jama.Matrix;
import Jama.CholeskyDecomposition;
import java.util.Random;

class basket {
	static double rho;			// Correlation
	static double T;			// Time to expiry
	static double K;			// Strike
	static double[][] C;		// Correlation Matrix
	static double[][] Cov;		// Covariance matrix (lognormal vols)
	static double[][] Covn;		// Covariance matrix (normal vols)
	static double[] s;			// Lognormal volatilities
	static double[] sn;			// Equivalent normal volatilities
	static double[] w;			// Weights
	static double[] A0;			// Starting values
	static double[] A;			// Asset values
	static int N;				// Number of assets
	static boolean neg;			// Basket has a negative weight 
	
	public basket(int num, double timeToExp, double cor, double weight){
		rho  = cor;
		N    = num;
		T    = timeToExp;
		neg	 = false;
		A    = new double[N];
		A0   = new double[N];
		s    = new double[N];
		sn   = new double[N];
		w    = new double[N];
		C	 = new double[N][N];
		Cov  = new double[N][N];
		Covn = new double[N][N];
		
		A0[N-1]= 100;
		s[N-1] = 0.4;
		w[N-1] = weight;
		if (w[N-1] < 0) neg = true;
		double B0 = w[N-1] * A0[N-1];
		for (int i=0; i<N-1; i++){
			A0[i]= A0[N-1];
			s[i] = 0.4;
			w[i] = (1-w[N-1])/(N-1);
			B0 += w[i] * A0[i];
		}
		double calibration = 0.9;
		K    = B0;
		
		double Ki = 0;
		for (int i=0; i<N; i++){
			Ki = calibration * A0[i];
			if (Math.abs((Ki-A0[i])/Ki)<0.001) {
				sn[i]  = s[i] * Math.sqrt(A0[i] * K) / A0[i];
				sn[i] *= (1 + Math.log(A0[i]/K)/24.)/(1+s[i]*s[i]*T/24.);
			}
			else {
				sn[i]  = s[i] * (A0[i] - Ki)/Math.log(A0[i]/Ki) / A0[i];
				sn[i] /= 1 +  s[i] * s[i] * T / 24.;
			}
			/*
			System.out.format("Option ");
			System.out.format("%d:", i+1);
			System.out.format("%n");
			System.out.println(calcBachelierPrice(i, Ki));
			System.out.println(calcBSMPrice(i, Ki));
			System.out.format("%n");
			*/
		}
		
		for(int i=0; i<N; i++){
			for (int j=0; j<N; j++){
				if (i == j) C[i][j] = 1;
				else C[i][j] = rho;
				Cov[i][j] = s[i] * s[j] * C[i][j];
				Covn[i][j] = sn[i] * sn[j] * C[i][j];
			}
		}
		/*
		Matrix matC = new Matrix(C);
		matC.print(3, 3);
		*/
	}
	public double calcBachelierPrice(int i, double Ki){
		double d = (A0[i] - Ki)/(A0[i] * sn[i] * Math.sqrt(T));
		
		double N1 = 0.5 * (1 + StatUtil.erf(d/Math.sqrt(2)));
		double N2 = 1/Math.sqrt(2 * Math.PI) * Math.exp(- 0.5 * d * d);
				
		return (A0[i] - Ki) * N1 + A0[i] * sn[i] * Math.sqrt(T) * N2;
	}
	public double calcBSMPrice(int i, double Ki){
		double d1 = (Math.log(A0[i]/Ki) + 0.5 * s[i] * s[i] * T)/(s[i] * Math.sqrt(T));
		double d2 = d1 - s[i] * Math.sqrt(T);
		double N1 = 0.5 * (1 + StatUtil.erf(d1/Math.sqrt(2)));
		double N2 = 0.5 * (1 + StatUtil.erf(d2/Math.sqrt(2)));

		return A0[i] * N1 - Ki * N2;
	}
	
	public static double Deng(){
		double[][] rho 	= new double[2][2];
		double[] mu 	= new double[N];
		double[] muH	= new double[2];
		double[] nu 	= new double[N];
		double[] nuH	= new double[2];
		double m  		= 0;
		double m0  		= 0;
		double v		= 0;
		double r		= 0;
		double vK  		= 0;
		double dK		= 0;
		double M		= N-1;
		
		double H0 	= 0;
		double H1 	= -w[N-1] * A0[N-1];
		for (int i=0; i<N-1; i++){
			H0 += w[i] * A0[i];
		}
		for (int i=0; i<N; i++){
			mu[i] = Math.log(Math.abs(w[i]) * A0[i]) - 0.5 * s[i] * s[i] * T;
			nu[i] = s[i] * Math.sqrt(T);
		}
		muH[0] = 0;
		nuH[0] = 0;
		for (int i=0; i<N-1; i++){
			muH[0] += Math.exp(mu[i] + 0.5 * nu[i] * nu[i]);
			for (int j=0; j<N-1; j++){
				nuH[0] += C[i][j] * nu[i] * nu[j];
			}
		}
		nuH[0] = Math.sqrt(nuH[0]) / M;
		muH[0] = Math.log(muH[0]) - 0.5 * nuH[0] * nuH[0];
		muH[1] = mu[N-1];
		nuH[1] = nu[N-1];
		
		rho[0][0] = 1;
		rho[1][0] = 0;
		for (int i=0; i<N-1; i++){
			rho[1][0] += C[i][N-1] * nu[1];
		}
		rho[1][0] /= (M * nuH[0]);
		rho[0][1]  = rho[1][0];
		rho[1][1]  = C[N-1][N-1];
		/*
		Matrix matr = new Matrix(rho);
		matr.print(3, 3);
		*/
		m   = Math.exp(muH[1] + 0.5 * nuH[1] * nuH[1]);
		m0  = Math.exp(muH[0] + 0.5 * nuH[0] * nuH[0]) / (m + K);
		m   = m / (m + K);
		v   = nuH[1];
		r	= rho[0][1];
		
		vK 	= Math.sqrt(nuH[0] * nuH[0] - 2 * r * nuH[0] * v * m + v * v * m * m);
		dK	= Math.log(m0) / vK;
		
		double N1 = 0.5 * (1 + StatUtil.erf((dK+vK/2)/Math.sqrt(2)));
		double N2 = 0.5 * (1 + StatUtil.erf((dK-vK/2)/Math.sqrt(2)));
				
		return H0 * N1 - (H1 + K) * N2;
	}
	
	public static double Ju(){
		double[] S 	= new double[N];
		double[][] CT = new double[N][N];
		double levy = 0;
		double ju 	= 0;
		
		double U1 		= 0;
		double U20 		= 0;
		double U21 		= 0;
		double dU20 	= 0;
		double ddU20 	= 0;
		double dddU20 	= 0;
		double a1 		= 0;
		double a2 		= 0;
		double a3 		= 0;
		double b1 		= 0;
		double b1E 		= 0;
		double b2 		= 0;
		double c1 		= 0;
		double c2 		= 0;
		double c2E1 	= 0;
		double c2E2 	= 0;
		double c3 		= 0;
		double c3E1 	= 0;
		double c3E2 	= 0;
		double c4 		= 0;
//		double d1 		= 0;
		double d2 		= 0;
		double d3 		= 0;
		double d4 		= 0;
		double m1 		= 0;
		double v1 		= 0;
		double y 		= 0;
		double y1 		= 0;
		double y2 		= 0;
		double z1 		= 0;
		double z2 		= 0;
		double z3 		= 0;
		
		for(int i=0; i<N; i++){
			for (int j=0; j<N; j++){
				CT[i][j] = Cov[i][j] * T;
			}
		}
		
		for (int i=0; i<N; i++) {
			S[i] = w[i] * A0[i];
		}
			
		for (int i=0; i<N; i++){
			U1 += S[i];
			for (int j=0; j<N; j++){
				U20    += S[i] * S[j];
				U21    += S[i] * S[j] * Math.exp(CT[i][j]);
				dU20   += S[i] * S[j] * CT[i][j];
				ddU20  += S[i] * S[j] * CT[i][j]* CT[i][j];
				dddU20 += S[i] * S[j] * CT[i][j]* CT[i][j]* CT[i][j];
				for (int k=0; k<N; k++){
					b1E  += S[i] * S[j] * S[k] * CT[i][k] * CT[j][k];
					c3E1 += S[i] * S[j] * S[k] * CT[i][k] * CT[j][k] * CT[j][k];
					c3E2 += S[i] * S[j] * S[k] * CT[i][j] * CT[i][k] * CT[j][k];
					for (int l=0; l<N; l++){
						c2E1 += S[i] * S[j] * S[k] * S[l] * CT[i][l] * CT[j][k] * CT[k][l];
						c2E2 += S[i] * S[j] * S[k] * S[l] * CT[i][l] * CT[j][l] * CT[k][l];
					}
				}
			}
		}
		
		b1E  = 2 * b1E; 
		c2E1 = 8 * c2E1 + 2 * dU20 * ddU20;
		c2E2 = 6 * c2E2;
		c3E1 = 6 * c3E1;
		c3E2 = 8 * c3E2;
		
		a1 = - dU20 / (2 * U20);
		a2 = 2 * a1 * a1 - ddU20 / (2 * U20);
		a3 = 6 * a1 * a2 - 4 * a1 *a1 *a1 - dddU20 / (2 * U20);
		
		b1 = b1E / (4 * U1 * U1 * U1);
		b2 = a1 * a1 - 0.5 * a2;
		
		c1 = - a1 * b1;
		c2 = (9 * c2E1 + 4 * c2E2) / (144 * U1 * U1 * U1 * U1);
		c3 = (4 * c3E1 + c3E2) / (48 * U1 * U1 * U1);
		c4 = a1 * a2 - 2./3 * a1 * a1 * a1 - 1./6 * a3;

		/*
		d1 = 1/2 * (6 * a1 * a1 + a2 - 4 *b1 + 2 * b1) -
				1./6 * (120 * a1 * a1 * a1 - a3 + 
				6 * (24 * c1 - 6 * c2 + 2 * c3 - c4));
		*/
		d2 = 1/2 * (10 * a1 * a1 + a2 - 6 * b1 + 2 * b2) - 
				(128./3 * a1 * a1 * a1 - a3 / 6 + 2 * a1 * b1 - 
				a1 * b2 + 50 * c1 - 11 * c2 + 3 * c3 - c4);
		d3 = (2 * a1 * a1 - b1) - 1/3 * (88 * a1 * a1 * a1 + 
				3 * a1 * (5 * b1 - 2 * b2) + 3 * (35 * c1 -
				6 * c2 + c3));
		d4 = (-20./3 * a1 * a1 * a1 + a1 * (- 4 * b1 + b2) -
				10 * c1 + c2);
		
		z1 = d2 - d3 + d4;
		z2 = d3 - d4;
		z3 = d4;
		
		m1 = 2 * Math.log(U1) - 0.5 * Math.log(U21);
		v1 = Math.log(U21) - 2 * Math.log(U1);
		
		y  = Math.log(K);
		y1 = (m1 - y) / Math.sqrt(v1) + Math.sqrt(v1);
		y2 = y1 - Math.sqrt(v1);
		
		double N1 = 0.5 * (1 + StatUtil.erf(y1/Math.sqrt(2)));
		double N2 = 0.5 * (1 + StatUtil.erf(y2/Math.sqrt(2)));
		
		levy =  U1 * N1 - K * N2;
		
		double p   = 1 / Math.sqrt(2 * Math.PI * v1) * Math.exp(- (y - m1) * (y - m1)/(2 * v1));
		double py  = - p * (y - m1) / v1;
		double pyy = - (py * (y - m1)/v1 + p / v1);
		
		ju   = levy + K *( z1 * p + z2 * py + z3 * pyy);
		/*
		System.out.format("%.3f%n", K);
		System.out.format("%.3f%n", py);
		System.out.format("%.3f%n", pyy);
		*/
		return ju;
	}
	
	public static double MClog(){
		int Nruns = 100000;
		double V = 0;
		Random gen = new Random();
		for(int run = 0; run < Nruns; run++){
			double[] W = new double[N];
			for (int i=0; i<N; i++){
				W[i] = gen.nextGaussian();
			}
			Matrix matW = new Matrix(W,1);
				
			Matrix matC = new Matrix(C);
		    CholeskyDecomposition matChol = new CholeskyDecomposition(matC);
		    Matrix matL = matChol.getL();
		    Matrix matWcor = matL.times(matW.transpose());
		    /*
		    matW.print(3, 3);
		    matL.print(3, 3);
		    matWcor.print(3, 3);
		    */
		    double expo = 0;
		    double sum = 0;
		    for (int i=0; i<N; i++){
		    	expo = s[i] * matWcor.get(i,0) * Math.sqrt(T) - 0.5 * s[i] * s[i] * T;
		    	A[i] = A0[i] * Math.exp(expo);
		    	sum += w[i] * A[i];
		    }
		    //System.out.println(sum);
		    V += Math.max(0, sum - K);
		}
		V /= Nruns;
	    return V;
	}
	
	public static double MCnormal(){
		int Nruns = 100000;
		double V = 0;
		Random gen = new Random();
		for(int run = 0; run < Nruns; run++){
			double[] W = new double[N];
			for (int i=0; i<N; i++){
				W[i] = gen.nextGaussian();
			}
			Matrix matW = new Matrix(W,1);
				
			Matrix matC = new Matrix(C);
		    CholeskyDecomposition matChol = new CholeskyDecomposition(matC);
		    Matrix matL = matChol.getL();
		    Matrix matWcor = matL.times(matW.transpose());
		    /*
		    matW.print(3, 3);
		    matL.print(3, 3);
		    matWcor.print(3, 3);
		    */
		    double sum = 0;
		    for (int i=0; i<N; i++){
		    	A[i] = A0[i] * (1 + sn[i] * matWcor.get(i,0) * Math.sqrt(T));
		    	sum += w[i] * A[i];
		    }
		    //System.out.println(sum);
		    V += Math.max(0, sum - K);
		}
		V /= Nruns;
	    return V;
	}

	public static double Bachelier(){
		double B0 = 0;
		double sumWeights = 0;
		for (int i=0; i<N; i++){
			B0 += w[i] * A0[i];
			sumWeights += w[i];
		}
		if (sumWeights == 0) sumWeights = 1;
		
		Matrix matCov = new Matrix(Covn);
		Matrix matw = new Matrix(w,N);
		matCov = matw.transpose().times(matCov);
		matCov = matCov.times(matw);
		double sig = B0 * Math.sqrt(T * matCov.get(0, 0)) / sumWeights;
		
		double d = (B0 - K) / sig;
		double N1 = 0.5 * (1 + StatUtil.erf(d/Math.sqrt(2)));
		double N2 = 1/Math.sqrt(2 * Math.PI) * Math.exp(- 0.5 * d * d);
				
		return (B0 - K) * N1 + sig * N2;
	}
	
	public double calcFV(boolean lognormal, boolean MC){
		if (lognormal){
			if (MC) return MClog();
			else 	if (w[N-1] < 0) return Deng();
					else return Ju();
		}
		else {
			if (MC) return MCnormal();
			else return Bachelier();			
		}
	}
	
	public static void main(String args[]){
		System.out.println("rho; Log;  MClog; Norm; MCnor");

		for (int cor=25; cor<26; cor++){
			// basket(number of stocks, time to expiry, correlation, weight of first stock)
			basket B = new basket(3, 5.0, cor/100., -1.);
			
			
			double logPrice = B.calcFV(true, false);
			double normalPrice = B.calcFV(false, false);
			double MCnPrice = B.calcFV(false, true);
			double MClPrice = B.calcFV(true, true);
			System.out.format("%2d;  ", cor);
			System.out.format("%.3f; ", logPrice);
			System.out.format("%.3f; ", MClPrice);
			System.out.format("%.3f; ", normalPrice);
			System.out.format("%.3f%n",MCnPrice);			
		}
	}
}
