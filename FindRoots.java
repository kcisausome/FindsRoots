import java.io.PrintStream;
import java.util.*;
import java.util.Scanner;
import java.io.*;
import java.io.FileReader;
import java.util.StringTokenizer;
import java.util.Arrays;
//solves a system of equations using partial pivoting gaussian elimination

public class FindRoots {
    private Scanner scn = new Scanner(System.in);
  
    public static void main(String[] args) throws Exception{
    	new FindRoots();
    }

    public FindRoots()throws Exception{
      System.out.println("find roots of _, type:\n'1' ~ f(x) = 2x^3 - 11.7x^2 + 17.7x - 5\n'2' ~ f(x) = x + 10 - xcosh(50/x)");
      int eq = scn.nextInt();
      System.out.println("type: \n'b' = bisection method\n'n' = newton-raphson method\n's' = secant method\n'f' = false position method\n'm' = modified secant ");
      String method = scn.next();
      
    	switch (method) {
	    case "b":
	      bisection(eq);
	      break;
	    case "n":
	      newtonRapson(eq);
	      break;
	    case "s":
	      secant(eq);
	      break;
	    case "f":
	       falsePosition(eq);
	       break;
	    case "m":
	      	modifiedSecant(eq);
	    default:
	      System.out.println("\nSo I poured my root beer in a square glass...\n\n\n\n\n\n\n\nnow I just have beer\n\n");
	    }
    }
   /*
    *
    *
    */
   	public void bisection(int eq){
   		System.out.println("enter 2 numbers:");//assumes numbers 'bracket' the root
   		double ptA = scn.nextDouble();
   		double ptB = scn.nextDouble();
   		double avg = 0;
      double fA;
      double fB;
      double fAvg = 1;
      double prevAvg = 1;
   		int iteration = 0;
   		double error = 1;
      System.out.println("Bisection method");
      System.out.printf("%2s |  %3s  |  %4s  |  %3s   |  %3s  |  %3s  |  %3s  |  %3s", "n","an","bn","cn","fA","fB","fc","error\n");
      System.out.println("---------------------------------------------------------------");
   		while(iteration < 100 && error > 0.01){
          avg = (ptA+ptB)/2;
          if(eq == 1){
            fA = equationA(ptA);
            fB = equationA(ptB);
            fAvg = equationA(avg);
          }else{
            fA = equationB(ptA);
            fB = equationB(ptB);
            fAvg = equationB(avg);
          }
          
          System.out.printf("%2d | %5.2f |  %5.2f |  %5.2f | %5.2f | %5.2f | %5.2f |  ",iteration,ptA,ptB,avg,fA,fB,fAvg);
          if(fA<0 && fB>0){
            if(fAvg > 0){
              ptB = avg;
            }else{
              ptA = avg;
            }
          }else if(fA>0 && fB<0){
            if(fAvg > 0){
              ptA = avg;
            }else{
              ptB = avg;
            }
          }else if(fAvg < 0.0001){
            System.out.println("root: "+avg);
          }
          
          if(iteration != 0){
            error = Math.abs(avg - prevAvg)/Math.abs(avg);
            System.out.printf("%.3f\n",error);
          }else{
            System.out.println(" - ");
          }
          prevAvg = avg;
          iteration++;  
   		}
      System.out.println("root: "+avg);
   	}

   	public void newtonRapson(int eq){
      System.out.print("starting point: ");
   		double xi = scn.nextDouble();
      double fXi;
      double fprime;
      double xi2 = 0;
      double error = 1;
      int iteration = 0;
      System.out.println("newtonRapson method");
      System.out.println(" n |   xn  |   xn+1 | f(xn) | f'(xn) |  error");
      System.out.println("---------------------------------------------");
      while(iteration<100 && error > 0.01){
        if(eq == 1){
          fXi = equationA(xi);
          fprime = equationAderivative(xi);
        }else{
          fXi = equationB(xi);
          fprime = equationBderivative(xi);
        }
        xi2 = xi - (fXi/fprime);
        error = (Math.abs(xi2-xi))/Math.abs(xi2);
        System.out.printf("%2d |%5.2f | %5.2f | %5.2f | %6.2f | %5.3f \n",iteration,xi,xi2,fXi,fprime,error);
        xi = xi2;
        iteration++;
      }
      System.out.println("root: "+xi2);
   	}

   	public void secant(int eq){
      System.out.println("2 starting points: ");
      double xi0 = scn.nextDouble();
      double xi1 = scn.nextDouble();
      double fxi0;
      double fxi1;
      double xi2 = 0;
      double fxi2;
      double prevXi2 = 1;
      double error = 1;
      int iteration = 0;
      System.out.println("Secant method");
      System.out.println(" n |  Xn-1 |   Xn   |  Xn+1  | f(Xn-1)| f(Xn) |f(Xn+1)| error");
      System.out.println("---------------------------------------------------------------");
      while(iteration < 100 && error > 0.01){
        if(eq == 1){
          fxi0 = equationA(xi0);
          fxi1 = equationA(xi1);
          xi2 = xi1 - ((xi1-xi0)/(fxi1-fxi0))*fxi1;
          fxi2 = equationA(xi2); 
        }else{
          fxi0 = equationB(xi0);
          fxi1 = equationB(xi1);
          xi2 = xi1 - ((xi1-xi0)/(fxi1-fxi0))*fxi1;
          fxi2 = equationB(xi2);
        }
        
        System.out.printf("%2d |%5.2f | %5.2f | %5.2f | %6.2f | %5.2f | %5.2f |  ",iteration,xi0,xi1,xi2,fxi0,fxi1,fxi2);
        if(iteration != 0){
          error = Math.abs(xi2 - prevXi2)/Math.abs(xi2);
          System.out.printf("%.3f\n",error);
        }else{
          System.out.println(" - ");
        }
        prevXi2 = xi2;
        
        if(Math.abs(fxi0)<Math.abs(fxi1)){
          xi0 = xi1;
          xi1 = xi2;
        }else{
          xi0 = xi2;
        }

        iteration++;
      }
      System.out.println("root: "+xi2);
   	}

   	public void falsePosition(int eq){
      System.out.println("2 starting points: ");
      double a = scn.nextDouble();
      double b = scn.nextDouble();
      double fa;
      double fb;
      double c = 0;
      double fc;
      double prevC = 1;
      double error = 1;
      int iteration = 0;
      System.out.println("falsePosition method");
      System.out.println(" n |   An  |    Bn  |    Cn  |  f(An) | f(Bn) | f(Cn) | error");
      System.out.println("----------------------------------------------------------------");
      while(error > 0.01 && iteration<100 && Math.abs(a-b)>0.001){
        if(eq == 1){
          fa = equationA(a);
          fb = equationA(b);
        }else{
          fa = equationB(a);
          fb = equationB(b);
        }

        c = a - (fa*(b-a))/(fb-fa);
        if(eq == 1){
          fc = equationA(c);
        }else{
          fc = equationB(c);
        }

        System.out.printf("%2d |%5.2f | %5.2f | %5.2f | %6.2f | %5.2f | %5.2f |  ",iteration,a,b,c,fa,fb,fc);
        if(iteration != 0){
          error = Math.abs(c - prevC)/Math.abs(c);
          System.out.printf("%.3f\n",error);
        }else{
          System.out.println(" - ");
        }
        prevC = c;

        if(fc > 0){
          b = c;
        }else{
          a = c;
        }

        iteration++;
      }
      System.out.println("root: "+ c);
   	}

   	public void modifiedSecant(int eq){
      System.out.print("starting point:");
      double xi1 = scn.nextDouble();
      double deltaXi = 0.01*xi1; 
      double xi2 = 1;
      double error = 1;
      double fXi;
      double fXi_DXi;
      double prevXi2 = 1;
      int iteration = 0;
      System.out.println("modifiedSecant method");
      System.out.println(" n |   xn  |   xn+1 | error");
      System.out.println("------------------------------");
      while(iteration < 100 && error > 0.01){
        if(eq == 1){
          fXi = equationA(xi1);
          fXi_DXi = equationA(xi1+(0.01*xi1));
        }else{
          fXi = equationB(xi1);
          fXi_DXi = equationB(xi1+(0.01*xi1));
        }
        xi2 = xi1 - (deltaXi/(fXi_DXi-fXi))*fXi;

        System.out.printf("%2d |%5.2f | %5.2f | ",iteration,xi1,xi2);
        if(iteration != 0){
          error = Math.abs(xi2 - prevXi2)/Math.abs(xi2);
          System.out.printf("%.3f\n",error);
        }else{
          System.out.println(" - ");
        }
        prevXi2 = xi2;
        xi1 = xi2;
        iteration++;
      }
      System.out.println("root: "+ xi2);     
   	}

   	public double equationA(double pt){
   		return 2*Math.pow(pt,3) - 11.7*Math.pow(pt,2) + 17.7*pt -5;
   	}

   	public double equationAderivative(double pt){
   		return 6*Math.pow(pt,2) - 23.4*pt + 17.7;
   	}

    public double equationB(double pt){
      return pt + 10 - pt*Math.cosh(50/pt);
    }

    public double equationBderivative(double pt){
      return (50*Math.sinh(50/pt))/pt - Math.cosh(50/pt) + 1;
    }

    public static void exit(){
    	System.out.println("Usage error: please type b,n,s,f, or m");
    	System.exit(0);
    }
}
