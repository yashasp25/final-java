//package final_year_project;
import java.util.Scanner;

public class Main {

    //solid segment area (h < t)
    static double areaSolid(double R, double h) {
        return (R * R * Math.acos((R - h) / R)) - ((R - h) * Math.sqrt((2 * R * h) - h * h));
    }

    //hollow segment area (h > t)
    static double areaHollow(double R, double t, double h) {
        double outer = areaSolid(R, h);
        double inner = areaSolid(R - t, h - t);
        return outer - inner;
    }

    // Function to get difference for bisection
    static double areaDifference(double R, double t, double h, double targetArea) {
        if (h < t) {
            return areaSolid(R, h) - targetArea;
        } else {
            return areaHollow(R, t, h) - targetArea;
        }
    }

    public static void main (String[] args) {
        Scanner sc = new Scanner(System.in);

        // Take inputs
        System.out.print("Enter Rt (top chord radius in mm): "); double Rt = sc.nextDouble();
        
        System.out.print("Enter Tt (top chord thickness in mm): "); double Tt = sc.nextDouble();
        
        System.out.print("Enter Rb (bottom chord radius in mm): "); double Rb = sc.nextDouble();
        
        System.out.print("Enter Tb (bottom chord thickness in mm): "); double Tb = sc.nextDouble();
        
        System.out.print("Enter d (center to center distance between the top and bottom chords in mm): "); double d = sc.nextDouble();
        
        System.out.print("Enter Fysc of Steel (yield stress of steel under compression in MPa): "); double fysc = sc.nextDouble();
        
        System.out.print("Enter Fyst of Steel (yield stress of steel under tension in MPa): "); double fyst = sc.nextDouble();
        
        //Steel Calculations
        double As1 = (2 * (Math.PI * ((2 * Rt*Tt) - (Tt*Tt))));System.out.printf("\n\nSteel Calculations:-\narea of steel of top chords = %.6f sqmm\n", As1);
        double As2 = ((Math.PI) * ((2 * Rb*Tb) - (Tb*Tb)));System.out.printf("area of steel of bottom chord = %.6f sqmm\n", As2);
        
        if(As2>As1){
            // Calculate target area for bottom chord segment
            double targetArea = (As2 - As1) / 2.0;

            System.out.printf("\nsegment area of steel = %.6f sqmm\n", targetArea);

            // Use bisection method to find h for bottom chord
            double low = 0.0;
            double high = 2 * Rb;  // Extended high limit for safety
            double mid = 0.0;
            double tol = 1e-6;

            while ((high - low) > tol) {
                mid = (low + high) / 2.0;
                double diff = areaDifference(Rb, Tb, mid, targetArea);
                if (diff > 0) {
                    high = mid;
                } else {
                    low = mid;
                }   
            }
            System.out.printf("EAA axis for steel will be at bottom chord (from topmost point) = %.2f mm\n", mid);
            if(mid<=Tb){
                double chord = 2*Math.sqrt((2*Rb*mid)-(mid *mid));
                double theta = 2*Math.asin(chord/(2*Rb));
                
                double centroid1 = (4 * Rb * Math.pow(Math.sin(theta / 2), 3))/(3 * (theta - Math.sin(theta)));
                
                System.out.printf("Centroid of minor segmented area of steel (from the center of bottom chord) = %.2f mm\n", centroid1);
                
                theta = 2*Math.PI - 2*Math.asin(chord/(2*Rb));
                double oarea = ((Rb*Rb)/2)*(theta - Math.sin(theta));
                double outerCentroid = (4 * Rb * Math.pow(Math.sin(theta / 2), 3)) / (3 * (theta - Math.sin(theta)));
                
                double numerator = (oarea * outerCentroid);
                double denominator = oarea - (Math.PI*(Rb-Tb)*(Rb-Tb));
                double centroid2 = numerator / denominator;
                
                System.out.printf("Centroid of major segmented area of steel (from the center of bottom chord) = %.2f mm\n", centroid2);
                
                double c1 = fysc*As1; System.out.printf("  c1 = %.3f N\n",c1);
                double c2 = fysc*targetArea; System.out.printf("  c2 = %.3f N\n",c2);
                double t = fyst*(As2 - targetArea); System.out.printf("  t = %.3f N\n\n",t);
                
                double Mr = (c1*(d-Rb+mid)) + (c2*(mid-Rb+centroid1)) + (t*(Rb-mid+centroid2));
                
                System.out.printf("Moment of resistence created by steel = "+Mr+" Nmm");
                
            } else{
                double ochord = 2*Math.sqrt((2*Rb*mid)-(mid *mid));
                double ichord = 2*Math.sqrt((2*(Rb-Tb)*(mid-Tb))-((mid-Tb) *(mid-Tb)));
                
                double otheta = 2*Math.asin(ochord/(2*Rb));
                double itheta = 2*Math.asin(ichord/(2*(Rb-Tb)));
                
                double oarea = ((Rb*Rb)/2)*(otheta - Math.sin(otheta));
                double iarea = (((Rb-Tb)*(Rb-Tb))/2)*(itheta - Math.sin(itheta));
                
                double outerCentroid = (4 * Rb * Math.pow(Math.sin(otheta / 2), 3)) / (3 * (otheta - Math.sin(otheta)));
                double innerCentroid = (4 * (Rb - Tb) * Math.pow(Math.sin(itheta / 2), 3)) / (3 * (itheta - Math.sin(itheta)));
                
                double numerator = (oarea * outerCentroid) - (iarea * innerCentroid);
                double denominator = oarea - iarea;
                double centroid1 = numerator / denominator;
                
                System.out.printf("Centroid of minor segmented area of steel (from the center of bottom chord) = %.2f mm\n", centroid1);
                
                otheta = 2*Math.PI - 2*Math.asin(ochord/(2*Rb));
                itheta = 2*Math.PI - 2*Math.asin(ichord/(2*(Rb-Tb)));
                
                oarea = ((Rb*Rb)/2)*(otheta - Math.sin(otheta));
                iarea = (((Rb-Tb)*(Rb-Tb))/2)*(itheta - Math.sin(itheta));
                
                outerCentroid = (4 * Rb * Math.pow(Math.sin(otheta / 2), 3)) / (3 * (otheta - Math.sin(otheta)));
                innerCentroid = (4 * (Rb - Tb) * Math.pow(Math.sin(itheta / 2), 3)) / (3 * (itheta - Math.sin(itheta)));
                
                numerator = (oarea * outerCentroid) - (iarea * innerCentroid);
                denominator = oarea - iarea;
                double centroid2 = numerator / denominator;
                
                System.out.printf("Centroid of major segmented area of steel (from the center of bottom chord) = %.2f mm\n", centroid2);
                
                double c1 = fysc*As1; System.out.printf("  c1 = %.3f N\n",c1);
                double c2 = fysc*targetArea; System.out.printf("  c2 = %.3f N\n",c2);
                double t = fyst*(As2 - targetArea);System.out.printf("  t = %.3f N\n",t);
                
                double Mr = (c1*(d-Rb+mid)) + (c2*(mid-Rb+centroid1)) + (t*(Rb-mid+centroid2));
                
                System.out.printf("Moment of resistence created by steel = "+Mr+(" Nmm"));
                
            }
        } 
        else if(As2<As1){
            // Calculate target area for bottom chord segment
            double targetArea = (As1 - As2 ) / 2.0;

            System.out.printf("\nsegment area for steel = %.6f sqmm\n", targetArea);
            
            targetArea = targetArea/2; 
             
            // Use bisection method to find h for bottom chord
            double low = 0.0;
            double high = 2 * Rt;  // Extended high limit for safety
            double mid = 0.0;
            double tol = 1e-6;

            while ((high - low) > tol) {
                mid = (low + high) / 2.0;
                double diff = areaDifference(Rt, Tt, mid, targetArea);

                if (diff > 0) {
                    high = mid;
                } else {
                    low = mid;
                }   
            }
            System.out.printf("EAA axis for steel will be at top chords (from bottomost point) = %.2f mm\n", mid);
            if(mid<=Tt){
                double chord = 2*Math.sqrt((2*Rt*mid)-(mid * mid));
                double theta = 2*Math.asin(chord/(2*Rt));
                
                double centroid1 = (4 * Rt * Math.pow(Math.sin(theta / 2), 3))/(3 * (theta - Math.sin(theta)));
                
                System.out.printf("Centroid of minor segmented area of steel (from the center of bottom chord) = %.2f mm\n", centroid1);
                
                theta = 2*Math.PI - 2*Math.asin(chord/(2*Rt));
                double oarea = ((Rt*Rt)/2)*(theta - Math.sin(theta));
                double outerCentroid = (4 * Rt * Math.pow(Math.sin(theta / 2), 3)) / (3 * (theta - Math.sin(theta)));
                
                double numerator = (oarea * outerCentroid);
                double denominator = oarea - (Math.PI*(Rt-Tt)*(Rt-Tt));
                
                double centroid2 = numerator / denominator;
                System.out.printf("Centroid of major segmented area of steel (from the center of bottom chord) = %.2f mm\n", centroid2);
                
                double c = fysc*(As1-(2*targetArea));System.out.printf("  c = %.3f N\n",c);
                double t1 = fyst*(2*targetArea);System.out.printf("  t1 = %.3f N\n",t1);
                double t2 = fyst*As2;System.out.printf("  t2 = %.3f N\n",t2);
                
                double Mr = (c*(Rt-mid+centroid2)) + (t1*(mid-Rt+centroid1)) + (t2*(d-Rt+mid));
                
                System.out.printf("Moment of resistence created by steel = "+Mr+(" Nmm"));
                
            } else{
                double ochord = 2*Math.sqrt((2*Rt*mid)-(mid *mid));
                double ichord = 2*Math.sqrt((2*(Rt-Tt)*(mid-Tt))-((mid-Tt) *(mid-Tt)));
                
                double otheta = 2*Math.asin(ochord/(2*Rt));
                double itheta = 2*Math.asin(ichord/(2*(Rt-Tt)));
                
                double oarea = ((Rt*Rt)/2)*(otheta - Math.sin(otheta));
                double iarea = (((Rt-Tt)*(Rt-Tt))/2)*(itheta - Math.sin(itheta));
                
                double outerCentroid = (4 * Rt * Math.pow(Math.sin(otheta / 2), 3)) / (3 * (otheta - Math.sin(otheta)));
                double innerCentroid = (4 * (Rt - Tt) * Math.pow(Math.sin(itheta / 2), 3)) / (3 * (itheta - Math.sin(itheta)));
                
                double numerator = (oarea * outerCentroid) - (iarea * innerCentroid);
                double denominator = oarea - iarea;
                double centroid1 = numerator / denominator;  
                
                System.out.printf("Centroid of minor segmented area of steel (from the center of top chord) = %.2f mm\n", centroid1);
                
                otheta = 2*Math.PI - 2*Math.asin(ochord/(2*Rt));
                itheta = 2*Math.PI - 2*Math.asin(ichord/(2*(Rt-Tt)));
                
                oarea = ((Rt*Rt)/2)*(otheta - Math.sin(otheta));
                iarea = (((Rt-Tt)*(Rt-Tt))/2)*(itheta - Math.sin(itheta));
                
                outerCentroid = (4 * Rt * Math.pow(Math.sin(otheta / 2), 3)) / (3 * (otheta - Math.sin(otheta)));
                innerCentroid = (4 * (Rt - Tt) * Math.pow(Math.sin(itheta / 2), 3)) / (3 * (itheta - Math.sin(itheta)));
                
                numerator = (oarea * outerCentroid) - (iarea * innerCentroid);
                denominator = oarea - iarea;
                double centroid2 = numerator / denominator;
                
                System.out.printf("Centroid of major segmented area of steel (from the center of top chord) = %.2f mm\n", centroid2);
                
                double c = fysc*(As1-(2*targetArea));System.out.printf("  c = %.3f N\n",c);
                double t1 = fyst*(2*targetArea);System.out.printf("  t1 = %.3f N\n",t1);
                double t2 = fyst*As2;System.out.printf("  t2 = %.3f N\n",t2);
                
                double Mr = (c*(Rt-mid+centroid2)) + (t1*(mid-Rt+centroid1)) + (t2*(d-Rt+mid));
                
                System.out.printf("Moment of resistence created by steel = "+Mr+(" Nmm"));
            }
        }
        else if(As1==As2){
            System.out.printf("EAA axis for steel will be at somewhere in between the Top and Bottom chords");
        }
        
    }
}