import java.io.IOException;
import java.math.BigDecimal;

public class Main {
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        // TODO code application logic here
		
		double sigma = 0.0; int b = 0;
        if(args.length < 7) {
            sigma = 0.5;
        } else
			sigma = Double.parseDouble(args[6]);
        
        if(args.length < 6) {
            b = 300;
        } else
			b = Integer.parseInt(args[5]);
        
        if(args.length > 4) {
            int k = Integer.parseInt(args[2]), m = Integer.parseInt(args[3]), num_seqs = Integer.parseInt(args[1]), min2mkminus1 = Math.min(2*m, k - 1), min2mk = Math.min(2*m, k);
            BigDecimal alphabet_size = new BigDecimal(args[4]);
            String dataFile = args[0];

            OperationsForKernelComputation kop = new OperationsForKernelComputation(dataFile, k, m, min2mkminus1, min2mk, alphabet_size, b, sigma, num_seqs);
            kop.pairsOfKmersAtDistanceMin2mk();
            kop.writeToFile("Kernel-k"+k+"-m"+m+".txt");
        } else {
			System.out.println("Some arguments are missing.");
        }
    }
}
