import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

public class OperationsForKernelComputation {

    int k, m, N, min2mkminus1, min2mk;
    BigDecimal alphabet_size;
    int[][] sequens;
    BigDecimal[][] K;
    BigDecimal[][] K1;
    BigDecimal[] intSizes;
    int kmers[][];
    int totalKmers;
    double sigma;
    int B;
    List<Integer> repititions;
    
    public OperationsForKernelComputation() {
        
    }

    /**
    * Constructor for new Approximation based Algorithm
    * 
    * Initialization of Parameters, and Load Sequence File
    */
    public OperationsForKernelComputation(String filename, int k, int m, int min2mkminus1, int min2mk, BigDecimal alphabet_size, int b, double sigma, int num_seqs) throws FileNotFoundException, IOException {
        BufferedReader dataBR = new BufferedReader(new FileReader(new File(filename)));
        this.k = k;
        this.m = m;
        this.N = num_seqs;
        this.min2mkminus1 = min2mkminus1;
        this.min2mk = min2mk;
        this.alphabet_size = alphabet_size;
        this.sigma = sigma;
        B = b;
        sequens = new int[this.N][];
        intSizes = new BigDecimal[this.min2mk + 1];
        K = new BigDecimal[this.N][this.N];
        
        String inputLine;
        int i = 0;
        while ((inputLine = dataBR.readLine()) != null) {
            sequens[i] = convertToInts(inputLine.split(" "));
            totalKmers += sequens[i].length - k + 1;
            i++;
        }
        kmers = new int[totalKmers + 1][this.k + 1];
        extractKmers();
    }
    
    /**
    * Constructor for exact Brute-Force Algorithm
    * 
    * Initialization of Parameters, and Load Sequence File
    */
    public OperationsForKernelComputation(String filename, int k, int m, int min2mkminus1, int min2mk, BigDecimal alphabet_size, int num_seqs) throws FileNotFoundException, IOException {
        BufferedReader dataBR = new BufferedReader(new FileReader(new File(filename)));
        this.k = k;
        this.m = m;
        this.N = num_seqs;
        this.min2mkminus1 = min2mkminus1;
        this.min2mk = min2mk;
        this.alphabet_size = alphabet_size;

        sequens = new int[this.N][];
        intSizes = new BigDecimal[this.min2mk + 1];
        
        // For Brute Force only
        intersecSizes();
        
        K = new BigDecimal[this.N][this.N];
        String inputLine;
        int i = 0;
        while ((inputLine = dataBR.readLine()) != null) {
            sequens[i] = convertToInts(inputLine.split(" "));
            i++;
        }
    }
    
    /**
    * This function extracts each k-mer from all the sequences
    * 
    * and store into kmers array
    */
    public void extractKmers() {
        int c = 0;
        for(int i = 0 ; i < sequens.length ; i++) {
            for(int ii = 0 ; ii < sequens[i].length - k + 1 ; ii++) {
                for(int jj = 0 ; jj < k ; jj++)
                    kmers[c][jj] = sequens[i][ii+jj];
                kmers[c][k] = i;
                c++;
            }
        }
        
        for(int jj = 0 ; jj <= k ; jj++)
            kmers[c][jj] = 99999;
    }
    
    /**
    * This function computes the final kernel matrix,
    * and stores in K
    */
    public void pairsOfKmersAtDistanceMin2mk() throws FileNotFoundException {
        double[] Mis = new double[min2mk + 1];
        
        double[][][] meanCis = new double[min2mkminus1 + 1][this.N][this.N];
        double[][][] M2 = new double[min2mkminus1 + 1][this.N][this.N];
        double[][] currCis;
        int dist = 0, iter, temp, seqi, seqj;
        
        double[][] nchoosek = new double[min2mkminus1 + 1][min2mkminus1];
        double[] kchoosed = new double[min2mkminus1 + 1];
        
        for(temp = 0 ; temp <= min2mkminus1 ; temp++)
            kchoosed[temp] = xchoosey_int(k, temp);

        repititions = new ArrayList<Integer>();
        int[] comb = getNextCombination(dist, (int) kchoosed[dist]);
        
        // computation of exact mismatches, dist == 0
        sort(comb, dist);
        
        currCis = new double[this.N][this.N];
        computeCis(currCis, comb, dist);
        
        double max_sd = online_variance(currCis, meanCis[dist], M2[dist], 0);

        // computation commulative mismatches, dist > 0
        for (dist = 1 ; dist <= min2mkminus1 ; dist++) {   // dist in {1....min}
            meanCis[dist] = new double[this.N][this.N];
            M2[dist] = new double[this.N][this.N];
            
            repititions = new ArrayList<Integer>();
            for (iter = 0 ; iter < B ; iter++) {
                comb = getNextCombination(dist, (int) kchoosed[dist]);
                sort(comb, dist);
                
                currCis = new double[this.N][this.N];
                computeCis(currCis, comb, dist);

                max_sd = online_variance(currCis, meanCis[dist], M2[dist], iter);
                max_sd = Math.sqrt(max_sd / (iter + 1) * (1 - (iter + 1) / kchoosed[dist]));
                if (iter < 3 || kchoosed[dist] < 70 && iter + 1 < kchoosed[dist]) continue;
                if(max_sd < sigma)
                    break;
            }
        } //end of dist loop
        
        for(dist = 1 ; dist <= min2mkminus1 ; dist++) {
            for(temp = 0 ; temp < dist ; temp++)
                nchoosek[dist][temp] = xchoosey_int(k - temp, dist - temp);
        }
        
        //Scaling Sample Ci's up
        for(seqi = 0 ; seqi < this.N ; seqi++) {
            for(seqj = seqi ; seqj < this.N ; seqj++) {
                for(dist = 1 ; dist <= min2mkminus1 ; dist++) {
                    meanCis[dist][seqi][seqj] = meanCis[dist][seqi][seqj] * kchoosed[dist];
                }
            }
        }
        
        // Computing Ii's
        intersecSizes();
        
        for(seqi = 0 ; seqi < this.N ; seqi++) {
            for(seqj = seqi ; seqj < this.N ; seqj++) {
                // Computing Mi's
                Mis[0] = meanCis[0][seqi][seqj];
                for(dist = 1 ; dist <= min2mkminus1 ; dist++) {
                    Mis[dist] = meanCis[dist][seqi][seqj];
                    for(temp = 0 ; temp < dist ; temp++) {
                        Mis[dist] -= nchoosek[dist][temp] * Mis[temp];
                    }
                    if(Mis[dist] < 0)
                        Mis[dist] = 0;
                }
                if(min2mkminus1 < min2mk) {
                    int sumMis = 0;
                    for(temp = 0 ; temp <= min2mkminus1 ; temp++)
                        sumMis += Mis[temp];
                    Mis[min2mk] = (sequens[seqi].length - k + 1) * (sequens[seqj].length - k + 1) - sumMis;
                    if(Mis[min2mk] < 0)
                        Mis[min2mk] = 0;
                }
                // Computing K(X,Y|k,m)
                K[seqi][seqj] = BigDecimal.ZERO;
                for(temp = 0 ; temp <= min2mkminus1 ; temp++)
                    K[seqi][seqj] = K[seqi][seqj].add(intSizes[temp].multiply(BigDecimal.valueOf(Mis[temp])));
                if(min2mkminus1 < min2mk)
                    K[seqi][seqj] = K[seqi][seqj].add(intSizes[min2mk].multiply(BigDecimal.valueOf(Mis[min2mk])));
                if(seqi != seqj)
                    K[seqj][seqi] = K[seqi][seqj];
            }
        }
    }
    
    /**
    * This function returns the variance calculated between currentCis and meanCis
    */
    public double online_variance(double[][] data, double[][] mean, double M2[][], int n) {
        double max_var = 0.0; n += 1;
        for(int i = 0 ; i < this.N ; i++) {
            for(int j = 0 ; j < this.N ; j++) {
                double delta = data[i][j] - mean[i][j];
                mean[i][j] += delta / n;
                double delta2 = data[i][j] - mean[i][j];
                M2[i][j] += delta * delta2;
                
                if(M2[i][j] > max_var)
                    max_var = M2[i][j];
            }
        }
        return n > 1 ? max_var / (n - 1) : 9999999;
    }
    
    /**
    * This function implements xchoosey for primitive data type int
    */
    public double xchoosey_int(int n, int k) {
        if (k < 0 || k > n)
            return 0;
        if (k == 0 || k == n)
            return 1;
        k = Math.min(k, n - k); // take advantage of symmetry
        int c = 1;
        for (int i = 0 ; i < k ; i++)
            c = c * (n - i) / (i + 1);
        return c;
    }
    
    /**
    * Returns the next (not already used) combination
    */
    public int[] getNextCombination(int dist, int total) {
        Random random = new Random();
        
        int rand_perm = random.nextInt(total);
        
        while(repititions.contains(rand_perm))
            rand_perm = random.nextInt(total);
        
        repititions.add(rand_perm);
        int[] comb = element(k, k - dist, rand_perm);
        
        return comb;
    }
    
    /**
    * Returns the next random combination
    */
    public int[] element(int n, int k, int m) {
        int[] ans = new int[k];

        int a = n;
        int b = k;
        int x = ((int) xchoosey_int(n, k) - 1) - m;  // x is the "dual" of m

        for (int i = 0; i < k; ++i) {
            a = largestV(a, b, x);          // largest value v, where v < a and vCb < x
            x = x - (int) xchoosey_int(a, b);
            b = b - 1;
            ans[i] = (n - 1) - a;
        }

        return ans;
    }
    
    /**
    * returns the largest value v where v < a and Choose(v,b) <= x
    */
    public int largestV(int a, int b, int x) {
        int v = a - 1;

        while ((int) xchoosey_int(v, b) > x)
            --v;

        return v;
    }
    
    /**
    * This function sorts all the k-mers stored in kmers array,
    * 
    * corresponding to a given combination in kminusi
    */
    public void sort(final int[] kminusi, final int d) {
        Arrays.sort(kmers, new Comparator<int[]>() {
            @Override
            public int compare(int[] o1, int[] o2) {
                int move = 0;
                for(int p = 0 ; p < k - d ; p++)
                {
                    if(o1[kminusi[p]] < o2[kminusi[p]]) {
                        move = -1;
                        break;
                    }
                    else if(o1[kminusi[p]] > o2[kminusi[p]]) {
                        move = 1;
                        break;
                    }
                }
                if(move == 0) {
                    if(o1[k] < o2[k])
                        move = -1;
                    else if(o1[k] > o2[k])
                        move = 1;
                    else
                        move = 0;
                }
                return move;
            }
        });
    }
    
    /**
    * This function counts the similar kmers for all pairs of sequences and,
    * 
    * returns the counts into currentCis
    */
    public void computeCis(double Ks[][], int[] comb, int d) {
        int i = 0, j = 0, st, en; boolean same;
        int[] sameKmers = new int[this.N];
        int[] sameKmersSupport = new int[this.N];
        int numKeys;
        
        while(i < totalKmers) {
            st = i; same = true;
            while(i < totalKmers && same) {
                for(j = 0 ; j < k - d ; j++) {
                    if(kmers[i][comb[j]] != kmers[i + 1][comb[j]]) {
                        same = false;
                        break;
                    }
                }
                i++;
            }
            en = i - 1;
            
            if (en - st + 1 > 2 * this.N) {
                for(j = 0 ; j < this.N ; j++)
                    sameKmers[j] = 0;

                for(j = st ; j <= en ; j++)
                    sameKmers[kmers[j][k]]++;

                numKeys = 0;
                for(j = 0 ; j < this.N ; j++){
                    if(sameKmers[j] > 0)
                    {
                        sameKmersSupport[numKeys] = j;
                        numKeys++;
                    }
                }

                for(j = 0 ; j < numKeys ; j++) {
                    for(int jj = 0 ; jj < numKeys ; jj++)
                        Ks[sameKmersSupport[j]][sameKmersSupport[jj]] += sameKmers[sameKmersSupport[j]] * sameKmers[sameKmersSupport[jj]]; 
                }
            } else {
                for (j = st ; j <= en ; j++) {
                    for (int jj = st ; jj <= en ; jj++)
                        Ks[kmers[j][k]][kmers[jj][k]]++;
                }
            }
        }
    }
    
    /**
    * This function computes the Intersection Sizes for d = 0,...,min(2*m, k)
    */
    public void intersecSizes() {

        // Intersection sizes
        for (int d = 0; d <= min2mk ; d++) {
            intSizes[d] = BigDecimal.ZERO;
            for (int i = 0; i <= m; i++) {
                for (int j = 0; j <= m; j++) {
                    for (int t = 0; t <= (i + j - d) / 2; t++)
                        intSizes[d] = intSizes[d].add(xchoosey(2 * d - i - j + 2 * t, d - i + t).multiply(xchoosey(d, i + j - 2 * t - d)).multiply(pow(alphabet_size.subtract(BigDecimal.valueOf(2.0)), i + j - 2 * t - d)).multiply(xchoosey(k - d, t)).multiply(pow(alphabet_size.subtract(BigDecimal.valueOf(1.0)), t)));
                }
            }
        }
    }
    
    
    /**
    * This function read an existing Kernel file and,
    * 
    * Normalize it for SVM
    */
    /** Kernel Matrix Normalization **/
    public void readFile(String filename, int num_seqs) throws FileNotFoundException, IOException {
        BufferedReader dataBR = new BufferedReader(new FileReader(new File(filename)));
        this.N = num_seqs;
        K1 = new BigDecimal[this.N][this.N];
        String inputLine;
        int i = 0;
        while ((inputLine = dataBR.readLine()) != null) {
            convertToBigInts(inputLine.split(" "), i);
            i++;
        }

        normalize();
    }

    /**
    * This function normalize the Kernel File
    */
    public void normalize() {
        SQRT obj = new SQRT();
        double frac = 0.000000000001;
        for (int i = 0; i < this.N; i++) {
            for (int j = i; j < this.N; j++) {
                if (i != j) {
                    K1[i][j] = K1[i][j].divide(obj.bigSqrt(K1[i][i], new MathContext(16, RoundingMode.HALF_EVEN)).add(BigDecimal.valueOf(frac)).multiply(obj.bigSqrt(K1[j][j], new MathContext(16, RoundingMode.HALF_EVEN)).add(BigDecimal.valueOf(frac))), 16, RoundingMode.HALF_EVEN);
                    K1[j][i] = K1[j][i].divide(obj.bigSqrt(K1[j][j], new MathContext(16, RoundingMode.HALF_EVEN)).add(BigDecimal.valueOf(frac)).multiply(obj.bigSqrt(K1[i][i], new MathContext(16, RoundingMode.HALF_EVEN)).add(BigDecimal.valueOf(frac))), 16, RoundingMode.HALF_EVEN);
                }
            }
            K1[i][i] = K1[i][i].divide(obj.bigSqrt(K1[i][i], new MathContext(16, RoundingMode.HALF_EVEN)).add(BigDecimal.valueOf(frac)).multiply(obj.bigSqrt(K1[i][i], new MathContext(16, RoundingMode.HALF_EVEN)).add(BigDecimal.valueOf(frac))), 16, RoundingMode.HALF_EVEN);
        }
    }

    /**
    * This function converts the numbers read from Kernel file into BigInts
    */
    private void convertToBigInts(String[] numbers, int row) {
        for (int i = 0; i < numbers.length; i++) {
            K1[row][i] = new BigDecimal(numbers[i]);
        }
    }

    /**
    * This function converts the numbers read from sequence file into integers
    */
    private int[] convertToInts(String[] numbers) {
        int[] ints = new int[numbers.length];
        for (int i = 0; i < ints.length; i++) {
            ints[i] = Integer.parseInt(numbers[i]);
        }
        return ints;
    }

    /**
    * This function writes the Kernel file into a given filename
    */
    public void writeToFile(String filename) throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(new FileOutputStream(new File(filename)));

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                pw.print(K[i][j].setScale(0, RoundingMode.CEILING).toString()+" ");
            }
            pw.println();
        }
        pw.close();
    }
    
    /**
    * This function writes the Normalized file into a given filename
    */
    public void writeToFileAfterNormalization(String filename) throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(new FileOutputStream(new File(filename)));

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                pw.print(K1[i][j].toString()+" ");
            }
            pw.println();
        }
        pw.close();
    }

    /**
    * This function implements a^b for BigDecimals
    */
    public static BigDecimal pow(BigDecimal a, int b) {
        BigDecimal result = BigDecimal.ONE;

        for (int i = 0; i < b; i++) {
            result = result.multiply(a);
        }

        return result;
    }

    /**
    * This function implements xchoosey for BigDecimals
    */
    public static BigDecimal xchoosey(int x, int y) {

        if (x < 0) {
            return BigDecimal.ZERO;
        }
        if (y < 0) {
            return BigDecimal.ZERO;
        }
        if (y == 0) {
            return BigDecimal.ONE;
        }
        if (x == 0 && y != 0) {
            return BigDecimal.ZERO;
        }
        if (x < y) {
            return BigDecimal.ZERO;
        }

        return factorial(x).divide(factorial(y).multiply(factorial(x - y)));
    }

    /**
    * This function implements x! for BigDecimals
    */
    public static BigDecimal factorial(int x) {

        if (x == 0) {
            return BigDecimal.ONE;
        } else {
            BigDecimal prod = BigDecimal.ONE;
            for (int i = x; i > 0; i--) {
                prod = prod.multiply(BigDecimal.valueOf(i));
            }

            return prod;
        }
    }
    
    /**
    * This function computes exact Kernel matrix using Brute Force approach
    */
    public Thread kernelUsing_BruteForce(final int st, final int en) {

        Thread t = new Thread(new Runnable() {
            public void run() {
                for (int i = st ; i < st + en ; i++) {
                    for (int j = i ; j < N ; j++) {
                        double[] nPairs = bruteForce(i, j);
                        K[i][j] = BigDecimal.ZERO;
                        for (int h = min2mk ; h >= 0 ; h--) {
                            K[i][j] = K[i][j].add(intSizes[h].multiply(BigDecimal.valueOf(nPairs[h])));
                        }
                        if (i != j) {
                            K[j][i] = K[i][j];
                        }
                    }
                }
            }
        });
        t.start();
        return t;
    }

    /**
    * This function computes Mis for Brute Force algorithms
    */
    public double[] bruteForce(int s1, int s2) {
        double M_dist[] = new double[k + 1];
        for (int i = 0; i < sequens[s1].length - (k - 1); i++) {
            for (int j = 0; j < sequens[s2].length - (k - 1); j++) {
                int dist = 0;
                for (int s = 0 ; s < k ; s++) {
                    if (sequens[s1][i + s] != sequens[s2][j + s]) {
                        dist++;
                    }
                }
                M_dist[dist] += 1.0;
            }
        }
        return M_dist;
    }
}
