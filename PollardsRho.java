/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pollardsrho;

import java.math.BigInteger;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author hanxinlei
 * implementation of Pollard's Rho algorithm
 * specifically for elliptic curve
 */
public class PollardsRho {

    private static final BigInteger THREE = new BigInteger("3");
    private static final BigInteger TWO = new BigInteger("2");
    private static final BigInteger ONE = new BigInteger("1");
    private static final BigInteger ZERO = new BigInteger("0");
    private static final BigInteger[] NEUTRAL = {BigInteger.ZERO, BigInteger.ONE};

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        System.err.close();
        
        int N=100;
        long k=0;
        long tmp=(2<<15)-17;
        BigInteger p = new BigInteger(tmp+"");
        BigInteger d = new BigInteger("154");
        BigInteger n = new BigInteger("16339");
        BigInteger[] a = {new BigInteger("12"), new BigInteger("61833")};
        for(int i=0;i<N;i++)
        {
            System.err.println("pollardsrho.PollardsRho.main()");
            try {
                k+=check(a,d,p,n);
            } catch (Exception ex) {
                Logger.getLogger(PollardsRho.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        k/=N;
        System.out.println("k.avg="+k+" #under the condition N="+N);
    }

    /**
     *
     * @param a1
     * @param a2
     * @param d
     * @param p
     * @return a3
     */
    public static BigInteger[] mul(BigInteger[] a1, BigInteger[] a2, BigInteger d, BigInteger p) {
        BigInteger x, y;
        //(a - b) mod p = ((a mod p - b mod p) + p) mod p
        //(a / b) mod p = ((a mod p) * (b^(-1) mod p)) mod p
        BigInteger x1 = a1[0], x2 = a2[0];
        BigInteger y1 = a1[1], y2 = a2[1];

        BigInteger x1y2 = x1.mod(p).multiply(y2.mod(p)).mod(p);
        BigInteger y1x2 = y1.mod(p).multiply(x2.mod(p)).mod(p);

        BigInteger x1x2 = x1.mod(p).multiply(x2.mod(p)).mod(p);
        BigInteger y1y2 = y1.mod(p).multiply(y2.mod(p)).mod(p);

        BigInteger x1x2y1y2 = x1y2.multiply(y1x2).mod(p);
        BigInteger dx = d.mod(p).multiply(x1x2y1y2).mod(p);
        
        BigInteger xu=x1y2.add(y1x2).mod(p);
        BigInteger xd=BigInteger.ONE.add(dx).mod(p);
                
        BigInteger yu=y1y2.subtract(x1x2).mod(p);
        BigInteger yd=BigInteger.ONE.subtract(dx).mod(p);
        x = xu.multiply(xd.modInverse(p)).mod(p);
        y = yu.multiply(yd.modInverse(p)).mod(p);

        BigInteger[] result = {x, y};
        return result;
    }

    /**
     *
     * @param a
     * @param m
     * @param d
     * @param p
     * @return a^m
     */
    public static BigInteger[] exp(BigInteger[] a, BigInteger m, BigInteger d, BigInteger p) {
        BigInteger[] b = NEUTRAL.clone();

        int k = m.bitLength();
        for (int i = k-1; i >=0; i--) {
            b = mul(b, b, d, p);
            if (m.testBit(i)) {
                b = mul(b, a, d, p);
            }

        }
        return b;
    }

    public static BigInteger[] rho(BigInteger[] a, BigInteger[] b, 
            BigInteger d, BigInteger p, BigInteger n) throws Exception {
        BigInteger m = null;
        BigInteger k = ZERO;
        BigInteger ai = ZERO, bi = ZERO, a2i = ZERO, b2i = ZERO;
        BigInteger[] zi = NEUTRAL.clone(), z2i = NEUTRAL.clone();

        while (true) {
            k = k.add(ONE);
            System.err.println("pollardsrho.PollardsRho.rho()");
            System.err.println(" k="+k);
            System.err.println("a="+"("+a[0]+","+a[1]+")");
            System.err.println("b="+"("+b[0]+","+b[1]+")");
            System.err.print("zi="+"("+zi[0]+","+zi[1]+")");
            System.err.println(" z2i="+"("+z2i[0]+","+z2i[1]+")");
            if (zi[0].mod(THREE).equals(ZERO)) {
                zi = mul(b, zi, d, p);
                ai = ai.add(ONE).mod(n);
            } else if (zi[0].mod(THREE).equals(ONE)) {
                zi = mul(zi, zi, d, p);
                ai = ai.add(ai).mod(n);
                bi = bi.add(bi).mod(n);

            } else if (zi[0].mod(THREE).equals(TWO)) {
                zi = mul(a, zi, d, p);
                bi = bi.add(ONE).mod(n);
            }else{
                    throw new Exception("zi move err");
                    }
            

            for (int i = 0; i < 2; i++) {
                if (z2i[0].mod(THREE).equals(ZERO)) {
                    z2i = mul(b, z2i, d, p);
                    a2i = a2i.add(ONE).mod(n);
                } else if (z2i[0].mod(THREE).equals(ONE)) {
                    z2i = mul(z2i, z2i, d, p);
                    a2i = a2i.add(a2i).mod(n);
                    b2i = b2i.add(b2i).mod(n);

                } else if (z2i[0].mod(THREE).equals(TWO)) {
                    z2i = mul(a, z2i, d, p);
                    b2i = b2i.add(ONE).mod(n);
                }
                else{
                    throw new Exception("z2i move err");
                    }
            }
            

            if(zi[0].equals(z2i[0])&&zi[1].equals(z2i[1]))
                break;
        }
        //TODO something is wrong here
        //cause collision do happens, but m and _m doesn't match
        
        m=b2i.subtract(bi).mod(n).multiply(ai.subtract(a2i).mod(n).modInverse(n)).mod(n);

        BigInteger[] result = {m, k};
        return result;
    }

    public static long check(BigInteger[] a, BigInteger d, BigInteger p, BigInteger n) throws Exception {
        BigInteger m=new BigInteger(n.bitLength(), new Random()).mod(n);
        BigInteger[] b=exp(a, m, d, p);
        System.err.println("pollardsrho.PollardsRho.check()");
        System.err.println("a="+"("+a[0]+","+a[1]+")");
        System.err.println("b="+"("+b[0]+","+b[1]+")");
        BigInteger[] testResult=rho(a, b, d, p, n);
        BigInteger _m=testResult[0];
        System.err.println("m="+m+", _m="+_m);
        if(!m.equals(_m))
            throw new RuntimeException();
        return testResult[1].longValue();
    }
}
