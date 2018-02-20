/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pollardsrho;

import java.math.BigInteger;
import java.util.Random;

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
        int N=1;
        long k=0;
        long tmp=(2<<16)-17;
        BigInteger p = new BigInteger(tmp+"");
        BigInteger d = new BigInteger("154");
        BigInteger n = new BigInteger("16339");
        BigInteger[] a = {new BigInteger("12"), new BigInteger("61833")};
        for(int i=0;i<N;i++)
        {
            k+=check(a,d,p,n);
        }
        k/=N;
        System.out.println(k);

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

        BigInteger x1 = a1[0], x2 = a2[0];
        BigInteger y1 = a1[1], y2 = a2[1];

        BigInteger x1y2 = x1.multiply(y2).mod(p);
        BigInteger y1x2 = y1.multiply(x2).mod(p);

        BigInteger x1x2 = x1.multiply(x2).mod(p);
        BigInteger y1y2 = y1.multiply(y2).mod(p);

        BigInteger x1x2y1y2 = x1y2.multiply(y1x2).mod(p);
        BigInteger dx = d.multiply(x1x2y1y2).mod(p);

        x = x1y2.add(y1x2).mod(p).divide(BigInteger.ONE.add(dx).mod(p)).mod(p);
        y = y1y2.subtract(x1x2).mod(p).divide(BigInteger.ONE.subtract(dx).mod(p)).mod(p);

        BigInteger[] result = {x, y};
//        System.out.println("("+result[0]+","+result[1]+")");
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
        for (int i = 0; i <= k - 1; i++) {
            b = mul(b, b, d, p);
            if (m.testBit(i)) {
                b = mul(b, a, d, p);
            }

        }
        return b;
    }

    public static BigInteger[] rho(BigInteger[] a, BigInteger[] b, BigInteger d, BigInteger p, BigInteger n) {
        BigInteger m = null;
        BigInteger k = ZERO;
        BigInteger ai = ZERO, bi = ZERO, a2i = ZERO, b2i = ZERO;
        BigInteger[] zi = NEUTRAL.clone(), z2i = NEUTRAL.clone();

        while (true) {
            k = k.add(ONE);
            if (zi[0].mod(THREE) == ZERO) {
                zi = mul(b, zi, d, p);
                ai = ai.add(ONE).mod(p);
            } else if (zi[0].mod(THREE) == ONE) {
                zi = mul(zi, zi, d, p);
                ai = ai.add(ai).mod(p);
                bi = bi.add(bi).mod(p);

            } else if (zi[0].mod(THREE) == TWO) {
                zi = mul(a, zi, d, p);
                bi = bi.add(ONE).mod(p);
            }

            for (int i = 0; i < 2; i++) {
                if (z2i[0].mod(THREE) == ZERO) {
                    z2i = mul(b, z2i, d, p);
                    a2i = a2i.add(ONE).mod(p);
                } else if (z2i[0].mod(THREE) == ONE) {
                    z2i = mul(z2i, z2i, d, p);
                    a2i = a2i.add(a2i).mod(p);
                    b2i = b2i.add(b2i).mod(p);

                } else if (z2i[0].mod(THREE) == TWO) {
                    z2i = mul(a, z2i, d, p);
                    b2i = b2i.add(ONE).mod(p);
                }
            }
            if(zi.equals(z2i))
                break;
        }
        
        m=b2i.subtract(bi).mod(n).divide(ai.subtract(a2i).mod(n)).mod(n);

        BigInteger[] result = {m, k};
        return result;
    }

    public static long check(BigInteger[] a, BigInteger d, BigInteger p, BigInteger n) {
        BigInteger m=new BigInteger(n.bitLength(), new Random());
        BigInteger[] b=exp(a, m, d, p);
        BigInteger[] testResult=rho(a, b, d, p, n);
        BigInteger _m=testResult[0];
        
        if(!m.equals(_m))
            throw new RuntimeException();
        return testResult[1].longValue();
    }
}
