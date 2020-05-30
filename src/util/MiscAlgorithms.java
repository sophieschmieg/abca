package util;

import java.math.BigInteger;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;

public class MiscAlgorithms {
	/**
	 * Finds the solutions to a set of modular equations.
	 * 
	 * @param moduloValues the value x the result is supposed to have modulo the
	 *                     corresponding moduli.
	 * @param moduli       the moduli m for the equations.
	 * @param balanced     the solution is unique modulo the lcm of the moduli. If
	 *                     balanced == true, the returned solution will be between
	 *                     -lcm/2 and lcm/2, otherwise between 0 and lcm.
	 * @return n with n = x mod m.
	 */
	public static BigInteger chineseRemainder(List<BigInteger> moduloValues, List<BigInteger> moduli,
			boolean balanced) {
		if (moduloValues.size() != moduli.size()) {
			throw new ArithmeticException("Lists need to be same size for chinese remainder theorem");
		}
		BigInteger modulus = BigInteger.ONE;
		for (BigInteger m : moduli) {
			modulus = modulus.multiply(m).divide(modulus.gcd(m));
		}
		BigInteger result = BigInteger.ZERO;
		for (int i = 0; i < moduli.size(); i++) {
			BigInteger m = moduli.get(i);
			BigInteger multiplier = modulus.divide(m);
			multiplier = multiplier.multiply(multiplier.modInverse(m)).mod(modulus);
			result = result.add(moduloValues.get(i).multiply(multiplier)).mod(modulus);
		}
		if (balanced && result.compareTo(modulus.shiftRight(1)) > 0) {
			result = result.subtract(modulus);
		}
		return result;
	}

	public static String toHexDigit(byte b) {
		String s = Integer.toHexString(Byte.toUnsignedInt(b));
		while (s.length() < 2) {
			s = "0" + s;
		}
		return s;
	}

	public static String toHex(byte[] bytes) {
		StringBuffer buf = new StringBuffer();
		buf.append("0x");
		for (byte b : bytes) {
			buf.append(toHexDigit(b));
		}
		return buf.toString();
	}

	public static BigInteger randomBigInteger(Random rand, BigInteger max) {
		BigInteger num;
		do {
			byte bytes[] = new byte[max.toByteArray().length];
			rand.nextBytes(bytes);
			bytes[0] = (byte)(Byte.toUnsignedInt(bytes[0]) & 0x7f);
			num = new BigInteger(bytes);
			for (int i = num.bitLength(); i > max.bitLength(); i--) {
				num = num.clearBit(i-1);
			}
		} while (num.compareTo(BigInteger.ZERO) < 0 || num.compareTo(max) >= 0);
		return num;
	}

	public static Map<BigInteger, Integer> primeDecomposition(BigInteger n) {
		return naivePrimeDecomposition(n, n);
	}

	public static Map<BigInteger, Integer> naivePrimeDecomposition(BigInteger n, BigInteger max) {
		Map<BigInteger, Integer> result = new TreeMap<>();
		BigInteger zero = BigInteger.ZERO;
		BigInteger one = BigInteger.ONE;
		BigInteger two = BigInteger.TWO;

		BigInteger prime = two;
		int power = 0;
		while (prime.compareTo(n.sqrt().add(one)) < 0 && prime.compareTo(max) < 0) {
			if (n.mod(prime).equals(zero)) {
				power++;
				n = n.divide(prime);
				result.put(prime, power);
			} else {
				prime = prime.nextProbablePrime();
				power = 0;
			}
		}
		return result;
	}
}
