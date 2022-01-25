package fields.tests;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.math.BigInteger;
import java.security.SecureRandom;

import org.junit.jupiter.api.Test;

import fields.finitefields.ModuloIntegerRing;
import fields.finitefields.ModuloIntegerRing.ModuloIntegerRingElement;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Ring.ExtendedEuclideanResult;

class RecoverRSAKeyTest {

	private IntE recoverPublicKey(int publicExponent, IntE message1, IntE signature1, IntE message2, IntE signature2) {
		Integers z = Integers.z();
		IntE recovered1 = z.power(signature1, publicExponent);
		System.out.println("signature1^exp done"/* + recovered1*/);
		IntE recovered2 = z.power(signature2, publicExponent);
		System.out.println("signature2^exp done"/* + recovered2*/);
			IntE modulus = z.gcd(z.subtract(recovered1, message1), z.subtract(recovered2, message2));
			System.out.println("gcd = " + modulus);
			
			for (IntE prime : z.setOfPrimes()) {
			if (prime.compareTo(z.getInteger(257)) > 0) {
				break;
			}
			while (z.isDivisible(modulus, prime)) {
				modulus = z.divideChecked(modulus, prime);
				System.out.println("Dividing off " + prime + " = " + modulus);
				}
		}
		return modulus;
	}

	@Test
	void test() {
		Integers z = Integers.z();
		int exponent = 17;
		int modulusSize = 32;
		IntE n;
		IntE d;
		ExtendedEuclideanResult<IntE> ee;
		ModuloIntegerRing zModN;
		do {
			IntE p = z.getInteger(BigInteger.probablePrime(modulusSize / 2, new SecureRandom()));
			IntE q = z.getInteger(BigInteger.probablePrime(modulusSize / 2, new SecureRandom()));
			n = z.multiply(p, q);
			zModN = new ModuloIntegerRing(n.getValue());
			IntE toitent = z.multiply(z.subtract(p, z.one()), z.subtract(q, z.one()));
			ee = z.extendedEuclidean(z.getInteger(exponent), toitent);
			d = ee.getCoeff1();
		} while (!ee.getGcd().equals(z.one()));
		System.out.println("Public key modulus is " + n);
		ModuloIntegerRingElement message1 = zModN.getRandomElement();
		System.out.println("Message one is " + message1);
		ModuloIntegerRingElement message2 = zModN.getRandomElement();
		System.out.println("Message two is " + message2);
		ModuloIntegerRingElement signature1 = zModN.power(message1, d);
		System.out.println("Signature one is " + signature1);
		ModuloIntegerRingElement signature2 = zModN.power(message2, d);
		System.out.println("Signature two is " + signature2);
		assertEquals(message1, zModN.power(signature1, exponent));
		assertEquals(message2, zModN.power(signature2, exponent));
		assertEquals(n, recoverPublicKey(exponent, zModN.lift(message1), zModN.lift(signature1), zModN.lift(message2),
				zModN.lift(signature2)));
	}

}
