package cryptography;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.math.BigInteger;

import org.junit.jupiter.api.Test;

import cryptography.Sidh.SidhPrivateKey;
import cryptography.Sidh.SidhPublicKey;
import cryptography.interfaces.DhScheme.Role;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import varieties.curves.elliptic.EllipticCurve;

class SidhTest {

	void doTest(IntE prime) {
		Sidh sidh = new Sidh(prime);
		EllipticCurve<FFE> supersingular = sidh.getPublicCurve();
		System.out.println(supersingular + "; j = " + supersingular.jInvariant());
		System.out.println("Number of points: " + supersingular.getNumberOfElements());
		System.out.println("Parameters: " + sidh);
		SidhPrivateKey alicePrivateKey = sidh.createPrivateKey(Role.Alice);
		System.out.println("Private key Alice: " + alicePrivateKey);
		SidhPublicKey alicePublicKey = sidh.createPublicKey(alicePrivateKey);
		System.out.println("Public key Alice: " + alicePublicKey);
		SidhPrivateKey bobPrivateKey = sidh.createPrivateKey(Role.Bob);
		System.out.println("Private key Bob: " + bobPrivateKey);
		SidhPublicKey bobPublicKey = sidh.createPublicKey(bobPrivateKey);
		System.out.println("Public key Bob: " + bobPublicKey);

		FFE sharedSecretAlice = sidh.agree(alicePrivateKey, bobPublicKey);
		FFE sharedSecretBob = sidh.agree(bobPrivateKey, alicePublicKey);

		System.out.println("shared secret Alice: " + sharedSecretAlice);
		System.out.println("shared secret Bob:   " + sharedSecretBob);
		assertEquals(sharedSecretAlice, sharedSecretBob);
	}

	@Test
	void test() {
		Integers z = Integers.z();
		BigInteger one = BigInteger.ONE;
		BigInteger two = BigInteger.TWO;
		BigInteger three = BigInteger.valueOf(3);
		PrimeField f54 = PrimeField.getPrimeField(two.pow(5).multiply(three.pow(4)).subtract(one));
		doTest(z.getInteger(f54.getNumberOfElements()));
		PrimeField f9754 = PrimeField.getPrimeField(two.pow(97).multiply(three.pow(54)).subtract(one));
		doTest(z.getInteger(f9754.getNumberOfElements()));
		PrimeField fp434 = PrimeField.getPrimeField(two.pow(0xd8).multiply(three.pow(0x89)).subtract(one));
		doTest(z.getInteger(fp434.getNumberOfElements()));
	}

}
