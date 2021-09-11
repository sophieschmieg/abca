package fields.local;

import java.math.BigInteger;

import org.junit.jupiter.api.Test;

import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.PAdicField.PAdicNumber;

class HenselLiftTest {

	@Test
	void test() {
		PAdicField z2 = new PAdicField(BigInteger.TWO, 32);
		UnivariatePolynomialRing<PAdicNumber> r = z2.getUnivariatePolynomialRing();
		UnivariatePolynomial<PAdicNumber> f = r.getPolynomial(z2.getInteger(2), z2.one(), z2.one());
		System.out.println(z2.factorization(f));
		UnivariatePolynomial<PAdicNumber> g1 = r.getPolynomial(z2.getInteger(3), z2.zero(), z2.one());
		System.out.println(z2.ringOfIntegers().theMontesAlgorithm(g1));
		UnivariatePolynomial<PAdicNumber> g2 = r.getPolynomial(z2.getInteger(1), z2.zero(), z2.one());
		System.out.println(z2.ringOfIntegers().theMontesAlgorithm(g2));
		UnivariatePolynomial<PAdicNumber> g3 = r.getPolynomial(z2.getInteger(15), z2.getInteger(8), z2.one());
		System.out.println(z2.ringOfIntegers().theMontesAlgorithm(g3));
	}

}
