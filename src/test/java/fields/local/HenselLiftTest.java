package fields.local;

import java.math.BigInteger;

import org.junit.jupiter.api.Test;

import fields.interfaces.Polynomial;
import fields.interfaces.Ring.FactorizationResult;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.PAdicField.PAdicNumber;

class HenselLiftTest {

	@Test
	void test() {
		PAdicField z2 = new PAdicField(BigInteger.TWO, 70);
		UnivariatePolynomialRing<PAdicNumber> r = z2.getUnivariatePolynomialRing();
		UnivariatePolynomial<PAdicNumber> f = r.getPolynomial(z2.getInteger(-17), z2.zero(), z2.one());
		FactorizationResult<Polynomial<PAdicNumber>, PAdicNumber> factorization = z2.factorization(f);
		System.out.println(factorization);
		for (Polynomial<PAdicNumber> factor : factorization.primeFactors()) {
			int root = z2.roundToInteger(z2.negative(r.toUnivariate(factor).univariateCoefficient(0)), 32)
					.intValueExact();
			System.out.println(root);
			System.out.println(root + "^2 = " + (root * root));
		}
		f = r.getPolynomial(z2.getInteger(2), z2.one(), z2.one());
		System.out.println(z2.factorization(f));
		UnivariatePolynomial<PAdicNumber> g1 = r.getPolynomial(z2.getInteger(3), z2.zero(), z2.one());
		System.out.println(z2.ringOfIntegers().theMontesAlgorithm(g1));
		UnivariatePolynomial<PAdicNumber> g2 = r.getPolynomial(z2.getInteger(1), z2.zero(), z2.one());
		System.out.println(z2.ringOfIntegers().theMontesAlgorithm(g2));
		UnivariatePolynomial<PAdicNumber> g3 = r.getPolynomial(z2.getInteger(15), z2.getInteger(8), z2.one());
		System.out.println(z2.ringOfIntegers().theMontesAlgorithm(g3));
	}

}
