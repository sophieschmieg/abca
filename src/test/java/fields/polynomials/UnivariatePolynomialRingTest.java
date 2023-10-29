package fields.polynomials;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.Set;
import java.util.TreeSet;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;

class UnivariatePolynomialRingTest {

	@Test
	void testIterator() {
		int[] primes = { 2, 3, 5, 7, 11 };
		for (int prime : primes) {
			int expected = 1;
			for (int degree = -1; degree < 4; degree++) {
				PrimeField f = PrimeField.getPrimeField(prime);
				UnivariatePolynomialRing<PFE> ring = f.getUnivariatePolynomialRing();
				Set<Polynomial<PFE>> set = new TreeSet<>();
				for (Polynomial<PFE> t : ring.polynomialSet(degree)) {
					set.add(t);
				}
				assertEquals(expected, set.size());
				expected *= prime;
			}
		}
	}

	@Test
	void testMultiplication() {
		Integers z = Integers.z();
		UnivariatePolynomialRing<IntE> polynomials = z.getUnivariatePolynomialRing();
		PolynomialRing<IntE> multivariate = AbstractPolynomialRing.getPolynomialRing(z, Monomial.GREVLEX,
				new String[] { "X", "Y" });
		for (int size1 = -1; size1 < 10; size1++) {
			for (int size2 = -1; size2 < 10; size2++) {
				UnivariatePolynomial<IntE> polynomial1 = polynomials.getRandomElement(size1);
				UnivariatePolynomial<IntE> polynomial2 = polynomials.getRandomElement(size2);
				Polynomial<IntE> embed1 = multivariate.getEmbedding(polynomial1);
				Polynomial<IntE> embed2 = multivariate.getEmbedding(polynomial2);
				assertEquals(polynomials.getEmbedding(multivariate.multiply(embed1, embed2)),
						polynomials.multiply(polynomial1, polynomial2));
			}
		}
	}
}
