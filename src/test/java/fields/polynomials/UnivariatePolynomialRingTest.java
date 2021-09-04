package fields.polynomials;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.Set;
import java.util.TreeSet;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.interfaces.Polynomial;
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

}
