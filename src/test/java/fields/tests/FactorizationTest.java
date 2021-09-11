package fields.tests;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.PrimeField;
import fields.integers.Integers;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Element;
import fields.interfaces.Polynomial;
import fields.interfaces.Ring;
import fields.interfaces.Ring.FactorizationResult;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.numberfields.NumberField;

class FactorizationTest {

	private <T extends Element<T>> void testFactorizationOfRing(Ring<T> ring) {
		UnivariatePolynomialRing<T> polynomialRing = ring.getUnivariatePolynomialRing();
		List<UnivariatePolynomial<T>> polynomials = new ArrayList<>();
		UnivariatePolynomial<T> p;
		do {
			p = polynomialRing.getRandomElement(3);
		} while (p.equals(polynomialRing.zero()));
		polynomials.add(p);
		do {
			p = polynomialRing.getRandomElement(3);
		} while (p.equals(polynomialRing.zero()));
		polynomials.add(p);
		do {
			p = polynomialRing.getRandomElement(1);
		} while (p.equals(polynomialRing.zero()));
		polynomials.add(p);
		polynomials.add(polynomialRing.toUnivariate(polynomialRing.subtract(polynomialRing.getVarPower(4), polynomialRing.one())));
		for (int i = 0; i < polynomials.size(); i++) {
			for (int j = 0; j < polynomials.size(); j++) {
				Map<Polynomial<T>, Integer> factors = new TreeMap<>();
				UnivariatePolynomial<T> product = polynomialRing.multiply(polynomials.get(i), polynomials.get(j));
				FactorizationResult<Polynomial<T>, T> fs = ring.factorization(polynomials.get(i));
				assertTrue(ring.isUnit(fs.getUnit()));
				for (Polynomial<T> f : fs.primeFactors()) {
					int multiplicity = fs.multiplicity(f);
					f = polynomialRing.upToUnit(f);
					if (factors.containsKey(f)) {
						factors.put(f, factors.get(f) + multiplicity);
					} else {
						factors.put(f, multiplicity);
					}
				}
				fs = ring.factorization(polynomials.get(j));
				assertTrue(ring.isUnit(fs.getUnit()));
				for (Polynomial<T> f : fs.primeFactors()) {
					int multiplicity = fs.multiplicity(f);
						f = polynomialRing.upToUnit(f);
						if (factors.containsKey(f)) {
						factors.put(f, factors.get(f) + multiplicity);
					} else {
						factors.put(f, multiplicity);
					}
				}
				FactorizationResult<Polynomial<T>, T> factorization = ring.factorization(product);
				Polynomial<T> test = polynomialRing.getEmbedding(factorization.getUnit());
				assertTrue(ring.isUnit(factorization.getUnit()));
				//assertEquals(unit, factorization.getUnit());
				for (Polynomial<T> factor : factorization.primeFactors()) {
					assertTrue(factors.containsKey(polynomialRing.upToUnit(factor)));
					assertEquals(factors.get(polynomialRing.upToUnit(factor)), factorization.multiplicity(factor));
					test = polynomialRing.multiply(test, polynomialRing.power(factor, factorization.multiplicity(factor)));
				}
				assertEquals(product, test);
			}
		}
	}

	@Test
	void testFactorizationF2() {
		PrimeField f2 = PrimeField.getPrimeField(2);
		testFactorizationOfRing(f2);
	}

	@Test
	void testFactorizationF4() {
		FiniteField f4 = FiniteField.getFiniteField(4);
		testFactorizationOfRing(f4);
	}

	@Test
	void testFactorizationF8() {
		FiniteField f8 = FiniteField.getFiniteField(8);
		testFactorizationOfRing(f8);
	}

	@Test
	void testFactorizationF3() {
		testFactorizationOfRing(PrimeField.getPrimeField(3));
		FiniteField f3 = FiniteField.getFiniteField(3);
		testFactorizationOfRing(f3);
	}

	@Test
	void testFactorizationF9() {
		FiniteField f9 = FiniteField.getFiniteField(9);
		testFactorizationOfRing(f9);
	}

	@Test
	void testFactorizationF27() {
		FiniteField f27 = FiniteField.getFiniteField(27);
		testFactorizationOfRing(f27);
	}

	@Test
	void testFactorizationF5() {
		testFactorizationOfRing(PrimeField.getPrimeField(5));
		FiniteField f5 = FiniteField.getFiniteField(5);
		testFactorizationOfRing(f5);
	}

	@Test
	void testFactorizationF25() {
		FiniteField f25 = FiniteField.getFiniteField(25);
		testFactorizationOfRing(f25);
	}

	@Test
	void testFactorizationQ() {
		testFactorizationOfRing(Rationals.q());
	}

	@Test
	void testFactorizationZ() {
		Integers z = Integers.z();
		z.factorization(z.getUnivariatePolynomialRing().getPolynomial(z.getInteger(-126), z.getInteger(-18), z.getInteger(-1)));
		testFactorizationOfRing(z);
	}

	@Test
	void testFactorizationGauss() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> ring = q.getUnivariatePolynomialRing();
		testFactorizationOfRing(new NumberField(ring.getPolynomial(q.one(), q.zero(), q.one())));
	}

	@Test
	void testFactorizationEisenstein() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> ring = q.getUnivariatePolynomialRing();
		testFactorizationOfRing(new NumberField(ring.getPolynomial(q.one(), q.one(), q.one())));
	}

	@Test
	void testFactorizationSqrt5() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> ring = q.getUnivariatePolynomialRing();
		testFactorizationOfRing(new NumberField(ring.getPolynomial(q.getInteger(5), q.zero(), q.one())));
	}

	@Test
	void testFactorizationF5Bivariate() {
		PrimeField fp = PrimeField.getPrimeField(5);
		testFactorizationOfRing(fp.getUnivariatePolynomialRing());
	}

}
