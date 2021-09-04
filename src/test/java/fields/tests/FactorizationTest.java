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
import fields.interfaces.Field;
import fields.interfaces.Polynomial;
import fields.interfaces.Ring.FactorizationResult;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.numberfields.NumberField;

class FactorizationTest {

	private <T extends Element<T>> void testFactorizationOfField(Field<T> field) {
		UnivariatePolynomialRing<T> ring = field.getUnivariatePolynomialRing();
		List<UnivariatePolynomial<T>> polynomials = new ArrayList<>();
		UnivariatePolynomial<T> p;
		do {
			p = ring.getRandomElement(3);
		} while (p.equals(ring.zero()));
		polynomials.add(p);
		do {
			p = ring.getRandomElement(3);
		} while (p.equals(ring.zero()));
		polynomials.add(p);
		do {
			p = ring.getRandomElement(1);
		} while (p.equals(ring.zero()));
		polynomials.add(p);
		polynomials.add(ring.toUnivariate(ring.subtract(ring.getVarPower(4), ring.one())));
		for (int i = 0; i < polynomials.size(); i++) {
			for (int j = 0; j < polynomials.size(); j++) {
				Map<Polynomial<T>, Integer> factors = new TreeMap<>();
				UnivariatePolynomial<T> product = ring.multiply(polynomials.get(i), polynomials.get(j));
				FactorizationResult<Polynomial<T>> fs = field.factorization(polynomials.get(i));
				Polynomial<T> unit = fs.getUnit();
				assertTrue(ring.isUnit(unit));
				for (Polynomial<T> f : fs.primeFactors()) {
					if (factors.containsKey(f)) {
						factors.put(f, factors.get(f) + fs.multiplicity(f));
					} else {
						factors.put(f, fs.multiplicity(f));
					}
				}
				fs = field.factorization(polynomials.get(j));
				assertTrue(ring.isUnit(fs.getUnit()));
				unit = ring.multiply(unit, fs.getUnit());
				for (Polynomial<T> f : fs.primeFactors()) {
					if (factors.containsKey(f)) {
						factors.put(f, factors.get(f) + fs.multiplicity(f));
					} else {
						factors.put(f, fs.multiplicity(f));
					}
				}
				FactorizationResult<Polynomial<T>> factorization = field.factorization(product);
				Polynomial<T> test = factorization.getUnit();
				assertTrue(ring.isUnit(factorization.getUnit()));
				assertEquals(unit, factorization.getUnit());
				for (Polynomial<T> factor : factorization.primeFactors()) {
					assertTrue(factors.containsKey(factor));
					assertEquals(factors.get(factor), factorization.multiplicity(factor));
					test = ring.multiply(test, ring.power(factor, factorization.multiplicity(factor)));
				}
				assertEquals(product, test);
			}
		}
	}

	@Test
	void testFactorizationF2() {
		PrimeField f2 = PrimeField.getPrimeField(2);
		testFactorizationOfField(f2);
	}

	@Test
	void testFactorizationF4() {
		FiniteField f4 = FiniteField.getFiniteField(4);
		testFactorizationOfField(f4);
	}

	@Test
	void testFactorizationF8() {
		FiniteField f8 = FiniteField.getFiniteField(8);
		testFactorizationOfField(f8);
	}

	@Test
	void testFactorizationF3() {
		testFactorizationOfField(PrimeField.getPrimeField(3));
		FiniteField f3 = FiniteField.getFiniteField(3);
		testFactorizationOfField(f3);
	}

	@Test
	void testFactorizationF9() {
		FiniteField f9 = FiniteField.getFiniteField(9);
		testFactorizationOfField(f9);
	}

	@Test
	void testFactorizationF27() {
		FiniteField f27 = FiniteField.getFiniteField(27);
		testFactorizationOfField(f27);
	}

	@Test
	void testFactorizationF5() {
		testFactorizationOfField(PrimeField.getPrimeField(5));
		FiniteField f5 = FiniteField.getFiniteField(5);
		testFactorizationOfField(f5);
	}

	@Test
	void testFactorizationF25() {
		FiniteField f25 = FiniteField.getFiniteField(25);
		testFactorizationOfField(f25);
	}

	@Test
	void testFactorizationQ() {
		Integers z = Integers.z();
		z.factorization(z.getUnivariatePolynomialRing().getPolynomial(z.getInteger(126), z.getInteger(18), z.one()));
		testFactorizationOfField(Rationals.q());
	}

	@Test
	void testFactorizationGauss() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> ring = q.getUnivariatePolynomialRing();
		testFactorizationOfField(new NumberField(ring.getPolynomial(q.one(), q.zero(), q.one())));
	}

	@Test
	void testFactorizationEisenstein() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> ring = q.getUnivariatePolynomialRing();
		testFactorizationOfField(new NumberField(ring.getPolynomial(q.one(), q.one(), q.one())));
	}

	@Test
	void testFactorizationSqrt5() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> ring = q.getUnivariatePolynomialRing();
		testFactorizationOfField(new NumberField(ring.getPolynomial(q.getInteger(5), q.zero(), q.one())));
	}

}
