package fields.polynomials;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring.FactorizationResult;

class MultivariateFactorizationTest {

	@Test
	void squareFreeTest() {
		PrimeField f5 = PrimeField.getPrimeField(5);
		PolynomialRing<PFE> f5x = AbstractPolynomialRing.getPolynomialRing(f5, 2, Monomial.GREVLEX);
		Polynomial<PFE> irrdX1 = f5x.subtract(f5x.getVarPower(1, 2), f5x.getInteger(2));
		Polynomial<PFE> irrdY1 = f5x.subtract(f5x.getVarPower(2, 2), f5x.getInteger(2));
		Polynomial<PFE> redX1 = f5x.add(f5x.getVarPower(1, 2), f5x.one());
		Polynomial<PFE> redY1 = f5x.add(f5x.getVarPower(2, 2), f5x.one());
		Polynomial<PFE> unit = f5x.getEmbedding(f5.getInteger(3));
		List<Polynomial<PFE>> tests = new ArrayList<>();
		tests.add(f5x.multiply(irrdX1, irrdX1));
		tests.add(f5x.multiply(unit, irrdX1, irrdX1));
		tests.add(f5x.multiply(unit, redX1, redY1));
		tests.add(f5x.multiply(irrdY1, irrdY1));
		tests.add(f5x.multiply(irrdY1, f5x.multiply(unit, irrdX1, irrdY1)));
		tests.add(f5x.multiply(unit, redX1, redY1));
		for (Polynomial<PFE> test : tests) {
			FactorizationResult<Polynomial<PFE>, PFE> squareFree = f5x.squareFreeFactorization(test);
			assertTrue(f5.isUnit(squareFree.getUnit()));
			Polynomial<PFE> product = f5x.getEmbedding(squareFree.getUnit());
			for (Polynomial<PFE> factor : squareFree.primeFactors()) {
				product = f5x.multiply(product, f5x.power(factor, squareFree.multiplicity(factor)));
			}
			assertEquals(test, product);
		}
	}

	@Test
	void standardBivariateTest() {
		//fail();
		PrimeField fp = PrimeField.getPrimeField(5);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fp, 2, Monomial.GREVLEX);
		List<Polynomial<PFE>> polynomials = new ArrayList<>();
		Polynomial<PFE> p;
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
		polynomials.add(polynomialRing.subtract(polynomialRing.getVarPower(1, 4), polynomialRing.one()));
		polynomials.add(polynomialRing.subtract(polynomialRing.getVarPower(2, 4), polynomialRing.one()));
		for (int i = 0; i < polynomials.size(); i++) {
			for (int j = 0; j < polynomials.size(); j++) {
				Map<Polynomial<PFE>, Integer> factors = new TreeMap<>();
				Polynomial<PFE> product = polynomialRing.multiply(polynomials.get(i), polynomials.get(j));
				FactorizationResult<Polynomial<PFE>, Polynomial<PFE>> fs = polynomialRing.uniqueFactorization(polynomials.get(i));
				Polynomial<PFE> unit = fs.getUnit();
				assertTrue(polynomialRing.isUnit(unit));
				for (Polynomial<PFE> f : fs.primeFactors()) {
					int multiplicity = fs.multiplicity(f);
					f = polynomialRing.upToUnit(f);
					if (factors.containsKey(f)) {
						factors.put(f, factors.get(f) + multiplicity);
					} else {
						factors.put(f, multiplicity);
					}
				}
				fs = polynomialRing.uniqueFactorization(polynomials.get(j));
				assertTrue(polynomialRing.isUnit(fs.getUnit()));
				unit = polynomialRing.multiply(unit, fs.getUnit());
				for (Polynomial<PFE> f : fs.primeFactors()) {
					int multiplicity = fs.multiplicity(f);
					f = polynomialRing.upToUnit(f);
					if (factors.containsKey(f)) {
						factors.put(f, factors.get(f) + multiplicity);
					} else {
						factors.put(f, multiplicity);
					}
				}
				System.out.println();
				System.out.println("Polynomial: " + product);
				FactorizationResult<Polynomial<PFE>, Polynomial<PFE>> factorization = polynomialRing.uniqueFactorization(product);
				Polynomial<PFE> test = factorization.getUnit();
				assertTrue(polynomialRing.isUnit(factorization.getUnit()));
				System.out.println("Actual Factors: " + factors);
				System.out.println("Computed Factors: " + factorization);
				for (Polynomial<PFE> factor : factorization.primeFactors()) {
					assertTrue(factors.containsKey(polynomialRing.upToUnit(factor)));
					assertEquals(factors.get(polynomialRing.upToUnit(factor)), factorization.multiplicity(factor));
					test = polynomialRing.multiply(test,
							polynomialRing.power(factor, factorization.multiplicity(factor)));
				}
				assertEquals(product, test);
				System.out.println("Test case done");
				System.out.println();
			}
		}
	}

	@Test
	void hangingTest() {
		PrimeField fp = PrimeField.getPrimeField(5);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fp, 2, Monomial.GREVLEX);
		// Y^2 + -2*Y + 1
		Polynomial<PFE> candidate1 = polynomialRing.getVarPower(2, 2);
		candidate1 = polynomialRing.add(candidate1, polynomialRing.multiply(-2, polynomialRing.getVar(2)));
		candidate1 = polynomialRing.add(candidate1, polynomialRing.one());
		FactorizationResult<Polynomial<PFE>, Polynomial<PFE>> factors = polynomialRing.uniqueFactorization(candidate1);
		assertTrue(polynomialRing.isUnit(factors.getUnit()));
		Polynomial<PFE> test = factors.getUnit();
		for (Polynomial<PFE> factor : factors.primeFactors()) {
			test = polynomialRing.multiply(test, polynomialRing.power(factor, factors.multiplicity(factor)));
		}
		assertEquals(candidate1, test);

		// -1*X^3Y + -2*X^2Y^2 + -2*XY^3 + Y^4 + -2*X^3 + -1*XY^2 + 2*Y^3 + X^2 + XY +
		// -2*Y^2
		Polynomial<PFE> candidate2 = polynomialRing.multiply(-1, polynomialRing.getVarPower(1, 3),
				polynomialRing.getVarPower(2, 1));
		candidate2 = polynomialRing.add(candidate2,
				polynomialRing.multiply(-2, polynomialRing.getVarPower(1, 2), polynomialRing.getVarPower(2, 2)));
		candidate2 = polynomialRing.add(candidate2,
				polynomialRing.multiply(-2, polynomialRing.getVarPower(1, 1), polynomialRing.getVarPower(2, 3)));
		candidate2 = polynomialRing.add(candidate2, polynomialRing.getVarPower(2, 4));
		candidate2 = polynomialRing.add(candidate2, polynomialRing.multiply(-2, polynomialRing.getVarPower(1, 3)));
		candidate2 = polynomialRing.add(candidate2,
				polynomialRing.multiply(-1, polynomialRing.getVarPower(1, 1), polynomialRing.getVarPower(2, 2)));
		candidate2 = polynomialRing.add(candidate2, polynomialRing.multiply(2, polynomialRing.getVarPower(2, 3)));
		candidate2 = polynomialRing.add(candidate2, polynomialRing.multiply(1, polynomialRing.getVarPower(1, 2)));
		candidate2 = polynomialRing.add(candidate2,
				polynomialRing.multiply(1, polynomialRing.getVarPower(1, 1), polynomialRing.getVarPower(2, 1)));
		candidate2 = polynomialRing.add(candidate2, polynomialRing.multiply(-2, polynomialRing.getVarPower(2, 2)));
		factors = polynomialRing.uniqueFactorization(candidate2);
		assertTrue(polynomialRing.isUnit(factors.getUnit()));
		test = factors.getUnit();
		for (Polynomial<PFE> factor : factors.primeFactors()) {
			test = polynomialRing.multiply(test, polynomialRing.power(factor, factors.multiplicity(factor)));
		}
		assertEquals(candidate2, test);

		// -2*X^2Y + 2*XY^2 + -2*XY + Y^2 + -1
		Polynomial<PFE> candidate3 = polynomialRing.multiply(-2, polynomialRing.getVarPower(1, 2),
				polynomialRing.getVarPower(2, 1));
		candidate3 = polynomialRing.add(candidate3,
				polynomialRing.multiply(2, polynomialRing.getVarPower(1, 1), polynomialRing.getVarPower(2, 2)));
		candidate3 = polynomialRing.add(candidate3,
				polynomialRing.multiply(-2, polynomialRing.getVarPower(1, 1), polynomialRing.getVarPower(2, 1)));
		candidate3 = polynomialRing.add(candidate3, polynomialRing.multiply(1, polynomialRing.getVarPower(2, 2)));
		candidate3 = polynomialRing.add(candidate3, polynomialRing.getInteger(-1));
		factors = polynomialRing.uniqueFactorization(candidate3);
		assertTrue(polynomialRing.isUnit(factors.getUnit()));
		test = factors.getUnit();
		for (Polynomial<PFE> factor : factors.primeFactors()) {
			test = polynomialRing.multiply(test, polynomialRing.power(factor, factors.multiplicity(factor)));
		}
		assertEquals(candidate3, test);
	}

	@Test
	void extensionRequiredTest() {
		PrimeField fp = PrimeField.getPrimeField(5);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fp, 2, Monomial.GREVLEX);
		Polynomial<PFE> frobenius = polynomialRing.subtract(polynomialRing.getVarPower(2, 5), polynomialRing.getVar(2));
		Polynomial<PFE> testPolynomial = polynomialRing.subtract(polynomialRing.getVarPower(1, 2),
				polynomialRing.power(frobenius, 2));
		List<Polynomial<PFE>> testPolynomials = new ArrayList<>();
		testPolynomials.add(testPolynomial);
		testPolynomials.add(polynomialRing.add(frobenius, polynomialRing.power(polynomialRing
				.subtract(polynomialRing.getVarPower(1, 2), polynomialRing.getEmbedding(fp.primitiveRoot())), 2)));
		for (Polynomial<PFE> toTest : testPolynomials) {
			FactorizationResult<Polynomial<PFE>, Polynomial<PFE>> factors = polynomialRing.uniqueFactorization(toTest);
			assertTrue(polynomialRing.isUnit(factors.getUnit()));
			Polynomial<PFE> test = factors.getUnit();
			for (Polynomial<PFE> factor : factors.primeFactors()) {
				test = polynomialRing.multiply(test, polynomialRing.power(factor, factors.multiplicity(factor)));
			}
			assertEquals(toTest, test);
		}
	}
}
