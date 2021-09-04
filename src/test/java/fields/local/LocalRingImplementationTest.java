package fields.local;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.Reals;
import fields.helper.GenericExtensionField;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Field.Extension;
import fields.interfaces.FieldExtension;
import fields.interfaces.LocalField;
import fields.interfaces.LocalRing;
import fields.interfaces.LocalRing.OkutsuType;
import fields.interfaces.LocalRing.TheMontesResult;
import fields.interfaces.Polynomial;
import fields.interfaces.Ring.FactorizationResult;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.PAdicField.PAdicNumber;

class LocalRingImplementationTest {

	@SuppressWarnings("unchecked")
	private <T extends Element<T>, S extends Element<S>> void testHenselLift(LocalField<T, S> field) {
		LocalRing<T, S> ints = field.ringOfIntegers();
		UnivariatePolynomialRing<T> ring = field.getUnivariatePolynomialRing();
		UnivariatePolynomialRing<S> reducedRing = field.reduction().getUnivariatePolynomialRing();
		UnivariatePolynomial<T> xsqp1 = ring.getPolynomial(field.one(), field.zero(), field.one());
		UnivariatePolynomial<S> xsqp1r = ints.reduceUnivariatePolynomial(xsqp1);
		Map<S, Integer> roots = field.reduction().roots(xsqp1r);
		if (roots.size() == 2) {
			Iterator<S> it = roots.keySet().iterator();
			S root1 = it.next();
			S root2 = it.next();
			T lift1 = ints.henselLift(xsqp1, root1, field.getAccuracy());
			T lift2 = ints.henselLift(xsqp1, root2, field.getAccuracy());
			assertEquals(field.zero(), ring.evaluate(xsqp1, lift1));
			assertEquals(field.zero(), ring.evaluate(xsqp1, lift2));
			Polynomial<T> factor1 = ring.subtract(ring.getVar(), ring.getEmbedding(lift1));
			Polynomial<T> factor2 = ring.subtract(ring.getVar(), ring.getEmbedding(lift2));
			assertEquals(xsqp1, ring.multiply(factor1, factor2));
			UnivariatePolynomial<S> factor1r = reducedRing
					.toUnivariate(reducedRing.subtract(reducedRing.getVar(), reducedRing.getEmbedding(root1)));
			UnivariatePolynomial<S> factor2r = reducedRing
					.toUnivariate(reducedRing.subtract(reducedRing.getVar(), reducedRing.getEmbedding(root2)));
			Polynomial<T> factor1lifted = ints.henselLiftFactor(xsqp1, factor1r, field.getAccuracy());
			Polynomial<T> factor2lifted = ints.henselLiftFactor(xsqp1, factor2r, field.getAccuracy());
			assertEquals(xsqp1, ring.multiply(factor1lifted, factor2lifted));
			assertEquals(factor1, factor1lifted);
			assertEquals(factor2, factor2lifted);
		}

	}

	private <T extends Element<T>, S extends Element<S>> void testFactorizationOfLocalField(LocalField<T, S> field) {
		UnivariatePolynomialRing<T> ring = field.ringOfIntegers().getUnivariatePolynomialRing();
		List<UnivariatePolynomial<T>> polynomials = new ArrayList<>();
		for (int i = 0; i < 8; i++) {
			polynomials.add(field.ringOfIntegers().roundUnivariatePolynomial(
					ring.add(ring.getRandomElement(2), ring.getVarPower(3)), field.getAccuracy() / 6));
		}
		int noninfinite = 0;
		for (int i = 0; i < polynomials.size(); i++) {
			for (int j = 0; j < polynomials.size(); j++) {
				if (polynomials.get(i).equals(polynomials.get(j))) {
					continue;
				}
				Map<Polynomial<T>, Integer> factors = new TreeMap<>();
				FactorizationResult<Polynomial<T>> fs = field.factorization(polynomials.get(i));
				for (Polynomial<T> f : fs.primeFactors()) {
					if (factors.containsKey(f)) {
						factors.put(f, factors.get(f) + fs.multiplicity(f));
					} else {
						factors.put(f, fs.multiplicity(f));
					}
				}
				fs = field.factorization(polynomials.get(j));
				for (Polynomial<T> f : fs.primeFactors()) {
					if (factors.containsKey(f)) {
						factors.put(f, factors.get(f) + fs.multiplicity(f));
					} else {
						factors.put(f, fs.multiplicity(f));
					}
				}
				UnivariatePolynomial<T> product = ring.multiply(polynomials.get(i), polynomials.get(j));
				FactorizationResult<Polynomial<T>> factorization = field.factorization(product);
				UnivariatePolynomial<T> test = ring.one();
				for (Polynomial<T> factor : factorization.primeFactors()) {
					Value maxValue = new Value(-10);
					Polynomial<T> max = null;
					for (Polynomial<T> knownFactor : factors.keySet()) {
						Value value = field.ringOfIntegers()
								.valuationOfUnivariatePolynomial(ring.subtract(factor, knownFactor));
						if (value.compareTo(maxValue) > 0) {
							maxValue = maxValue.max(value);
							max = knownFactor;
						}
					}
					if (!maxValue.isInfinite()) {
						noninfinite++;
					}
					if (maxValue.compareTo(new Value(field.getAccuracy() * 3 / 4)) <= 0) {
						System.out.println("Factor " + factor + " not found. Best approximation was " + max
								+ " with value " + maxValue);
					}
					assertTrue(maxValue.compareTo(new Value(field.getAccuracy() * 3 / 4)) > 0);
					int multiplicity = factors.remove(max);
					if (multiplicity > 0) {
						factors.put(max, multiplicity - 1);
					}
//					if (!factors.containsKey(factor)) {
//						System.out.println("Factor not found!");
//						System.out.println("P1:      " + polynomials.get(i));
//						System.out.println("P2:      " + polynomials.get(j));
//						System.out.println("P1 * P2: " + product);
//						System.out.println("Factor:  " + factor);
//						System.out.println("All:     " + factorization.toString());
//						System.out.println(factors.toString());
//					} else if (factors.get(factor) != factorization.getFactors().get(factor)) {
//						System.err.println("Factor multiplicity too low!");
//						System.out.println("P1:      " + polynomials.get(i));
//						System.out.println("P2:      " + polynomials.get(j));
//						System.out.println("P1 * P2: " + product);
//						System.out.println("Factor:  " + factor);
//						System.out.println("All:     " + factorization.toString());
//						System.out.println(factors.toString());
//					}
					// factors.containsKey(factor));
					// assertEquals(factors.get(factor), factorization.getFactors().get(factor));
					for (int k = 0; k < factorization.multiplicity(factor); k++) {
						test = ring.multiply(test, factor);
					}
				}
				test = ring.multiply(test, factorization.getUnit());
				if (!test.equals(product)) {
					System.err.println("i: " + i + " j = " + j);
					System.err.println("Multiplied: " + test);
					System.err.println("Original:   " + product);
					System.err.println();
				}
				assertEquals(test, product);
			}
		}
		if (noninfinite > 0) {
			System.err.println("Non infinite values: " + noninfinite);
		}
	}

	private <T extends Element<T>, S extends Element<S>> void testIntegralBasis(LocalRing<T, S> ring,
			UnivariatePolynomial<T> minimalPolynomial) {
		Field<S> reduction = ring.reduction();
		testIntegralBasis(ring, minimalPolynomial,
				reduction.getExtension(reduction.getUnivariatePolynomialRing().getVar()));
	}

	private <T extends Element<T>, S extends Element<S>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> void testIntegralBasis(
			LocalRing<T, S> ring, UnivariatePolynomial<T> minimalPolynomial,
			Extension<S, R, RE, RFE> reductionExtension) {
		TheMontesResult<T, S, R, RE, RFE> theMontes = ring.theMontesAlgorithm(minimalPolynomial, reductionExtension);
		List<UnivariatePolynomial<T>> integralBasis = ring.triagonalizeIntegralBasis(minimalPolynomial, ring.integralBasis(minimalPolynomial, theMontes, false));
		GenericExtensionField<T> extension = new GenericExtensionField<>(minimalPolynomial, ring.fieldOfFractions());
		for (UnivariatePolynomial<T> b : integralBasis) {
			UnivariatePolynomial<T> mipo = extension.minimalPolynomial(extension.fromPolynomial(b));
			assertEquals(ring.one(), mipo.leadingCoefficient());
			for (int i = 0; i < mipo.degree(); i++) {
				assertTrue(ring.isElement(mipo.univariateCoefficient(i)));
			}
			boolean found = false;
			for (OkutsuType<T, S, R, RE, RFE> type : theMontes.getTypes()) {
				Value value = type.valuation(b);
				if (value.isInfinite()) {
					continue;
				}
				assertTrue(value.value() >= 0);
				if (value.value() < type.ramificationIndex()) {
					found = true;
				}
			}
			assertTrue(found);
		}

	}

	@Test
	void testIntegralBasisUnramified() {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		UnivariatePolynomial<Fraction> minimalPolynomial = polynomials.getPolynomial(q.getInteger(93), q.getInteger(63),
				q.getInteger(58), q.getInteger(39), q.getInteger(14), q.getInteger(3), q.one());
		FactorizationResult<IntE> primes = z
				.uniqueFactorization(polynomials.discriminant(minimalPolynomial).getNumerator());
		for (IntE prime : primes.primeFactors()) {
			LocalRing<Fraction, PFE> zp = z.localize(prime);
			testIntegralBasis(zp, minimalPolynomial);
		}
	}

	@Test
	void testIntegralBasisLarge() {
		Rationals q = Rationals.q();
		LocalRing<Fraction, PFE> z2 = Integers.z().localize(new IntE(2));
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		UnivariatePolynomial<Fraction> minimalPolynomial = polynomials.getPolynomial(
				q.getInteger(BigInteger.valueOf(1125899906842816L)), q.getInteger(1600), q.getInteger(1728),
				q.getInteger(2080), q.getInteger(1328), q.getInteger(992), q.getInteger(544), q.getInteger(240),
				q.getInteger(100), q.getInteger(36), q.getInteger(12), q.getInteger(2), q.one());
		testIntegralBasis(z2, minimalPolynomial);
	}

	@Test
	void testHenselLiftZ2() {
		PAdicField z2 = new PAdicField(BigInteger.valueOf(2), 10);
		testHenselLift(z2);
	}

	@Test
	void testHenselLiftZ2i() {
		PAdicField z2 = new PAdicField(BigInteger.valueOf(2), 10);
		testHenselLift(z2
				.getExtension(
						z2.getUnivariatePolynomialRing().getPolynomial(z2.getInteger(2), z2.getInteger(2), z2.one()))
				.extension());
	}

	@Test
	void testHenselLiftZ2sqrt2() {
		PAdicField z2 = new PAdicField(BigInteger.valueOf(2), 10);
		testHenselLift(
				z2.getExtension(z2.getUnivariatePolynomialRing().getPolynomial(z2.getInteger(2), z2.zero(), z2.one()))
						.extension());
	}

	@Test
	void testHenselLiftZ2F4() {
		PAdicField z2 = new PAdicField(BigInteger.valueOf(2), 10);
		testHenselLift(z2.getExtension(z2.getUnivariatePolynomialRing().getPolynomial(z2.one(), z2.one(), z2.one()))
				.extension());
	}

	@Test
	void testHenselLiftZ3() {
		PAdicField z3 = new PAdicField(BigInteger.valueOf(3), 10);
		testHenselLift(z3);
		testHenselLift(new PAdicField(BigInteger.valueOf(3), 20));
		testHenselLift(z3
				.getExtension(
						z3.getUnivariatePolynomialRing().getPolynomial(z3.getInteger(3), z3.getInteger(3), z3.one()))
				.extension());
		testHenselLift(
				z3.getExtension(z3.getUnivariatePolynomialRing().getPolynomial(z3.getInteger(3), z3.zero(), z3.one()))
						.extension());
		testHenselLift(z3.getExtension(z3.getUnivariatePolynomialRing().getPolynomial(z3.one(), z3.zero(), z3.one()))
				.extension());
	}

	@Test
	void testHenselLiftZ57() {
		PAdicField z5 = new PAdicField(BigInteger.valueOf(5), 10);
		testHenselLift(z5);
		testHenselLift(new PAdicField(BigInteger.valueOf(7), 5));
		testHenselLift(
				z5.getExtension(z5.getUnivariatePolynomialRing().getPolynomial(z5.getInteger(2), z5.zero(), z5.one()))
						.extension());
		testHenselLift(
				z5.getExtension(z5.getUnivariatePolynomialRing().getPolynomial(z5.getInteger(5), z5.zero(), z5.one()))
						.extension());
	}

	@Test
	void testHenselLiftPowerSeries2() {
		testHenselLift(new FormalPowerSeries<>(PrimeField.getPrimeField(2), 10));
	}

	@Test
	void testHenselLiftPowerSeries3() {
		testHenselLift(new FormalPowerSeries<>(PrimeField.getPrimeField(3), 10));
	}

	@Test
	void testHenselLiftPowerSeries5() {
		testHenselLift(new FormalPowerSeries<>(PrimeField.getPrimeField(5), 10));
	}

	@Test
	void testHenselLiftPowerSeries8() {
		testHenselLift(new FormalPowerSeries<>(FiniteField.getFiniteField(4), 10));
		testHenselLift(new FormalPowerSeries<>(FiniteField.getFiniteField(8), 10));
	}

	@Test
	void testHenselLiftPowerSeries9() {
		testHenselLift(new FormalPowerSeries<>(FiniteField.getFiniteField(9), 10));
	}

	@Test
	void testHenselLiftPowerSeries27() {
		testHenselLift(new FormalPowerSeries<>(FiniteField.getFiniteField(27), 10));
	}

	@Test
	void testHenselLiftPowerSeries25() {
		testHenselLift(new FormalPowerSeries<>(FiniteField.getFiniteField(25), 5));
	}

	@Test
	void testHenselLiftPowerSeriesQ() {
		testHenselLift(new FormalPowerSeries<>(Rationals.q(), 10));
	}

	@Test
	void testHenselLiftPowerSeriesR() {
		testHenselLift(new FormalPowerSeries<>(Reals.r(1024), 10));
	}

	@Test
	void testFactorizationZ2() {
		fail();
		PAdicField z2 = new PAdicField(BigInteger.valueOf(2), 40);
		UnivariatePolynomialRing<PAdicNumber> ring = z2.getUnivariatePolynomialRing();
		z2.ringOfIntegers()
				.theMontesAlgorithm(ring.getPolynomial(z2.getInteger(0), z2.getInteger(1431 * 2), z2.getInteger(2461),
						z2.getInteger(645 * 4), z2.getInteger(29 * 32), z2.getInteger(15 * 4), z2.getInteger(1)));
		z2.ringOfIntegers().theMontesAlgorithm(ring.getPolynomial(z2.getInteger(305), z2.getInteger(3375),
				z2.getInteger(3493), z2.getInteger(359 * 2), z2.getInteger(483), z2.getInteger(61), z2.getInteger(1)));
		assertEquals(1,
				z2.factorization(ring.getPolynomial(z2.one(), z2.getInteger(6), z2.one())).primeFactors().size());
		assertEquals(1, z2.factorization(ring.getPolynomial(z2.getInteger(-1), z2.getInteger(4), z2.one()))
				.primeFactors().size());
		assertEquals(2, z2.factorization(ring.getPolynomial(z2.getInteger(37), z2.getInteger(6), z2.one()))
				.primeFactors().size());
		assertEquals(2,
				z2.factorization(ring.getPolynomial(z2.getInteger(4), z2.one(), z2.one())).primeFactors().size());
		FactorizationResult<Polynomial<PAdicNumber>> factors = z2.factorization(ring.getPolynomial(z2.getInteger(4),
				z2.getInteger(2), z2.getInteger(5), z2.getInteger(2), z2.getInteger(1)));
		assertNotEquals(0, factors.multiplicity(ring.getPolynomial(z2.one(), z2.zero(), z2.one())));
		assertNotEquals(0, factors.multiplicity(ring.getPolynomial(z2.getInteger(4), z2.getInteger(2), z2.one())));
		assertEquals(1, z2.factorization(
				ring.getPolynomial(z2.getInteger(5), z2.getInteger(4), z2.getInteger(3), z2.getInteger(2), z2.one()))
				.primeFactors().size());
//		LocalFieldExtension<PAdicNumber, PFE, FFE> z2sqrt2 = z2.getExtension(ring.getPolynomial(z2.getInteger(-2), z2.zero(), z2.one())).extension();
//		UnivariatePolynomial<PAdicNumber> f1 = ring.getPolynomial(z2.getInteger(16), z2.getInteger(6), z2.getInteger(4), z2.one());
//		UnivariatePolynomial<Ext<PAdicNumber, PFE, FFE>> p1 = z2sqrt2.getUnivariatePolynomialRing().getEmbedding(f1, z2sqrt2.getPrimeEmbeddingMap());
//		List<Polynomial<Ext<PAdicNumber, PFE, FFE>>> factors2 = z2sqrt2.factorization(p1);
//		Polynomial<Ext<PAdicNumber, PFE, FFE>> product = z2sqrt2.getUnivariatePolynomialRing().one();
//		for (Polynomial<Ext<PAdicNumber, PFE, FFE>> factor : factors2) {
//			System.out.println(factor);
//			product = z2sqrt2.getUnivariatePolynomialRing().multiply(product, factor);
//		}
//		System.out.println(product);
		assertEquals(2,
				z2.factorization(ring.getPolynomial(z2.getInteger(16), z2.getInteger(6), z2.getInteger(4), z2.one()))
						.primeFactors().size());
		factors = z2.factorization(ring.getPolynomial(z2.getInteger(4), z2.getInteger(4), z2.one()));
		assertEquals(1, factors.primeFactors().size());
		assertEquals(2, factors.multiplicity(factors.firstPrimeFactor()));
//		factors = z2.withAccuracy(200).factorization(ring.getPolynomial(z2.getInteger(220), z2.getInteger(1776), z2.getInteger(3792),
//				z2.getInteger(3188), z2.getInteger(940), z2.getInteger(58), z2.getInteger(1)));
//		assertEquals(2, factors.size());
//		assertTrue(factors.get(0)
//				.equals(ring.getPolynomial(z2.getInteger(10), z2.getInteger(58), z2.getInteger(26), z2.one()))
//				|| factors.get(1)
//						.equals(ring.getPolynomial(z2.getInteger(22), z2.getInteger(50), z2.getInteger(32), z2.one())));
//		assertTrue(factors.get(1)
//				.equals(ring.getPolynomial(z2.getInteger(10), z2.getInteger(58), z2.getInteger(26), z2.one()))
//				|| factors.get(0)
//						.equals(ring.getPolynomial(z2.getInteger(22), z2.getInteger(50), z2.getInteger(32), z2.one())));
		testFactorizationOfLocalField(z2);
	}

	@Test
	void testFactorizationZ2i() {
		PAdicField z2 = new PAdicField(BigInteger.valueOf(2), 30);
		testFactorizationOfLocalField(z2
				.getExtension(
						z2.getUnivariatePolynomialRing().getPolynomial(z2.getInteger(2), z2.getInteger(2), z2.one()))
				.extension());
	}

	@Test
	void testFactorizationZ2sqrt2() {
		PAdicField z2 = new PAdicField(BigInteger.valueOf(2), 30);
		testFactorizationOfLocalField(
				z2.getExtension(z2.getUnivariatePolynomialRing().getPolynomial(z2.getInteger(2), z2.zero(), z2.one()))
						.extension());
	}

	@Test
	void testFactorizationZ2F4() {
		PAdicField z2 = new PAdicField(BigInteger.valueOf(2), 30);
		testFactorizationOfLocalField(
				z2.getExtension(z2.getUnivariatePolynomialRing().getPolynomial(z2.one(), z2.one(), z2.one()))
						.extension());
	}

	@Test
	void testFactorizationZ3() {
		fail();

		PAdicField z3 = new PAdicField(BigInteger.valueOf(3), 20);
		UnivariatePolynomialRing<PAdicNumber> ring = z3.getUnivariatePolynomialRing();
		z3.ringOfIntegers()
				.theMontesAlgorithm(ring.getPolynomial(z3.getInteger(657655148), z3.getInteger(1308416681),
						z3.getInteger(1744039174 * 81), z3.getInteger(184016488 * 3), z3.getInteger(128910175),
						z3.getInteger(28979416 * 3), z3.getInteger(1)));
		testFactorizationOfLocalField(z3);
		testFactorizationOfLocalField(new PAdicField(BigInteger.valueOf(3), 40));
	}

	@Test
	void testFactorizationZ3Extensions() {
		PAdicField z3 = new PAdicField(BigInteger.valueOf(3), 30);
		testFactorizationOfLocalField(z3
				.getExtension(
						z3.getUnivariatePolynomialRing().getPolynomial(z3.getInteger(3), z3.getInteger(3), z3.one()))
				.extension());
		testFactorizationOfLocalField(
				z3.getExtension(z3.getUnivariatePolynomialRing().getPolynomial(z3.getInteger(3), z3.zero(), z3.one()))
						.extension());
		testFactorizationOfLocalField(
				z3.getExtension(z3.getUnivariatePolynomialRing().getPolynomial(z3.one(), z3.zero(), z3.one()))
						.extension());
	}

	@Test
	void testFactorizationZ57() {
		fail();

		PAdicField z5 = new PAdicField(BigInteger.valueOf(5), 30);
		testFactorizationOfLocalField(z5);
		testFactorizationOfLocalField(new PAdicField(BigInteger.valueOf(7), 20));
		testFactorizationOfLocalField(
				z5.getExtension(z5.getUnivariatePolynomialRing().getPolynomial(z5.getInteger(2), z5.zero(), z5.one()))
						.extension());
		testFactorizationOfLocalField(
				z5.getExtension(z5.getUnivariatePolynomialRing().getPolynomial(z5.getInteger(5), z5.zero(), z5.one()))
						.extension());
	}

	@Test
	void testFactorizationPowerSeries2() {
		testFactorizationOfLocalField(new FormalPowerSeries<>(PrimeField.getPrimeField(2), 10));
	}

	@Test
	void testFactorizationPowerSeries3() {
		testFactorizationOfLocalField(new FormalPowerSeries<>(PrimeField.getPrimeField(3), 10));
	}

	@Test
	void testFactorizationPowerSeries5() {
		testFactorizationOfLocalField(new FormalPowerSeries<>(PrimeField.getPrimeField(5), 10));
	}

	@Test
	void testFactorizationPowerSeries8() {
		testFactorizationOfLocalField(new FormalPowerSeries<>(FiniteField.getFiniteField(4), 10));
		testFactorizationOfLocalField(new FormalPowerSeries<>(FiniteField.getFiniteField(8), 10));
	}

	@Test
	void testFactorizationPowerSeries9() {
		testFactorizationOfLocalField(new FormalPowerSeries<>(FiniteField.getFiniteField(9), 10));
	}

	@Test
	void testFactorizationPowerSeries27() {
		testFactorizationOfLocalField(new FormalPowerSeries<>(FiniteField.getFiniteField(27), 20));
	}

	@Test
	void testFactorizationPowerSeries25() {
		testFactorizationOfLocalField(new FormalPowerSeries<>(FiniteField.getFiniteField(25), 20));
	}

	@Test
	void testFactorizationPowerSeriesQ() {
		testFactorizationOfLocalField(new FormalPowerSeries<>(Rationals.q(), 10));
	}

	@Test
	void testFactorizationPowerSeriesR() {
		testFactorizationOfLocalField(new FormalPowerSeries<>(Reals.r(1024), 10));
	}

}
