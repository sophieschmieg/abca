package fields.numberfields;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.floatingpoint.Complex;
import fields.floatingpoint.Complex.ComplexNumber;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.AlgebraicExtensionElement;
import fields.interfaces.DiscreteValuationRing;
import fields.interfaces.DiscreteValuationRing.OkutsuType;
import fields.interfaces.DiscreteValuationRing.TheMontesResult;
import fields.interfaces.Element;
import fields.interfaces.Field.Extension;
import fields.interfaces.FieldExtension;
import fields.interfaces.Ideal;
import fields.interfaces.Ring.FactorizationResult;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.PAdicField;
import fields.local.Value;
import fields.numberfields.CompletedNumberField.Ext;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;

class NumberFieldTest {

	@Test
	void findIntegerMinimalPolynomialTest() {
		Rationals q = Rationals.q();
		Complex c = Complex.c(128);
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		NumberField gauss = NumberField.getNumberField(polynomials.getPolynomial(q.one(), q.zero(), q.one()));
		EmbeddedNumberField<ComplexNumber, Complex> embedded = gauss.complexEmbeddings().get(0);
		for (int tc = 0; tc < 20; tc++) {
			int degree = new Random().nextInt(10) + 1;
			List<NFE> rng = new ArrayList<>();
			for (int i = 0; i <= degree; i++) {
				rng.add(gauss.maximalOrder().getRandomElement());
			}
			rng.add(gauss.one());
			UnivariatePolynomial<NFE> polynomial = gauss.getUnivariatePolynomialRing().getPolynomial(rng);
			List<ComplexNumber> conjugates = new ArrayList<>();
			conjugates.addAll(
					c.roots(c.getUnivariatePolynomialRing().getEmbedding(polynomial, embedded.embeddingMap()))
							.keySet());
			System.out.println(polynomial);
			assertEquals(polynomial, gauss.integerMinimalPolynomial(c, conjugates));
		}
	}

	@Test
	void fixTheStupidValuationBugTest() {
		Integers z = Integers.z();
		DiscreteValuationRing<Fraction, PFE> z2 = z.localize(z.getInteger(2));
		PrimeField f2 = PrimeField.getPrimeField(2);
		Extension<PFE, PFE, FFE, FiniteField> trivial = f2.getExtension(f2.getUnivariatePolynomialRing().getVar());
		UnivariatePolynomialRing<Fraction> polynomials = z2.getUnivariatePolynomialRing();
		UnivariatePolynomial<Fraction> polynomial1 = polynomials.getPolynomial(z2.getInteger(23), z2.zero(), z2.one());
		UnivariatePolynomial<Fraction> polynomial2 = polynomials.getPolynomial(z2.getInteger(6), z2.one(), z2.one());
		TheMontesResult<Fraction, PFE, PFE, FFE, FiniteField> theMontes1 = z2.theMontesAlgorithm(polynomial1, trivial);
		TheMontesResult<Fraction, PFE, PFE, FFE, FiniteField> theMontes2 = z2.theMontesAlgorithm(polynomial2, trivial);
		UnivariatePolynomial<Fraction> test1 = polynomials.getPolynomial(z2.getInteger(-3), z2.one());
		UnivariatePolynomial<Fraction> test2 = polynomials.getPolynomial(z2.getInteger(3), z2.one());
		for (OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> type : theMontes1.getTypes()) {
			type.valuation(test1);
			type.valuation(test2);
		}
		for (OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> type : theMontes2.getTypes()) {
			type.valuation(test1);
			type.valuation(test2);
		}
	}

	@Test
	void padicEmbeddingTest() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		NumberField gauss = NumberField.getNumberField(polynomials.getPolynomial(q.one(), q.zero(), q.one()));
		PAdicField base = new PAdicField(5, 20);
		List<EmbeddedNumberField<Ext, CompletedNumberField>> padicEmbeddings = gauss.padicEmbeddings(base);
		NFE random = gauss.maximalOrder().getRandomElement();
		System.out.println(random);
		Ext embedded = padicEmbeddings.get(0).embedding(random);
		System.out.println(embedded);
		NFE recovered = padicEmbeddings.get(0).roundToInteger(embedded);
		System.out.println(recovered);

		NumberField eisenstein = NumberField
				.getNumberField(polynomials.getPolynomial(q.getInteger(3), q.zero(), q.one()));
		List<EmbeddedNumberField<Ext, CompletedNumberField>> padicEmbeddingsEisenstein = eisenstein
				.padicEmbeddings(base);
		NFE randomEisenstein = eisenstein.maximalOrder().getRandomElement();
		System.out.println(randomEisenstein);
		Ext embeddedEisenstein = padicEmbeddingsEisenstein.get(0).embedding(randomEisenstein);
		System.out.println(embeddedEisenstein);
		NFE recoveredEisenstein = padicEmbeddingsEisenstein.get(0).roundToInteger(embeddedEisenstein);
		System.out.println(recoveredEisenstein);
	}

	@Test
	void numberFieldNoPowerBasisTest() {
		Rationals q = Rationals.q();
		Integers z = Integers.z();
		UnivariatePolynomialRing<Fraction> rationalPolynomialRing = q.getUnivariatePolynomialRing();
		NumberField nf = NumberField.getNumberField(
				rationalPolynomialRing.getPolynomial(q.getInteger(-8), q.getInteger(-2), q.getInteger(-1), q.one()));
		NumberFieldIntegers order = nf.maximalOrder();
		List<NumberFieldIdeal> ideals = order.idealsOver(z.getIdeal(z.getInteger(2)));
		for (NumberFieldIdeal primeIdeal : ideals) {
			NumberFieldIdeal newlyGenerated = order.getIdeal(primeIdeal.generators());
			assertEquals(newlyGenerated, primeIdeal);
			NumberFieldIdeal primePower = order.power(primeIdeal, 2);
			FactorizationResult<Ideal<NFE>, Ideal<NFE>> factorization = order.idealFactorization(primePower);
			assertEquals(1, factorization.primeFactors().size());
			assertEquals(2, factorization.multiplicity(primeIdeal));
			NumberFieldIdeal newlyGeneratedPower = order.getIdeal(primePower.generators());
			assertEquals(newlyGeneratedPower, primePower);
		}
	}

	@Test
	void testNoPowerBasisIdeals() {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		NumberField nf = NumberField.getNumberField(
				polynomials.getPolynomial(q.getInteger(-8), q.getInteger(-2), q.getInteger(-1), q.one()));
		NumberFieldIntegers order = nf.maximalOrder();
		List<IntE> primes = new ArrayList<>();
		primes.add(new IntE(2));
		primes.add(new IntE(3));
		primes.add(new IntE(5));
		primes.add(new IntE(7));
		primes.add(new IntE(31));
		primes.add(new IntE(1741));
		primes.add(new IntE(257));
		primes.add(new IntE(503));
		primes.add(new IntE(65537));
		List<NFE> toGenerate = new ArrayList<>();
		toGenerate.add(nf.one());
		toGenerate.add(nf.alpha());
		toGenerate.add(nf.getInteger(2));
		toGenerate.add(nf.getInteger(3));
		toGenerate.add(nf.getInteger(4));
		toGenerate.add(nf.getInteger(5));
		toGenerate.add(nf.getInteger(6));
		toGenerate.add(nf.add(nf.alpha(), nf.one()));
		toGenerate.add(nf.add(nf.alpha(), nf.getInteger(2)));
		toGenerate.add(nf.add(nf.alpha(), nf.getInteger(3)));
		toGenerate.add(nf.add(nf.alpha(), nf.getInteger(4)));
		toGenerate.add(nf.add(nf.alpha(), nf.getInteger(5)));
		toGenerate.add(nf.add(nf.multiply(-3, nf.alpha()), nf.getInteger(3)));
		for (IntE prime : primes) {
			List<NumberFieldIdeal> ideals = order.idealsOver(z.getIdeal(Collections.singletonList(prime)));
			for (Ideal<NFE> ideal : ideals) {
				for (NFE g : toGenerate) {
					List<NFE> coefficients = ideal.generate(g);
					NFE reconstructed = ideal.residue(g);
					if (ideal.contains(g)) {
						assertEquals(order.zero(), reconstructed);
					}
					for (int i = 0; i < coefficients.size(); i++) {
						assertTrue(nf.isInteger(coefficients.get(i)));
						reconstructed = nf.add(reconstructed,
								nf.multiply(coefficients.get(i), ideal.generators().get(i)));
					}
					assertEquals(g, reconstructed);
				}
			}
		}
	}

	@Test
	void testNoPowerBasisIdealFactorization() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		NumberField nf = NumberField.getNumberField(
				polynomials.getPolynomial(q.getInteger(-8), q.getInteger(-2), q.getInteger(-1), q.one()));
		NumberFieldIntegers order = nf.maximalOrder();
		List<NFE> toFactorize = new ArrayList<>();
		toFactorize.add(nf.one());
		toFactorize.add(nf.getInteger(2));
		toFactorize.add(nf.getInteger(3));
		toFactorize.add(nf.getInteger(4));
		toFactorize.add(nf.getInteger(5));
		toFactorize.add(nf.getInteger(6));
		toFactorize.add(nf.alpha());
		toFactorize.add(nf.add(nf.alpha(), nf.one()));
		toFactorize.add(nf.add(nf.alpha(), nf.getInteger(2)));
		toFactorize.add(nf.add(nf.alpha(), nf.getInteger(3)));
		toFactorize.add(nf.add(nf.alpha(), nf.getInteger(4)));
		toFactorize.add(nf.add(nf.alpha(), nf.getInteger(5)));
		toFactorize.add(nf.add(nf.multiply(-3, nf.alpha()), nf.getInteger(3)));
		for (NFE factorize : toFactorize) {
			Ideal<NFE> ideal = order.getIdeal(factorize);
			System.out.println(factorize + ": " + ideal);
			for (NFE g : toFactorize) {
				System.out.println(g);
				List<NFE> coefficients = ideal.generate(g);
				NFE reconstructed = ideal.residue(g);
				if (ideal.contains(g)) {
					assertEquals(order.zero(), reconstructed);
				}
				for (int i = 0; i < coefficients.size(); i++) {
					assertTrue(nf.isInteger(coefficients.get(i)));
					reconstructed = nf.add(reconstructed, nf.multiply(coefficients.get(i), ideal.generators().get(i)));
				}
				assertEquals(g, reconstructed);
			}
		}
	}

	private <T extends Element<T>, S extends Element<S>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> boolean checkUniformizer(
			DiscreteValuationRing<T, S> ring, OkutsuType<T, S, R, RE, RFE> type) {
		int[] logUniformizer = type.uniformizer();
		OkutsuType<T, S, R, RE, RFE> it = type;
		int value = 0;
		while (it.level() > 0) {
			value += type.valuation(it.phi()).value() * logUniformizer[it.level()];
			it = it.previousLevel();
		}
		value += type.valuation(ring.getUnivariatePolynomialRing().getEmbedding(ring.uniformizer())).value()
				* logUniformizer[0];
		return value == 1;
	}

	@Test
	void testTestvector() {
		Rationals q = Rationals.q();
		Integers z = Integers.z();
		UnivariatePolynomialRing<Fraction> ring = q.getUnivariatePolynomialRing();
		UnivariatePolynomial<Fraction> minimalPolynomial = ring.getPolynomial(
				q.getInteger(BigInteger.valueOf(1125899906842816L)), q.getInteger(1600), q.getInteger(1728),
				q.getInteger(2080), q.getInteger(1328), q.getInteger(992), q.getInteger(544), q.getInteger(240),
				q.getInteger(100), q.getInteger(36), q.getInteger(12), q.getInteger(2), q.one());
		DiscreteValuationRing<Fraction, PFE> z2 = Integers.z().localize(BigInteger.TWO);
		PrimeField f2 = PrimeField.getPrimeField(2);
		UnivariatePolynomialRing<PFE> reductionRing = f2.getUnivariatePolynomialRing();
		FiniteField f2e = f2.getExtension(reductionRing.getVar()).extension();
		UnivariatePolynomialRing<FFE> reductionRingE = f2e.getUnivariatePolynomialRing();
		TheMontesResult<Fraction, PFE, PFE, FFE, FiniteField> result = z2.theMontesAlgorithm(minimalPolynomial,
				f2.getExtension(reductionRing.getVar()));
		assertEquals(2, result.getTypes().size());
		OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> firstType = result.getTypes().get(0);
		OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> secondType = result.getTypes().get(1);

		if (firstType.level() == 4) {
			OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> tmp = secondType;
			secondType = firstType;
			firstType = tmp;
		}
		firstType = firstType.previousLevel();
		secondType = secondType.previousLevel();

		assertEquals(2, firstType.level());
		assertEquals(ring.getPolynomial(q.getInteger(12), q.getInteger(4), q.getInteger(4), q.getInteger(2), q.one()),
				firstType.representative());
		assertEquals(ring.getPolynomial(q.getInteger(2), q.zero(), q.one()), firstType.phi());
		assertEquals(q.getFraction(z.one(), z.one()), firstType.lambda());
		assertEquals(reductionRingE.getPolynomial(f2e.one(), f2e.one(), f2e.one()), firstType.psi());
		assertEquals(new Value(2), firstType.value());
		assertEquals(2, firstType.ramificationIndex());
		assertEquals(2, firstType.residueDegree());
		assertTrue(firstType.complete());
		assertTrue(checkUniformizer(z2, firstType));
		assertEquals(BigInteger.valueOf(4), firstType.reduction().extension().getNumberOfElements());

		firstType = firstType.previousLevel();
		assertEquals(1, firstType.level());
		assertEquals(ring.getPolynomial(q.getInteger(2), q.zero(), q.one()), firstType.representative());
		assertEquals(ring.getVar(), firstType.phi());
		assertEquals(q.getFraction(z.one(), z.getInteger(2)), firstType.lambda());
		assertEquals(reductionRingE.getPolynomial(f2e.one(), f2e.one()), firstType.psi());
		assertEquals(Value.ZERO, firstType.value());
		assertEquals(2, firstType.ramificationIndex());
		assertEquals(1, firstType.residueDegree());
		assertFalse(firstType.complete());
		assertTrue(checkUniformizer(z2, firstType));
		assertEquals(BigInteger.valueOf(2), firstType.reduction().extension().getNumberOfElements());

		firstType = firstType.previousLevel();
		assertEquals(0, firstType.level());
		assertEquals(ring.getVar(), firstType.representative());
		assertEquals(reductionRingE.getVar(), firstType.psi());
		assertEquals(1, firstType.ramificationIndex());
		assertEquals(1, firstType.residueDegree());
		assertFalse(firstType.complete());
		assertTrue(checkUniformizer(z2, firstType));
		assertEquals(BigInteger.valueOf(2), firstType.reduction().extension().getNumberOfElements());

		assertEquals(3, secondType.level());
		assertEquals(
				ring.getPolynomial(q.getInteger(16), q.getInteger(128), q.getInteger(96), q.getInteger(96),
						q.getInteger(24), q.getInteger(16), q.getInteger(8), q.getInteger(0), q.one()),
				secondType.representative());
		assertEquals(ring.getPolynomial(q.getInteger(4), q.getInteger(8), q.getInteger(4), q.getInteger(0), q.one()),
				secondType.phi());
		assertEquals(q.getFraction(z.one(), z.getInteger(2)), secondType.lambda());
		assertEquals(reductionRingE.getPolynomial(f2e.one(), f2e.one()), secondType.psi());
		assertEquals(new Value(14), secondType.value());
		assertEquals(8, secondType.ramificationIndex());
		assertEquals(1, secondType.residueDegree());
		assertTrue(secondType.complete());
		assertTrue(checkUniformizer(z2, secondType));
		assertEquals(BigInteger.valueOf(2), secondType.reduction().extension().getNumberOfElements());

		secondType = secondType.previousLevel();
		assertEquals(2, secondType.level());
		assertEquals(ring.getPolynomial(q.getInteger(4), q.getInteger(8), q.getInteger(4), q.getInteger(0), q.one()),
				secondType.representative());
		assertEquals(ring.getPolynomial(q.getInteger(2), q.getInteger(0), q.one()), secondType.phi());
		assertEquals(q.getFraction(z.getInteger(3), z.getInteger(2)), secondType.lambda());
		assertEquals(reductionRingE.getPolynomial(f2e.one(), f2e.one()), secondType.psi());
		assertEquals(new Value(2), secondType.value());
		assertEquals(4, secondType.ramificationIndex());
		assertEquals(1, secondType.residueDegree());
		assertFalse(secondType.complete());
		assertTrue(checkUniformizer(z2, secondType));
		assertEquals(BigInteger.valueOf(2), secondType.reduction().extension().getNumberOfElements());

		secondType = secondType.previousLevel();
		assertEquals(1, secondType.level());
		assertEquals(ring.getPolynomial(q.getInteger(2), q.getInteger(0), q.one()), secondType.representative());
		assertEquals(ring.getVar(), secondType.phi());
		assertEquals(q.getFraction(z.one(), z.getInteger(2)), secondType.lambda());
		assertEquals(reductionRingE.getPolynomial(f2e.one(), f2e.one()), secondType.psi());
		assertEquals(Value.ZERO, secondType.value());
		assertEquals(2, secondType.ramificationIndex());
		assertEquals(1, secondType.residueDegree());
		assertFalse(secondType.complete());
		assertTrue(checkUniformizer(z2, secondType));
		assertEquals(BigInteger.valueOf(2), secondType.reduction().extension().getNumberOfElements());

		secondType = secondType.previousLevel();
		assertEquals(0, secondType.level());
		assertEquals(ring.getVar(), secondType.representative());
		assertEquals(reductionRingE.getVar(), secondType.psi());
		assertEquals(1, secondType.ramificationIndex());
		assertEquals(1, secondType.residueDegree());
		assertFalse(secondType.complete());
		assertTrue(checkUniformizer(z2, secondType));
		assertEquals(BigInteger.valueOf(2), secondType.reduction().extension().getNumberOfElements());
	}

	private <T extends Element<T>, S extends Element<S>, R extends Element<R>, RE extends AlgebraicExtensionElement<R, RE>, RFE extends FieldExtension<R, RE, RFE>> boolean testType(
			OkutsuType<T, S, R, RE, RFE> type, DiscreteValuationRing<T, S> ring) {
		for (int tc = 0; tc < 10; tc++) {
			RFE reduction = type.reduction().extension();
			RE rng;
			do {
				rng = reduction.getRandomElement();
			} while (rng.equals(reduction.zero()));
			int value = new Random().nextInt(6) - 3 + 6;
			assertEquals(type.ramificationIndex(),
					type.valuation(ring.getUnivariatePolynomialRing().getEmbedding(ring.uniformizer())).value());
			UnivariatePolynomial<T> lifted = type.lift(rng, value);
			assertEquals(value, type.valuation(lifted).value());
			assertEquals(rng, type.reduce(lifted));
		}
		return true;
	}

	// @Test
	void testInvertInType() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> ring = q.getUnivariatePolynomialRing();
		UnivariatePolynomial<Fraction> minimalPolynomial = ring.getPolynomial(
				q.getInteger(BigInteger.valueOf(1125899906842816L)), q.getInteger(1600), q.getInteger(1728),
				q.getInteger(2080), q.getInteger(1328), q.getInteger(992), q.getInteger(544), q.getInteger(240),
				q.getInteger(100), q.getInteger(36), q.getInteger(12), q.getInteger(2), q.one());
		DiscreteValuationRing<Fraction, PFE> z2 = Integers.z().localize(BigInteger.TWO);
		PrimeField f2 = PrimeField.getPrimeField(2);
		UnivariatePolynomialRing<PFE> reductionRing = f2.getUnivariatePolynomialRing();
		TheMontesResult<Fraction, PFE, PFE, FFE, FiniteField> result = z2.theMontesAlgorithm(minimalPolynomial,
				f2.getExtension(reductionRing.getVar()));
		assertEquals(2, result.getTypes().size());
		OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> firstType = result.getTypes().get(0);
		OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> secondType = result.getTypes().get(1);

		if (firstType.level() == 2) {
			OkutsuType<Fraction, PFE, PFE, FFE, FiniteField> tmp = secondType;
			secondType = firstType;
			firstType = tmp;
		}
		z2.singleFactorLifting(secondType, 100);
		// secondType.lift(secondType.reduction().extension().alpha(), 4);
		testType(firstType, z2);
		testType(secondType, z2);
		for (int tc = 0; tc < 10; tc++) {
			UnivariatePolynomial<Fraction> test;
			do {
				List<Fraction> c = new ArrayList<>();
				for (int i = 0; i < firstType.representative().degree(); i++) {
					c.add(z2.getRandomElement());
				}
				test = ring.getPolynomial(c);
			} while (firstType.valuation(test).compareTo(Value.ZERO) > 0);
			UnivariatePolynomial<Fraction> testInverse = z2.invertInType(test, firstType, 20);
			assertTrue(
					firstType.valuation(ring.toUnivariate(ring.subtract(ring.multiply(test, testInverse), ring.one())))
							.compareTo(new Value(20)) > 0);
		}
		for (int tc = 0; tc < 10; tc++) {
			UnivariatePolynomial<Fraction> test;
			do {
				List<Fraction> c = new ArrayList<>();
				for (int i = 0; i < secondType.representative().degree(); i++) {
					c.add(z2.getRandomElement());
				}
				test = ring.getPolynomial(c);
			} while (secondType.valuation(test).compareTo(Value.ZERO) > 0);
			UnivariatePolynomial<Fraction> testInverse = z2.invertInType(test, secondType, 20);
			assertTrue(
					secondType.valuation(ring.toUnivariate(ring.subtract(ring.multiply(test, testInverse), ring.one())))
							.compareTo(new Value(20)) > 0);
		}
	}

	private long lastTime;

	private void printTime(String name, long current, long start) {
		long diff = current - start;
		double seconds = (double) diff / 1000.0;
		long diffLast = current - lastTime;
		double secondsLast = (double) diffLast / 1000.0;
		lastTime = current;
		System.out.printf("%1.3f, %1.3f: %s\n", secondsLast, seconds, name);
	}

	@Test
	void testHiddenFieldIdeals() {
		long start = System.currentTimeMillis();
		lastTime = start;
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		List<UnivariatePolynomial<Fraction>> minimalPolynomials = new ArrayList<>();
		// minimalPolynomials.add(
		// polynomials.getPolynomial(q.getInteger(7), q.getInteger(0), q.getInteger(5),
		// q.getInteger(0), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(93), q.getInteger(63), q.getInteger(58),
				q.getInteger(39), q.getInteger(14), q.getInteger(3), q.one()));
		printTime("Setup", System.currentTimeMillis(), start);
		for (UnivariatePolynomial<Fraction> minimalPolynomial : minimalPolynomials) {
			NumberField nf = NumberField.getNumberField(minimalPolynomial);
			// System.out.println(nf);
			NumberFieldIntegers order = nf.maximalOrder();
			List<IntE> primes = new ArrayList<>();
			primes.add(new IntE(2));
//			primes.add(new IntE(3));
//			primes.add(new IntE(5));
//			primes.add(new IntE(7));
			primes.add(new IntE(31));
			primes.add(new IntE(1741));
			primes.add(new IntE(257));
//			primes.add(new IntE(65537));
			List<NFE> toGenerate = new ArrayList<>();
			toGenerate.add(nf.one());
			toGenerate.add(nf.alpha());
			toGenerate.add(nf.add(nf.alpha(), nf.one()));
			toGenerate.add(nf.add(nf.alpha(), nf.getInteger(2)));
			toGenerate.add(nf.add(nf.alpha(), nf.getInteger(3)));
			toGenerate.add(nf.add(nf.alpha(), nf.getInteger(4)));
			toGenerate.add(nf.add(nf.alpha(), nf.getInteger(5)));
			toGenerate.add(nf.add(nf.multiply(-3, nf.alpha()), nf.getInteger(3)));
			printTime("Setup " + nf, System.currentTimeMillis(), start);
			for (IntE prime : primes) {
				printTime("Start prime " + prime, System.currentTimeMillis(), start);
				List<NumberFieldIdeal> ideals = order.idealsOver(z.getIdeal(prime));
				printTime("Ideals prime " + prime, System.currentTimeMillis(), start);
				for (Ideal<NFE> ideal : ideals) {
					// printTime("Start ideal " + ideal, System.currentTimeMillis(), start);
					// System.out.println(ideal);
					for (NFE g : toGenerate) {
						List<NFE> coefficients = ideal.generate(g);
						NFE reconstructed = ideal.residue(g);
						// System.out.print(g + " = " + reconstructed);
						for (int i = 0; i < coefficients.size(); i++) {
							// System.out.print(" + (" + coefficients.get(i) + ")*(" +
							// ideal.generators().get(i) + ")");
							assertTrue(nf.isInteger(coefficients.get(i)));
							reconstructed = nf.add(reconstructed,
									nf.multiply(coefficients.get(i), ideal.generators().get(i)));
						}
						// System.out.println();
						assertEquals(g, reconstructed);
					}
					// printTime("End ideal " + ideal, System.currentTimeMillis(), start);
				}
			}
			printTime("End " + nf, System.currentTimeMillis(), start);
		}
		printTime("End", System.currentTimeMillis(), start);
	}

	@Test
	void testQuadraticFieldIdeals() {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		List<UnivariatePolynomial<Fraction>> minimalPolynomials = new ArrayList<>();
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(-6), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(-3), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(-2), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(1), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(2), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(3), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(5), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(6), q.zero(), q.one()));
		for (UnivariatePolynomial<Fraction> minimalPolynomial : minimalPolynomials) {
			NumberField nf = NumberField.getNumberField(minimalPolynomial);
			NumberFieldIntegers order = nf.maximalOrder();
			List<IntE> primes = new ArrayList<>();
			primes.add(new IntE(2));
			primes.add(new IntE(3));
			primes.add(new IntE(5));
			primes.add(new IntE(7));
			List<NFE> toGenerate = new ArrayList<>();
			toGenerate.add(nf.one());
			toGenerate.add(nf.alpha());
			toGenerate.add(nf.getInteger(2));
			toGenerate.add(nf.getInteger(3));
			toGenerate.add(nf.getInteger(4));
			toGenerate.add(nf.getInteger(5));
			toGenerate.add(nf.getInteger(6));
			toGenerate.add(nf.add(nf.alpha(), nf.one()));
			toGenerate.add(nf.add(nf.alpha(), nf.getInteger(2)));
			toGenerate.add(nf.add(nf.alpha(), nf.getInteger(3)));
			toGenerate.add(nf.add(nf.alpha(), nf.getInteger(4)));
			toGenerate.add(nf.add(nf.alpha(), nf.getInteger(5)));
			toGenerate.add(nf.add(nf.multiply(-3, nf.alpha()), nf.getInteger(3)));
			for (IntE prime : primes) {
				List<NumberFieldIdeal> ideals = order.idealsOver(z.getIdeal(Collections.singletonList(prime)));
				for (Ideal<NFE> ideal : ideals) {
					for (NFE g : toGenerate) {
						List<NFE> coefficients = ideal.generate(g);
						NFE reconstructed = ideal.residue(g);
						if (ideal.contains(g)) {
							assertEquals(order.zero(), reconstructed);
						}
						for (int i = 0; i < coefficients.size(); i++) {
							assertTrue(nf.isInteger(coefficients.get(i)));
							reconstructed = nf.add(reconstructed,
									nf.multiply(coefficients.get(i), ideal.generators().get(i)));
						}
						assertEquals(g, reconstructed);
					}
				}
			}
		}
	}

	@Test
	void testQuadraticFieldIdealFactorization() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		List<UnivariatePolynomial<Fraction>> minimalPolynomials = new ArrayList<>();
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(-6), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(-3), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(-2), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(1), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(2), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(3), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(5), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(6), q.zero(), q.one()));
		for (UnivariatePolynomial<Fraction> minimalPolynomial : minimalPolynomials) {
			NumberField nf = NumberField.getNumberField(minimalPolynomial);
			NumberFieldIntegers order = nf.maximalOrder();
			List<NFE> toFactorize = new ArrayList<>();
			toFactorize.add(nf.one());
			toFactorize.add(nf.getInteger(2));
			toFactorize.add(nf.getInteger(3));
			toFactorize.add(nf.getInteger(4));
			toFactorize.add(nf.getInteger(5));
			toFactorize.add(nf.getInteger(6));
			toFactorize.add(nf.alpha());
			toFactorize.add(nf.add(nf.alpha(), nf.one()));
			toFactorize.add(nf.add(nf.alpha(), nf.getInteger(2)));
			toFactorize.add(nf.add(nf.alpha(), nf.getInteger(3)));
			toFactorize.add(nf.add(nf.alpha(), nf.getInteger(4)));
			toFactorize.add(nf.add(nf.alpha(), nf.getInteger(5)));
			toFactorize.add(nf.add(nf.multiply(-3, nf.alpha()), nf.getInteger(3)));
			for (NFE factorize : toFactorize) {
				Ideal<NFE> ideal = order.getIdeal(factorize);
				for (NFE g : toFactorize) {
					List<NFE> coefficients = ideal.generate(g);
					NFE reconstructed = ideal.residue(g);
					if (ideal.contains(g)) {
						assertEquals(order.zero(), reconstructed);
					}
					for (int i = 0; i < coefficients.size(); i++) {
						assertTrue(nf.isInteger(coefficients.get(i)));
						reconstructed = nf.add(reconstructed,
								nf.multiply(coefficients.get(i), ideal.generators().get(i)));
					}
					assertEquals(g, reconstructed);
				}
			}
		}
	}

	@Test
	void testNumberFieldIndexCalculus() {
		PrimeField base = PrimeField.getPrimeField(11);
		FiniteField field = FiniteField.getFiniteField(
				base.getUnivariatePolynomialRing().getPolynomial(base.getInteger(5), base.zero(), base.one()), base);
		FFE primitive = field.add(field.one(), field.alpha()); // field.primitiveRoot();
		FFE power;
//		do {
//			power = field.getRandomElement();
//		} while (power.equals(field.zero()));
		power = field.add(field.getInteger(4), field.multiply(-3, field.alpha()));
		System.out.println(primitive);
		System.out.println(power);
		IntE discreteLog = field.discreteLogarithm(primitive, power);
		assertEquals(power, field.power(primitive, discreteLog));
//		Rationals q = Rationals.q();
//		Integers z = Integers.z();
//		FiniteField field = FiniteField.getFiniteField(125);
//		DiscreteValuationRing<Fraction, PFE> z5 = z.localize(field.characteristic());
//		NumberField nf = NumberField.getNumberField(z5.liftUnivariatePolynomial(field.minimalPolynomial()));
//		NumberFieldIntegers order = nf.maximalOrder();
//		DiscreteValuationRing<NFE, FFE> localized = order
//				.localize(order.idealsOver(z.getIdeal(new IntE(field.characteristic()))).get(0));
//		Set<IntE> primes = new TreeSet<>();
//		Iterator<IntE> primeIt = z.primes();
//		while (primes.size() < Math.ceil(Math.log(field.getNumberOfElements().doubleValue()))) {
//			IntE prime = primeIt.next();
//			if (prime.getValue().equals(field.characteristic())) {
//				continue;
//			}
//			primes.add(prime);
//		}
//		List<Ideal<NFE>> ideals = new ArrayList<>();
//		for (IntE prime : primes) {
//			ideals.addAll(order.idealsOver(z.getIdeal(prime)));
//		}
//		FFE base = field.primitiveRoot();
//		base = field.power(base, 4);
//		for (int i = 0; i < 64; i++) {
//			NFE lifted = localized.lift(field.power(base, i));
//			System.out.println("i: " + i + " result: " + lifted + " norm: " + nf.norm(lifted));
//		}
//		List<List<Fraction>> sieved = new ArrayList<>();
//		List<Fraction> rhs = new ArrayList<>();
//		Vector<Fraction> exponents;
//		while (true) {
//			BigInteger rng = MiscAlgorithms.randomBigInteger(new Random(), field.getNumberOfUnits().shiftRight(2));
//			FFE power = field.power(base, rng);
//			NFE lift = localized.lift(power);
//			Optional<NumberFieldIdeal> ideal = order.getIdealIfSmoothOver(Collections.singletonList(lift), primes);
//			if (ideal.isEmpty()) {
//				continue;
//			}
//			FactorizationResult<Ideal<NFE>, Ideal<NFE>> factorization = order.idealFactorization(ideal.get());
//			List<Fraction> powers = new ArrayList<>();
//			for (Ideal<NFE> primeIdeal : ideals) {
//				powers.add(q.getInteger(factorization.multiplicity(primeIdeal)));
//			}
//			sieved.add(powers);
//			rhs.add(q.getInteger(rng));
//			Matrix<Fraction> m = new Matrix<>(sieved);
//			Vector<Fraction> b = new Vector<>(rhs);
//			MatrixModule<Fraction> module = m.getModule(q);
//			if (module.kernelBasis(m).size() != 0) {
//				continue;
//			}
//			exponents = module.solve(m, b);
//			break;
//		}
//		for (FFE power : field.getMultiplicativeGroup()) {
//			Fraction exponent;
//			while (true) {
//				BigInteger rng = MiscAlgorithms.randomBigInteger(new Random(), field.getNumberOfUnits());
//				FFE result = field.divide(power, field.power(base, rng));
//				NFE lift = localized.lift(result);
//				Optional<NumberFieldIdeal> ideal = order.getIdealIfSmoothOver(Collections.singletonList(lift), primes);
//				if (ideal.isEmpty()) {
//					continue;
//				}
//				FactorizationResult<Ideal<NFE>, Ideal<NFE>> factorization = order.idealFactorization(ideal.get());
//				exponent = q.getInteger(rng);
//				for (int i = 0; i < ideals.size(); i++) {
//					Ideal<NFE> primeIdeal = ideals.get(i);
//					exponent = q.add(exponent,
//							q.multiply(factorization.multiplicity(primeIdeal), exponents.get(i + 1)));
//				}
//				break;
//			}
//			assertEquals(power, field.power(base, exponent.asInteger().getValue()));
//		}
	}

	@Test
	void testIntegralBasisQuadraticField() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		List<UnivariatePolynomial<Fraction>> minimalPolynomials = new ArrayList<>();
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(-6), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(-3), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(-2), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(1), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(2), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(3), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(5), q.zero(), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(6), q.zero(), q.one()));
		for (UnivariatePolynomial<Fraction> minimalPolynomial : minimalPolynomials) {
			NumberField nf = NumberField.getNumberField(minimalPolynomial);
			NumberFieldIntegers order = nf.maximalOrder();
			for (NFE b : order.getModuleGenerators()) {
				UnivariatePolynomial<Fraction> mipo = nf.minimalPolynomial(b);
				mipo = polynomials.normalize(mipo);
				for (int i = 0; i < mipo.degree(); i++) {
					assertTrue(q.isInteger(mipo.univariateCoefficient(i)));
				}
			}
		}
	}

	@Test
	void testHiddenIntegralBasis() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		List<UnivariatePolynomial<Fraction>> minimalPolynomials = new ArrayList<>();
		minimalPolynomials.add(
				polynomials.getPolynomial(q.getInteger(7), q.getInteger(0), q.getInteger(5), q.getInteger(0), q.one()));
		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(93), q.getInteger(63), q.getInteger(58),
				q.getInteger(39), q.getInteger(14), q.getInteger(3), q.one()));
		for (UnivariatePolynomial<Fraction> minimalPolynomial : minimalPolynomials) {
			NumberField nf = NumberField.getNumberField(minimalPolynomial);
			NumberFieldIntegers order = nf.maximalOrder();
			for (NFE b : order.getModuleGenerators()) {
				UnivariatePolynomial<Fraction> mipo = nf.minimalPolynomial(b);
				mipo = polynomials.normalize(mipo);
				for (int i = 0; i < mipo.degree(); i++) {
					assertTrue(q.isInteger(mipo.univariateCoefficient(i)));
				}
			}
		}
	}
}
