package fields.numberfields;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.math.BigInteger;

import org.junit.jupiter.api.Test;

import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Ring.FactorizationResult;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;

class ClassNumberTest {

//	@Test
//	void testGeneratorExpressionQuadraticFields() {
//		Integers z = Integers.z();
//		Rationals q = Rationals.q();
//		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
//		List<UnivariatePolynomial<Fraction>> minimalPolynomials = new ArrayList<>();
//		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(-6), q.zero(), q.one()));
//		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(-3), q.zero(), q.one()));
//		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(-2), q.zero(), q.one()));
//		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(1), q.zero(), q.one()));
//		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(2), q.zero(), q.one()));
//		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(3), q.zero(), q.one()));
//		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(5), q.zero(), q.one()));
//		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(6), q.zero(), q.one()));
//		minimalPolynomials.add(
//				polynomials.getPolynomial(q.getInteger(7), q.getInteger(0), q.getInteger(5), q.getInteger(0), q.one()));
//		for (UnivariatePolynomial<Fraction> minimalPolynomial : minimalPolynomials) {
//			NumberField nf = new NumberField(minimalPolynomial);
//			System.out.println(nf);
//			NumberFieldIntegers order = nf.maximalOrder();
//			for (int tc = 0; tc < 10; tc++) {
//				NFE rng;
//				do {
//					rng = order.getRandomElement();
//				} while (nf.norm(rng).equals(q.zero()) || z.isUnit(nf.norm(rng).asInteger()));
//				IntE norm = nf.norm(rng).asInteger();
//				List<IntE> factors = z.factors(norm);
//				IntE factor;
//				do {
//					factor = factors.get(new Random().nextInt(factors.size()));
//				} while (z.isUnit(factor));
//				Ideal<NFE> id = order.getIdeal(rng, order.getInteger(factor));
//				List<NFE> generators = new ArrayList<>();
//				for (int i = 0; i < 10; i++) {
//					generators.add(id.getRandomElement());
//				}
//				IdealResult<NFE> ideal = order.getIdealWithTransforms(generators);
//				for (int i = 0; i < ideal.getIdeal().generators().size(); i++) {
//					NFE generator = ideal.getIdeal().generators().get(i);
//					NFE test = order.zero();
//					for (int j = 0; j < generators.size(); j++) {
//						NFE coefficient = ideal.getGeneratorExpressions().get(i).get(j);
//						test = order.add(test, order.multiply(coefficient, generators.get(j)));
//					}
//					assertEquals(generator, test);
//				}
//			}
//		}
//	}
//
//	@Test
//	void testGeneratorExpressionHighDegreeFields() {
//		Integers z = Integers.z();
//		Rationals q = Rationals.q();
//		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
//		List<UnivariatePolynomial<Fraction>> minimalPolynomials = new ArrayList<>();
//		minimalPolynomials.add(polynomials.getPolynomial(q.getInteger(93), q.getInteger(63), q.getInteger(58),
//				q.getInteger(39), q.getInteger(14), q.getInteger(3), q.one()));
//		for (UnivariatePolynomial<Fraction> minimalPolynomial : minimalPolynomials) {
//			NumberField nf = new NumberField(minimalPolynomial);
//			System.out.println(nf);
//			NumberFieldIntegers order = nf.maximalOrder();
//			for (int tc = 0; tc < 3; tc++) {
//				NFE rng;
//				do {
//					rng = order.getRandomElement();
//				} while (nf.norm(rng).equals(q.zero()) || z.isUnit(nf.norm(rng).asInteger()));
//				IntE norm = nf.norm(rng).asInteger();
//				List<IntE> factors = z.factors(norm);
//				IntE factor;
//				do {
//					factor = factors.get(new Random().nextInt(factors.size()));
//				} while (z.isUnit(factor));
//				Ideal<NFE> id = order.getIdeal(rng, order.getInteger(factor));
//				List<NFE> generators = new ArrayList<>();
//				for (int i = 0; i < 10; i++) {
//					generators.add(id.getRandomElement());
//				}
//				IdealResult<NFE> ideal = order.getIdealWithTransforms(generators);
//				for (int i = 0; i < ideal.getIdeal().generators().size(); i++) {
//					NFE generator = ideal.getIdeal().generators().get(i);
//					NFE test = order.zero();
//					for (int j = 0; j < generators.size(); j++) {
//						NFE coefficient = ideal.getGeneratorExpressions().get(i).get(j);
//						test = order.add(test, order.multiply(coefficient, generators.get(j)));
//					}
//					assertEquals(generator, test);
//				}
//			}
//		}
//	}

	@Test
	void testClassImaginaryQuadraticNumberOne() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		int[] d = new int[] { 1, 2, 3, 7, 11, 19, 43, 67, 163 };
		for (int i : d) {
			NumberField nf = new NumberField(polynomials.getPolynomial(q.getInteger(i), q.zero(), q.one()));
			assertEquals(BigInteger.ONE, nf.classNumber());
		}
	}

	@Test
	void testClassImaginaryQuadraticNumberTwo() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		int[] d = new int[] { 15, 5, 6, 35, 10, 51, 13, 22, 91, 115, 123, 37, 187, 58, 235, 267, 403, 427 };
		for (int i : d) {
			NumberField nf = new NumberField(polynomials.getPolynomial(q.getInteger(i), q.zero(), q.one()));
			assertEquals(BigInteger.TWO, nf.classNumber());
		}
	}

	@Test
	void testClassImaginaryQuadraticNumberThree() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		int[] d = new int[] { 23, 31, 59, 83, 107, 139, 211, 283, 307, 331, 379, 499, 547, 643, 883, 907 };
		System.out.println("Checking class number 3");
		for (int i : d) {
			NumberField nf = new NumberField(polynomials.getPolynomial(q.getInteger(i), q.zero(), q.one()));
			System.out.println(nf);
			assertEquals(BigInteger.valueOf(3), nf.classNumber());
		}
	}

	@Test
	void testClassImaginaryQuadraticNumberFour() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		int[] d = new int[] { 39, 55, 14, 17, 21, 30, 33, 34, 155, 42, 46, 195, 203, 219, 57, 259, 70, 291, 73, 78, 323,
				82, 85, 355, 93, 97, 102, 435, 483, 130, 133, 555, 142, 595, 627, 667, 177, 715, 723, 190, 763, 193,
				795, 955, 1003, 253, 1027, 1227, 1243, 1387, 1411, 1435, 1507, 1555 };
		System.out.println("Checking class number 4");
		for (int i : d) {
			NumberField nf = new NumberField(polynomials.getPolynomial(q.getInteger(i), q.zero(), q.one()));
			System.out.println(nf);
			assertEquals(BigInteger.valueOf(4), nf.classNumber());
		}
	}

	@Test
	void testClassImaginaryQuadraticNumberFive() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		int[] d = new int[] { 47, 79, 103, 127, 131, 179, 227, 347, 443, 523, 571, 619, 683, 691, 739, 787, 947, 1051,
				1123, 1723, 1747, 1867, 2203, 2347, 2683 };
		System.out.println("Checking class number 5");
		for (int i : d) {
			NumberField nf = new NumberField(polynomials.getPolynomial(q.getInteger(i), q.zero(), q.one()));
			System.out.println(nf);
			assertEquals(BigInteger.valueOf(5), nf.classNumber());
		}
	}

	@Test
	void testClassRealQuadraticNumberOne() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		int[] d = new int[] { 2, 5, 13, 29, 53, 173, 293, 151, 199, 367, 622, 863, 1151 };
		System.out.println("Checking real quadratic class number 1");
		for (int i : d) {
			NumberField nf = new NumberField(polynomials.getPolynomial(q.getInteger(-i), q.zero(), q.one()));
			System.out.println(nf);
			assertEquals(BigInteger.ONE, nf.classNumber());
		}
	}

	@Test
	void testClassRealQuadraticNumberTwo() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		int[] d = new int[] { 10, 26, 85, 122, 362, 365, 533, 629, 965, 1685, 1853, 2813, 15 };
		System.out.println("Checking real quadratic class number 2");
		for (int i : d) {
			NumberField nf = new NumberField(polynomials.getPolynomial(q.getInteger(-i), q.zero(), q.one()));
			System.out.println(nf);
			assertEquals(BigInteger.TWO, nf.classNumber());
		}
	}

	@Test
	void testClassRealQuadraticNumberThree() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		int[] d = new int[] {229, 257, 79, 321, 469, 473, 142, 733, 761, 223, 993, 254 };
		System.out.println("Checking real quadratic class number 3");
		for (int i : d) {
			NumberField nf = new NumberField(polynomials.getPolynomial(q.getInteger(-i), q.zero(), q.one()));
			System.out.println(nf);
			assertEquals(BigInteger.valueOf(3), nf.classNumber());
		}
	}

	@Test
	void testClassNumberDegreeThree() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		UnivariatePolynomial<Fraction> f = polynomials.getPolynomial(q.getInteger(4), q.one(), q.zero(), q.one());
		NumberField field = new NumberField(f);
		assertEquals(BigInteger.valueOf(2), field.classNumber());
	}

	@Test
	void testDiscriminant() {
		Rationals q = Rationals.q();
		Integers z = Integers.z();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		testCaseLoop: for (int i = 1; i < 20; i++) {
			FactorizationResult<IntE> f = z.uniqueFactorization(z.getInteger(i));
			for (IntE prime : f.primeFactors()) {
				if (f.multiplicity(prime) > 1) {
					continue testCaseLoop;
				}
			}
			int d = -4 * i;
			if (i % 4 == 3) {
				d = -i;
			}
			assertEquals(d, new NumberField(polynomials.getPolynomial(q.getInteger(i), q.zero(), q.one()))
					.discriminant().intValueExact());
		}
	}

}
