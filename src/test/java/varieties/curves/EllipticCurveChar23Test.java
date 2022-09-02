package varieties.curves;

import static org.junit.Assert.assertTrue;
import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.IOException;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.vectors.FiniteVectorSpace;
import fields.vectors.Vector;
import varieties.curves.elliptic.EllipticCurve;
import varieties.curves.elliptic.EllipticCurve.Isomorphism;
import varieties.curves.elliptic.Isogeny;
import varieties.projective.ProjectivePoint;

class EllipticCurveChar23Test {

	@Test
	void testChar2() throws IOException {
		FiniteField f16 = FiniteField.getFiniteField(16);
		FFE one = f16.one();
		EllipticCurve<FFE> curve = new EllipticCurve<>(f16, f16.alpha(), one, one, one, one);
		List<Isomorphism<FFE>> automorphism = curve.getAutomorphisms();
		System.out.println(automorphism);
		for (Isogeny<FFE> morphism : automorphism) {
			for (ProjectivePoint<FFE> point : curve) {
				assertTrue(curve.hasRationalPoint(morphism.evaluate(point)));
			}
		}
		Isogeny<FFE> weierstrass = curve.getWeierstrassForm();
		EllipticCurve<FFE> range = weierstrass.getRange();
		for (ProjectivePoint<FFE> point : curve) {
			assertTrue(range.hasRationalPoint(weierstrass.evaluate(point)));
		}
		int counter = 0;
		for (ProjectivePoint<FFE> point : curve) {
			counter++;
			System.out.println(point);
		}
		System.out.println(counter);
		System.out.println(curve.getNumberOfElements());
		assertEquals(counter, curve.getNumberOfElements().intValueExact());
	}

	@Test
	void testChar2Supersingular() {
		FiniteField f16 = FiniteField.getFiniteField(16);
		FFE one = f16.one();
		FFE zero = f16.zero();
		EllipticCurve<FFE> curve = new EllipticCurve<>(f16, zero, zero, one, zero, zero/*f16.alpha(), one, one, one*/);
		assertTrue(curve.isSupersingular());
		List<Isomorphism<FFE>> automorphism = curve.getAutomorphisms();
		assertEquals(24, automorphism.size());
		System.out.println(automorphism);
		for (Isogeny<FFE> morphism : automorphism) {
			for (ProjectivePoint<FFE> point : curve) {
				assertTrue(curve.hasRationalPoint(morphism.evaluate(point)));
			}
		}
		Isogeny<FFE> weierstrass = curve.getWeierstrassForm();
		EllipticCurve<FFE> range = weierstrass.getRange();
		for (ProjectivePoint<FFE> point : curve) {
			assertTrue(range.hasRationalPoint(weierstrass.evaluate(point)));
		}
		int counter = 0;
		for (ProjectivePoint<FFE> point : curve) {
			counter++;
			System.out.println(point);
		}
		System.out.println(counter);
		System.out.println(curve.getNumberOfElements());
		assertEquals(counter, curve.getNumberOfElements().intValueExact());
	}

	@Test
	void testChar3() {
		FiniteField f27 = FiniteField.getFiniteField(27);
		EllipticCurve<FFE> fromJInvariant = EllipticCurve.fromJInvariant(f27, f27.alpha());
		assertEquals(f27.alpha(), fromJInvariant.jInvariant());
		FFE one = f27.one();
		EllipticCurve<FFE> curve = new EllipticCurve<>(f27, f27.alpha(), one, one, one, one);
		List<Isomorphism<FFE>> automorphism = curve.getAutomorphisms();
		assertEquals(12, automorphism.size());
		System.out.println(automorphism);
		for (Isogeny<FFE> morphism : automorphism) {
			for (ProjectivePoint<FFE> point : curve) {
				assertTrue(curve.hasRationalPoint(morphism.evaluate(point)));
			}
		}
		Isogeny<FFE> weierstrass = curve.getWeierstrassForm();
		EllipticCurve<FFE> range = weierstrass.getRange();
		for (ProjectivePoint<FFE> point : curve) {
			assertTrue(range.hasRationalPoint(weierstrass.evaluate(point)));
		}
		int counter = 0;
		for (ProjectivePoint<FFE> point : curve) {
			counter++;
			System.out.println(point);
		}
		System.out.println(counter);
		System.out.println(curve.getNumberOfElements());
		assertEquals(counter, curve.getNumberOfElements().intValueExact());
	}

	@Test
	void testChar3Base() {
		PrimeField f3 = PrimeField.getPrimeField(3);
		FiniteVectorSpace<PFE> f3Sqr = new FiniteVectorSpace<>(f3, 2);
		UnivariatePolynomialRing<PFE> r = f3.getUnivariatePolynomialRing();
		for (Vector<PFE> vector : f3Sqr) {
			PFE a = vector.get(1);
			PFE b = vector.get(2);
			UnivariatePolynomial<PFE> rhs = r.getPolynomial(b, a, f3.one(), f3.one());
			if (!r.squareFreeFactorization(rhs).squareFree()) {
				continue;
			}
			EllipticCurve<PFE> curve = new EllipticCurve<>(f3, f3.zero(), f3.one(), f3.zero(), a, b);
			int counter = 0;
			for (ProjectivePoint<PFE> point : curve) {
				counter++;
				assertTrue(curve.hasRationalPoint(point));
			}
			assertEquals(counter, curve.getNumberOfElements().intValueExact());
		}
	}

	@Test
	void testChar3Base2() {
		FiniteField f3 = FiniteField.getFiniteField(3);
		FiniteVectorSpace<FFE> f3Sqr = new FiniteVectorSpace<>(f3, 2);
		UnivariatePolynomialRing<FFE> r = f3.getUnivariatePolynomialRing();
		for (Vector<FFE> vector : f3Sqr) {
			FFE a = vector.get(1);
			FFE b = vector.get(2);
			UnivariatePolynomial<FFE> rhs = r.getPolynomial(b, a, f3.one(), f3.one());
			if (!r.squareFreeFactorization(rhs).squareFree()) {
				continue;
			}
			EllipticCurve<FFE> curve = new EllipticCurve<>(f3, f3.zero(), f3.one(), f3.zero(), a, b);
			int counter = 0;
			for (ProjectivePoint<FFE> point : curve) {
				counter++;
				assertTrue(curve.hasRationalPoint(point));
			}
			assertEquals(counter, curve.getNumberOfElements().intValueExact());
		}
	}

	@Test
	void testF81() {
		FiniteField f81 = FiniteField.getFiniteField(81);
		FiniteVectorSpace<FFE> f81Fifth = new FiniteVectorSpace<>(f81, 5);
		for (int tc = 0; tc < 100; tc++) {
			Vector<FFE> vector = f81Fifth.getRandomElement();
			FFE a1 = vector.get(1);
			FFE a2 = vector.get(2);
			FFE a3 = vector.get(3);
			FFE a4 = vector.get(4);
			FFE a6 = vector.get(5);
			try {
				EllipticCurve<FFE> curve = new EllipticCurve<>(f81, a1, a2, a3, a4, a6);
				System.out.println(curve);
				int counter = 0;
				for (ProjectivePoint<FFE> point : curve) {
					counter++;
					assertTrue(curve.hasRationalPoint(point));
				}
				System.out.println(counter);
				assertEquals(counter, curve.getNumberOfElements().intValueExact());
			} catch (Exception e) {
				continue;
			}
		}
	}

	@Test
	void testF64() {
		FiniteField f64 = FiniteField.getFiniteField(64);
		FiniteVectorSpace<FFE> f64Fifth = new FiniteVectorSpace<>(f64, 5);
		for (int tc = 0; tc < 100; tc++) {
			Vector<FFE> vector = f64Fifth.getRandomElement();
			FFE a1 = vector.get(1);
			FFE a2 = vector.get(2);
			FFE a3 = vector.get(3);
			FFE a4 = vector.get(4);
			FFE a6 = vector.get(5);
			try {
				EllipticCurve<FFE> curve = new EllipticCurve<>(f64, a1, a2, a3, a4, a6);
				System.out.println(curve);
				int counter = 0;
				for (ProjectivePoint<FFE> point : curve) {
					counter++;
					assertTrue(curve.hasRationalPoint(point));
				}
				System.out.println(counter);
				assertEquals(counter, curve.getNumberOfElements().intValueExact());
			} catch (Exception e) {
				continue;
			}
		}
	}

	@Test
	void testF65536() {
		FiniteField f65536 = FiniteField.getFiniteField(65536);
		FiniteVectorSpace<FFE> f64Fifth = new FiniteVectorSpace<>(f65536, 5);
		for (int tc = 0; tc < 10; tc++) {
			Vector<FFE> vector = f64Fifth.getRandomElement();
			FFE a1 = vector.get(1);
			FFE a2 = vector.get(2);
			FFE a3 = vector.get(3);
			FFE a4 = vector.get(4);
			FFE a6 = vector.get(5);
			try {
				EllipticCurve<FFE> curve = new EllipticCurve<>(f65536, a1, a2, a3, a4, a6);
				System.out.println(curve);
				int count = curve.getNumberOfElements().intValueExact();
				System.out.println(count);
				for (int i = 0; i < 20; i++) {
					assertEquals(curve.neutral(), curve.multiply(count, curve.getRandomElement()));
				}
			} catch (Exception e) {
				tc--;
				System.err.println(e.getMessage());
				e.printStackTrace(System.err);
				continue;
			}
		}
	}

	@Test
	void testF2to128() {
		PrimeField f2 = PrimeField.getPrimeField(2);
		UnivariatePolynomialRing<PFE> polynomials = f2.getUnivariatePolynomialRing();
		UnivariatePolynomial<PFE> minimalPolynomial = polynomials.add(polynomials.getVarPower(128), polynomials
				.add(polynomials.getVarPower(7), polynomials.getVarPower(2), polynomials.getVar(), polynomials.one()));
		FiniteField f2to128 = FiniteField.getFiniteField(minimalPolynomial, f2);
		FiniteVectorSpace<FFE> f64Fifth = new FiniteVectorSpace<>(f2to128, 5);
		for (int tc = 0; tc < 10; tc++) {
			Vector<FFE> vector = f64Fifth.getRandomElement();
			FFE a1 = vector.get(1);
			FFE a2 = vector.get(2);
			FFE a3 = vector.get(3);
			FFE a4 = vector.get(4);
			FFE a6 = vector.get(5);
			try {
				EllipticCurve<FFE> curve = new EllipticCurve<>(f2to128, a1, a2, a3, a4, a6);
				System.out.println(curve);
				int count = curve.getNumberOfElements().intValueExact();
				System.out.println(count);
				for (int i = 0; i < 20; i++) {
					assertEquals(curve.neutral(), curve.multiply(count, curve.getRandomElement()));
				}
			} catch (Exception e) {
				tc--;
				System.err.println(e.getMessage());
				e.printStackTrace(System.err);
				continue;
			}
		}
	}

	@Test
	void testChar3Supersingular() {
		FiniteField f27 = FiniteField.getFiniteField(27);
		FFE one = f27.one();
		FFE zero = f27.zero();
		EllipticCurve<FFE> curve = new EllipticCurve<>(f27, zero, zero, zero, one, one);
		assertTrue(curve.isSupersingular());
		List<Isomorphism<FFE>> automorphism = curve.getAutomorphisms();
		System.out.println(automorphism);
		for (Isogeny<FFE> morphism : automorphism) {
			for (ProjectivePoint<FFE> point : curve) {
				assertTrue(curve.hasRationalPoint(morphism.evaluate(point)));
			}
		}
		Isogeny<FFE> weierstrass = curve.getWeierstrassForm();
		EllipticCurve<FFE> range = weierstrass.getRange();
		for (ProjectivePoint<FFE> point : curve) {
			assertTrue(range.hasRationalPoint(weierstrass.evaluate(point)));
		}
		int counter = 0;
		for (ProjectivePoint<FFE> point : curve) {
			counter++;
			System.out.println(point);
		}
		System.out.println(counter);
		System.out.println(curve.getNumberOfElements());
		assertEquals(counter, curve.getNumberOfElements().intValueExact());
	}

	@Test
	void testChar5() {
		FiniteField f125 = FiniteField.getFiniteField(125);
		FFE one = f125.one();
		EllipticCurve<FFE> curve = new EllipticCurve<>(f125, f125.one(), f125.zero(), f125.zero(), one, one);
		List<Isomorphism<FFE>> automorphism = curve.getAutomorphisms();
		System.out.println(automorphism);
		for (Isogeny<FFE> morphism : automorphism) {
			for (ProjectivePoint<FFE> point : curve) {
				assertTrue(curve.hasRationalPoint(morphism.evaluate(point)));
			}
		}
		int counter = 0;
		for (ProjectivePoint<FFE> point : curve) {
			counter++;
			System.out.println(point);
		}
		System.out.println(counter);
		System.out.println(curve.getNumberOfElements());
		assertEquals(counter, curve.getNumberOfElements().intValueExact());
	}

	@Test
	void testChar5Supersingular() {
		FiniteField f125 = FiniteField.getFiniteField(125);
		FFE one = f125.one();
		FFE zero = f125.zero();
		EllipticCurve<FFE> curve = new EllipticCurve<>(f125, zero, zero, zero, zero, one);
		assertTrue(curve.isSupersingular());
		List<Isomorphism<FFE>> automorphism = curve.getAutomorphisms();
		assertEquals(6, automorphism.size());
		System.out.println(automorphism);
		for (Isogeny<FFE> morphism : automorphism) {
			for (ProjectivePoint<FFE> point : curve) {
				assertTrue(curve.hasRationalPoint(morphism.evaluate(point)));
			}
		}
		Isogeny<FFE> weierstrass = curve.getWeierstrassForm();
		EllipticCurve<FFE> range = weierstrass.getRange();
		for (ProjectivePoint<FFE> point : curve) {
			assertTrue(range.hasRationalPoint(weierstrass.evaluate(point)));
		}
		int counter = 0;
		for (ProjectivePoint<FFE> point : curve) {
			counter++;
			System.out.println(point);
		}
		System.out.println(counter);
		System.out.println(curve.getNumberOfElements());
		assertEquals(counter, curve.getNumberOfElements().intValueExact());
	}

	@Test
	void countPointsChar2() throws IOException {
		PrimeField f2 = PrimeField.getPrimeField(2);
		FiniteField f128 = FiniteField.getFiniteField(f2.getUnivariatePolynomialRing().parse("X^7 + X + 1"), f2);
		EllipticCurve<FFE> curve = new EllipticCurve<FiniteField.FFE>(f128, f128.one(), f128.zero(), f128.zero(),
				f128.zero(), f128.inverse(f128.parse("x^5 + x + 1")));
		int counter = 0;
		for (ProjectivePoint<FFE> point : curve) {
			counter++;
			System.out.println(point);
		}
		System.out.println(counter);
		System.out.println(curve.getNumberOfElements());
		assertEquals(counter, curve.getNumberOfElements().intValueExact());
		curve = curve.getQuadraticTwist();
		counter = 0;
		for (ProjectivePoint<FFE> point : curve) {
			counter++;
			System.out.println(point);
		}
		System.out.println(counter);
		System.out.println(curve.getNumberOfElements());
		assertEquals(counter, curve.getNumberOfElements().intValueExact());
	}

}
