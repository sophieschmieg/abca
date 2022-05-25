package varieties.curves;

import static org.junit.Assert.assertTrue;

import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import varieties.curves.elliptic.EllipticCurve;
import varieties.curves.elliptic.Isogeny;
import varieties.projective.ProjectivePoint;

class EllipticCurveChar23Test {

	@Test
	void testChar2() {
		FiniteField f16 = FiniteField.getFiniteField(16);
		FFE one = f16.one();
		EllipticCurve<FFE> curve = new EllipticCurve<>(f16, f16.alpha(), one, one, one, one);
		List<Isogeny<FFE>> automorphism = curve.getAutomorphisms();
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
	}

	@Test
	void testChar2Supersingular() {
		FiniteField f16 = FiniteField.getFiniteField(16);
		FFE one = f16.one();
		EllipticCurve<FFE> curve = new EllipticCurve<>(f16, f16.zero(), one, one, one, one);
		List<Isogeny<FFE>> automorphism = curve.getAutomorphisms();
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
	}


	@Test
	void testChar3() {
		FiniteField f27 = FiniteField.getFiniteField(27);
		FFE one = f27.one();
		EllipticCurve<FFE> curve = new EllipticCurve<>(f27, f27.alpha(), one, one, one, one);
		List<Isogeny<FFE>> automorphism = curve.getAutomorphisms();
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
	}

	@Test
	void testChar3Supersingular() {
		FiniteField f27 = FiniteField.getFiniteField(27);
		FFE one = f27.one();
		EllipticCurve<FFE> curve = new EllipticCurve<>(f27, f27.zero(), one, one, one, one);
		List<Isogeny<FFE>> automorphism = curve.getAutomorphisms();
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
	}


	@Test
	void testChar5() {
		FiniteField f125 = FiniteField.getFiniteField(125);
		FFE one = f125.one();
		EllipticCurve<FFE> curve = new EllipticCurve<>(f125, f125.alpha(), one, one, one, one);
		List<Isogeny<FFE>> automorphism = curve.getAutomorphisms();
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
	}

	@Test
	void testChar5Supersingular() {
		FiniteField f125 = FiniteField.getFiniteField(125);
		FFE one = f125.one();
		EllipticCurve<FFE> curve = new EllipticCurve<>(f125, f125.zero(), one, one, one, one);
		List<Isogeny<FFE>> automorphism = curve.getAutomorphisms();
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
	}

}
