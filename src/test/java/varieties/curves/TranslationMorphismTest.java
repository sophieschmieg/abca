package varieties.curves;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import varieties.curves.elliptic.EllipticCurve;
import varieties.curves.elliptic.Isogeny;
import varieties.curves.elliptic.KernelPointIsogeny;
import varieties.projective.ProjectiveMorphism;
import varieties.projective.ProjectivePoint;

class TranslationMorphismTest {

	// @Test
	void test() {
		FiniteField f125 = FiniteField.getFiniteField(125);
		EllipticCurve<FFE> curve = new EllipticCurve<>(f125, f125.one(), f125.alpha(), f125.one(), f125.one(),
				f125.one());
		for (ProjectivePoint<FFE> point1 : curve) {
			ProjectiveMorphism<FFE> translation = curve.translationMorphism(point1);
			assertEquals(curve.asGenericProjectiveScheme(), translation.getDomain());
			assertEquals(curve.asGenericProjectiveScheme(), translation.getRange());
			for (ProjectivePoint<FFE> point2 : curve) {
				ProjectivePoint<FFE> added = curve.add(point1, point2);
				System.out.println(point1 + " + " + point2 + " = " + added);
				assertTrue(curve.hasRationalPoint(added));
				assertEquals(added, translation.evaluate(point2));
			}
		}
	}

	@Test
	void testVelu() {
		FiniteField f125 = FiniteField.getFiniteField(125);
		EllipticCurve<FFE> curve = new EllipticCurve<>(f125, f125.one(), f125.alpha(), f125.one(), f125.one(),
				f125.one());
		List<ProjectivePoint<FFE>> torsionPointBasis = curve.getTorsionPointBasis(4);
		System.out.println(torsionPointBasis);
		Isogeny<FFE> isogeny = new KernelPointIsogeny<>(curve, torsionPointBasis.get(0), 4);
		for (ProjectivePoint<FFE> point : curve) {
			System.out.println("f(" + point + ") = " + isogeny.evaluate(point));
			assertTrue(isogeny.getRange().hasRationalPoint(isogeny.evaluate(point)));
		}
	}

}
