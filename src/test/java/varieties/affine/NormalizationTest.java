package varieties.affine;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Optional;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.integers.Rationals;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import varieties.GluedMorphism;
import varieties.GluedScheme;
import varieties.GluedScheme.GluedPoint;
import varieties.projective.GenericProjectiveScheme;

class NormalizationTest {

	private <T extends Element<T>> void doTest(Field<T> field, String equation) throws IOException {
		PolynomialRing<T> polynomials = AbstractPolynomialRing.getPolynomialRing(field, Monomial.GREVLEX,
				new String[] { "X", "Y" });
		Polynomial<T> elliptic = polynomials.parse(equation);
		System.out.println(elliptic);
		AffineScheme<T> scheme = new AffineScheme<>(field,
				polynomials.getIdeal(Collections.singletonList(elliptic)).divideOut());
		System.out.println(scheme.singularLocus());
		AffineMorphism<T> normalization = scheme.normalization().getNormalizationMorphism();
		System.out.println(normalization);
		Optional<AffineMorphism<T>> normalizationSingularLocus = normalization.getDomain().singularLocus();
		System.out.println(normalizationSingularLocus);
		assertTrue(scheme.dimension() != 1 || normalizationSingularLocus.isEmpty());
		System.out.println(normalization.getDomain().irreducibleComponents());
		assertTrue(!scheme.isIntegral() || normalization.getDomain().isIntegral());
		List<AffinePoint<T>> singularPoints = scheme.singularPoints();
		for (AffinePoint<T> singular : singularPoints) {
			System.out.println(normalization.preimageList(singular));
		}
		if (scheme.isFinite()) {
			for (AffinePoint<T> point : scheme) {
				List<AffinePoint<T>> preimage = normalization.preimageList(point);
				System.out.println("f^-1(" + point + ") = " + preimage);
				assertTrue(singularPoints.contains(point) || preimage.size() == 1);
			}
		}
	}

	// @Test
	void testAffine() throws IOException {
		doTest(PrimeField.getPrimeField(5), "X^2 + X^3 + -1*Y^2");
		doTest(PrimeField.getPrimeField(5), "X^3 + -1*Y^2");
		doTest(PrimeField.getPrimeField(5), "X^4 + Y^4 + -1*Y^2");
		doTest(PrimeField.getPrimeField(7), "X^2 + X^3 + -1*Y^2");
		doTest(PrimeField.getPrimeField(7), "X^3 + -1*Y^2");
		doTest(PrimeField.getPrimeField(7), "X^4 + Y^4 + -1*Y^2");
		doTest(Rationals.q(), "X^2 + X^3 + -1*Y^2");
		doTest(Rationals.q(), "X^3 + -1*Y^2");
		doTest(Rationals.q(), "X^4 + Y^4 + -1*Y^2");
	}

	@Test
	void testGlued() throws IOException {
		PrimeField f5 = PrimeField.getPrimeField(5);
		PolynomialRing<PFE> proj = AbstractPolynomialRing.getPolynomialRing(f5, Monomial.GREVLEX,
				new String[] { "X", "Y", "Z" });
		Polynomial<PFE> nodalPolynomial = proj.parse("-1*Z^4 + X^4 + -1*Y^2*Z^2");
		GluedScheme<PFE> nodal = new GluedScheme<>(f5,
				new GenericProjectiveScheme<>(f5, proj, Collections.singletonList(nodalPolynomial)).getAffineCover());
		nodal = nodal.simplify().getSimplification().getRange();
//		for (GluedPoint<PFE> point : nodal) {
//			System.out.println(point);
//			for (int index : nodal.affineCoverIndex(point)) {
//				System.out.println(nodal.asAffinePoint(point, index) + "@" + index);
//			}
//			System.out.println("----");
//		}
//		System.out.println();
		GluedMorphism<PFE> normalization = nodal.normalization();
//		for (GluedPoint<PFE> point : normalization.getDomain()) {
//			System.out.println(point);
//			for (int index : normalization.getDomain().affineCoverIndex(point)) {
//				System.out.println(normalization.getDomain().asAffinePoint(point, index) + "@" + index);
//			}
//			System.out.println("----");
//		}
		for (GluedPoint<PFE> point : nodal) {
			System.out.println(point + ": " + normalization.preimageList(point));
		}
	}
}
