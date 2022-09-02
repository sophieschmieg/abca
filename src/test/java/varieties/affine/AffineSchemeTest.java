package varieties.affine;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import varieties.curves.elliptic.EllipticCurve;
import varieties.projective.ProjectivePoint;

class AffineSchemeTest {

	@Test
	void testPoints() {
		PrimeField f5 = PrimeField.getPrimeField(5);
		EllipticCurve<PFE> j1728 = EllipticCurve.fromJInvariant(f5, f5.getElement(1728));
		Set<AffinePoint<PFE>> points = new TreeSet<>();
		for (ProjectivePoint<PFE> point : j1728) {
			if (!point.equals(j1728.neutral())) {
				points.add(point.getDehomogenous(3));
			}
		}
		assertEquals(j1728.getNumberOfElements().intValueExact() - 1, points.size());
		AffineScheme<PFE> affine = j1728.getAffineCover().getCover().get(2);
		int counter = 0;
		for (AffinePoint<PFE> point : affine) {
			assertTrue(points.contains(point));
			counter++;
		}
		assertEquals(j1728.getNumberOfElements().intValueExact() - 1, counter);
	}

	@Test
	void testSimplify() throws IOException {
		PrimeField f5 = PrimeField.getPrimeField(5);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(f5, Monomial.GREVLEX,
				new String[] { "X_1", "X_2", "X_3", "X_4" });
		List<Polynomial<PFE>> generators = new ArrayList<>();
		String[] polynomials = "X_1^3*X_4 + -1*X_4^2 + X_1*X_4 + X_3^2 + -2*X_1^2 + X_3, X_3^3 + -1*X_1^2*X_3 + -1*X_1^2, X_2*X_3^2 + -1*X_1*X_3 + -1*X_1, X_1*X_3^2 + -1*X_1^3 + -1*X_2*X_3 + -1*X_1*X_3 + X_4 + 2*X_1, X_2^2*X_3 + -1*X_3 + -1, X_3*X_4 + -1*X_2*X_3 + X_4 + X_1, X_2*X_4 + -1*X_3^2 + X_1^2 + X_3 + -1, X_1*X_2 + -1*X_3"
				.split(",");
		for (String polynomial : polynomials) {
			generators.add(polynomialRing.parse(polynomial));
		}
		AffineScheme<PFE> scheme = new AffineScheme<>(f5, polynomialRing.getIdeal(generators).divideOut());
		scheme.simplify();
	}
}
