package varieties;

import java.io.IOException;
import java.util.Collections;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;
import varieties.affine.AffinePoint;
import varieties.projective.GenericProjectiveScheme;

class GluedSchemeTest {

	//@Test
	void test1() throws IOException {
		PrimeField f5 = PrimeField.getPrimeField(5);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(f5, 2, Monomial.GREVLEX);
		Polynomial<PFE> nodalPolynomial = polynomialRing.parse("-1*Y^2 + X^3 + X^2");
		GenericProjectiveScheme<PFE> nodalCurve = GenericProjectiveScheme.fromAffineCoordinateRing(
				polynomialRing.getIdeal(Collections.singletonList(nodalPolynomial)).divideOut());
		GluedScheme<PFE> glued = new GluedScheme<>(f5, nodalCurve.getAffineCover());
		glued.fromAffinePoint(new AffinePoint<>(f5, f5.getInteger(2), f5.getInteger(2)), 0);
		glued.fromAffinePoint(new AffinePoint<>(f5, f5.getInteger(3), f5.getInteger(1)), 2);
	}

	@Test
	void test2() throws IOException {
		PrimeField f5 = PrimeField.getPrimeField(5);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(f5, 2, Monomial.GREVLEX);
		Polynomial<PFE> intersection = polynomialRing.parse("X*Y");
		GenericProjectiveScheme<PFE> nodalCurve = GenericProjectiveScheme.fromAffineCoordinateRing(
				polynomialRing.getIdeal(Collections.singletonList(intersection)).divideOut());
		GluedScheme<PFE> glued = new GluedScheme<>(f5, nodalCurve.getAffineCover());
		glued.fromAffinePoint(new AffinePoint<>(f5, f5.getInteger(0), f5.getInteger(2)), 0);
		glued.fromAffinePoint(new AffinePoint<>(f5, f5.getInteger(3), f5.getInteger(0)), 2);
	glued.irreducibleComponents();
	}

}
