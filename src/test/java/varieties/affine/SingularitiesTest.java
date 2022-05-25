package varieties.affine;

import java.io.IOException;
import java.util.Collections;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;

class SingularitiesTest {

	@Test
	void test() throws IOException {
		PrimeField f257 = PrimeField.getPrimeField(257);
		PolynomialRing<PFE> polynomials2 = AbstractPolynomialRing.getPolynomialRing(f257, 2, Monomial.REVLEX);
		Polynomial<PFE> elliptic = polynomials2.parse("X^3 + -1*X + -1*Y^2");
		AffineScheme<PFE> ellCurve = new AffineScheme<>(f257,
				polynomials2.getIdeal(Collections.singletonList(elliptic)).divideOut());
		ellCurve.singularPoints();
		Polynomial<PFE> nodal = polynomials2.parse("X^3 + X^2 + -1*Y^2");
		AffineScheme<PFE> nodalCurve = new AffineScheme<>(f257,
				polynomials2.getIdeal(Collections.singletonList(nodal)).divideOut());
		nodalCurve.singularPoints();
		}

}
