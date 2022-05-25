package varieties.affine;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;

class SingularitiesTest {

	@Test
	void testAffine2() throws IOException {
		PrimeField f257 = PrimeField.getPrimeField(257);
		PolynomialRing<PFE> polynomials2 = AbstractPolynomialRing.getPolynomialRing(f257, 2, Monomial.GREVLEX);
		Polynomial<PFE> elliptic = polynomials2.parse("X^3 + -1*X + -1*Y^2");
		AffineScheme<PFE> ellCurve = new AffineScheme<>(f257,
				polynomials2.getIdeal(Collections.singletonList(elliptic)).divideOut());
		System.out.println(ellCurve.singularPoints());
		Polynomial<PFE> nodal = polynomials2.parse("X^3 + X^2 + -1*Y^2");
		AffineScheme<PFE> nodalCurve = new AffineScheme<>(f257,
				polynomials2.getIdeal(Collections.singletonList(nodal)).divideOut());
		System.out.println(nodalCurve.singularPoints());
	}

	@Test
	void testTwistedCubic() throws IOException {
		PrimeField f257 = PrimeField.getPrimeField(257);
		PolynomialRing<PFE> polynomials3 = AbstractPolynomialRing.getPolynomialRing(f257, 3, Monomial.GREVLEX);
		List<Polynomial<PFE>> generators = new ArrayList<>();
		generators.add(polynomials3.parse("X*Z + -1*Y^2"));
		generators.add(polynomials3.parse("Y + -1*Z^2"));
		generators.add(polynomials3.parse("X + -1*Y*Z"));
		AffineScheme<PFE> twistedCubic = new AffineScheme<>(f257, polynomials3.getIdeal(generators).divideOut());
		System.out.println(twistedCubic.singularPoints());
	}

}
