package varieties;

import java.io.IOException;
import java.util.Collections;
import java.util.List;

import org.junit.jupiter.api.Test;

import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.Monomial;
import varieties.affine.AffinePoint;
import varieties.affine.AffineScheme;
import varieties.affine.AffineScheme.BlowUpResult;

class BlowupTest {

	@Test
	void test() throws IOException {
		Rationals q = Rationals.q();
		PolynomialRing<Fraction> polynomials = AbstractPolynomialRing.getPolynomialRing(q, Monomial.GREVLEX,
				new String[] { "X", "Y" });
		Polynomial<Fraction> nodalPolynomial = polynomials.parse("-1*Y^2 + X^3 + -1*X^2");
		CoordinateRing<Fraction> cr = polynomials.getIdeal(Collections.singletonList(nodalPolynomial)).divideOut();
		AffineScheme<Fraction> scheme = new AffineScheme<>(q, cr);
		List<AffinePoint<Fraction>> singularities = scheme.singularPoints();
		System.out.println(singularities);
		for (AffinePoint<Fraction> point : singularities) {
			BlowUpResult<Fraction> blowup = scheme.blowupPoint(point);
			System.out.println(blowup.getBlowUp().getDomain().getAffineCover().get(0));
			for (int i = 0; i < blowup.getBlowUp().getDomain().getAffineCover().size(); i++) {
				System.out.println(
						blowup.getBlowUp().getDomain().getAffineCover().get(i).singularLocus().get().simplifyDomain());
			}
			System.out.println(blowup.getInverse().getMorphism().closure().simplifyDomain());
		}

	}

}
