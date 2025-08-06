package fields.integers;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.IOException;
import java.util.Collections;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.interfaces.PolynomialRing;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.interfaces.UnivariatePolynomialRing.ExtendedResultantResult;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;

class ModGcdTest {

	@Test
	void test() throws IOException {
		PrimeField field = PrimeField.getPrimeField(7);
		PolynomialRing<PFE> ring = AbstractPolynomialRing.getPolynomialRing(field, Monomial.REVLEX,
				new String[] { "z" });
		PolynomialIdeal<PFE> ideal = ring.getIdeal(Collections.singletonList(ring.parse("z^2 + -2")));
		CoordinateRing<PFE> cr = ideal.divideOut();
		UnivariatePolynomialRing<CoordinateRingElement<PFE>> poly = cr.getUnivariatePolynomialRing();
		UnivariatePolynomial<CoordinateRingElement<PFE>> t1 = poly.getPolynomial(cr.parse("z"), cr.parse("2"),
				cr.parse("z + -3"));
		UnivariatePolynomial<CoordinateRingElement<PFE>> t2 = poly.getPolynomial(cr.parse("2"), cr.parse("z + -2"),
				cr.parse("z"));
		System.out.println(t1);
		System.out.println(t2);
		ExtendedResultantResult<CoordinateRingElement<PFE>> er = poly.extendedResultant(t1, t2);
		assertEquals(er.getGcd(), poly.add(poly.multiply(er.getCoeff1(), t1), poly.multiply(er.getCoeff2(), t2)));
	}

}
