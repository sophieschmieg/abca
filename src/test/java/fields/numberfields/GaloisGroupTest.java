package fields.numberfields;

import java.io.IOException;

import org.junit.jupiter.api.Test;

import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.UnivariatePolynomialRing;

class GaloisGroupTest {

	@Test
	void testNormal() throws IOException {
		Integers z = Integers.z();
		UnivariatePolynomialRing<IntE> polynomials = z.getUnivariatePolynomialRing();
		NumberField nf1 = NumberField.getNumberFieldFromIntegerPolynomial(polynomials.parse("X^3 + -3*X + -1"));
		System.out.println(nf1.galoisGroup());
		NumberField nf2 = NumberField.getNumberFieldFromIntegerPolynomial(polynomials.parse("X^6 + 108"));
		System.out.println(nf2.galoisGroup());
	}

}
