package fields.numberfields;

import java.io.IOException;

import org.junit.jupiter.api.Test;

import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.UnivariatePolynomialRing;
import fields.numberfields.NumberField.NFE;

class NumberFieldIntegerTest {

	@Test
	void test() throws IOException {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomials = q.getUnivariatePolynomialRing();
		NumberField nf = NumberField.getNumberField(polynomials.parse("X^4 + 1"));
		NFE onetwothreefour = nf.fromPolynomial(polynomials.parse("1 + 2*X + 3*X^2 + 4*X^3"));
		System.out.println(onetwothreefour);
		nf.maximalOrder().idealsOver(2);
		nf.maximalOrder().idealsOver(3);
			System.out.println(nf.maximalOrder().projectToUnit(onetwothreefour));
		System.out.println(nf.maximalOrder().upToUnit(onetwothreefour));
	}

}
