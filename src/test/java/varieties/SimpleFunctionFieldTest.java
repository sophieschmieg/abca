package varieties;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.UnivariatePolynomial;
import fields.numberfields.NumberField;
import fields.numberfields.NumberField.NFE;
import varieties.curves.elliptic.EllipticCurve;

class SimpleFunctionFieldTest {

	@Test
	void testFunctionField() {
		Rationals q = Rationals.q();
		UnivariatePolynomial<Fraction> eisenstein = q.getUnivariatePolynomialRing().getPolynomial(q.one(), q.one(),
				q.one());
		NumberField nf = NumberField.getNumberField(eisenstein);
		EllipticCurve<NFE> curve = new EllipticCurve<>(nf, nf.zero(), nf.one());
		assertEquals(3, curve.degree());
	//	SimpleFunctionFieldFromCoordinateRingOverExtension<Fraction, NFE> sff = SimpleFunctionField.forProjectiveVarietyOverExtensionField(nf, curve);
	}

}
