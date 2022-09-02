package varieties.curves;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import varieties.projective.GenericProjectiveScheme;

class HyperEllipticCurvesTest {

	@Test
	void test() {
		PrimeField fp = PrimeField.getPrimeField(7);
		UnivariatePolynomialRing<PFE> r = fp.getUnivariatePolynomialRing();
		UnivariatePolynomial<PFE> rhs = r.toUnivariate(r.subtract(r.getVarPower(6), r.one()));
		HyperEllipticCurve<PFE> curve = new HyperEllipticCurve<>(fp, rhs, r.zero());
		GenericProjectiveScheme<PFE> generic = curve.asGenericProjectiveScheme();
		assertEquals(1, generic.dimension());
		assertTrue(generic.singularPoints().isEmpty());
	}

}
