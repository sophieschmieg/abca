package varieties.curves;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;
import java.util.Collections;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.interfaces.DiscreteValuationRing.TheMontesResult;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.UnivariatePolynomial;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.LocalizedCoordinateRing;
import fields.polynomials.LocalizedCoordinateRing.LocalizedElement;
import fields.polynomials.Monomial;
import varieties.RationalFunction;
import varieties.curves.elliptic.EllipticCurve;
import varieties.projective.ProjectivePoint;

class RationalFunctionTest {

	// @Test
	void test() {
		PrimeField f29 = PrimeField.getPrimeField(29);
		EllipticCurve<PFE> curve = EllipticCurve.fromJInvariant(f29, f29.getInteger(1728));
		ProjectivePoint<PFE> zero = curve.getRandomElement();
		ProjectivePoint<PFE> pole1;
		ProjectivePoint<PFE> pole2;
		do {
			pole1 = curve.getRandomElement();
			pole2 = curve.getRandomElement();
		} while (pole1.equals(zero) || pole2.equals(zero));
		RationalFunction<PFE> function = curve.getRationalFunction(zero, pole1, pole2);
		System.out.println(zero);
		System.out.println(pole1);
		System.out.println(pole2);
		System.out.println(function);
		System.out.println(function.getZeroes());
		System.out.println(function.getPoles());
		assertTrue(function.getZeroes().contains(zero));
		assertTrue(function.getPoles().contains(pole1));
		assertTrue(function.getPoles().contains(pole2));
	}

	@Test
	void normalizationTest() throws IOException {
		PrimeField f5 = PrimeField.getPrimeField(5);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(f5, Monomial.GREVLEX,
				new String[] { "X", "Z" });
		Polynomial<PFE> def = polynomialRing.parse("Z^4 + X^4 + -1*Z^2");
		CoordinateRing<PFE> cr = f5.getUnivariatePolynomialRing()
				.getZeroIdeal().divideOut();
		LocalizedCoordinateRing<PFE> line = new LocalizedCoordinateRing<>(f5, cr,
				cr.getIdeal(Collections.singletonList(cr.getVar(1))));
	UnivariatePolynomial<LocalizedElement<PFE>> mipo = //line.getUnivariatePolynomialRing().getPolynomial(line.getEmbedding(cr.add(cr.getEmbedding(f5.getUnivariatePolynomialRing().getVarPower(3)), cr.getEmbedding(f5.getUnivariatePolynomialRing().getVarPower(2)))), line.zero(), line.one());
 line.getUnivariatePolynomialRing().getPolynomial(line.getEmbedding(f5.getUnivariatePolynomialRing().getVarPower(4)), line.zero(), line.negative(line.one()), line.zero(), line.one());
TheMontesResult<LocalizedElement<PFE>, PFE, PFE, FFE, FiniteField> theMontes=	line.ringOfIntegers().theMontesAlgorithm(mipo, f5.getExtension(f5.getUnivariatePolynomialRing().getVar()));
	System.out.println(line.ringOfIntegers().integralBasis(mipo, theMontes, false));
	}
}
