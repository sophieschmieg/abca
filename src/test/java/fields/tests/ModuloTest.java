package fields.tests;

import org.junit.jupiter.api.Test;

import fields.integers.Integers;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.UnivariatePolynomialRing;
import fields.numberfields.NumberField;
import fields.numberfields.NumberFieldIntegers;

class ModuloTest {

	@Test
	void numberFieldModuloTest() {
		Rationals q = Rationals.q();
		Integers z = Integers.z();
		UnivariatePolynomialRing<Fraction> rationalPolynomialRing = q.getUnivariatePolynomialRing();
	NumberField nf = new NumberField(rationalPolynomialRing.getPolynomial(q.getInteger(3), q.zero(), q.one()));
	NumberFieldIntegers order = nf.maximalOrder();
	System.out.println(order.idealsOver(z.getIdeal(z.getInteger(2))));
	}

	@Test
	void numberFieldModuloTest2() {
		Rationals q = Rationals.q();
	//	Integers z = Integers.z();
		UnivariatePolynomialRing<Fraction> rationalPolynomialRing = q.getUnivariatePolynomialRing();
	NumberField nf = new NumberField(rationalPolynomialRing.getPolynomial(q.getInteger(-8), q.getInteger(-2), q.getInteger(-1), q.one()));
	NumberFieldIntegers order = nf.maximalOrder();
	System.out.println(order.getModuleGenerators());
	}

}
