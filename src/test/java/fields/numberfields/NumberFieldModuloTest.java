package fields.numberfields;

import org.junit.jupiter.api.Test;

import fields.integers.Integers;
import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.UnivariatePolynomialRing;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;

class NumberFieldModuloTest {

	@Test
	void testSqrtMinus5() {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomialRing = q.getUnivariatePolynomialRing();
		NumberField nf = NumberField.getNumberField(polynomialRing.getPolynomial(q.getInteger(5), q.zero(), q.one()));
		System.out.println(nf);
		NumberFieldIntegers order = nf.maximalOrder();
		System.out.println(nf.classNumber());
		NumberFieldIdeal over2 = order.idealsOver(z.getIdeal(z.getInteger(2))).get(0);
		System.out.println(over2);
		System.out.println(order.multiply(over2, over2));
		System.out.println(order.multiply(over2, order.multiply(over2, over2)));
		NumberFieldIdeal over3 = order.idealsOver(z.getIdeal(z.getInteger(3))).get(0);
		System.out.println(over3);
		System.out.println(order.multiply(over3, over3));
		System.out.println(order.multiply(over3, order.multiply(over3, over3)));
	}

	@Test
	void testSqrtMinus3() {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomialRing = q.getUnivariatePolynomialRing();
		NumberField nf = NumberField.getNumberField(polynomialRing.getPolynomial(q.getInteger(3), q.zero(), q.one()));
		System.out.println(nf);
		NumberFieldIntegers order = nf.maximalOrder();
		System.out.println(nf.classNumber());
		NumberFieldIdeal over2 = order.idealsOver(z.getIdeal(z.getInteger(2))).get(0);
		System.out.println(over2);
		System.out.println(order.multiply(over2, over2));
		System.out.println(order.multiply(over2, order.multiply(over2, over2)));
		NumberFieldIdeal over3 = order.idealsOver(z.getIdeal(z.getInteger(3))).get(0);
		System.out.println(over3);
		System.out.println(order.multiply(over3, over3));
		System.out.println(order.multiply(over3, order.multiply(over3, over3)));
	}

	@Test
	void testNoPowerBasis() {
		Integers z = Integers.z();
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> polynomialRing = q.getUnivariatePolynomialRing();
		NumberField nf = NumberField.getNumberField(
				polynomialRing.getPolynomial(q.getInteger(-8), q.getInteger(-2), q.getInteger(-1), q.one()));
		System.out.println(nf);
		NumberFieldIntegers order = nf.maximalOrder();
		System.out.println(nf.classNumber());
		NumberFieldIdeal over2 = order.idealsOver(z.getIdeal(z.getInteger(2))).get(0);
		System.out.println(over2);
		System.out.println(order.multiply(over2, over2));
		System.out.println(order.multiply(over2, order.multiply(over2, over2)));
		NumberFieldIdeal over3 = order.idealsOver(z.getIdeal(z.getInteger(3))).get(0);
		System.out.println(over3);
		System.out.println(order.multiply(over3, over3));
		System.out.println(order.multiply(over3, order.multiply(over3, over3)));
	}

}
