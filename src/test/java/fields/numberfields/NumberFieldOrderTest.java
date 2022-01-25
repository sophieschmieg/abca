package fields.numberfields;

import org.junit.jupiter.api.Test;

import fields.integers.Rationals;

class NumberFieldOrderTest {

	@Test
	void test() {
		Rationals q = Rationals.q();
		NumberField nf = new NumberField(
				q.getUnivariatePolynomialRing().getPolynomial(q.getInteger(3), q.zero(), q.one()));
		NumberFieldOrder order = nf.getOrder(nf.alpha());
		System.out.println(order.conductor());
		System.out.println(order.getModuleGenerators());
		System.out.println(order.asVector(nf.alpha()));
		System.out.println(nf.maximalOrder().asOrder().conductor());
	}

}
