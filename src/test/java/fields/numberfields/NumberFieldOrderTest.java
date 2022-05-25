package fields.numberfields;

import java.util.Collections;

import org.junit.jupiter.api.Test;

import fields.integers.Rationals;
import fields.numberfields.NumberField.NFE;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;
import fields.numberfields.NumberFieldOrder.NumberFieldOrderIdeal;
import fields.numberfields.PicardGroup.OrderIdealClass;

class NumberFieldOrderTest {

	@Test
	void testEisenstein() {
		Rationals q = Rationals.q();
		NumberField nf = NumberField
				.getNumberField(q.getUnivariatePolynomialRing().getPolynomial(q.getInteger(3), q.zero(), q.one()));
		NumberFieldOrder order = nf.getOrder(nf.alpha());
		System.out.println(order.conductor());
		System.out.println(order.conductorInMaximalOrder());
		System.out.println(order.getModuleGenerators());
		System.out.println(order.asVector(nf.alpha()));
		System.out.println(nf.maximalOrder().asOrder().conductor());
	}

	@Test
	void testGauss() {
		Rationals q = Rationals.q();
		NumberField nf = NumberField
				.getNumberField(q.getUnivariatePolynomialRing().getPolynomial(q.one(), q.zero(), q.one()));
		NumberFieldOrder order = nf.getOrder(nf.multiply(nf.getInteger(6), nf.alpha()));
		System.out.println(order.conductor());
		System.out.println(order.conductorInMaximalOrder());
		System.out.println(order.getModuleGenerators());
		for (OrderIdealClass ic : order.picardGroup()) {
			System.out.println(ic);
		}

		NumberFieldOrderIdeal ideal = order.getIdeal(nf.getInteger(4), nf.multiply(nf.getInteger(6), nf.alpha()));
		System.out.println(ideal);
		System.out.println(order.extensibleToMaximalOrder(ideal));
		System.out.println(order.extendToMaximalOrder(ideal));
		System.out.println(order.restrictFromMaximalOrder(order.extendToMaximalOrder(ideal)));
		OrderIdealGroup idealGroup = order.idealGroup();
		System.out.println(ideal.isPrincipal());
		FractionalOrderIdeal fractional = new FractionalOrderIdeal(order, ideal);
		System.out.println(idealGroup.operate(fractional, fractional).getNumerator().isPrincipal());
		FractionalOrderIdeal inverse = new FractionalOrderIdeal(order, order.getUnitIdeal(), ideal);
		System.out.println(idealGroup.operate(fractional, inverse));
		NFE t1 = nf.add(nf.getInteger(6), nf.multiply(nf.getInteger(6), nf.alpha()));
		NFE t2 = nf.add(nf.getInteger(18), nf.multiply(nf.getInteger(-18), nf.alpha()));
		NumberFieldOrderIdeal ideal2 = order.getIdeal(t1, t2);
		FractionalOrderIdeal fractional2 = new FractionalOrderIdeal(order, ideal2);
		System.out.println(fractional2);
		System.out.println(fractional2.getNumerator().isPrincipal());
		System.out.println(idealGroup.power(2, fractional2));
		System.out.println(idealGroup.power(2, fractional2).getNumerator().isPrincipal());
		System.out.println(idealGroup.power(4, fractional2));
		System.out.println(idealGroup.power(4, fractional2).getNumerator().isPrincipal());

		order = nf.getOrder(nf.multiply(nf.getInteger(2), nf.alpha()));
		System.out.println(order.conductor());
		System.out.println(order.conductorInMaximalOrder());
		System.out.println(order.getModuleGenerators());

		NumberFieldIdeal ideal1 = nf.maximalOrder()
				.getIdeal(Collections.singletonList(nf.add(nf.getInteger(4), nf.alpha())));
		System.out.println(ideal1);
		System.out.println(order.restrictableFromMaximalOrder(ideal1));
		System.out.println(order.restrictFromMaximalOrder(ideal1));
		System.out.println(order.extendToMaximalOrder(order.restrictFromMaximalOrder(ideal1)));
	}

}
