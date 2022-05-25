package fields.tests;

import java.util.Collections;

import org.junit.jupiter.api.Test;

import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.UnivariatePolynomialRing;
import fields.numberfields.ModuloNumberFieldIdeal;
import fields.numberfields.ModuloNumberFieldIdeal.ModNFE;
import fields.numberfields.NumberField;
import fields.numberfields.NumberFieldIntegers;
import fields.numberfields.NumberFieldIntegers.NumberFieldIdeal;

class ModuloTest {

	@Test
	void numberFieldModuloTest() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> rationalPolynomialRing = q.getUnivariatePolynomialRing();
		NumberField nf = NumberField
				.getNumberField(rationalPolynomialRing.getPolynomial(q.getInteger(5), q.zero(), q.one()));
		NumberFieldIntegers order = nf.maximalOrder();
		NumberFieldIdeal ideal = order.getIdeal(Collections.singletonList(order.getInteger(18)));
		// NumberFieldIdeal idealPower = order.power(ideal, 4);
		ModuloNumberFieldIdeal mod = ideal.modOut();
		System.out.println(mod.getNumberOfElements());
		System.out.println(mod.getNumberOfUnits());
		System.out.println(mod.getUnitGenerators());
		int count = 0;
		for (ModNFE t : mod) {
			count++;
			System.out.println(count + ": " + t);
		}
		count = 0;
		System.out.println();
		for (ModNFE t : mod.getMultiplicativeGroup()) {
			count++;
			System.out.println(count + ": " + t);
		}
	}

	@Test
	void numberFieldModuloTest2() {
		NumberField nf = NumberField.getNumberField();
		NumberFieldIntegers order = nf.maximalOrder();
		NumberFieldIdeal ideal = order.getIdeal(Collections.singletonList(nf.getInteger(3500)));
		ModuloNumberFieldIdeal mod = ideal.modOut();
		System.out.println(mod.getUnitGenerators());
	}

	@Test
	void numberFieldModuloTest3() {
		Rationals q = Rationals.q();
		UnivariatePolynomialRing<Fraction> rationalPolynomialRing = q.getUnivariatePolynomialRing();
		NumberField nf = NumberField.getNumberField(rationalPolynomialRing.getPolynomial(q.getInteger(507),
				q.getInteger(48), q.getInteger(49), q.getInteger(2), q.one()));
		NumberFieldIntegers order = nf.maximalOrder();
		NumberFieldIdeal ideal = order.getIdeal(Collections.singletonList(order.getInteger(4)));//order.power(order.idealsOver(2).get(0), 2);
		ModuloNumberFieldIdeal mod = ideal.modOut();
		System.out.println(mod.getNumberOfElements());
		System.out.println(mod.getNumberOfUnits());
		int count = 0;
		int unitCount = 0;
		for (ModNFE t : mod) {
			count++;
			System.out.println(count + ": " + t);
			boolean unit = mod.isUnit(t);
			System.out.println(unit ? "Unit" : "Not a unit");
			unitCount += unit ? 1 : 0;
		}
		System.out.println("Unit count: " + unitCount);
		count = 0;
		System.out.println();
		System.out.println(mod.getUnitGenerators());
		System.out.println();
		for (ModNFE t : mod.getMultiplicativeGroup()) {
			count++;
			System.out.println(count + ": " + t);
		}
	}

}
