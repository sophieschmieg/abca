package fields.local;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.math.BigInteger;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField.PFE;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.PAdicField.PAdicNumber;

class CompleteLocalFieldExtensionTest {

	@Test
	void testZ2Extensions() {
		PAdicField z2 = new PAdicField(BigInteger.TWO, 40);
		UnivariatePolynomialRing<PAdicNumber> ring = z2.getUnivariatePolynomialRing();
		CompleteLocalFieldExtension<PAdicNumber, PFE, FFE, FiniteField> z2sqrt2 = z2
				.getExtension(ring.getPolynomial(z2.uniformizer(), z2.zero(), z2.one())).extension();
		assertEquals(2, z2sqrt2.ramificationIndex());
		assertEquals(1, z2sqrt2.residueDegree());
		assertEquals(1, z2sqrt2.valuation(z2sqrt2.uniformizer()).value());
		CompleteLocalFieldExtension<PAdicNumber, PFE, FFE, FiniteField> z2zeta3 = z2
				.getExtension(ring.getPolynomial(z2.one(), z2.one(), z2.one())).extension();
		assertEquals(z2zeta3.ramificationIndex(), 1);
		assertEquals(z2zeta3.residueDegree(), 2);
		/*Ext<PAdicNumber, PFE, FFE, FiniteField> dividend = z2zeta3.uniformizer();
		Ext<PAdicNumber, PFE, FFE, FiniteField> divisor = z2zeta3.multiply(3, z2zeta3.power(z2zeta3.uniformizer(), -2));
		assertEquals(dividend, z2zeta3.multiply(z2zeta3.divide(dividend, divisor), divisor));
		dividend = z2zeta3.uniformizer();
		divisor = z2zeta3.divide(z2zeta3.negative(z2zeta3.one()), z2zeta3.power(z2zeta3.uniformizer(), 2));
		assertEquals(dividend, z2zeta3.multiply(z2zeta3.divide(dividend, divisor), divisor));*/
		CompleteLocalFieldExtension<PAdicNumber, PFE, FFE, FiniteField> z2i = z2
				.getExtension(ring.getPolynomial(z2.one(), z2.zero(), z2.one())).extension();
		assertEquals(z2i.ramificationIndex(), 2);
		assertEquals(z2i.residueDegree(), 1);
		CompleteLocalFieldExtension<PAdicNumber, PFE, FFE, FiniteField> z2zeta3i = z2.getExtension(
				ring.getPolynomial(z2.getInteger(5), z2.getInteger(-2), z2.getInteger(1), z2.getInteger(-2), z2.one()))
				.extension();
		assertEquals(2, z2zeta3i.ramificationIndex());
		assertEquals(2, z2zeta3i.residueDegree());
		assertEquals(1, z2zeta3i.valuation(z2zeta3i.uniformizer()).value());
		assertEquals(0, z2zeta3i.valuation(z2zeta3i.alpha()).value());
		CompleteLocalFieldExtension<PAdicNumber, PFE, FFE, FiniteField> z2sqrt2var = z2
				.getExtension(ring.getPolynomial(z2.one(), z2.getInteger(6), z2.one())).extension();
		assertEquals(z2sqrt2var.ramificationIndex(), 2);
		assertEquals(z2sqrt2var.residueDegree(), 1);
		CompleteLocalFieldExtension<PAdicNumber, PFE, FFE, FiniteField> z2unramified = z2.getExtension(ring.getPolynomial(z2.getInteger(93),
				z2.getInteger(63), z2.getInteger(58), z2.getInteger(39), z2.getInteger(14), z2.getInteger(3), z2.one()))
				.extension();
		assertEquals(z2unramified.ramificationIndex(), 1);
		assertEquals(z2unramified.residueDegree(), 6);
	}

	@Test
	void testZ3Extensions() {
		PAdicField z3 = new PAdicField(BigInteger.valueOf(3), 10);
		UnivariatePolynomialRing<PAdicNumber> ring = z3.getUnivariatePolynomialRing();
		CompleteLocalFieldExtension<PAdicNumber, PFE, FFE, FiniteField> z3sqrt2 = z3
				.getExtension(ring.getPolynomial(z3.getInteger(-2), z3.zero(), z3.one())).extension();
		assertEquals(1, z3sqrt2.ramificationIndex());
		assertEquals(2, z3sqrt2.residueDegree());
		assertEquals(1, z3sqrt2.valuation(z3sqrt2.uniformizer()).value());
		CompleteLocalFieldExtension<PAdicNumber, PFE, FFE, FiniteField> z3zeta3 = z3
				.getExtension(ring.getPolynomial(z3.one(), z3.one(), z3.one())).extension();
		assertEquals(2, z3zeta3.ramificationIndex());
		assertEquals(1, z3zeta3.residueDegree());
		CompleteLocalFieldExtension<PAdicNumber, PFE, FFE, FiniteField> z3i = z3
				.getExtension(ring.getPolynomial(z3.one(), z3.zero(), z3.one())).extension();
		assertEquals(1, z3i.ramificationIndex());
		assertEquals(2, z3i.residueDegree());
		CompleteLocalFieldExtension<PAdicNumber, PFE, FFE, FiniteField> z3zeta3i = z3.getExtension(
				ring.getPolynomial(z3.getInteger(1), z3.getInteger(4), z3.getInteger(5), z3.getInteger(2), z3.one()))
				.extension();
		assertEquals(2, z3zeta3i.ramificationIndex());
		assertEquals(2, z3zeta3i.residueDegree());
		assertEquals(1, z3zeta3i.valuation(z3zeta3i.uniformizer()).value());
		assertEquals(0, z3zeta3i.valuation(z3zeta3i.alpha()).value());
		CompleteLocalFieldExtension<PAdicNumber, PFE, FFE, FiniteField> z3sqrt3 = z3
				.getExtension(ring.getPolynomial(z3.getInteger(-3), z3.zero(), z3.one())).extension();
		assertEquals(2, z3sqrt3.ramificationIndex());
		assertEquals(1, z3sqrt3.residueDegree());
		// CompleteLocalFieldExtension<PAdicNumber, PFE, FFE, FiniteField> z3zeta3unramified3 = z3
		// .getExtension(ring.getPolynomial(z3.getInteger(93), z3.getInteger(63),
		// z3.getInteger(58), z3.getInteger(39),z3.getInteger(14),z3.getInteger(3),
		// z3.one()));
		// assertEquals(2, z3zeta3unramified3.ramificationIndex());
		// assertEquals(3, z3zeta3unramified3.residueDegree());
	}

}
