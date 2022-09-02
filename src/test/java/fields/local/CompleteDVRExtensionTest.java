//package fields.local;
//
//import static org.junit.jupiter.api.Assertions.assertEquals;
//
//import java.math.BigInteger;
//
//import org.junit.jupiter.api.Test;
//
//import fields.finitefields.FiniteField;
//import fields.finitefields.FiniteField.FFE;
//import fields.finitefields.PrimeField.PFE;
//import fields.interfaces.UnivariatePolynomialRing;
//import fields.local.CompleteDVRExtension.Ext;
//import fields.local.PAdicField.PAdicNumber;
//
//class CompleteDVRExtensionTest {
//
//	@Test
//	void testZ2Extensions() {
//		PAdicField z2 = new PAdicField(BigInteger.TWO, 40);
//		UnivariatePolynomialRing<PAdicNumber> ring = z2.getUnivariatePolynomialRing();
//		CompleteDVRExtension<PAdicNumber, PFE, PFE, FFE, FiniteField> z2sqrt2 = z2
//				.getExtension(ring.getPolynomial(z2.uniformizer(), z2.zero(), z2.one())).extension();
//		assertEquals(2, z2sqrt2.ramificationIndex());
//		assertEquals(1, z2sqrt2.residueDegree());
//		assertEquals(1, z2sqrt2.valuation(z2sqrt2.uniformizer()).value());
//		@SuppressWarnings("unchecked")
//		CompleteDVRExtension<PAdicNumber, PFE, PFE, FFE, FiniteField> z2sqrt2i = z2sqrt2
//				.getExtension(z2sqrt2.getUnivariatePolynomialRing().getPolynomial(z2sqrt2.one(), z2sqrt2.zero(), z2sqrt2.one())).extension();
//		assertEquals(z2sqrt2i.ramificationIndex(), 4);
//		assertEquals(z2sqrt2i.residueDegree(), 1);
//		System.out.println(z2sqrt2i.uniformizer());
//		CompleteDVRExtension<PAdicNumber, PFE, PFE, FFE, FiniteField> z2zeta3 = z2
//				.getExtension(ring.getPolynomial(z2.one(), z2.one(), z2.one())).extension();
//		assertEquals(z2zeta3.ramificationIndex(), 1);
//		assertEquals(z2zeta3.residueDegree(), 2);
//		Ext<PAdicNumber> dividend = z2zeta3.uniformizer();
//		Ext<PAdicNumber> divisor = z2zeta3.multiply(3, z2zeta3.power(z2zeta3.uniformizer(), -2));
//		System.out.println(dividend);
//		System.out.println(divisor);
//		System.out.println(z2zeta3.inverse(z2zeta3.alpha()));
//		assertEquals(dividend, z2zeta3.multiply(z2zeta3.divide(dividend, divisor), divisor));
//		dividend = z2zeta3.uniformizer();
//		divisor = z2zeta3.divide(z2zeta3.negative(z2zeta3.one()), z2zeta3.power(z2zeta3.uniformizer(), 2));
//		System.out.println(dividend);
//		System.out.println(divisor);
//		assertEquals(dividend, z2zeta3.multiply(z2zeta3.divide(dividend, divisor), divisor));
//		CompleteDVRExtension<PAdicNumber, PFE, PFE, FFE, FiniteField> z2i = z2
//				.getExtension(ring.getPolynomial(z2.one(), z2.zero(), z2.one())).extension();
//		assertEquals(z2i.ramificationIndex(), 2);
//		assertEquals(z2i.residueDegree(), 1);
//		CompleteDVRExtension<PAdicNumber, PFE, PFE, FFE, FiniteField> z2zeta3i = z2.getExtension(
//				ring.getPolynomial(z2.getInteger(5), z2.getInteger(-2), z2.getInteger(1), z2.getInteger(-2), z2.one()))
//				.extension();
//		assertEquals(2, z2zeta3i.ramificationIndex());
//		assertEquals(2, z2zeta3i.residueDegree());
//		assertEquals(1, z2zeta3i.valuation(z2zeta3i.uniformizer()).value());
//		assertEquals(0, z2zeta3i.valuation(z2zeta3i.alpha()).value());
//		CompleteDVRExtension<PAdicNumber, PFE, PFE, FFE, FiniteField> z2sqrt2var = z2
//				.getExtension(ring.getPolynomial(z2.one(), z2.getInteger(6), z2.one())).extension();
//		assertEquals(z2sqrt2var.ramificationIndex(), 2);
//		assertEquals(z2sqrt2var.residueDegree(), 1);
//		CompleteDVRExtension<PAdicNumber, PFE, PFE, FFE, FiniteField> z2unramified = z2
//				.getExtension(ring.getPolynomial(z2.getInteger(93), z2.getInteger(63), z2.getInteger(58),
//						z2.getInteger(39), z2.getInteger(14), z2.getInteger(3), z2.one()))
//				.extension();
//		assertEquals(z2unramified.ramificationIndex(), 1);
//		assertEquals(z2unramified.residueDegree(), 6);
//		CompleteDVRExtension<PAdicNumber, PFE, PFE, FFE, FiniteField> z2zeta3b = z2
//				.getExtension(ring.getPolynomial(z2.getInteger(3), z2.zero(), z2.one())).extension();
//		assertEquals(z2zeta3b.ramificationIndex(), 1);
//		assertEquals(z2zeta3b.residueDegree(), 2);
//		System.out.println(z2zeta3b.uniformizer());
//		System.out.println(z2zeta3b.alpha());
//		System.out.println(z2zeta3b.inverse(z2zeta3b.alpha()));
//			}
//
//	@Test
//	void testZ3Extensions() {
//		PAdicField z3 = new PAdicField(BigInteger.valueOf(3), 10);
//		UnivariatePolynomialRing<PAdicNumber> ring = z3.getUnivariatePolynomialRing();
//		CompleteDVRExtension<PAdicNumber, PFE, PFE, FFE, FiniteField> z3sqrt2 = z3
//				.getExtension(ring.getPolynomial(z3.getInteger(-2), z3.zero(), z3.one())).extension();
//		assertEquals(1, z3sqrt2.ramificationIndex());
//		assertEquals(2, z3sqrt2.residueDegree());
//		assertEquals(1, z3sqrt2.valuation(z3sqrt2.uniformizer()).value());
//		CompleteDVRExtension<PAdicNumber, PFE, PFE, FFE, FiniteField> z3zeta3 = z3
//				.getExtension(ring.getPolynomial(z3.one(), z3.one(), z3.one())).extension();
//		assertEquals(2, z3zeta3.ramificationIndex());
//		assertEquals(1, z3zeta3.residueDegree());
//		CompleteDVRExtension<PAdicNumber, PFE, PFE, FFE, FiniteField> z3i = z3
//				.getExtension(ring.getPolynomial(z3.one(), z3.zero(), z3.one())).extension();
//		assertEquals(1, z3i.ramificationIndex());
//		assertEquals(2, z3i.residueDegree());
//		CompleteDVRExtension<PAdicNumber, PFE, PFE, FFE, FiniteField> z3zeta3i = z3.getExtension(
//				ring.getPolynomial(z3.getInteger(1), z3.getInteger(4), z3.getInteger(5), z3.getInteger(2), z3.one()))
//				.extension();
//		assertEquals(2, z3zeta3i.ramificationIndex());
//		assertEquals(2, z3zeta3i.residueDegree());
//		assertEquals(1, z3zeta3i.valuation(z3zeta3i.uniformizer()).value());
//		assertEquals(0, z3zeta3i.valuation(z3zeta3i.alpha()).value());
//		CompleteDVRExtension<PAdicNumber, PFE, PFE, FFE, FiniteField> z3sqrt3 = z3
//				.getExtension(ring.getPolynomial(z3.getInteger(-3), z3.zero(), z3.one())).extension();
//		assertEquals(2, z3sqrt3.ramificationIndex());
//		assertEquals(1, z3sqrt3.residueDegree());
//	}
//
//}
