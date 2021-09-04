package fields.tests;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.helper.FieldEmbedding;
import fields.integers.Rationals;
import fields.numberfields.NumberField;

class FieldExtensionTest {

	@Test
	void testF5() {
		PrimeField f5 = PrimeField.getPrimeField(5);
		FiniteField f25 = f5.getExtension(f5.getUnivariatePolynomialRing().getPolynomial(f5.getInteger(2), f5.zero(), f5.one())).extension();
		System.out.println(f25.genericNorm());
	}
	
	@Test
	void testF2() {
		PrimeField f2 = PrimeField.getPrimeField(2);
		FiniteField f4 = f2.getExtension(f2.getUnivariatePolynomialRing().getPolynomial(f2.one(), f2.one(), f2.one())).extension();
		FieldEmbedding<PFE, FFE, FiniteField> f16overf4 = f4.getEmbeddedExtension(f4.getUnivariatePolynomialRing().getPolynomial(f4.alpha(), f4.one(), f4.one()));
		System.out.println(f16overf4.minimalPolynomial());
		System.out.println(f16overf4.getField().minimalPolynomial());
	}

	@Test
	void testGauss() {
		Rationals q  = Rationals.q();
		NumberField gauss = new NumberField(q.getUnivariatePolynomialRing().getPolynomial(q.one(), q.zero(), q.one()));
		System.out.println(gauss.genericNorm());
	}

	@Test
	void testSqrt2() {
		Rationals q  = Rationals.q();
		NumberField sqrt2 = new NumberField(q.getUnivariatePolynomialRing().getPolynomial(q.getInteger(-2), q.zero(), q.one()));
		System.out.println(sqrt2.genericNorm());
	}

	@Test
	void testEisenstein() {
		Rationals q  = Rationals.q();
		NumberField eisenstein = new NumberField(q.getUnivariatePolynomialRing().getPolynomial(q.one(), q.one(), q.one()));
		NumberField eisensteini = eisenstein.getExtension(eisenstein.getUnivariatePolynomialRing().getPolynomial(eisenstein.one(), eisenstein.zero(), eisenstein.one())).extension();
		System.out.println(eisenstein.genericNorm());
		System.out.println(eisensteini.genericNorm());
		}

}
