package fields.tests;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;

import java.io.IOException;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.integers.Integers;
import fields.integers.Rationals;
import fields.interfaces.PolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.Monomial;

class ParserTest {

	@Test
	void testInteger() throws IOException {
		Integers z = Integers.z();
		assertEquals(z.getInteger(15), z.parse("15"));
		assertEquals(z.getInteger(-7), z.parse("-7"));
		assertEquals(z.zero(), z.parse("0"));
		assertThrows(IOException.class, () -> {
			z.parse("");
		});
		assertThrows(IOException.class, () -> {
			z.parse("+1");
		});
		assertThrows(IOException.class, () -> {
			z.parse("/1");
		});
		assertThrows(IOException.class, () -> {
			z.parse("1/");
		});
	}

	@Test
	void testFraction() throws IOException {
		Rationals q = Rationals.q();
		assertEquals(q.getInteger(15), q.parse("15"));
		assertEquals(q.getInteger(-7), q.parse("-7"));
		assertEquals(q.zero(), q.parse("0"));
		assertEquals(q.zero(), q.parse("0/4"));
		assertEquals(q.getFraction(1, 5), q.parse("1/5"));
		assertEquals(q.getFraction(1, 5), q.parse("3/15"));
		assertEquals(q.getFraction(-1, 5), q.parse("-3/15"));
		assertEquals(q.getFraction(1, 5), q.parse("-3/-15"));
		assertEquals(q.getFraction(-1, 5), q.parse("3/-15"));
		assertThrows(IOException.class, () -> {
			q.parse("");
		});
		assertThrows(IOException.class, () -> {
			q.parse("+1");
		});
		assertThrows(IOException.class, () -> {
			q.parse("/1");
		});
		assertThrows(IOException.class, () -> {
			q.parse("1/");
		});
	}

	@Test
	void testPrimeField() throws IOException {
		PrimeField fp = PrimeField.getPrimeField(31);
		assertEquals(fp.getInteger(15), fp.parse("15"));
		assertEquals(fp.getInteger(-7), fp.parse("-7"));
		assertEquals(fp.zero(), fp.parse("0"));
		assertThrows(IOException.class, () -> {
			fp.parse("");
		});
		assertThrows(IOException.class, () -> {
			fp.parse("+1");
		});
		assertThrows(IOException.class, () -> {
			fp.parse("/1");
		});
		assertThrows(IOException.class, () -> {
			fp.parse("1/");
		});
	}

	@Test
	void testPolynomialPrimeField() throws IOException {
		PrimeField fp = PrimeField.getPrimeField(31);
		PolynomialRing<PFE> polynomialRing = AbstractPolynomialRing.getPolynomialRing(fp, Monomial.GREVLEX,
				new String[] { "X", "Y", "a", "b" });
		assertEquals(polynomialRing.add(polynomialRing.negative(polynomialRing.getVarPower(2, 2)),
				polynomialRing.getVarPower(1, 3),
				polynomialRing.multiply(polynomialRing.getVar(3), polynomialRing.getVar(1)), polynomialRing.getVar(4)),
				polynomialRing.parse("-1*Y^2 + X^3 + a*X + b"));
		assertEquals(
				polynomialRing.multiply(polynomialRing.subtract(polynomialRing.getVar(1), polynomialRing.getVar(3)),
						polynomialRing.subtract(polynomialRing.getVar(2), polynomialRing.getVar(4))),
				polynomialRing.parse("X*Y + -1*b*X + -1*a*Y + a*b"));
	}
}
