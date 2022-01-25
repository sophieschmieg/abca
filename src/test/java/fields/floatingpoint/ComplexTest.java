package fields.floatingpoint;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Random;

import org.junit.jupiter.api.Test;

import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.UnivariatePolynomial;

class ComplexTest {

	@Test
	void testAddition() {
		Complex c = Complex.c(1024);
		assertEquals(c.one(), c.add(c.one(), c.zero()));
		assertEquals(c.getInteger(2), c.add(c.one(), c.one()));
		assertEquals(c.getDouble(1.0, 1.0), c.add(c.one(), c.i()));
		assertEquals(c.getDouble(2.0, 2.0), c.add(c.getDouble(1.0, 1.0), c.getDouble(1.0, 1.0)));
		assertEquals(c.getDouble(0.0, 2.0), c.add(c.getDouble(-1.0, 1.0), c.getDouble(1.0, 1.0)));
		assertEquals(c.getDouble(1.0, 2.0), c.add(c.getDouble(.5, 1.0), c.getDouble(.5, 1.0)));
		assertEquals(c.zero(), c.add(c.getDouble(0.0, -1.0), c.getDouble(0.0, 1.0)));
	}

	@Test
	void testMultiplication() {
		Complex c = Complex.c(1024);
		assertEquals(c.zero(), c.multiply(c.one(), c.zero()));
		assertEquals(c.one(), c.multiply(c.one(), c.one()));
		assertEquals(c.i(), c.multiply(c.one(), c.i()));
		assertEquals(c.negative(c.one()), c.multiply(c.i(), c.i()));
		assertEquals(c.getDouble(-2.0, 0.0),
				c.multiply(c.getDouble(-1.0, 1.0), c.getDouble(1.0, 1.0)));
		assertEquals(c.getDouble(0.0, 2.0), c.multiply(c.getDouble(1.0, 1.0), c.getDouble(1.0, 1.0)));
		assertEquals(c.getDouble(-5.0, 5.0),
				c.multiply(c.getDouble(2.0, -1.0), c.getDouble(-3.0, 1.0)));
	}

	@Test
	void testDivision() {
		Complex c =  Complex.c(1024);
		assertEquals(c.one(), c.divide(c.one(), c.one()));
		assertEquals(c.negative(c.i()), c.divide(c.one(), c.i()));
		assertEquals(c.one(), c.divide(c.i(), c.i()));
		assertEquals(c.getDouble(-1.0, 1.0), c.divide(c.getDouble(-2.0, 0.0), c.getDouble(1.0, 1.0)));
		assertEquals(c.getDouble(1.0, 1.0), c.divide(c.getDouble(0.0, 2.0), c.getDouble(1.0, 1.0)));
		assertEquals(c.getDouble(-3.0, 1.0),
				c.divide(c.getDouble(-5.0, 5.0), c.getDouble(2.0, -1.0)));
		assertEquals(c.getDouble(2.0, -1.0),
				c.divide(c.getDouble(-5.0, 5.0), c.getDouble(-3.0, 1.0)));
	}

	@Test
	void testValue() {
		Complex c =  Complex.c(1024);
		Reals r = Reals.r(1024);
		assertEquals(r.zero(), c.value(c.zero()));
		assertEquals(r.one(), c.value(c.one()));
		assertEquals(r.getInteger(2), c.value(c.getInteger(2)));
		assertEquals(r.one(), c.value(c.i()));
		assertEquals(r.positiveSqrt(r.getInteger(2)), c.value(c.getDouble(-1.0, 1.0)));
		assertTrue(r.close(r.getInteger(5), c.value(c.getDouble(3.0, 4.0))));
		assertTrue(r.close(r.getInteger(13), c.value(c.getDouble(12.0, -5.0))));
	}

}
