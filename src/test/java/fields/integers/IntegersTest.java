package fields.integers;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.SortedSet;
import java.util.TreeSet;

import org.junit.jupiter.api.Test;

import fields.integers.Integers.IntE;

class IntegersTest {

	@Test
	void testFactors() {
		Integers z = Integers.z();
		SortedSet<IntE> factors = new TreeSet<>();
		factors.addAll(z.factors(z.getInteger(60)));
		assertEquals(12, factors.size());
		assertTrue(factors.contains(z.getInteger(1)));
		assertTrue(factors.contains(z.getInteger(2)));
		assertTrue(factors.contains(z.getInteger(3)));
		assertTrue(factors.contains(z.getInteger(4)));
		assertTrue(factors.contains(z.getInteger(5)));
		assertTrue(factors.contains(z.getInteger(6)));
		assertTrue(factors.contains(z.getInteger(10)));
		assertTrue(factors.contains(z.getInteger(12)));
		assertTrue(factors.contains(z.getInteger(15)));
		assertTrue(factors.contains(z.getInteger(20)));
		assertTrue(factors.contains(z.getInteger(30)));
		assertTrue(factors.contains(z.getInteger(60)));
	}

}
