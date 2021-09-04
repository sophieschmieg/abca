package fields.quaternions;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import fields.integers.Rationals;

class QuaternionOrderTest {
	
	@Test
	void testLipschitzOrder() {
		Rationals q = Rationals.q();
		RationalQuaternions h = RationalQuaternions.quaternions(q.getInteger(-1), q.getInteger(-1));
		RationalQuaternionOrder lipschitz = new RationalQuaternionOrder(h, h.i(), h.j(), h.k());
		assertEquals(q.getInteger(4), q.getInteger(lipschitz.reducedDiscriminant()));
		assertFalse(lipschitz.isMaximal());
	}

	@Test
	void testHurwitzOrder() {
		Rationals q = Rationals.q();
		RationalQuaternions h = RationalQuaternions.quaternions(q.getInteger(-1), q.getInteger(-1));
		RationalQuaternionOrder hurwitz = new RationalQuaternionOrder(h, h.i(), h.j(),
				h.divide(h.getElement(q.getInteger(-1), q.one(), q.one(), q.one()), h.getInteger(2)));
		assertEquals(h.discriminant(), q.getInteger(hurwitz.reducedDiscriminant()));
		assertTrue(hurwitz.isMaximal());
	}

}
