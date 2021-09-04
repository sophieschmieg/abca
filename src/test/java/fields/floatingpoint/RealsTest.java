package fields.floatingpoint;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import fields.floatingpoint.Reals.Real;

class RealsTest {

	private boolean closeEnough(Reals r, Real t1, Real t2) {
		return r.abs(r.subtract(t1, t2)).compareTo(r.getDouble(1e-100)) < 0;
	}
	
	@Test
	void testRound() {
		Reals r = Reals.r(1024);
		assertEquals(r.getDouble(0.0625), r.roundToFixedPoint(r.getDouble(0.09375), 4));
	}
	
	@Test
	void testToString() {
		Reals r = Reals.r(1024);
		assertEquals("1.5", r.getDouble(1.5).toString());
		assertEquals("0.5", r.getDouble(.5).toString());
		assertEquals("0.0625", r.getDouble(.0625).toString());
		assertEquals("0.09375", r.getDouble(.09375).toString());
		assertTrue(r.close(r.getDouble(1.5), r.fromString("1.5")));
		assertTrue(r.close(r.getDouble(0.09375), r.fromString("0.09375")));
		assertTrue(r.close(r.getDouble(-1.5), r.fromString("-1.5")));
		assertTrue(r.close(r.getDouble(-0.09375), r.fromString("-0.09375")));
		assertTrue(r.close(r.getDouble(-0.09375), r.fromString("-.09375")));
		assertTrue(r.close(r.getDouble(0.09375), r.fromString(".09375")));
	}

	@Test
	void testAddition() {
		Reals r = Reals.r(1024);
		assertEquals(r.one(), r.add(r.zero(), r.one()));
		assertEquals(r.getDouble(-3.0), r.add(r.getInteger(-4), r.one()));
		assertEquals(r.one(), r.add(r.getDouble(.5), r.getDouble(.5)));
		assertEquals(r.getInteger(16), r.add(r.getInteger(7), r.getInteger(9)));
		assertEquals(r.getInteger(4), r.subtract(r.getInteger(12), r.getInteger(8)));
	}

	@Test
	void testMultiplication() {
		Reals r = Reals.r(1024);
		assertEquals(r.getDouble(6.0), r.multiply(r.getInteger(2), r.getInteger(3)));
		assertEquals(r.one(), r.multiply(r.getInteger(2), r.getDouble(.5)));
	}

	@Test
	void testInversion() {
		Reals r = Reals.r(1024);
		for (int tc = 0; tc < 10; tc++) {
			Real n;
			do {
				n = r.getRandomElement();
			} while (n.equals(r.zero()));
			assertTrue(r.close(r.one(), r.multiply(n, r.inverse(n))));
		}
	}

	@Test
	void testPi() {
		Reals r = Reals.r(1024);
		Real pi = r.fromString("3.14159265358979323846264338327950288419716939937510582097494459"
				+ "23078164062862089986280348253421170679821480865132823066470938446095505822317"
				+ "25359408128481117450284102701938521105559644622948954930381964428810975665933"
				+ "44612847564823378678316527120190914564856692346034861045432664821339360726024"
				+ "91412737245870066063155881748815209209628292540917153643678925903600113305305"
				+ "48820466521384146951941511609433057270365759591953092186117381932611793105118"
				+ "54807446237996274956735188575272489122793818301194912");
		assertTrue(r.close(pi, r.pi()));
	}

	@Test
	void testE() {
		Reals r = Reals.r(1024);
		Real e = r.fromString("2.718281828459045235360287471352662497757247093699959574966967627"
				+ "72407663035354759457138217852516642742746639193200305992181741359662904357290"
				+ "03342952605956307381323286279434907632338298807531952510190115738341879307021"
				+ "54089149934884167509244761460668082264800168477411853742345442437107539077744"
				+ "99206955170276183860626133138458300075204493382656029760673711320070932870912"
				+ "74437470472306969772093101416928368190255151086574637721112523897844250569536"
				+ "9677078544996996794686445490598793163688923009879312");
		assertTrue(r.close(e, r.e()));
	}
	
	@Test
	void testExpLog() {
		Reals r = Reals.r(1024);
		for (int tc = 0; tc < 10; tc++) {
			Real test = r.getRandomElement();
			assertTrue(closeEnough(r, test, r.log(r.exp(test))));
		}
	}
	
	@Test
	void testExpLogBinarySearch() {
		Reals r = Reals.r(1024);
		for (int tc = 0; tc < 10; tc++) {
			Real test = r.getRandomElement();
			assertTrue(closeEnough(r, test, r.invertMonotonic(r.exp(), r.exp(test))));
		}
	}
	
	
	@Test
	void testArctan() {
		Reals r = Reals.r(1024);
		assertTrue(r.close(r.zero(), r.arctan(r.getInteger(0))));
		assertTrue(closeEnough(r, r.divide(r.pi(), r.getInteger(12)), r.arctan(r.subtract(r.getInteger(2), r.positiveSqrt(r.getInteger(3))))));
		assertTrue(closeEnough(r, r.divide(r.pi(), r.getInteger(8)), r.arctan(r.subtract(r.positiveSqrt(r.getInteger(2)), r.one()))));
		assertTrue(closeEnough(r, r.divide(r.pi(), r.getInteger(4)), r.arctan(r.one())));
		assertTrue(closeEnough(r, r.divide(r.pi(), r.getInteger(3)), r.arctan(r.positiveSqrt(r.getInteger(3)))));
		assertTrue(closeEnough(r, r.divide(r.multiply(3, r.pi()), r.getInteger(8)), r.arctan(r.add(r.positiveSqrt(r.getInteger(2)), r.one()))));
		assertTrue(closeEnough(r, r.divide(r.multiply(-3, r.pi()), r.getInteger(8)), r.arctan(r.negative(r.add(r.positiveSqrt(r.getInteger(2)), r.one())))));
		assertTrue(closeEnough(r, r.divide(r.multiply(-1, r.pi()), r.getInteger(3)), r.arctan(r.negative(r.positiveSqrt(r.getInteger(3))))));
		assertTrue(closeEnough(r, r.divide(r.multiply(-1, r.pi()), r.getInteger(8)), r.arctan(r.subtract(r.one(), r.positiveSqrt(r.getInteger(2))))));
		}
}
