package fields.local;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.MathMap;
import fields.local.PAdicField.PAdicNumber;

class PAdicFieldTest {

	@Test
	void testLogarithm() {
		PAdicField z2 = new PAdicField(2, 128);
		MathMap<PAdicNumber, PAdicNumber> log = new MathMap<>() {

			@Override
			public PAdicNumber evaluate(PAdicNumber t) {
				PAdicNumber result = z2.zero();
				PAdicNumber x = z2.subtract(t, z2.one());
				PAdicNumber power = z2.one();
				for (int i = 0; i <= 2 * z2.getAccuracy(); i++) {
					int sign = i % 2 == 0 ? 1 : -1;
					power = z2.multiply(x, power);
					result = z2.add(z2.divide(z2.multiply(sign, power), z2.getInteger(i + 1)), result);
				}
				return result;
			}
		};
		PAdicNumber log211 = log.evaluate(z2.getInteger(211));
		PAdicNumber log3 = log.evaluate(z2.getInteger(3));
		System.out.println(log211);
		System.out.println(log3);
		System.out.println(z2.divide(log211, log3));
		System.out.println(z2.power(z2.getInteger(3), 13));
		System.out.println(z2.getInteger(211));
	}

	@Test
	void testRationalApproximation() {
		Rationals q = Rationals.q();
		int[] testIntegers = new int[] { 3, 2, 1, -1, 5, 16, 15, 29, 101, 1024, 11, 66, 45, -65536, -3517 };
		int[] testPrimes = new int[] { 2, 3, 5, 7, 11, 65537, 257 };
		for (int prime : testPrimes) {
			PAdicField field = new PAdicField(prime, (int) Math.ceil(40 * Math.log(2) / Math.log(prime)));
			System.out.println(field);
			for (int i : testIntegers) {
				for (int j : testIntegers) {
					Fraction fraction = q.getFraction(i, j);
					System.out.println(fraction);
					PAdicNumber asPadic = field.getFraction(fraction);
					System.out.println(asPadic);
					assertEquals(fraction, field.toRational(asPadic));
				}
			}
		}
	}

}
