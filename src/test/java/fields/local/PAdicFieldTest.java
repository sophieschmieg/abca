package fields.local;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.local.PAdicField.PAdicNumber;

class PAdicFieldTest {

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
