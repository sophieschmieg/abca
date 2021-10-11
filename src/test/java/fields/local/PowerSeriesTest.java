package fields.local;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.util.Random;

import org.junit.jupiter.api.Test;

import fields.integers.Rationals;
import fields.integers.Rationals.Fraction;
import fields.interfaces.Polynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.local.FormalPowerSeries.PowerSeries;
import fields.polynomials.LocalizedCoordinateRing;
import fields.polynomials.LocalizedCoordinateRing.LocalizedElement;

class PowerSeriesTest {

	@Test
	void testFromRationalFunction() {
		Rationals q = Rationals.q();
		int degree = 10;
		for (int i = 0; i < 2; i++) {
			FormalPowerSeries<Fraction> powerSeries = new FormalPowerSeries<>(q, 2 * degree + i);
			UnivariatePolynomialRing<Fraction> polynomialRing = q.getUnivariatePolynomialRing();
			LocalizedCoordinateRing<Fraction> localizedRing = powerSeries.exact().getField();
			for (int tc = 0; tc < 100; tc++) {
				Polynomial<Fraction> numerator = polynomialRing.getRandomElement(new Random().nextInt(degree));
				Polynomial<Fraction> denominator;
				do {
					denominator = polynomialRing.getRandomElement(new Random().nextInt(degree + 1));
				} while (denominator.equals(polynomialRing.zero()));
				denominator = polynomialRing.normalize(denominator);
				LocalizedElement<Fraction> rationalFunction = localizedRing.getEmbedding(numerator, denominator);
				PowerSeries<Fraction> asPowerSeries = powerSeries.fromRationalFunction(rationalFunction);
				LocalizedElement<Fraction> fromPowerSeries = powerSeries.roundToRationalFunction(asPowerSeries);
				assertEquals(rationalFunction, fromPowerSeries);
			}
		}
	}

}
