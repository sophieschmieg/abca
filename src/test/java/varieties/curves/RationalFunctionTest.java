package varieties.curves;

import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import varieties.RationalFunction;
import varieties.curves.elliptic.EllipticCurve;
import varieties.projective.ProjectivePoint;

class RationalFunctionTest {

	@Test
	void test() {
		PrimeField f29 = PrimeField.getPrimeField(29);
		EllipticCurve<PFE> curve = EllipticCurve.fromJInvariant(f29, f29.getInteger(1728));
		ProjectivePoint<PFE> zero = curve.getRandomElement();
		ProjectivePoint<PFE> pole1;
		ProjectivePoint<PFE> pole2;
		do {
			pole1 = curve.getRandomElement();
			pole2 = curve.getRandomElement();
		} while (pole1.equals(zero) || pole2.equals(zero));
		RationalFunction<PFE> function = curve.getRationalFunction(zero, pole1, pole2);
		System.out.println(zero);
		System.out.println(pole1);
		System.out.println(pole2);
		System.out.println(function);
		System.out.println(function.getZeroes());
		System.out.println(function.getPoles());
		assertTrue(function.getZeroes().contains(zero));
		assertTrue(function.getPoles().contains(pole1));
		assertTrue(function.getPoles().contains(pole2));
	}

}
