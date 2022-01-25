package varieties.curves;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.math.BigInteger;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import varieties.curves.elliptic.EllipticCurve;

class TorsionPointTest {

	@Test
	void testSupersingular() {
		FiniteField fq = FiniteField.getFiniteField(23 * 23);
		EllipticCurve<FFE> curve = EllipticCurve.fromJInvariant(fq, fq.getInteger(1728));
		assertEquals(BigInteger.valueOf(24 * 24), curve.getNumberOfElements());
		assertEquals(2, curve.getTorsionPointBasis(2).size());
		assertEquals(2, curve.getTorsionPointBasis(4).size());
		assertEquals(2, curve.getTorsionPointBasis(8).size());
		assertEquals(0, curve.getTorsionPointBasis(16).size());
		assertEquals(2, curve.getTorsionPointBasis(3).size());
		assertEquals(0, curve.getTorsionPointBasis(9).size());
		assertEquals(2, curve.getTorsionPointBasis(6).size());
		assertEquals(2, curve.getTorsionPointBasis(12).size());
		assertEquals(2, curve.getTorsionPointBasis(24).size());
		assertEquals(0, curve.getTorsionPointBasis(23).size());
		assertEquals(0, curve.getTorsionPointBasis(15).size());
	}

	@Test
	void testOrdinary() {
		FiniteField fq = FiniteField.getFiniteField(625);
		EllipticCurve<FFE> curve = EllipticCurve.fromJInvariant(fq, fq.getInteger(1728));
		assertEquals(BigInteger.valueOf(640), curve.getNumberOfElements());
		assertEquals(1, curve.getTorsionPointBasis(5).size());
		assertEquals(2, curve.getTorsionPointBasis(2).size());
		assertEquals(2, curve.getTorsionPointBasis(4).size());
		assertEquals(2, curve.getTorsionPointBasis(8).size());
		assertEquals(1, curve.getTorsionPointBasis(16).size());
		assertEquals(0, curve.getTorsionPointBasis(32).size());
		assertEquals(0, curve.getTorsionPointBasis(64).size());
		assertEquals(0, curve.getTorsionPointBasis(128).size());
		assertEquals(0, curve.getTorsionPointBasis(256).size());
		assertEquals(0, curve.getTorsionPointBasis(512).size());
		assertEquals(1, curve.getTorsionPointBasis(10).size());
		assertEquals(0, curve.getTorsionPointBasis(6).size());
	}

}
