package varieties.affine;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.util.Set;
import java.util.TreeSet;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import varieties.curves.elliptic.EllipticCurve;
import varieties.projective.ProjectivePoint;

class AffineSchemeTest {

	@Test
	void test() {
		PrimeField f5 = PrimeField.getPrimeField(5);
		EllipticCurve<PFE> j1728 = EllipticCurve.fromJInvariant(f5, f5.getElement(1728));
		Set<AffinePoint<PFE>> points = new TreeSet<>();
		for (ProjectivePoint<PFE> point : j1728) {
			if (!point.equals(j1728.neutral())) {
				points.add(point.getDehomogenous(3));
			}
		}
		assertEquals(j1728.getNumberOfElements().intValueExact() - 1, points.size());
		AffineScheme<PFE> affine = j1728.getAffineCover().getCover().get(2);
		int counter = 0;
		for (AffinePoint<PFE> point : affine) {
			assertTrue(points.contains(point));
			counter++;
		}
		assertEquals(j1728.getNumberOfElements().intValueExact() - 1, counter);
	}

}
