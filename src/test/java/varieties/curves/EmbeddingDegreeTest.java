package varieties.curves;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.interfaces.Field.Extension;
import fields.vectors.FiniteVectorSpace;
import fields.vectors.Vector;
import varieties.curves.elliptic.EllipticCurve;
import varieties.curves.elliptic.EllipticCurve.EllipticCurveExtensionResult;
import varieties.projective.ProjectivePoint;

class EmbeddingDegreeTest {

	@Test
	void test() {
		PrimeField f5 = PrimeField.getPrimeField(7);
		Extension<PFE, PFE, FFE, FiniteField> extension = f5
				.getExtension(FiniteField.getFiniteField(7, 3).minimalPolynomial());
		for (Vector<PFE> coeff : new FiniteVectorSpace<>(f5, 2)) {
			if (f5.add(f5.multiply(4, f5.power(coeff.get(1), 3)), f5.multiply(27, f5.power(coeff.get(2), 2)))
					.equals(f5.zero())) {
				continue;
			}
			//coeff = new Vector<>(f5.one(), f5.getInteger(3));
			EllipticCurve<PFE> curve = new EllipticCurve<>(f5, coeff.get(1), coeff.get(2));
			if (curve.getTorsionPoints(3).size() == 1) {
				continue;
			}
			System.out.println();
			System.out.println(curve);
			int counter = 0;
			for (ProjectivePoint<PFE> point : curve) {
				counter++;
			}
			System.out.println(counter);
			System.out.println(curve.embeddingDegree());
			System.out.println("Number of 3 torison: " + curve.getTorsionPoints(3).size());
		//	System.out.println("Number of 9 torison: " + curve.getTorsionPoints(9).size());
			System.out.println("Number of Elements: " + curve.getNumberOfElements());
			// System.out.println("Number of 27 torison: " +
			// curve.getTorsionPoints(27).size());
			EllipticCurveExtensionResult<PFE, PFE, FFE, FiniteField> extended = curve.extendBaseField(extension);
			System.out.println("Number of extension Elements: " + extended.getCurve().getNumberOfElements());
			System.out.println("Number of extension 3 torison: " + extended.getCurve().getTorsionPoints(3).size());
		//	System.out.println("Number of extension 9 torison: " + extended.getCurve().getTorsionPoints(9).size());
			// System.out.println("Number of extension 27 torison: " +
			// extended.getCurve().getTorsionPoints(27).size());
		}
	}

}
